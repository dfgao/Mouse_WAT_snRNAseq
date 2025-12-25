##' clusterData_new: Improve the "clusterData" function to simplify its input
##'
##' 支持三种输入类型：Monocle3 的 cell_data_set, matrix, data.frame。
##' 通过 mfuzz、TCseq、kmeans、WGCNA 等方法聚类.
##'
##' @param data 一个 cell_data_set 对象或表达矩阵（matrix 或 data.frame）
##' @param scaleData 是否对数据进行标准化（mfuzz/TCseq）
##' @param method 聚类方法: "mfuzz", "TCseq", "kmeans", "wgcna"
##' @param TCseq_paramsList TCseq 额外参数列表
##' @param min.std mfuzz/TCseq 前对基因过滤的最小标准差
##' @param cluster.num 聚类簇数
##' @param subcluster 子簇列表，用于结果过滤
##' @param seed 随机数种子
##' @param object WGCNA 网络对象
##' @param ... 其它参数（传递给 pre_pseudotime_matrix 等）
##' @return 列表，包括 wide.res, long.res, cluster.list, method
clusterData_new <- function(
    data,
    scaleData      = TRUE,
    method         = c("mfuzz", "TCseq", "kmeans", "wgcna"),
    TCseq_paramsList = list(),
    min.std        = 0,
    cluster.num    = NULL,
    subcluster     = NULL,
    seed           = 5201314,
    object         = NULL,
    ...
) {
  # 统一选项
  method <- match.arg(method)
  
  # 准备表达矩阵
  if (inherits(data, "cell_data_set")) {
    # Monocle3 转矩阵
    input_mat <- pre_pseudotime_matrix(cds_obj = data, assays = "counts", ...)
  } else if (is.matrix(data) || is.data.frame(data)) {
    input_mat <- as.matrix(data)
  } else {
    stop("data must be a cell_data_set, matrix, or data.frame")
  }
  
  # 根据方法执行不同聚类
  if (method %in% c("mfuzz", "TCseq")) {
    ## mfuzz or TCseq
    eset <- Biobase::ExpressionSet(assayData = input_mat)
    ## 过滤低方差基因
    eset <- Mfuzz::filter.std(eset, min.std = min.std, visu = FALSE)
    if (scaleData) eset <- Mfuzz::standardise(eset)
    
    if (method == "mfuzz") {
      set.seed(seed)
      m_val <- Mfuzz::mestimate(eset)
      fuzz_res <- Mfuzz::mfuzz(eset, c = cluster.num, m = m_val)
      exprs_mat <- Biobase::exprs(eset)
      raw_anno <- data.frame(exprs_mat, cluster = fuzz_res$cluster)
      mem_df <- data.frame(fuzz_res$membership, cluster = fuzz_res$cluster)
      # 合并 membership 和 exprs
      long_df <- merge(
        reshape2::melt(raw_anno, varnames = c("gene","cell")),
        reshape2::melt(mem_df, varnames = c("gene","cluster")),
        by = c("gene"), suffixes = c(".expr",".mem")
      )
      final_wide <- raw_anno
    } else {
      # TCseq
      tca <- do.call(
        TCseq::timeclust,
        modifyList(list(x = input_mat, algo = "cm", k = cluster.num, standardize = scaleData), TCseq_paramsList)
      )
      dt <- data.frame(tca@data, gene = rownames(tca@data))
      clu <- data.frame(cluster = tca@cluster, gene = names(tca@cluster))
      mem <- reshape2::melt(data.frame(tca@membership), varnames = c("gene","cluster"), value.name = "membership")
      final_wide <- merge(dt, merge(clu, mem, by = c("gene","cluster")), by = "gene")
      long_df <- reshape2::melt(final_wide, id.vars = c("gene","cluster","membership"), variable.name = "cell_type", value.name = "expr")
    }
    
  } else if (method == "kmeans") {
    # kmeans
    mat_filt <- input_mat[rowSds(input_mat) >= min.std, , drop = FALSE]
    mat_scale <- if (scaleData) t(scale(t(mat_filt))) else mat_filt
    set.seed(seed)
    km <- kmeans(mat_scale, centers = cluster.num, nstart = 10)
    final_wide <- data.frame(mat_filt, gene = rownames(mat_filt), cluster = km$cluster)
    long_df <- reshape2::melt(final_wide, id.vars = c("gene","cluster"), variable.name = "cell_type", value.name = "expr")
    
  } else if (method == "wgcna") {
    # WGCNA
    if (is.null(object)) stop("WGCNA method requires 'object' input")
    colors <- object$colors + 1
    module_colors <- WGCNA::labels2colors(object$colors)
    expr_scaled <- t(scale(t(input_mat)))
    final_wide <- data.frame(expr_scaled, gene = rownames(input_mat), cluster = colors, modulecol = module_colors)
    long_df <- reshape2::melt(final_wide, id.vars = c("gene","cluster","modulecol"), variable.name = "cell_type", value.name = "expr")
    
  } else {
    stop("Unsupported method: ", method)
  }
  
  # 子簇过滤
  if (!is.null(subcluster)) {
    final_wide <- final_wide[final_wide$cluster %in% subcluster, , drop = FALSE]
    long_df <- long_df[long_df$cluster %in% subcluster, , drop = FALSE]
  }
  
  # 构建输出
  cluster_list <- split(final_wide$gene, final_wide$cluster)
  return(list(
    wide.res   = final_wide,
    long.res   = long_df,
    cluster.list = cluster_list,
    type       = method,         
    geneMode   = "none",        
    geneType   = "none",      
    pseudotime = NULL        
  ))
}

# 08. sankey plot ----
plot_cell_fraction_sankey <- function(obj, condition, cell.type, cols) {
  library(dplyr); library(tidyr); library(ggplot2); library(ggalluvial)
  
  # 抓数据
  df <- FetchData(obj, vars = c(condition, cell.type))
  colnames(df) <- c("condition", "celltype")
  
  # 补全 + 计算
  plot.data <- df %>%
    dplyr::count(condition, celltype, name = "n") %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(Prop = n / sum(n) * 100) %>%
    ungroup() %>%
    tidyr::complete(condition, celltype, fill = list(n = 0, Prop = 0))
  
  # 生成 alluvium id
  n_celltype <- length(unique(plot.data$celltype))
  n_condition <- length(unique(plot.data$condition))
  plot.data$alluvium <- rep(seq_len(n_celltype), times = n_condition)
  
  # 画图
  ggplot(plot.data,
         aes(x = condition, stratum = celltype, alluvium = alluvium,
             y = Prop, fill = celltype, label = celltype)) +
    scale_fill_manual(values = cols) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              width = 0.3, color = "darkgray") +
    geom_stratum(alpha = .8) +
    guides(fill = guide_legend(title = 'Cell type')) +
    ylab('Proportion (%)') +
    theme_bw() +
    theme(legend.position = "right",
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15))
}


# 08. sankey plot ----
plot_cell_fraction_sankey <- function(obj, condition, cell.type, cols) {
  library(dplyr); library(tidyr); library(ggplot2); library(ggalluvial)
  
  # 抓数据
  df <- FetchData(obj, vars = c(condition, cell.type))
  colnames(df) <- c("condition", "celltype")
  
  # 补全 + 计算
  plot.data <- df %>%
    dplyr::count(condition, celltype, name = "n") %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(Prop = n / sum(n) * 100) %>%
    ungroup() %>%
    tidyr::complete(condition, celltype, fill = list(n = 0, Prop = 0))
  
  # 生成 alluvium id
  n_celltype <- length(unique(plot.data$celltype))
  n_condition <- length(unique(plot.data$condition))
  plot.data$alluvium <- rep(seq_len(n_celltype), times = n_condition)
  
  # 画图
  ggplot(plot.data,
         aes(x = condition, stratum = celltype, alluvium = alluvium,
             y = Prop, fill = celltype, label = celltype)) +
    scale_fill_manual(values = cols) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              width = 0.3, color = "darkgray") +
    geom_stratum(alpha = .8) +
    guides(fill = guide_legend(title = 'Cell type')) +
    ylab('Proportion (%)') +
    theme_bw() +
    theme(legend.position = "right",
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15))
}


# 02. LTSR ----

plot_ranef.new <- function(ranef_tbl, vars, celltypes = NULL, celltype_order = "hclust", references = NULL,
                           maxFC = 3, LTSR2p = F, highlightLtsr = 0.0, filterLtsr = 0.0, swap_axes = F) {
  ranef_tbl <- .getCondValLtsr(ranef_tbl, vars, celltypes = celltypes, references = references)
  ranef_tbl <- ranef_tbl %>%
    dplyr::mutate(
      grpval = factor(grpval, levels = vars[[1]])
    )
  save(ranef_tbl, file = "ranef_tbl.RData")
  condval_mat <- ranef_tbl %>%
    dplyr::select(
      Celltype, grpval, condval
    ) %>%
    spread(
      "grpval", "condval"
    ) %>%
    column_to_rownames(
      var = "Celltype"
    ) %>%
    as.matrix()
  if (length(celltype_order) == 1 && celltype_order == "hclust") {
    dendy <- hclust(dist(condval_mat))
    ordered_celltype <- rownames(condval_mat)[dendy$ord]
  } else if (!is.null(celltype_order) && length(celltype_order) == dim(condval_mat)[1]) {
    ordered_celltype <- celltype_order
  }
  
  ranef_tbl <- ranef_tbl %>% mutate(
    Celltype = factor(Celltype, levels = ordered_celltype),
    grpval   = factor(grpval, levels = levels(ranef_tbl$grpval)),
    condval = condval %>% pmin(log(maxFC)) %>% pmax(log(1 / maxFC)),
    ltsr = ltsr %>% pmin(0.9999) %>% pmax(0.5)
  )
  
  # ranef_tbl <- ranef_tbl %>% mutate(
  #   Celltype = factor(Celltype, levels = ordered_celltype),
  #      grpval   = ifelse(grpvar=="Treatment",
  #                        factor(grpval, levels = vars$Treatment),
  #                        factor(grpval, levels = unique(grpval))),
  #   condval   = condval %>% pmin(log(maxFC)) %>% pmax(log(1/maxFC)),
  #   ltsr      = ltsr %>% pmin(0.9999) %>% pmax(0.5)
  # )
  
  if (swap_axes) {
    ranef_tbl$Celltype <- factor(ranef_tbl$Celltype, levels = rev(levels(ranef_tbl$Celltype)))
    ranef_tbl$grpval <- factor(ranef_tbl$grpval, levels = rev(levels(ranef_tbl$grpval)))
  }
  
  if (filterLtsr > 0) {
    filtered_celltypes <- ranef_tbl %>%
      group_by(Celltype) %>%
      summarise(maxLtsr = max(ltsr)) %>%
      dplyr::filter(maxLtsr >= filterLtsr) %>%
      dplyr::select(Celltype) %>%
      unlist(use.names = F)
    ranef_tbl <- ranef_tbl %>% dplyr::filter(Celltype %in% filtered_celltypes)
  }
  
  geom_dots <- geom_point(
    aes(
      fill = log2(exp(condval)),
      size = -log10(1 - ltsr)
    ),
    color = "white",
    shape = 21
  )
  
  if (swap_axes) {
    p <- (
      ggplot(ranef_tbl, aes(y = grpval, x = Celltype)) +
        facet_grid(grpvar ~ ., scales = "free_y", space = "free_y", switch = "x") +
        geom_dots
    )
  } else {
    p <- (
      ggplot(ranef_tbl, aes(x = grpval, y = Celltype)) +
        facet_grid(. ~ grpvar, scales = "free_x", space = "free_x", switch = "x") +
        geom_dots
    )
  }
  
  p <- (
    p + scale_fill_distiller(
      palette = "RdBu",
      limits = log2(c(1 / maxFC, maxFC)),
      breaks = log2(c(1 / maxFC, maxFC)),
      labels = c(paste0("1/", maxFC), maxFC),
      oob = squish,
      guide = guide_colorbar(
        title = "Fold change", title.position = "top", direction = "horizontal",
        barwidth = 5, barheight = 0.75, raster = F, order = 1)
    )
    + scale_size(
      limits = -log10(1 - c(0.5, 0.9999)),
      breaks = -log10(1 - c(0.5, 0.7, 0.9, 0.99)),
      range = c(0.5, 9),
      labels = ifelse(
        rep(LTSR2p, 4),
        c("0.5", "0.3", "0.1", "<0.01"),
        c("0.5", "0.7", "0.9", ">0.99")
      ),
      guide = guide_legend(
        title = ifelse(LTSR2p, "p", "LTSR"), reverse = T, order = 2,
        override.aes = list(fill = "black", color = "white")
      )
    )
    + theme_bw()
    + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      legend.spacing.y = unit(0.5, "line")
    )
  )
  
  if (highlightLtsr > 0) {
    p <- (
      p + geom_point(
        aes(
          color = ltsr > highlightLtsr,
          alpha = ltsr > highlightLtsr
        ),
        shape = 21, size = 9
      )
      + scale_color_manual(
        label = c(
          "", ifelse(LTSR2p, paste0("p < ", 1 - highlightLtsr), paste0("LTSR > ", highlightLtsr))
        ),
        values = c("white", "red"),
        guide = guide_legend(title = NULL, override.aes = list(size = 9), reverse = T, order = 3)
      )
      + scale_alpha_manual(values = c(0, 1), guide = F)
    )
  }
  
  p
}

run_celltype_composition <- function(
    seu,
    sample_col        = "sample",
    celltype_col      = "cell.type.percise.new",
    treatment_levels  = c('A','B'),
    reps              = c("one","two"),
    ltsr_vars         = c("Treatment","Rep"),
    FC                = 2,
    var_num1_name     = "Var_Num1"
) {
  
  levels(seu[[sample_col]]) <- paste0(
    rep(treatment_levels, each = length(reps)), "_", reps
  )
  
  cell_number <- FetchData(seu, vars = c(sample_col, celltype_col)) %>%
    as_tibble() %>%
    dplyr::count(!!sym(sample_col), !!sym(celltype_col)) %>%
    pivot_wider(
      names_from  = !!sym(sample_col),
      values_from = n,
      values_fill = 0
    )
  
  sample_ids         <- colnames(cell_number)[-1]
  cell_types         <- cell_number[[celltype_col]]
  n_cells_per_sample <- colSums(cell_number[,-1])
  
  sample_meta <- tibble(
    Sample_ID = sample_ids,
    Treatment = rep(treatment_levels, each = length(reps)),
    Rep       = rep(reps, times = length(treatment_levels)),
    cell.num  = n_cells_per_sample
  )
  
  obs_tbl <- data.frame(
    Sample_ID = rep(sample_ids, times = n_cells_per_sample),
    Treatment = rep(sample_meta$Treatment, times = n_cells_per_sample),
    Rep       = rep(sample_meta$Rep, times = n_cells_per_sample),
    Var_Num1  = rep(1, sum(n_cells_per_sample))
  )
  
  obs_tbl$Cell_type <- unlist(
    lapply(sample_ids, function(sid) {
      rep(cell_number[[celltype_col]], times = cell_number[[sid]])
    })
  )
  
  results <- CellTypeCompositionAnalysis(
    obs_tbl,
    colSample   = "Sample_ID",
    colCelltype = "Cell_type",
    colVarCats  = c(unlist(ltsr_vars)),
    colVarNums  = var_num1_name
  )
  ranef_tbl <- results$ranef
  sdse_tbl  <- results$sdse
  
  ranef_plot <- plot_ranef.new(
    ranef_tbl,
    vars           = list(Treatment = treatment_levels),
    celltypes      = cell_types,
    celltype_order = rev(cell_types),
    maxFC          = FC,
    LTSR2p         = FALSE
  ) +
    xlab("Condition") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  sdse_plot <- plot_sdse(
    sdse_tbl,
    colSample = "Sample_ID",
    ci                = 0.95,
    xlim              = c(0,1)
  )
  
  list(
    sample_meta = sample_meta,
    obs_tbl     = obs_tbl,
    ranef_tbl   = ranef_tbl,
    sdse_tbl    = sdse_tbl,
    ranef_plot  = ranef_plot,
    sdse_plot   = sdse_plot,
    combined    = ranef_plot + sdse_plot
  )
}


