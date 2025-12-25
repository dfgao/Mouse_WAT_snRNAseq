
############################################################################## scRNA-seq section ########################################################################################

# 01.umap plot ----
Clean_sct.inte.rm.edc$condition <- factor(Clean_sct.inte.rm.edc$condition, levels = c('LP','HP'))
## remove LpTs
p1 <- DimPlot(Clean_sct.inte.rm.edc, group.by = 'cell.type.percise.new', cols = c( use.cols,npg.cols), label = T) + ggtitle("") + NoAxes()
p2 <- DimPlot(Clean_sct.inte.rm.edc, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols), split.by = 'condition') + ggtitle("") + NoAxes() + NoLegend()
p1 + p2
Clean_sct.inte.rm.edc$cell.type.percise.new <- factor(Clean_sct.inte.rm.edc$cell.type.percise.new, levels = c('Adipocytes',
                                                                                                              'FAPs',
                                                                                                              'MCs',
                                                                                                              'ECs',
                                                                                                              'LECs',
                                                                                                              'Pericytes',
                                                                                                              'Schwann cells',
                                                                                                              'Epididymal cells',
                                                                                                              'Macrophages',
                                                                                                              'DCs',
                                                                                                              'Monocytes',
                                                                                                              'B cells',
                                                                                                              'T cells'))
table(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'LP',]$cell.type.percise.new)/nrow(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'LP',])
table(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'HP',]$cell.type.percise.new)/nrow(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'HP',])

## cirlize umap
library(plot1cell)

plot.cir.test <- function (data_plot, do.label = T, contour.levels = c(0.2, 0.4, 0.6), 
                           pt.size = 0.5, kde2d.n = 1000, contour.nlevels = 100, bg.color = "#F9F2E4", 
                           col.use = NULL, label.cex = 0.5, repel = FALSE) 
{
  centers <- data_plot %>% dplyr::group_by(Cluster) %>% summarise(x = median(x = x), 
                                                                  y = median(x = y))
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  celltypes <- names(table(data_plot$Cluster))
  cell_colors <- (scales::hue_pal())(length(celltypes))
  if (!is.null(col.use)) {
    cell_colors = col.use
    col_df <- data.frame(Cluster = celltypes, color2 = col.use)
    cells_order <- rownames(data_plot)
    data_plot <- merge(data_plot, col_df, by = "Cluster")
    rownames(data_plot) <- data_plot$cells
    data_plot <- data_plot[cells_order, ]
    data_plot$Colors <- data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0), track.height = 0.01, gap.degree = c(rep(2, (length(celltypes) -  1)), 12), points.overflow.warning = FALSE)
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + 
                               mm_y(4), CELL_META$sector.index, cex = 0.5, 
                             col = "black", facing = "bending.inside", niceFacing = T)
                 breaks = seq(0, 100, by = 50)
                 circos.axis(labels.cex = 0.3, col = "black", labels.col = "black",major.at = breaks,
                             labels = paste0(breaks, "%"))
               })
  for (i in 1:length(celltypes)) {
    dd <- data_plot[data_plot$Cluster == celltypes[i], ]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), 
                    y1 = 0, col = cell_colors[i], lwd = 3, sector.index = celltypes[i])
  }
  text(x = 1, y = 0.1, labels = "Cluster", cex = 0.4, col = "black", 
       srt = -90)
  points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$Colors, 
                                                         0.2), cex = pt.size)
  contour(z, drawlabels = F, nlevels = 100, levels = contour.levels, 
          col = "#ae9c76", add = TRUE)
  if (do.label) {
    if (repel) {
      textplot(x = centers$x, y = centers$y, words = centers$Cluster, 
               cex = label.cex, new = F, show.lines = F)
    }
    else {
      text(centers$x, centers$y, labels = centers$Cluster, 
           cex = label.cex, col = "black")
    }
  }
}

cell_order_test <- function (dat) 
{
  celltypes <- names(table(dat$Cluster))
  new_dat <- list()
  for (i in 1:length(celltypes)) {
    dat$Cluster <- as.character(dat$Cluster)
    dat1 <- dat[dat$Cluster == celltypes[i], ]
    dat1$x_polar <- (1:nrow(dat1))/nrow(dat1)
    new_dat[[i]] <- dat1
  }
  new_dat <- do.call("rbind", new_dat)
  new_dat
}

prepare_circlize_data_test <- function (seu_obj, scale = 1) 
{
  celltypes <- levels(seu_obj)
  cell_colors <- (scales::hue_pal())(length(celltypes))
  data_plot <- get_metadata(seu_obj, color = cell_colors, 
                            coord_scale = scale)
  data_plot <- cell_order_test(data_plot)
  data_plot$x_polar2 <- log10(data_plot$x_polar)
  data_plot
}

Idents(Clean_sct.inte.rm.edc) <- Clean_sct.inte.rm.edc$cell.type.percise.new
circ_data <- prepare_circlize_data_test(Clean_sct.inte.rm.edc, scale = 0.7 )
circ_data$Cluster <- factor(circ_data$Cluster, levels = c(levels(Clean_sct.inte.rm.edc$cell.type.percise.new)))
set.seed(1234)
circ_data$x_polar2 <- circ_data$x_polar*100

cluster_colors<- c(use.cols,npg.cols)[1:13]
con_colors <- c('#0B996F','#D6570D') # untreated & treated
circ_data$sample <- factor(circ_data$sample, levels = c('LP1','LP2','HP1','HP2'))
samp_colors <-  c('#fb8500','#f47068','#8ecae6', "#219ebc")

plot.cir.test(circ_data,do.label = T, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 1)

add_track(circ_data, group = "condition", colors = con_colors, track_num = 2) 
add_track(circ_data, group = "sample",colors = samp_colors, track_num = 4) 

# 02.miloR:composition of neighbourhood of cells ----
library(miloR)
library(scales)
library(SingleCellExperiment)
library(ggbeeswarm)

milo.all.seu <- subset(Clean_sct.inte.rm.edc, cells = c(sample( 
  rownames(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'LP',]), 16471),
  rownames(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'HP',] )))
milo.all.seu <- as.SingleCellExperiment(milo.all.seu,assay = 'SCT') %>%
  Milo() %>%
  miloR::buildGraph(k = 30, d = 50) %>% 
  makeNhoods(
    prop = 0.2,                 
    k = 30,         
    d=50,                   
    refined = T)

milo.all.seu <- countCells(milo.all.seu, 
                           meta.data = data.frame(colData(milo.all.seu)),                     
                           sample="sample")
milo.traj_design <- data.frame(colData(milo.all.seu))[,c("sample", "condition")]#分别是重复样本ID和分组
milo.traj_design$sample <- as.factor(milo.traj_design$sample)
milo.traj_design <- distinct(milo.traj_design)
rownames(milo.traj_design) <- milo.traj_design$sample

milo.all.seu <- calcNhoodDistance(milo.all.seu, d=50)
milo.da_results <- testNhoods(milo.all.seu,                 
                              design = ~ condition,            
                              design.df = milo.traj_design)
milo.all.seu <- buildNhoodGraph(milo.all.seu)

ggplot(milo.da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(milo.da_results, aes(logFC, -log10(SpatialFDR))) +  
  geom_point() + 
  geom_hline(yintercept = 1)
scater::plotReducedDim(milo.all.seu, dimred = "UMAP", colour_by="cell.type.percise.new", text_size = 3, point_size=0.1)

milo.da_results$logFC <- -milo.da_results$logFC
miloR::plotNhoodGraphDA(milo.all.seu, milo.da_results, alpha=0.1,) +  
  scale_fill_gradient2(low="#F2F2F2",#修改颜色            
                       mid="#F2F2F2",
                       high="#FF5E5B",                  
                       name="log2FC",                   
                       limits=c(-5,5),                  
                       oob=squish)

milo.da_results <- annotateNhoods(milo.all.seu, milo.da_results, coldata_col = "cell.type.percise.new")
milo.da_results$cell.type.percise.new <- factor(milo.da_results$cell.type.percise.new, levels = rev(c('Adipocytes',
                                                                                                              'FAPs',
                                                                                                              'MCs',
                                                                                                              'ECs',
                                                                                                              'LECs',
                                                                                                              'Pericytes',
                                                                                                              'Schwann cells',
                                                                                                              'Epididymal cells',
                                                                                                              'Macrophages',
                                                                                                              'DCs',
                                                                                                              'Monocytes',
                                                                                                              'B cells',
                                                                                                              'T cells')))

plotbee <- function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL) 
{
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
           group.by, ")?")
    }
    if (is.numeric(da.res[, group.by])) {
    }
    da.res <- mutate(da.res, group_by = da.res[, group.by])
  }
  else {
    da.res <- mutate(da.res, group_by = "g1")
  }
  if (!is.factor(da.res[, "group_by"])) {
    messHP("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, 
                                               levels = unique(group_by)))
  }
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }
  beeswarm_pos <- ggplot_build(da.res %>% mutate(is_signif = ifelse(SpatialFDR < 
                                                                      alpha, 1, 0)) %>% arrange(group_by) %>% ggplot(aes(group_by, 
                                                                                                                         logFC)) + geom_quasirandom())
  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y
  n_groups <- unique(da.res$group_by) %>% length()
  da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 
                                       1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 
                                                                                1, logFC, NA)) %>% arrange(group_by) %>% mutate(Nhood = factor(Nhood, 
                                                                                                                                               levels = unique(Nhood))) %>% mutate(pos_x = pos_x, pos_y = pos_y) %>% 
    ggplot(aes(pos_x, pos_y, color = logFC_color)) + scale_color_gradient2() + 
    guides(color = "none") + xlab(group.by) + ylab("Log Fold Change") + 
    scale_x_continuous(breaks = seq(1, n_groups), labels = setNames(levels(da.res$group_by), 
                                                                    seq(1, n_groups))) + geom_point(size = .5) + coord_flip() + 
    theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0))
}

plotbee (milo.da_results, alpha = 0.1,
         group.by = "cell.type.percise.new") +  
  scale_color_gradient2(low="#4cc9f0",           
                        mid="#F2F2F2",                  
                        high="#FF5E5B",                  
                        limits=c(-5,5),                  
                        oob=squish) + 
  labs(x="", y="Log2 Fold Change") +  
  theme_bw(base_size=10)+  
  theme(axis.text = element_text(colour = 'black')) + 
  scale_size_continuous(range = .1)

# 03.ltsr composition of cell types ----
library(magrittr)
library(lme4)
library(numDeriv)
## data perpare
Clean_sct.inte.rm.edc$sample <- factor(Clean_sct.inte.rm.edc$sample, levels = c(paste0('LP',seq(1,2)), paste0('HP',seq(1,2))))
# Clean_sct.inte.rm.edc$condition <- factor(Clean_sct.inte.rm.edc$condition, levels = c('LP','HP'))
cell.number <- FetchData(Clean_sct.inte.rm.edc, 
                         vars = c("sample", "cell.type.percise.new")) %>%
  dplyr::count(sample, cell.type.percise.new) %>% 
  tidyr::spread(sample, n) 
rio::export(cell.number, file = '25.8.add/cell number.xlsx')

sample_ids <- colnames(cell.number)[-1]
cell_types <- cell.number$cell.type.percise.new
n_cells_per_sample <- colSums(cell.number[,-1])
n_var_cats <- 2 

sample_cats <- tibble(
  Sample_ID = sample_ids,
  Treatment = c(rep('LP',2),rep('HP',2)),
  Rep = c(rep(c('one','two'),2)),
  cell.num = n_cells_per_sample
)

sample_num1_values<- rep(1,4)
obs_tbl <- data.frame(
  Sample_ID = rep(sample_ids, c(n_cells_per_sample)),
  Treatment = rep(sample_cats$Treatment, c(n_cells_per_sample)),
  Rep = rep(sample_cats$Rep, c(n_cells_per_sample)),
  Var_Num1 = rep(sample_num1_values, c(n_cells_per_sample))
)

obs_tbl$Cell_type <- c(rep(cell.number$cell.type.percise.new,c(cell.number$LP1)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$LP2)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$HP1)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$HP2)))


## RUN LTSR 
source('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/LTSR.raw.code.R')
results <- CellTypeCompositionAnalysis(obs_tbl, "Sample_ID", "Cell_type", c("Treatment",'Rep'), "Var_Num1")
ranef_tbl <- results$ranef
sdse_tbl <- results$sdse

vars1 <- list(Treatment = c('LP','HP'))

plot_ranef.new <- function(ranef_tbl, vars, celltypes = NULL, celltype_order = "hclust", references = NULL,
                           maxFC = 3, LTSR2p = F, highlightLtsr = 0.0, filterLtsr = 0.0, swap_axes = F) {
  ranef_tbl <- .getCondValLtsr(ranef_tbl, vars, celltypes = celltypes, references = references)
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
    condval = condval %>% pmin(log(maxFC)) %>% pmax(log(1 / maxFC)),
    ltsr = ltsr %>% pmin(0.9999) %>% pmax(0.5)
  )
  
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


ranef_plot <- plot_ranef.new(ranef_tbl, vars = vars1, celltypes = cell_types, celltype_order = rev(cell_types),
                             maxFC = 1.50, LTSR2p = FALSE) + xlab('Condition')
sdse_plot <- plot_sdse(sdse_tbl, "Sample_ID", ci = 0.95, xlim = c(0, 1))

ranef_plot + sdse_plot

plot_cell_fraction_sankey(Clean_sct.inte.rm.edc,condition = 'condition',cell.type = 'cell.type.percise.new',cols = cluster_colors)


# 04.Augur:prioritize the cell types most responsive to biological perturbations DISCARD-----
library(Augur)
augur <- calculate_auc(Clean_sct.inte.rm.edc, cell_type_col = 'cell.type.percise.new', label_col = 'condition', n_threads = 40)
head(augur$AUC,5)
p2 <- plot_lollipop(augur)

plot_umap(augur,
          Clean_sct.inte.rm.edc,
          mode = "default",
          reduction = "umap",
          palette = "YlGnBu", #  "viridis", "plasma", "magma", "inferno"
          # augur_mode = "default",
          cell_type_col = "cell.type.percise.new")
plot_loli <- function (augur) {
  aucs = augur$AUC
  size_sm = 10
  size_lg = 10
  range = range(aucs$auc)
  expand = abs(diff(range)) * 0.1
  p = aucs %>% ggplot(aes(x = reorder(cell_type, auc), y = auc)) + 
    geom_hline(aes(yintercept = 0.5), linetype = "dotted", 
               size = 0.8) +
    geom_point(size = 2) + 
    geom_text(aes(label = format(auc, digits = 3), 
                  y = ifelse(auc < 0.5, 0.5, auc)), 
              size = 4, nudge_y = expand, hjust = 0.5) + 
    geom_segment(aes(xend = cell_type, yend = 0.5)) + 
    scale_y_continuous("AUC", limits = c(min(range[1] - expand, 0.5), range[2] + expand * 1.5)) + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = size_sm), 
          axis.text.y = element_text(size = size_sm + 2),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_blank(), panel.grid = element_blank(), 
          strip.text = element_text(size = size_lg), strip.background = element_blank(),
          axis.line.y = element_blank(), axis.line.x = element_blank(), 
          legend.position = "top", legend.text = element_text(size = size_sm), 
          legend.title = element_text(size = size_sm), 
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          legend.background = element_blank(), 
          plot.title = element_text(size = size_lg, hjust = 0.5))
  p
}

p1 <- plot_loli(augur)+
  geom_segment(aes(xend=cell_type,yend=0.5),size=1)+
  geom_point(size=3,aes(color=cell_type))+
  scale_color_manual(values = c('Adipocytes' = use.cols[1],
                                'FAPs' = use.cols[2],
                                'MCs' = use.cols[3], 
                                'ECs' = use.cols[4],
                                'LECs' = use.cols[5],
                                'Pericytes' = use.cols[6],
                                'Schwann cells' = use.cols[7],
                                'Epididymal cells' = use.cols[8],
                                'Macrophages' = use.cols[9],
                                'DCs' = use.cols[10],
                                'Monocytes' = use.cols[11],
                                'B cells' = use.cols[12],
                                'T cells' = use.cols[13])) + 
  theme(legend.position = "") + ggtitle('Augur')
p1
### 08.2.scDist 
library(scDist)
sim <- list(Y = Clean_sct.inte.rm.edc@assays$SCT@data ,
            meta.data = Clean_sct.inte.rm.edc@meta.data %>% as.data.frame())
scdist <- scDist(normalized_counts = sim$Y,
                   meta.data = sim$meta.data,
                   d = 20,
                   fixed.effects = "condition",
                   # random.effects = 'orig.ident',
                   clusters="cell.type.percise.new"
)
DistPlot(scdist) + theme_bw() + ggtitle('scDist')


# 05.DEGs in cell types ----
library(ClusterGVis)
Clean_sct.inte.rm.edc <- PrepSCTFindMarkers(Clean_sct.inte.rm.edc)
rm.edc.ct.deg <- FindAllMarkers(Clean_sct.inte.rm.edc, only.pos = T,min.pct = 0.1,logfc.threshold = 0.25) %>% 
  dplyr::filter(p_val_adj < 0.05)
rm.edc.ct.deg.top10 <- rm.edc.ct.deg %>% dplyr::group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
ht.data <- prepareDataFromscRNA(Clean_sct.inte.rm.edc,
                                diffData = rm.edc.ct.deg,
                                showAverage = T,
                                assays = 'SCT',slot = 'data',
                                group.by = 'cell.type.percise.new',keep.uniqGene = F,
                                scale.data = T)

enrich.go <- enrichCluster(object = ht.data,
                           OrgDb = org.Mm.eg.db,
                           type = "BP",
                           organism = "mmu",
                           pvalueCutoff = 0.01,
                           topn = 5,
                           seed = 1234)
enrich.go$ratio <- -log(enrich.go$pvalue)
rio::export(rm.edc.ct.deg, file = 'table1.scRNAseq cell types DEGs.xlsx')

pdf('cell.type.degs.heatmap.pdf',height = 12,width = 14,onefile = F)

visCluster(object = ht.data,
           ht.col.list = list(col_range = c(-4, 0, 4)),
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           # markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich.go,
           line.side = "left",
           cluster.order = c(1:13),
           go.col = rep(c(use.cols),each = 5),
           sample.col = c(use.cols),
           ctAnno.col = c(use.cols),
           go.size = 8,
           add.bar = T)
dev.off()

rm.edc.ct.deg_pos_neg <- FindAllMarkers(Clean_sct.inte.rm.edc, only.pos = F,min.pct = 0.1,logfc.threshold = 0.25) %>% 
  dplyr::filter(p_val_adj < 0.05)
rio::export(rm.edc.ct.deg_pos_neg, file = '25.6.add/rm.edc.ct.deg_pos_neg.xlsx')

# 06.dotplot of key genes -----
library(plot1cell)

complex_dotplot_single <- function (seu_obj, feature, celltypes = NULL, groups, splitby = NULL, 
                                    color.palette = NULL, font.size = 12, strip.color = NULL, 
                                    do.scale = T, scale.by = "radius") 
{
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1", 
                                        "indianred1", "darkred"))(255)
  }
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (is.null(celltypes)) {
    celltypes <- levels(seu_obj)
  }
  if (length(groups) == 1) {
    groups_level <- levels(seu_obj@meta.data[, groups])
    if (is.null(groups_level)) {
      seu_obj@meta.data[, groups] <- factor(seu_obj@meta.data[, 
                                                              groups], levels = names(table(seu_obj@meta.data[, 
                                                                                                              groups])))
      groups_level <- levels(seu_obj@meta.data[, groups])
    }
    if (!is.null(splitby)) {
      if (is.null(levels(seu_obj@meta.data[, splitby]))) {
        seu_obj@meta.data[, splitby] <- factor(seu_obj@meta.data[, 
                                                                 splitby], levels = names(table(seu_obj@meta.data[, 
                                                                                                                  splitby])))
      }
      splitby_level <- levels(seu_obj@meta.data[, splitby])
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = c(groups, 
                                                                             splitby))
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], count_df[, splitby], 
                                  sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$splitby <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                            split = "___"), FUN = function(x) {
                                                              x[[3]]
                                                            }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$splitby <- factor(data_plot$splitby, levels = splitby_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    else {
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = groups)
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
      geom_tile(fill = "white", color = "white") + geom_point(aes(colour = avg.exp, 
                                                                  size = pct.exp)) + scale_color_gradientn(colours = color.palette) + 
      theme(panel.background = element_rect(fill = "white", 
                                            colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                          hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                            2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
            legend.text = element_text(size = (font.size - 
                                                 2)), legend.title = element_text(size = (font.size)), 
            strip.text = element_text(size = font.size), 
            legend.position = "right") + ylab("") + xlab("") + 
      ggtitle(feature)
    if (do.scale) {
      p = p + scale_size(range = c(0, 10))
    }
    else {
      if (max(data_plot$pct.exp) >= 20) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                        20))
      }
    }
    if (!is.null(splitby)) {
      p <- p + facet_wrap(~splitby, scales = "free_x")
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid.draw(g))
    }
    else {
      p
    }
  }
  else {
    gene_count <- extract_gene_count(seu_obj = seu_obj, 
                                     features = feature, cell.types = celltypes, meta.groups = c(groups, 
                                                                                                 splitby))
    allgroups <- c(groups, splitby)
    for (i in 1:length(allgroups)) {
      if (is.null(levels(seu_obj@meta.data[, allgroups[i]]))) {
        seu_obj@meta.data[, allgroups[i]] <- factor(seu_obj@meta.data[, 
                                                                      allgroups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                            allgroups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, allgroups[i]])
      gene_count[, allgroups[i]] <- factor(gene_count[, 
                                                      allgroups[i]], levels = group_level)
    }
    gene_count$celltype <- factor(gene_count$celltype, levels = celltypes)
    all_levels <- list()
    for (i in 1:length(groups)) {
      if (is.null(levels(seu_obj@meta.data[, groups[i]]))) {
        seu_obj@meta.data[, groups[i]] <- factor(seu_obj@meta.data[, 
                                                                   groups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                      groups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, groups[i]])
      all_levels[[i]] <- group_level
    }
    all_levels <- as.character(unlist(all_levels))
    data_plot <- list()
    for (i in 1:length(groups)) {
      count_df <- gene_count
      count_df$new_group <- paste(gene_count[, groups[i]], 
                                  gene_count[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      df1 <- merge(exp_df, pct_df, by = "new_group")
      df1$groupID <- groups[i]
      data_plot[[i]] <- df1
    }
    data_plot <- do.call("rbind", data_plot)
    data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                         split = "___"), FUN = function(x) {
                                                           x[[1]]
                                                         }))
    data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[2]]
                                                           }))
    data_plot$groups <- factor(data_plot$groups, levels = all_levels)
    data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    data_plot$groupID <- factor(data_plot$groupID, levels = groups)
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    if (is.null(splitby)) {
      p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_wrap(~groupID, scales = "free_x")
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid::grid.draw(g))
    }
    else {
      df2 <- reshape2::melt(gene_count[, c(groups, splitby)], 
                            measure.vars = groups)
      df2 <- df2[!duplicated(df2$value), ]
      colnames(df2)[colnames(df2) == "value"] <- "groups"
      data_plot2 <- list()
      for (i in 1:length(groups)) {
        df3 <- data_plot[data_plot$groupID == groups[i], 
        ]
        df4 <- df2[df2$variable == groups[i], c("groups", 
                                                splitby[i])]
        colnames(df4)[2] <- "split"
        df5 <- merge(df3, df4, by = "groups")
        data_plot2[[i]] <- df5
      }
      data_plot2 <- do.call("rbind", data_plot2)
      fill_x1 <- grDevices::rainbow(length(groups), alpha = 0.5)
      fill_x2 <- list()
      for (i in 1:length(splitby)) {
        n_col <- unique(gene_count[, splitby[i]])
        fill_x2[[i]] <- (scales::hue_pal(l = 90))(length(n_col))
      }
      fill_x2 <- as.character(unlist(fill_x2))
      fill_x <- c(fill_x1, fill_x2)
      p <- ggplot(data_plot2, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_nested(~groupID + split, 
                                        scales = "free_x", strip = strip_nested(background_x = elem_list_rect(fill = fill_x)))
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      p
    }
  }
}
com.dot.new <- com.dot.new <- function (seu_obj, features, celltypes = NULL, groups, color.palette = NULL, 
                                        strip.color = NULL) 
{
  pb <- progress_bar$new(format = "  Ploting [:bar] :percent eta: :eta", 
                         clear = FALSE, total = length(features), width = 100)
  plot_list <- list()
  for (i in 1:length(features)) {
    pp <- invisible(complex_dotplot_single(seu_obj = seu_obj, 
                                           feature = features[i], groups = groups, celltypes = celltypes))
    pp <- pp$data
    pp$gene <- features[i]
    plot_list[[i]] <- pp
    pb$tick()
    Sys.sleep(1/length(features))
  }
  all_data <- do.call("rbind", plot_list)
  all_data$gene <- factor(all_data$gene, levels = rev(features))
  all_data$celltype <- factor(all_data$celltype, levels = levels(seu_obj))
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1",
                                        "indianred1", "darkred"))(255)
  }
  p <- invisible(ggplot(all_data, aes(x = groups, y = gene)) + 
                   geom_tile(fill = "white", color = "white") + 
                   geom_point(aes(colour = avg.exp, size = pct.exp), alpha = 0.9) + 
                   scale_color_gradientn(colours = color.palette) + 
                   scale_size(range = c(0, 5)) +
                   theme(
                     panel.background = element_rect(fill = "white", colour = "black"), 
                     axis.text.x = element_text(angle = 45,hjust = 1),
                     axis.text.y = element_text(face = 'italic'),
                     plot.title = element_text(size = 10, hjust = 0.5,face = "bold"), 
                     axis.text = element_text(size = 12), 
                     axis.title = element_text(size = 8), 
                     legend.text = element_text(size = 8), 
                     legend.title = element_text(size = 12),
                     legend.position = "right", 
                     strip.text = element_text(size = 8, colour = "black",face = "bold")) + 
                   ylab("") + xlab("") + ggtitle("") + 
                   facet_wrap(~celltype, ncol = length(levels(seu_obj))))
  g <- change_strip_background(p, type = "top", strip.color = strip.color)
  print(grid.draw(g))
}
extract_gene_count <- function (seu_obj, features, cell.types = NULL, data.type = "data", 
                                meta.groups = NULL) 
{
  if (is.null(cell.types)) {
    cell.types = levels(seu_obj)
  }
  seu_obj@meta.data$celltype <- as.character(seu_obj@active.ident)
  if (is.null(meta.groups)) {
    meta.groups = colnames(seu_obj@meta.data)
  }
  if (!is.null(cell.types)) {
    new_seu <- subset(seu_obj, idents = cell.types)
  }
  feature_count <- Seurat::FetchData(new_seu, slot = data.type, 
                                     vars = c(features, meta.groups, "celltype"))
  umap_data <- data.frame(new_seu[["umap"]]@cell.embeddings)
  feature_count$UMAP1 <- umap_data$UMAP_1
  feature_count$UMAP2 <- umap_data$UMAP_2
  feature_count
}
change_strip_background <- function (ggplt_obj, type = "top", strip.color = NULL) 
{
  g <- ggplot_gtable(ggplot_build(ggplt_obj))
  if (type == "top") {
    strip_both <- which(grepl("strip-t", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else if (type == "right") {
    strip_both <- which(grepl("strip-r", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else {
    strip_t <- which(grepl("strip-t", g$layout$name))
    strip_r <- which(grepl("strip-r", g$layout$name))
    strip_both <- c(strip_t, strip_r)
    fills <- strip.color
    if (is.null(fills)) {
      fills <- c((scales::hue_pal(l = 90))(length(strip_t)), 
                 (scales::hue_pal(l = 90))(length(strip_r)))
    }
  }
  k <- 1
  for (i in strip_both) {
    j <- which(grepl("rect", g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  g
}

Idents(Clean_sct.inte.rm.edc) <- Clean_sct.inte.rm.edc$cell.type.percise.new
com.dot.new(Clean_sct.inte.rm.edc, feature = c('Ghr','Tshr','Adipor2',
                                               'Pdgfra','Fbn1','Egfr',
                                               'Upk3b','Upk1b','Msln',
                                               'Pecam1','Ptprb','Flt1',
                                               'Prox1','Mmrn1','Reln',
                                               'Abcc9','Notch3','Kcnq5',
                                               'Chl1','Grid2','Lrp1b',
                                               'Ank3','Abcb5','Pax2',
                                               'C1qa','Adgre1','Mrc1',
                                               'H2-Aa','H2-Eb1','Cd74',
                                               'Treml4','Pou2f2','Adgre4',
                                               'Cd79a','Cd79b','Ms4a1',
                                               'Cd3d','Cd3g','Skap1')
            ,groups = "condition",strip.color = c(use.cols)[1:13])

nad.enzyme <- c('Nampt','Nmnat1','Nmnat2','Nmnat3','Ido1','Cd38','Bst1','Sirt1','Sirt3','Sirt5','Parp2')
com.dot.new(Clean_sct.inte.rm.edc, feature = nad.enzyme,
            groups = "condition",strip.color = c(use.cols)[1:13])

# 07.trans noise ----

## cell levels
library(MASS)
library(ggpubr)
celltypes <- unique(Clean_sct.inte.rm.edc@meta.data$cell.type.percise.new)
celltypes <- celltypes[which(!is.na(celltypes))]
Idents(Clean_sct.inte.rm.edc) <- Clean_sct.inte.rm.edc$cell.type.percise.new
set.seed(12345)

getEuclideanDistance <- function(celltype, obj, assay, slot, ident1, ident2, group.by, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- subset(obj, cells = WhichCells(obj, idents = celltype))
  
  counts <- GetAssayData(object = tmp, slot = slot, assay = assay)
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) > 0
  expr <- counts[keep_genes, ]
  
  ifelse(min(table(tmp@meta.data[[group.by]])) > 300,
         expr <- expr[,c(rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,]),300)],
                         rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,]),300)])
         ],
         expr <- expr)
  tmp <- subset(tmp,cells = colnames(expr))
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data[[group.by]])[c(ident1,ident2)])
  
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  print(nsample)
  ident2_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident2)], nsample)
  ident1_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident1)], nsample)
  ds_expr_r <- ds_expr[, c(ident1_r, ident2_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, ident1, ident2){
    tmp <- data.matrix(sqrt(matr[genes, ident1]))
    mean <- rowMeans(sqrt(matr[genes, ident1]))
    d_ident1 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident1) <- ident1
    
    tmp <- data.matrix(sqrt(matr[genes, ident2]))
    mean <- rowMeans(sqrt(matr[genes, ident2]))
    d_ident2 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident2) <- ident2
    
    list(ident1 = d_ident1, ident2 = d_ident2)
  }
  ds <- calcEuclDist(matr = ds_expr_r, ident2 = ident2_r, ident1 = ident1_r)
  ds
}

res <- lapply(celltypes, function(x) getEuclideanDistance(x, 
                                                          obj = Clean_sct.inte.rm.edc,
                                                          assay = 'SCT',
                                                          slot = 'counts',
                                                          group.by = 'condition',
                                                          ident1 = 'LP',
                                                          ident2 = 'HP',
                                                          lowcv = F))
names(res) <- celltypes
res.df <- data.frame(TN.value = unlist(do.call(c, res))) %>% 
  rownames_to_column(var = 'info') %>% 
  separate(col = info,
           into = c('Celltype','Condition','Sample_cells'),
           sep = '\\.',
           remove = T,
           extra = "merge")
res.df$Condition <- factor(ifelse(res.df$Condition == 'ident1','LP','HP'), levels = c('LP','HP'))
res.df$Celltype <- factor(res.df$Celltype, levels = c('Adipocytes',
                                                      'FAPs',
                                                      'MCs',
                                                      'ECs',
                                                      'LECs',
                                                      'Pericytes',
                                                      'Schwann cells',
                                                      'Epididymal cells',
                                                      'Macrophages',
                                                      'DCs',
                                                      'Monocytes',
                                                      'B cells',
                                                      'T cells'))

my_comparisons <- list(c("LP","HP"))
ggplot(res.df, aes(x=Condition, y=TN.value, fill=Condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw() +
  labs(x = '', y = 'Transcriptional \nheterogeneity') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  facet_wrap(~Celltype, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

ggdensity(res.df, 
          x = "TN.value",
          add = "median", rug = F,
          color = "Condition", fill = "Condition",
          palette = c('#0B996F',"#D6570D")) +
  ylab("Density") + 
  xlab('Transcriptional \nheterogeneity') + 
  # ggtitle('infm') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) 

TN.cts <- data.frame(div = 0, pvalue = 0)
for (ct in celltypes) {
  tmp <- res.df[res.df$Celltype == ct,]
  div <- sum(tmp[tmp$Condition == 'HP',]$TN.value) / sum(tmp[tmp$Condition == 'LP',]$TN.value)
  pvalue <- wilcox.test(tmp[tmp$Condition == 'HP',]$TN.value, tmp[tmp$Condition == 'LP',]$TN.value)
  pvalue <- pvalue[['p.value']]
  tmp.df <- data.frame(div = div, pvalue = pvalue)
  TN.cts <- rbind(TN.cts,tmp.df)
}
TN.cts <- TN.cts[-1,]
rownames(TN.cts) <- celltypes

res.df.sample <- res.df %>% 
  dplyr::group_by(Celltype,Samples) %>% 
  summarise(tn.mean = mean(TN.value, na.rm = T)) %>% ungroup()

res.df.condition <- res.df %>% 
  dplyr::group_by(Celltype,Condition) %>% 
  summarise(tn.mean = mean(TN.value, na.rm = T)) %>% ungroup()

res.df.ratio <- res.df.condition[seq(2,26,2),]/res.df.condition[seq(1,25,2),]
res.df.ratio$Celltype <- levels(res.df.condition$Celltype)
res.df.ratio$log2 <- log2(res.df.ratio$tn.mean)
res.df.ratio$Celltype <- factor(res.df.ratio$Celltype, levels = rev(c('Adipocytes',
                                                                      'FAPs',
                                                                      'MCs',
                                                                      'ECs',
                                                                      'LECs',
                                                                      'Pericytes',
                                                                      'Schwann cells',
                                                                      'Epididymal cells',
                                                                      'Macrophages',
                                                                      'DCs',
                                                                      'Monocytes',
                                                                      'B cells',
                                                                      'T cells')))

ggplot(res.df.ratio, aes(x=log2, y=Celltype ,color = Celltype)) + 
  geom_point(size = 4) + 
  geom_vline(xintercept = 0) + 
  theme_bw() +
  labs(x = 'Log2(HP/LP)', y = 'Cell types') +
  theme(legend.position = "none") +
  scale_color_manual(values = rev(use.cols)) 

## sample levels discard
Clean_sct.inte.rm.edc$ct_con <- factor(paste0(Clean_sct.inte.rm.edc$condition, "_", Clean_sct.inte.rm.edc$cell.type.percise.new),
                                       levels = c(paste0("LP", "_", levels(Clean_sct.inte.rm.edc$cell.type.percise.new)),
                                                  paste0("HP", "_", levels(Clean_sct.inte.rm.edc$cell.type.percise.new))))
Clean_sct.inte.rm.edc$samp_con <- factor(paste0(Clean_sct.inte.rm.edc$sample, "_", Clean_sct.inte.rm.edc$cell.type.percise.new),
                                         levels = paste0(c("LP1","LP2","LP3","HP1",'HP2','HP3'), "_", 
                                                         rep(levels(Clean_sct.inte.rm.edc$cell.type.percise.new),each = 6)))
cell.type.counts.aggr <- AggregateExpression(Clean_sct.inte.rm.edc,assays = 'SCT',group.by = 'samp_con') %>% as.data.frame()

getEuclideanDistance.samp <- function(celltype, obj, ident1, ident2, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- obj[, grep(pattern = celltype,x = colnames(obj),fixed = T)]
  expr <- tmp[rowSums(tmp) > 0, ]
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  ident2_r <- grep(x = colnames(ds_expr), pattern = ident2,fixed = T)
  ident1_r <- grep(x = colnames(ds_expr), pattern = ident1,fixed = T)
  ds_expr_r <- ds_expr[, c(ident1_r, ident2_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, ident1, ident2){
    tmp <- data.matrix(sqrt(matr[genes, ident1]))
    mean <- rowMeans(sqrt(matr[genes, ident1]))
    d_ident1 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident1) <- ident1
    
    tmp <- data.matrix(sqrt(matr[genes, ident2]))
    mean <- rowMeans(sqrt(matr[genes, ident2]))
    d_ident2 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident2) <- ident2
    
    list(ident1 = d_ident1, ident2 = d_ident2)
  }
  ds <- calcEuclDist(matr = ds_expr_r, ident2 = ident2_r, ident1 = ident1_r)
  ds
}

res.samp <- lapply(paste0('.',sub("^(?:[^.]*\\.){2}","",colnames(cell.type.counts.aggr)) %>% unique()), 
                   function(x) getEuclideanDistance.samp(x, 
                                                         obj = cell.type.counts.aggr,
                                                         ident1 = 'LP',
                                                         ident2 = 'HP',
                                                         lowcv = F))
names(res.samp) <- sub("^(?:[^.]*\\.){2}","",colnames(cell.type.counts.aggr)) %>% unique()

res.samp.df <- data.frame(TN.value = unlist(do.call(c, res.samp))) %>% 
  rownames_to_column(var = 'info') %>% 
  separate(col = info,
           into = c('Celltype','Condition'),
           sep = '.i',
           remove = T,
           extra = "merge")
res.samp.df$Condition <- factor(rep(c(rep(c('LP','HP'),each = 3)), 15), levels = c('LP','HP'))
ggplot(res.samp.df, aes(x=Condition, y=TN.value, fill=Condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw() +
  labs(x = '', y = 'Transcriptional \nheterogeneity') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  facet_wrap(~Celltype, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

# 08.umap of key genes ----
p1 <- FeaturePlot(Clean_sct.inte.rm.edc, features = c('Cd38'),split.by = 'condition',order = T ,combine = F, cols = c('gray90','red3'))
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

## macrophages
macrop.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'Macrophages'))
VlnPlot(macrop.seu,features =  c('Cd38'),
        split.by = 'condition',
        split.plot = F,
        # group.by = 'condition',
        assay = 'SCT',layer = 'data',
        pt.size = 0,cols = c('#0B996F',"#D6570D"))
DimPlot(macrop.seu, group.by = 'condition',split.by = 'condition')
macrop.seu.ct <- macrop.seu@meta.data[,c('cells','condition')]
write.csv(macrop.seu.ct, file = '25.6.add/compass/macrop/macrop.seu.ct.csv',quote = F,row.names = T,col.names = T)
## adipocytes
adipo.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'Adipocytes'))
VlnPlot(adipo.seu,features =  c('Cd38'),
        split.by = 'condition',
        split.plot = F,
        # group.by = 'condition',
        assay = 'SCT',layer = 'data',
        pt.size = 0,cols = c('#0B996F',"#D6570D"))

DimPlot(adipo.seu, group.by = 'condition',split.by = 'condition')
FeaturePlot(adipo.seu,features = 'Cd38', split.by = 'condition')
adipo.seu.ct <- adipo.seu@meta.data[,c('cells','condition')]
write.csv(adipo.seu.ct, file = '25.6.add/compass/adipo/adipo.seu.ct.csv',quote = F,row.names = T,col.names = T)

## macrophages pseudobulk
library(DESeq2)

macrop.seu.bulk <- AggregateExpression(macrop.seu, assays = 'SCT',group.by = 'sample') %>% as.data.frame()
colData = data.frame(row.names = colnames(macrop.seu.bulk),group_list = c(rep('LP',2),rep('HP',2)))
dds <- DESeqDataSetFromMatrix(countData = macrop.seu.bulk,colData = colData,design = ~group_list)
dds <- DESeq(dds)
dds_ashr <- lfcShrink(dds, contrast = c('group_list','HP','LP'),type = 'normal') # adjust LFC
summary(dds_ashr)
macrop_deseq<- dds_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  na.omit() %>% 
  arrange(log2FoldChange)
macrop_deseq$change <- factor(ifelse(macrop_deseq$padj < 0.05 & macrop_deseq$log2FoldChange > 0.75,'HP up',
                                  ifelse(macrop_deseq$padj < 0.05 & macrop_deseq$log2FoldChange < -0.75, 'LP up','not')
))
table(macrop_deseq$change)
rio::export(macrop_deseq, file = '25.6.add/macrop_deseq.lfc0.75.xlsx')

# heatmap
macrop.seu.bulk.cpm <- cpm(macrop.seu.bulk)
plot.data <- na.omit(macrop.seu.bulk.cpm[macrop_deseq[macrop_deseq$change !='not',]$geneid,])
ann <- data.frame(row.names = colnames(macrop.seu.bulk),Treatment = factor(c(rep('LP',2),rep('HP',2)),levels = c('LP','HP')))
ha_left.col <- list(Treatment = c('#0B996F',"#D6570D"))
names(ha_left.col$Treatment) <- factor(unique(ann$Treatment))

cd38_pos <- match('Cd38', rownames(plot.data))hao

ha <- rowAnnotation(
  highlight = anno_mark(
    at        = cd38_pos,
    labels    = 'Cd38',
    labels_gp = gpar(col = "black", fontsize = 9),
    link_gp   = gpar(col = "darkgray", lty = 2)
  )
)
ComplexHeatmap::pheatmap(plot.data, 
                         name = "Z-score",
                         border_color = NA, 
                         clustering_method = 'ward.D',cluster_cols = F,
                         show_colnames = F,
                         show_rownames = F,
                         scale = 'row',
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann,
                         annotation_colors = ha_left.col,
                         right_annotation = ha,
                         main = 'Macrophages.DEGs'
)

## adipocytes pseudobulk
library(DESeq2)

adipo.seu.bulk <- AggregateExpression(adipo.seu, assays = 'SCT',group.by = 'sample') %>% as.data.frame()
colData = data.frame(row.names = colnames(adipo.seu.bulk),group_list = c(rep('LP',2),rep('HP',2)))
dds <- DESeqDataSetFromMatrix(countData = adipo.seu.bulk,colData = colData,design = ~group_list)
dds <- DESeq(dds)
dds_ashr <- lfcShrink(dds, contrast = c('group_list','HP','LP'),type = 'normal') # adjust LFC
summary(dds_ashr)
adipo_deseq<- dds_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  na.omit() %>% 
  arrange(log2FoldChange)
adipo_deseq$change <- factor(ifelse(adipo_deseq$padj < 0.05 & adipo_deseq$log2FoldChange > 0.75,'HP up',
                                     ifelse(adipo_deseq$padj < 0.05 & adipo_deseq$log2FoldChange < -0.75, 'LP up','not')
))
table(adipo_deseq$change)
rio::export(adipo_deseq, file = '25.8.add/adipo_deseq.lfc0.75.xlsx')

# heatmap
adipo.seu.bulk.cpm <- cpm(adipo.seu.bulk)
plot.data <- na.omit(adipo.seu.bulk.cpm[adipo_deseq[adipo_deseq$change !='not',]$geneid,])
ann <- data.frame(row.names = colnames(adipo.seu.bulk),Treatment = factor(c(rep('LP',2),rep('HP',2)),levels = c('LP','HP')))
ha_left.col <- list(Treatment = c('#0B996F',"#D6570D"))
names(ha_left.col$Treatment) <- factor(unique(ann$Treatment))

ComplexHeatmap::pheatmap(plot.data, 
                         name = "Z-score",
                         border_color = NA, 
                         clustering_method = 'ward.D',cluster_cols = F,
                         show_colnames = F,
                         show_rownames = F,
                         scale = 'row',
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann,
                         annotation_colors = ha_left.col,
                         # right_annotation = rowAnnotation(
                         #   mark = anno_mark(
                         #     at        = match(c('Adipoq','Cd38'), rownames(plot.data)),
                         #     labels    = c('Adipoq','Cd38'),
                         #     labels_gp = gpar(col = "black", fontsize = 9),
                         #     link_gp   = gpar(col = "black")
                         #   )
                         # ),
                         main = 'Adipocytes DEGs'
)

pick.term <- import('adipo.pick.term.xlsx',header = F)
# colnames(pick.term) <- c('term','P')
# library(forcats)
# pick.term$group <- c(rep('HP',10),rep('LP',10))
# pick.term$term <- factor(pick.term$term, levels = pick.term$term)
# pick.term$logp <- -pick.term$P

go.vis <- pick.term
colnames(go.vis) <- c('Terms','LogP')
go.vis$lp <- -go.vis$LogP
go.vis <- go.vis[order(go.vis$lp),]
go.vis$Terms <- factor(go.vis$Terms,levels = go.vis$Terms)
go.vis$group <- c(rep('HP',10),rep('LP',10))

ggplot() +
  geom_bar(data = go.vis[go.vis$group == 'LP',],
           aes(x = lp, y = Terms,fill = 'red3'),
           width = 0.5,
           alpha = 0.5,
           stat = 'identity') +
  theme_classic() + 
  scale_fill_manual(values = '#D6570D') +
  
  scale_x_continuous(expand = c(0,0)) + 
  theme(axis.text.y = element_blank()) +
  # scale_fill_manual(values = ct.cols) +
  geom_text(data = go.vis[go.vis$group == 'LP',],
            aes(x = 0.1, 
                y = Terms, 
                label = Terms),
            size = 4.5,
            hjust = 0) +
  labs(x = bquote(~-Log[10]~ '(Pvalue)'), 
       y = 'HP pathways') 

# 09.adipocytes -----

### integrated
adipo.seu[['RNA']] <- split(adipo.seu[['RNA']], f = adipo.seu$sample)
adipo.seu.inte <- adipo.seu %>%
  # SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = .5) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)
DimPlot(adipo.seu.inte, group.by = 'condition',split.by = 'condition',cols = c(npg.cols[c(2,1)])) + 
  ggtitle('adipocytes integrated') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),
        title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 

adipo.seu.inte <- adipo.seu.inte %>% 
  FindNeighbors( reduction = "integrated.dr", dims = 1:30) %>% 
  FindClusters(resolution = .3) 

DimPlot(adipo.seu.inte, split.by = 'condition',group.by = 'SCT_snn_res.0.3',label = T) + 
  scale_color_simpsons() +
  ggtitle('adipohHPs integrated') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 


### merged 
adipo.seu <- RunUMAP(adipo.seu,reduction = "integrated.dr", 
                      dims = 1:30)
DimPlot(adipo.seu, group.by = 'condition',split.by = 'condition',cols = c(npg.cols[c(2,1)]))
adipo.seu.split <- SplitObject(adipo.seu,split.by = 'sample')
adipo.seu.merge <- merge(x = adipo.seu.split$LP1, y = c(adipo.seu.split$LP2, 
                                                             adipo.seu.split$HP1,
                                                             adipo.seu.split$HP2))
adipo.seu.merge <- adipo.seu.merge %>% 
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>%
  RunPCA() %>% 
  RunUMAP(reduction = "pca", 
          dims = 1:30)
adipo.seu.merge$condition <- factor(adipo.seu.merge$condition, levels = c('LP','HP'))
DimPlot(adipo.seu.merge, group.by = 'condition',split.by = 'condition') + 
  scale_color_simpsons() +
  ggtitle('adipocytes only merged') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 


adipo.seu.merge <- adipo.seu.merge %>% FindNeighbors( reduction = "pca", 
                                                        dims = 1:30) %>% 
  FindClusters(resolution = .3) 

DimPlot(adipo.seu.merge, split.by = 'condition',group.by = 'SCT_snn_res.0.3',label = T) + 
  scale_color_simpsons() +
  ggtitle('adipohHPs merged') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 

Idents(adipo.seu.merge) <- 'SCT_snn_res.0.3'
macro.ct.deg <- FindAllMarkers(adipo.seu.merge,only.pos = T,logfc.threshold = 0.25,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
rio::export(macro.ct.deg, file = 'R.hd/adipohHPs.merged.cluster.deg.xlsx')

adipo.seu.merge <- adipo.seu.merge %>% 
  FindClusters(resolution = .8) 

p1 <- DimPlot(adipo.seu.merge, group.by = 'SCT_snn_res.0.8',split.by = 'condition',label = T) +
  scale_color_simpsons() +
  ggtitle('adipohHPs only merged') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p2 <- FeaturePlot(adipo.seu.merge, features = 'Cd38',split.by = 'condition', order = T,
                  # min.cutoff = 'q1',
                  cols =c("grey90", "red3"),) &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p1 / p2

Idents(adipo.seu.merge) <- 'SCT_snn_res.0.8'
adipo.seu.merge <- RenameIdents(adipo.seu.merge,
                                 '0' = 'Cd38 +',
                                 '1' = 'Cd38 +',
                                 '2' = 'Cd38 +',
                                 '3' = 'Cd38 +',
                                 '4' = 'Cd38 +',
                                 '5' = 'Cd38 +'
)
adipo.seu.merge$adipo.type <- as.character(Idents(adipo.seu.merge))
adipo.seu.merge$adipo.type[adipo.seu.merge$adipo.type != 'Cd38 +' ] <- 'Cd38 -'
adipo.seu.merge$adipo.type <- factor(adipo.seu.merge$adipo.type,levels = c('Cd38 +','Cd38 -'))

p1 <- DimPlot(adipo.seu.merge, group.by = 'adipo.type',label = T) +
  scale_color_simpsons() +
  ggtitle('adipohHPs subtypes') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p2 <- FeaturePlot(adipo.seu.merge, features = 'Cd38', order = T,
                  # min.cutoff = '0.2',
                  cols =c("grey90", "red3"),) &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p1 + p2

adipo.subtype.deg <- FindMarkers(adipo.seu.merge, group.by = 'adipo.type',min.pct = 0.1,logfc.threshold = 0.25,
                                  ident.1 = 'Cd38 +', ident.2 = 'Cd38 -', recorrect_umi = F) %>% 
  dplyr::filter( p_val_adj < 0.05) %>% rownames_to_column(var = 'genes')
adipo.subtype.deg$type <- factor(ifelse(adipo.subtype.deg$avg_log2FC > 0,'Cd38 + up', 'Cd38 - up'))
rio::export(adipo.subtype.deg, file = 'R.hd/adipo.subtype.deg.xlsx')

# FeaturePlot(adipo.seu.merge, features = c('Cd163','Cd86','Cxcl10','Arg1','Mrc1'), order = T,
#             min.cutoff = 'q50',split.by = 'condition',raster = F,
#             # cols =c("#440154FF",'#31688EFF','#35B779FF', "#FDE725FF")
#             cols = c('grey90','darkred')
#             ) &
#   NoAxes() & 
#   theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
# 
# FeaturePlot(adipo.seu.merge, features = c('Il1b','Spp1','Apoe','Hes1','Mki67'), order = T,
#             min.cutoff = 'q50',split.by = 'condition',raster = F,
#             # cols =c("#440154FF",'#31688EFF','#35B779FF', "#FDE725FF")
#             cols = c('grey90','darkred')
# ) &
#   NoAxes() & 
#   theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
# 
# FeaturePlot(adipo.seu.merge, features = c('Il1b','Cxcl1','Cxcl2','Cxcl3','Il1rn'), order = T,
#             min.cutoff = 'q20',
#             split.by = 'condition',raster = F,
#             pt.size = .3,
#             # cols =c("#440154FF",'#31688EFF','#35B779FF', "#FDE725FF")
#             cols = c('grey90','darkred')
# ) &
#   NoAxes() & 
#   theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
# 
# 
# FeaturePlot(adipo.seu.merge, features = c('Cd68','Cd86','Cd163','Cd274','Mrc1'), order = T,
#             min.cutoff = 'q20',
#             split.by = 'condition',raster = F,
#             pt.size = .3,
#             # cols =c("#440154FF",'#31688EFF','#35B779FF', "#FDE725FF")
#             cols = c('grey90','darkred')
# ) &
#   NoAxes() & 
#   theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
# 
# # DEGs in adipohHPs merged
# comlist <- t(combn(unique(adipo.seu.merge$condition), 2))
# macro.deg.list <- list()
# for (ct in as.character(c(0,1,2,3,5,6))) {
#   print(ct)
#   Idents(adipo.seu.merge) <- 'SCT_snn_res.0.3'
#   sub.seu.1 <- subset(adipo.seu.merge, idents = ct)
#   Idents(sub.seu.1) <- 'condition'
#   for (gp in nrow(comlist)) {
#     HP1 <- comlist[gp,1]
#     HP2 <- comlist[gp,2]
#     DEGs <- FindMarkers(sub.seu.1,
#                         ident.1 = HP1,
#                         ident.2 = HP2,
#                         logfc.threshold = 0.2,
#                         group.by = 'condition',
#                         recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
#     DEGs$genes <- rownames(DEGs)
#     DEGs.name <- paste('macro.deg',ct,HP1,HP2,sep = '_')
#     macro.deg.list[[DEGs.name]] <- DEGs
#   }
# }
# class(macro.deg.list$macro.deg_0_LP_HP)
# class(macro.deg.list)
# rio::export(macro.deg.list, file = 'R.hd/adipohHPs.deg.LP_vs_HP.xlsx')

library(ensembldb)
library(org.Mm.eg.db)
hub<-AnnotationHub::AnnotationHub()
annolist <- AnnotationHub::query(hub, "Mus musculus")
ensdb110 <- hub[["AH113713"]]
plot.list <- list()
adipo.term <- list()

library(clusterProfiler)
for (ct in levels(adipo.subtype.deg$type)) {
  tmp <- adipo.subtype.deg[adipo.subtype.deg$type == ct,]
  gene2entrzid <- bitr(tmp$genes, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,' up'))
  erich.name <- paste0('macro_',ct,'_plot')
  adipo.term[[erich.name]] <- erich.go.BP
  # assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}
patchwork::wrap_plots(c(plot.list),nrow = 1)
rio::export(adipo.term, file = 'R.hd/adipo.term.xlsx')
write.csv(adipo.term$`macro_Cd38 + up_plot`, file = 'R.hd/adipo.cd38+_up.term.csv')
write.csv(adipo.term$`macro_Cd38 - up_plot`, file = 'R.hd/adipo.cd38-_up.term.csv')

# 09.SASP ----
sasp.genenset <- rio::import('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/SASP.geneset.xlsx')
sasp.geneuse <- sasp.genenset$`Gene ID` %>% tolower() %>% stringr::str_to_title() 
sasp.geneuse <- list(c(sasp.geneuse))
Clean_sct.inte.rm.edc <- AddModuleScore(Clean_sct.inte.rm.edc, features = sasp.geneuse,name = 'SASP')

p1 <- FeaturePlot(Clean_sct.inte.rm.edc, features = 'SASP1', split.by = 'condition', order = T, min.cutoff = 'q50', cols = c('grey90', "red3"), combine = F)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

sasp.score <- Clean_sct.inte.rm.edc@meta.data[,c('cell.type.percise.new','SASP1','condition')] %>% na.omit()
sasp.score$condition <- factor(sasp.score$condition, levels = c('LP','HP'))

ggplot(sasp.score, aes(x=condition, y=SASP1, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'SASP score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

p1 <- ggpubr::ggviolin(sasp.score, x="condition", y="SASP1", width = 0.6, 
                       add = 'boxplot',
                       add.params = list(fill = 'gray90',alpha = 0.6),
                       color = "black",#轮廓颜色 
                       fill="condition",#填充
                       # palette = "npg",
                       xlab = F, #不显示x轴的标签
                       # bxp.errorbar=T,#显示误差条
                       # bxp.errorbar.width=0.5, #误差条大小
                       size=.2, #箱型图边线的粗细
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.9) + 
  ylab('SASP score')  + 
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  theme_bw() + 
  # scale_color_npg() +
  facet_wrap(~ cell.type.percise.new,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

p2 <- ggpubr::ggviolin(sasp.score, x="condition", y="SASP1", width = 0.6, 
                       # add = 'boxplot',
                       # add.params = list(fill = 'gray90',alpha = 0.6),
                       color = "black",#轮廓颜色 
                       fill="condition",#填充
                       # palette = "npg",
                       xlab = F, #不显示x轴的标签
                       # bxp.errorbar=T,#显示误差条
                       # bxp.errorbar.width=0.5, #误差条大小
                       size=.2, #箱型图边线的粗细
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.9) + 
  ylab('SASP score')  + 
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  theme_bw() + 
  # scale_color_npg() +
  facet_wrap(~ cell.type.percise.new,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p1 + p2

ggplot(sasp.score, aes(x = SASP1, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  xlim(-0.05,0.10) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

ggdensity(sasp.score, 
          x = "SASP1",
          add = "median", rug = F,
          color = "condition", fill = "condition",
          palette = c('#0B996F',"#D6570D")) +
  ylab("Density") + 
  xlab('SASP score') + 
  ggtitle('SASP') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) 
wilcox.test(sasp.score[sasp.score$condition == 'LP',]$SASP1, sasp.score[sasp.score$condition == 'HP',]$SASP1)
median(sasp.score[sasp.score$condition == 'LP' & sasp.score$cell.type.percise.new == 'DSCs',]$SASP1)
median(sasp.score[sasp.score$condition == 'HP' & sasp.score$cell.type.percise.new == 'DSCs',]$SASP1)

ggpubr::ggdensity(sasp.score, x = 'SASP1',add = 'median',rug = T,color = 'condition', fill = 'condition')

## AUC
library(AUCell)
expr_mat <- GetAssayData(Clean_sct.inte.rm.edc, assay="SCT", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('SASP' = c(sasp.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings)))
Clean_sct.inte.rm.edc$SASP_AUC <- as.numeric(getAUC(cells_AUC)['SASP',])

sasp.auc <- Clean_sct.inte.rm.edc@meta.data[,c('SASP_AUC','condition')] %>% na.omit()
table(sasp.auc$condition)

ggplot(sasp.auc, aes(x = SASP_AUC, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  ggtitle("Macrophages SASP") + xlab('AUC') + ylab("") + 
  xlim(0.0,0.15) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 


# 10.inflammatory reponse----
infm.genenset <- rio::import('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/inflammatory_response.xlsx')
infm.geneuse <- infm.genenset$inflammatory.response
infm.geneuse <- list(c(infm.geneuse))
# infm.go <- GO_DATA$PATHID2EXTID$`GO:0006954`
Clean_sct.inte.rm.edc <- AddModuleScore(Clean_sct.inte.rm.edc, features = infm.geneuse,name = 'infm',assay = 'SCT')

p1 <- FeaturePlot(Clean_sct.inte.rm.edc, features = 'infm1', split.by = 'condition', order = T, min.cutoff = 'q50', cols = c('grey90', "red3"), combine = F)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

infm.score <- Clean_sct.inte.rm.edc@meta.data[,c('cell.type.percise.new','infm1','condition')] %>% na.omit()
infm.score$condition <- factor(infm.score$condition, levels = c('LP','HP'))

ggplot(infm.score, aes(x=condition, y=infm1, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'infm score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

p1 <- ggpubr::ggviolin(infm.score, x="condition", y="infm1", width = 0.6, 
                       add = 'boxplot',
                       add.params = list(fill = 'gray90',alpha = 0.6),
                       color = "black",#轮廓颜色 
                       fill="condition",#填充
                       # palette = "npg",
                       xlab = F, #不显示x轴的标签
                       # bxp.errorbar=T,#显示误差条
                       # bxp.errorbar.width=0.5, #误差条大小
                       size=.2, #箱型图边线的粗细
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.9) + 
  ylab('Inflammation reponse score')  + 
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  theme_bw() + 
  # scale_color_npg() +
  facet_wrap(~ cell.type.percise.new,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p2 <- ggpubr::ggviolin(infm.score, x="condition", y="infm1", width = 0.6, 
                       # add = 'boxplot',
                       # add.params = list(fill = 'gray90',alpha = 0.6),
                       color = "black",#轮廓颜色 
                       fill="condition",#填充
                       # palette = "npg",
                       xlab = F, #不显示x轴的标签
                       # bxp.errorbar=T,#显示误差条
                       # bxp.errorbar.width=0.5, #误差条大小
                       size=.2, #箱型图边线的粗细
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.9) + 
  ylab('Inflammation reponse score')  + 
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  theme_bw() + 
  # scale_color_npg() +
  facet_wrap(~ cell.type.percise.new,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p1 + p2

ggdensity(infm.score, 
          x = "infm1",
          add = "median", rug = F,
          color = "condition", fill = "condition",
          palette = c('#0B996F',"#D6570D")) +
  ylab("Density") + 
  xlab('infm score') + 
  ggtitle('infm') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) 
wilcox.test(infm.score[infm.score$condition == 'LP',]$infm1, infm.score[infm.score$condition == 'HP',]$infm1)
median(infm.score[infm.score$condition == 'LP' & infm.score$cell.type.percise.new == 'DSCs',]$infm1)
median(infm.score[infm.score$condition == 'HP' & infm.score$cell.type.percise.new == 'DSCs',]$infm1)

ggplot(infm.score, aes(x = infm1, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  ggtitle("INFM") + xlab('AUC') + ylab("") + 
  xlim(-0.02,0.10) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

## AUC
library(AUCell)
expr_mat <- GetAssayData(macrop.seu, assay="SCT", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('INFM' = c(infm.use[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings)))
macrop.seu$INFM_AUC <- as.numeric(getAUC(cells_AUC)['INFM',])

infm.auc <- macrop.seu@meta.data[,c('INFM_AUC','condition')] %>% na.omit()
table(infm.auc$condition)

ggplot(infm.auc, aes(x = INFM_AUC, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  ggtitle("Macrophages INFM") + xlab('AUC') + ylab("") + 
  xlim(0.0,0.20) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 


# 11.DEGs cross condition -----
comlist <- t(combn(unique(Clean_sct.inte.rm.edc$condition), 2))
deg_corss_condition_age_up <- list()
for (ct in unique(Clean_sct.inte.rm.edc$cell.type.percise.new)) {
  print(ct)
  Idents(Clean_sct.inte.rm.edc) <- 'cell.type.percise.new'
  sub.seu.1 <- subset(Clean_sct.inte.rm.edc, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    HP1 <- comlist[gp,1]
    HP2 <- comlist[gp,2]
    DEGs <- FindMarkers(sub.seu.1,
                        ident.1 = HP2,
                        ident.2 = HP1,
                        min.pct = 0.1,
                        logfc.threshold = 0.2,
                        group.by = 'condition',
                        recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    # DEGs.name <- paste('sc.deg',ct,HP1,HP2,sep = '_')
    # assign(DEGs.name,DEGs)
    deg_corss_condition_age_up[[ct]] <- rownames(DEGs)
  }
}

deg_corss_condition_age_up <- lapply(deg_corss_condition_age_up, function(vec) {
  data.frame(gene = vec, stringsAsFactors = FALSE)
})
rio::export(deg_corss_condition_age_up, file = '25.6.add/deg_corss_condition_age_up.xlsx')

library(purrr)
max_len <- max(map_int(deg_corss_condition_age_up, length))
deg_corss_condition_age_up.df <- map_dfr(deg_corss_condition_age_up, ~ { length(.x) <- max_len; .x })

deg_corss_condition_age_up.df <- list2DF(lapply(deg_corss_condition_age_up, `length<-`, max(lengths(deg_corss_condition_age_up))))
colnames(deg_corss_condition_age_up.df)
rio::export(deg_corss_condition_age_up.df, file = '25.6.add/deg_corss_condition_age_up.df.xlsx')

deg_corss_condition_age_down <- list()
for (ct in unique(Clean_sct.inte.rm.edc$cell.type.percise.new)) {
  print(ct)
  Idents(Clean_sct.inte.rm.edc) <- 'cell.type.percise.new'
  sub.seu.1 <- subset(Clean_sct.inte.rm.edc, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    HP1 <- comlist[gp,1]
    HP2 <- comlist[gp,2]
    DEGs <- FindMarkers(sub.seu.1,
                        ident.1 = HP1,
                        ident.2 = HP2,
                        min.pct = 0.1,
                        logfc.threshold = 0.2,
                        group.by = 'condition',
                        recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    DEGs.name <- paste('sc.deg',ct,HP1,HP2,sep = '_')
    assign(DEGs.name,DEGs)
    deg_corss_condition_age_down[[ct]] <- rownames(DEGs)
  }
}

deg_corss_condition_age_down <- lapply(deg_corss_condition_age_down, function(vec) {
  data.frame(gene = vec, stringsAsFactors = FALSE)
})
rio::export(deg_corss_condition_age_down, file = '25.6.add/deg_corss_condition_age_down.xlsx')

hub<-AnnotationHub::AnnotationHub()
annolist <- AnnotationHub::query(hub, "Mus musculus")
ensdb110 <- hub[["AH113713"]]

plot.list <- list()
library(clusterProfiler)
for (ct in levels(Clean_sct.inte.rm.edc$cell.type.percise.new)[-c(5,6)]) {
  tmp <- base::get(paste('sc.deg',ct,HP1,HP2,sep = '_')) %>% rownames_to_column(var = 'gene')
  gene2entrzid <- bitr(tmp[tmp$avg_log2FC < 0,]$gene, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,' HP up'))
  erich.name <- paste0(ct,'_plot')
  # assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}
patchwork::wrap_plots(c(plot.list),nrow = 3)

plot.list <- list()
for (ct in levels(Clean_sct.inte.rm.edc$cell.type.percise.new)[-c(5,6)]) {
  tmp <- base::get(paste('sc.deg',ct,HP1,HP2,sep = '_')) %>% rownames_to_column(var = 'gene')
  gene2entrzid <- bitr(tmp[tmp$avg_log2FC > 0,]$gene, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,' LP up'))
  erich.name <- paste0(ct,'_plot')
  # assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}
patchwork::wrap_plots(c(plot.list),nrow = 3)

# 12.pseodutime -----
library(monocle3)
library(ggpmisc)

## run mono
Idents(Clean_sct.inte.rm.edc)
getMonocds <- function(obj, assay, slot) {
  cell_types <- levels(Idents(obj))
  
  cds_list <- lapply(cell_types, function(celltype) {
    print(paste("Working on", celltype))
    
    tmp <- subset(obj, cells = WhichCells(obj, idents = celltype))
    mono <- GetAssayData(tmp, assay = assay, slot = slot)
    mono.cell_meta <- tmp@meta.data
    mono.gene_annotation <- data.frame(gene_short_name = rownames(mono))
    rownames(mono.gene_annotation) <- rownames(mono)
    
    cds <- new_cell_data_set(mono,
                             cell_metadata = mono.cell_meta,
                             gene_metadata = mono.gene_annotation) %>% 
      preprocess_cds( num_dim = 30) %>%
      reduce_dimension(preprocess_method = "PCA", cores = 80) %>% 
      cluster_cells(reduction_method = 'UMAP',
                    resolution = 0.001) %>% 
      learn_graph(use_partition = F)
    
    return(cds)
  })
  names(cds_list) <- cell_types
  return(cds_list)
}
sc.cds_list <- getMonocds(obj = Clean_sct.inte.rm.edc,assay = 'RNA',slot = 'counts')
save.image(compress = F)

## then analysis at each cell types
plot.ct <- 'Adipocytes'
cds <- sc.cds_list[[plot.ct]]
p1 <- plot_cells(cds,
                 reduction_method="UMAP", 
                 color_cells_by="condition",
                 label_cell_groups = F,
                 show_trajectory_graph = F,
                 cell_stroke = .2,
                 group_label_size = 5,
                 cell_size = .8) +
  ggtitle(plot.ct) +
  scale_color_manual(values = scales::alpha(colour = c('#48cae4',"#e9c46a"),alpha = 0.7)) + 
  NoAxes()
p1
cds <- order_cells(cds)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", 
                 label_cell_groups = F, 
                 label_leaves = FALSE,  
                 label_branch_points = F,
                 group_label_size = 1,
                 cell_size = .8,
                 label_roots = F,
                 cell_stroke = .2,
                 trajectory_graph_segment_size = 1
) + NoAxes() 

p1 + p2

## bin plot 
cdsmeta <- do.call(cbind, lapply(cds@colData@listData, as.data.frame))
colnames(cdsmeta) <- names(cds@colData@listData)
pseudotime.df <- data.frame(pseudotime = pseudotime(cds)) %>% 
  rownames_to_column(var = 'cells') %>% 
  left_join(y = cdsmeta[,c(21,22)], by = 'cells')

pseudotime.df$bin <- cut(pseudotime.df$pseudotime, 
                         breaks = seq(0, ceiling(max(pseudotime.df$pseudotime)), by = 1), 
                         include.lowest = TRUE, right = FALSE)
result <- pseudotime.df %>%
  group_by(bin, condition) %>%
  summarise(count = n()) %>%
  spread(key = condition, value = count, fill = 0) %>%
  mutate(treat_ratio = HP / (LP + HP)) 

result$bin_numeric <- as.numeric(result$bin)

model <- lm(treat_ratio ~ bin_numeric, data = result)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
R <- cor(result$bin_numeric,result$treat_ratio, method = 'spearman')
p_value <- summary_model$coefficients[2, 4] 

ggplot(result, aes(x = bin_numeric, y = treat_ratio)) +
  geom_point(size = 2.5,color = cluster_colors[1]) + 
  geom_smooth(method = "lm", se = TRUE, color = "gray70",alpha = .2) +  
  annotate("text", x = 5, y = 0.8, 
           label = paste("R = ", round(R, 3), "\np = ", format(p_value, digits = 3)),
           color = "black", size = 5) +  
  labs(x = "Pseudotime", y = "HP cell ratio") + 
  ggtitle(plot.ct) +
  theme_minimal() + 
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18))

## heatmap
library(Mfuzz)
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 40)
genes <- row.names(subset(modulated_genes, q_value <= 0.0001 & morans_I > 0.1))
mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

ck <- clusterData_new(data = mat,cluster.num = 5,method = 'kmeans')
# names(ck)[4] <- 'type'
sasp.geneuse[[1]][match(sasp.geneuse[[1]], rownames(mat))]
pdf(paste0('./25.6.add/monocle3.heatmap.',plot.ct,'.pdf'),height = 10,width = 8,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           markGenes = c('Inha','Cxcl8','Fas','Cxcl12','Cxcl1','Il12b','Il6r','Igf1','Cd38'))
dev.off()
cluster.list <- lapply(ck$cluster.list, function(vec) {
  data.frame(gene = vec, stringsAsFactors = FALSE)
})
names(cluster.list) <- paste0("Cluster_", names(cluster.list))
rio::export(cluster.list, file = paste0('25.6.add/',plot.ct,'.heatmap.cluster.genes.xlsx'))


# 14.NAD pathway -----
library(KEGGREST)
nad.pathway <- keggGet('mmu00760')
nad.pathway <- sapply(seq(2,86,2), function(x){
  gns <- unlist(strsplit(nad.pathway[[1]][['GENE']][x],";"))[1]
})

nad.geneuse <- list(c(nad.pathway))
Clean_sct.inte.rm.edc <- AddModuleScore(Clean_sct.inte.rm.edc, features = nad.geneuse,name = 'nad')

p1 <- FeaturePlot(Clean_sct.inte.rm.edc, features = 'nad1', split.by = 'condition', order = T, min.cutoff = 'q10', cols = c('grey90', "red3"), combine = F)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

nad.score <- Clean_sct.inte.rm.edc@meta.data[,c('cell.type.percise.new','nad1','condition')] %>% na.omit()
nad.score$condition <- factor(nad.score$condition, levels = c('LP','HP'))

ggplot(nad.score, aes(x=condition, y=nad1, fill=condition)) + 
  geom_boxplot(alpha = .7,outlier.size = .5) + 
  theme_bw()+
  labs(x = '', y = 'nad score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")


p1 <- ggpubr::ggviolin(nad.score, x="condition", y="nad1", width = 0.6,
                       add = 'boxplot',
                       add.params = list(fill = 'gray90',alpha = 0.6),
                       color = "black",#轮廓颜色 
                       fill="condition",#填充
                       # palette = "npg",
                       xlab = F, #不显示x轴的标签
                       # bxp.errorbar=T,#显示误差条
                       # bxp.errorbar.width=0.5, #误差条大小
                       size=.2, #箱型图边线的粗细
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.9) + 
  ylab('NAD metabolism score')  + 
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  theme_bw() + 
  # scale_color_npg() +
  facet_wrap(~ cell.type.percise.new,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p2 <- ggpubr::ggviolin(nad.score, x="condition", y="nad1", width = 0.6, 
                       # add = 'boxplot',
                       # add.params = list(fill = 'gray90',alpha = 0.6),
                       color = "black",#轮廓颜色 
                       fill="condition",#填充
                       # palette = "npg",
                       xlab = F, #不显示x轴的标签
                       # bxp.errorbar=T,#显示误差条
                       # bxp.errorbar.width=0.5, #误差条大小
                       size=.2, #箱型图边线的粗细
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.9) + 
  ylab('NAD metabolism score')  + 
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  theme_bw() + 
  # scale_color_npg() +
  facet_wrap(~ cell.type.percise.new,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p1 + p2

nad.score.adi <- nad.score[nad.score$cell.type.percise.new == 'Adipocytes',]
nad.score.adi$nadlog <- log2(nad.score.adi$nad1 + 1)
ggpubr::ggviolin(nad.score.adi, x="condition", y="nad1", width = .5, 
                 add = 'nad.score.adi',
                 add.params = list(fill = 'gray90',alpha = 0.6),
                 color = "black",#轮廓颜色 
                 fill="condition",#填充
                 # scale = 'area',
                 # palette = "npg",
                 xlab = F, #不显示x轴的标签
                 # bxp.errorbar=T,#显示误差条
                 # bxp.errorbar.width=0.5, #误差条大小
                 size=.2, #箱型图边线的粗细
                 outlier.shape=NA, #不显示outlier
                 legend = "right",
                 alpha = 0.9) + 
  ylab('NAD metabolism score')  + 
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  theme_bw() +
  ylim(-0.06,0.08) +
  # scale_color_npg() +
  facet_wrap(~ cell.type.percise.new,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggpubr::ggdensity(nad.score.adi, 
          x = "nad1",
          add = "median", rug = F,
          color = "condition", fill = "condition",
          palette = c('#0B996F',"#D6570D")) +
  ylab("Density") + 
  # scale_x_log10() +
  xlab('NAD metabolism score') + 
  ggtitle('NAD metabolism') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) 
median(nad.score[nad.score$condition == 'LP' & nad.score$cell.type.percise.new == 'Adipocytes',]$NAD1)
median(nad.score[nad.score$condition == 'HP' & nad.score$cell.type.percise.new == 'Adipocytes',]$NAD1)

## AUC
library(AUCell)
expr_mat <- GetAssayData(adipo.seu, assay="SCT", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('NAD' = c(nad.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings)))
adipo.seu$NAD_AUC <- as.numeric(getAUC(cells_AUC)['NAD',])

nad.auc <- adipo.seu@meta.data[,c('NAD_AUC','condition')] %>% na.omit()
table(sasp.auc$condition)

ggplot(nad.auc, aes(x = NAD_AUC, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  ggtitle("Adipocytes NAD") + xlab('AUC') + ylab("") + 
  xlim(-0.03,0.18) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

# 15.metabolism-----

## NC random permutation test
library(pbapply)
library(future)
library(furrr)
plan(multisession, workers=50)

### dataset
kegg.meta <- fgsea::gmtPathways('/data/01.database/05.pathway.geneset/KEGG_metabolism_nc.gmt')
kegg.meta <- lapply(kegg.meta, function(strings) {
  stringr::str_to_title(strings)
})
### function not filter
calculate_Eij <- function(expr_matrix, cell_types) {
  split_cells <- split(1:ncol(expr_matrix), cell_types)
  Eij <- lapply(split_cells, function(cell_idx) {
    rowMeans(expr_matrix[, cell_idx, drop = FALSE])
  })
  Eij <- do.call(cbind, Eij) # 合并为矩阵
  colnames(Eij) <- names(split_cells)
  return(Eij)
} # its can used AverHPExpression
calculate_rij <- function(Eij) {
  rij <- sweep(Eij, 1, rowMeans(Eij), FUN = "/")
  return(rij)
}
calculate_Ptj <- function(rij, gene_sets) {
  Ptj <- sapply(names(gene_sets), function(pathway) {
    genes <- gene_sets[[pathway]]
    valid_genes <- intersect(genes, rownames(rij)) # 过滤存在的基因
    if (length(valid_genes) > 0) {
      apply(rij[valid_genes, , drop = FALSE], 2, mean) # 平均值作为活性分数
    } else {
      rep(NA, ncol(rij)) # 如果没有基因，返回NA
    }
  })
  return(t(Ptj)) # 转置方便查看
}
permutation_test <- function(expr_matrix, cell_types, gene_sets, n_permutations = 5000) {
  observed_Ptj <- calculate_Ptj(calculate_rij(calculate_Eij(expr_matrix, cell_types)), gene_sets)
  
  random_Ptj <- pbapply::pbreplicate(n_permutations, {
    shuffled_types <- sample(cell_types)
    shuffled_Eij <- calculate_Eij(expr_matrix, shuffled_types)
    shuffled_rij <- calculate_rij(shuffled_Eij)
    calculate_Ptj(shuffled_rij, gene_sets)
  })
  
  # 计算p值
  p_values <- sapply(1:nrow(observed_Ptj), function(i) {
    sapply(1:ncol(observed_Ptj), function(j) {
      mean(random_Ptj[i, j, ] >= observed_Ptj[i, j], na.rm = TRUE)
    })
  })
  
  return(list(observed_Ptj = observed_Ptj, p_values = p_values))
}
permutation_test_furrr <- function(expr_matrix, cell_types, gene_sets, n_permutations = 1000) {
  observed_Ptj <- calculate_Ptj(calculate_rij(calculate_Eij(expr_matrix, cell_types)), gene_sets)
  
  # 使用 furrr 进行并行置换
  random_Ptj <- future_map(1:n_permutations, ~ {
    shuffled_types <- sample(cell_types)
    shuffled_Eij <- calculate_Eij(expr_matrix, shuffled_types)
    shuffled_rij <- calculate_rij(shuffled_Eij)
    calculate_Ptj(shuffled_rij, gene_sets)
  }, .progress = TRUE)
  
  random_Ptj <- simplify2array(random_Ptj) # 转换为数组
  
  # 计算p值
  p_values <- sapply(1:nrow(observed_Ptj), function(i) {
    sapply(1:ncol(observed_Ptj), function(j) {
      mean(random_Ptj[i, j, ] >= observed_Ptj[i, j], na.rm = TRUE)
    })
  })
  
  return(list(observed_Ptj = observed_Ptj, p_values = p_values))
}

### function filter
calculate_mean_expression <- function(expr_matrix, cell_type) {
  
  cell_type <- factor(cell_type)
  E_ij <- sapply(levels(cell_type), function(ct) {
    rowMeans(expr_matrix[, cell_type == ct, drop = FALSE], na.rm = TRUE)
  })
  
  colnames(E_ij) <- levels(cell_type) 
  return(E_ij)
}
calculate_rij_filtered <- function(E_ij) {
  
  r_ij <- sweep(E_ij, 1, rowMeans(E_ij, na.rm = TRUE), FUN = "/")
  
  q75 <- apply(r_ij, 1, quantile, probs = 0.75, na.rm = TRUE)
  q25 <- apply(r_ij, 1, quantile, probs = 0.25, na.rm = TRUE)
  upper_limit <- 3 * q75
  lower_limit <- (1 / 3) * q25
  
  r_ij_filtered <- r_ij
  for (i in seq_len(nrow(r_ij))) {
    r_ij_filtered[i, r_ij[i, ] > upper_limit[i] | r_ij[i, ] < lower_limit[i]] <- NA
  }
  
  return(r_ij_filtered)
}
calculate_pathway_score <- function(r_ij_filtered, pathway_genes) {
  
  pathway_scores <- future_map(pathway_genes, function(genes) {
    valid_genes <- intersect(genes, rownames(r_ij_filtered))
    if (length(valid_genes) == 0) return(rep(NA, ncol(r_ij_filtered)))
    pathway_rij <- r_ij_filtered[valid_genes, , drop = FALSE]
    P_tj <- colMeans(pathway_rij, na.rm = TRUE) 
    return(P_tj)
  })
  
  pathway_scores <- do.call(rbind, pathway_scores)
  rownames(pathway_scores) <- names(pathway_genes)
  return(pathway_scores)
}
perform_permutation_test <- function(expr_matrix, cell_type, pathway_genes, n_permutations = 1000) {
  
  E_ij <- calculate_mean_expression(expr_matrix, cell_type)
  r_ij_filtered <- calculate_rij_filtered(E_ij)
  original_scores <- calculate_pathway_score(r_ij_filtered, pathway_genes)
  
  permuted_p_values <- future_map(1:n_permutations, function(i) {
    permuted_cell_type <- sample(cell_type)  # 随机置换细胞类型
    E_ij_perm <- calculate_mean_expression(expr_matrix, permuted_cell_type)
    r_ij_perm <- calculate_rij_filtered(E_ij_perm)
    permuted_scores <- calculate_pathway_score(r_ij_perm, pathway_genes)
    return(permuted_scores)
  })
  
  permuted_p_values <- do.call(abind::abind, list(permuted_p_values, along = 3)) 
  p_values <- apply(original_scores, 1:2, function(orig, perm) {
    mean(perm >= orig, na.rm = TRUE)
  }, perm = permuted_p_values)
  
  return(list(original_scores = original_scores, p_values = p_values))
}

### not use impute
single.celltype <- Clean_sct.inte.rm.edc$cell.type.percise.new
expr_matrix <- GetAssayData(Clean_sct.inte.rm.edc,assay = 'SCT',layer = 'data') %>% as.data.frame()
valid_gene_sets <- lapply(kegg.meta, function(genes) {
  intersect(genes, rownames(expr_matrix))
})
valid_gene_sets <- valid_gene_sets[sapply(valid_gene_sets, length) > 0]
expr_matrix <- expr_matrix[unique(as.vector(unlist(valid_gene_sets))),] %>% na.omit()
expr_matrix.LP <- expr_matrix[,c(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'LP',]$cells)]
expr_matrix.HP <- expr_matrix[,c(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'HP',]$cells)]

#### no filter result
no.impute.Eij <- calculate_Eij(expr_matrix, single.celltype)
no.impute.rij <- calculate_rij(no.impute.Eij)
no.impute.Ptj <- calculate_Ptj(no.impute.rij, valid_gene_sets)
no.impute.nofilter.result <- permutation_test_furrr(expr_matrix, single.celltype, valid_gene_sets)
observed_Ptj <- no.impute.result$observed_Ptj
p_values <- no.impute.result$p_values

#### filter result
no.impute.Eij <- calculate_mean_expression(expr_matrix, single.celltype)
no.impute.rij <- calculate_rij_filtered(no.impute.Eij)
pathway_scores <- calculate_pathway_score(no.impute.rij, valid_gene_sets)
no.impute.result <- perform_permutation_test(expr_matrix, single.celltype, valid_gene_sets, n_permutations = 100)
##### LP
LP.celltype <- Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'LP',]$cell.type.percise.new
names(LP.celltype) <- Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'LP',]$cells
no.impute.LP.result <- perform_permutation_test(expr_matrix.LP, LP.celltype, valid_gene_sets, n_permutations = 100)
view(no.impute.LP.result$original_scores)
no.impute.LP.result$original_scores <- as.data.frame(no.impute.LP.result$original_scores)
no.impute.LP.result$p_values <- as.data.frame(no.impute.LP.result$p_values)
no.impute.LP.result$original_scores[is.na(no.impute.LP.result$original_scores)] <- 0

p1 <- ComplexHeatmap::pheatmap(no.impute.LP.result$original_scores,
                               name = "Activity",
                               border_color = NA, 
                               clustering_method = 'ward.D',
                               show_colnames = T,
                               cluster_cols = F,
                               color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.5)(256),alpha = 1),
                               angle_col = '45',
                               fontsize_row = 8,
                               breaks = (seq(0,3,length.out	= 256))
)
##### HP
HP.celltype <- Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'HP',]$cell.type.percise.new
names(HP.celltype) <- Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'HP',]$cells
no.impute.HP.result <- perform_permutation_test(expr_matrix.HP, HP.celltype, valid_gene_sets, n_permutations = 100)
no.impute.HP.result$original_scores <- as.data.frame(no.impute.HP.result$original_scores)
no.impute.HP.result$p_values <- as.data.frame(no.impute.HP.result$p_values)
no.impute.HP.result$original_scores[is.na(no.impute.HP.result$original_scores)] <- 0

p2 <- ComplexHeatmap::pheatmap(no.impute.HP.result$original_scores,
                               name = "Activity",
                               border_color = NA, 
                               clustering_method = 'ward.D',
                               show_colnames = T,
                               cluster_cols = F,
                               color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.5)(256),alpha = 1),
                               angle_col = '45',
                               fontsize_row = 8,
                               breaks = (seq(0,3,length.out	= 256))
)
p1 + p2

no.impute.all <- cbind(no.impute.LP.result$original_scores, no.impute.HP.result$original_scores)[,c(rbind(1:15,(1:15) + 15))]
colnames(no.impute.all) <- c(paste0('LP_', levels(Idents(Clean_sct.inte.rm.edc))),
                             paste0('HP_', levels(Idents(Clean_sct.inte.rm.edc))))[c(rbind(1:15,(1:15) + 15))]

ComplexHeatmap::pheatmap(no.impute.all,
                         name = "Activity",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         show_colnames = T,
                         cluster_cols = F,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1.5)(256),alpha = 1),
                         angle_col = '45',
                         fontsize_row = 8,
                         breaks = (seq(0,3,length.out	= 256))
)

### use impute
library(scImpute)
HP.seu <- subset(Clean_sct.inte.rm.edc, subset = (condition == 'HP'))
HP.counts <- GetAssayData(HP.seu, assay = 'RNA',slot = 'counts') %>% as.data.frame()
HP.counts <- HP.counts[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()

#### HP
gene_len <- read.table('../04.rnaseq/03.align/A10_TP_R.counts.txt',header = T,sep = '\t')[,c(1,6)]
gene_len <- gene_len %>% 
  left_join(y = mus.gene.info, by = c('Geneid' = 'geneid')) %>% 
  dplyr::filter(biotype == 'protein_coding') 
gene_len <- gene_len[!duplicated(gene_len$symbol),] 
rownames(gene_len) <- gene_len$symbol
gene_len <- gene_len[rownames(HP.counts),] %>% na.omit()
gene_len.use <- as.numeric(as.vector(gene_len$Length))
table(gene_len$symbol == rownames(HP.counts))

HP.counts <- HP.counts[rownames(gene_len),]
labels.HP <- HP.seu$cell.type.percise.new 
HP.counts <- saveRDS(HP.counts, file = './add.analysis/rm.lpts/24.20.newplot/15.scimpute/HP.counts.rds',compress = F)
scimpute('./add.analysis/rm.lpts/24.20.newplot/15.scimpute/HP.counts.rds',
         infile = 'rds',outfile = 'rds',
         out_dir = './add.analysis/rm.lpts/24.20.newplot/15.scimpute/HP_',
         labeled = T, labels = as.vector(labels.HP),
         type = 'count',genelen = gene_len.use,
         drop_thre = 0.5,ncores = 30)

#### run NC
HP.counts <- readRDS('./add.analysis/rm.lpts/24.20.newplot/15.scimpute/HP_scimpute_count.rds')
HP.celltype <- HP.seu@meta.data[colnames(HP.counts),]$cells
names(HP.celltype) <- HP.seu@meta.data[colnames(HP.counts),]$cell.type.percise.new
HP.counts.cpm <- edgeR::cpm(HP.counts) %>% round(2)
HP.counts.cpm[1:10,1:5]

impute.HP.result <- perform_permutation_test(HP.counts.cpm, HP.celltype, valid_gene_sets, n_permutations = 100)
impute.HP.result$original_scores <- as.data.frame(impute.HP.result$original_scores)
impute.HP.result$p_values <- as.data.frame(impute.HP.result$p_values)
impute.HP.result$original_scores[is.na(impute.HP.result$original_scores)] <- 0


#### LP
LP.seu <- subset(Clean_sct.inte.rm.edc, subset = (condition == 'LP'))
LP.counts <- GetAssayData(LP.seu, assay = 'RNA',slot = 'counts') %>% as.data.frame()
LP.counts <- LP.counts[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()
table(gene_len$symbol == rownames(LP.counts))

LP.counts <- LP.counts[rownames(gene_len),]
labels.LP <- LP.seu$cell.type.percise.new 
LP.counts <- saveRDS(LP.counts, file = './add.analysis/rm.lpts/24.20.newplot/15.scimpute/LP.counts.rds',compress = F)
scimpute('./add.analysis/rm.lpts/24.20.newplot/15.scimpute/LP.counts.rds',
         infile = 'rds',outfile = 'csv',
         out_dir = './add.analysis/rm.lpts/24.20.newplot/15.scimpute/LP_',
         labeled = T, labels = as.vector(labels.LP),
         type = 'count',genelen = gene_len.use,
         drop_thre = 0.5,ncores = 30)


## compass:run program in linux
GCs.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'GCs'))
GCs.seu.data <- GetAssayData(GCs.seu, assay = 'RNA',slot = 'data') %>% as.data.frame()
GCs.seu.data <- GCs.seu.data[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()
write.table(GCs.seu.data, file = '/data/02.project/00.other.lab/01.hualun/02.mus.HP.placenta/03.analysis/add.analysis/rm.lpts/24.20.newplot/15.meta.compass/GCs/GCs.seu.data.txt',sep = '\t',row.names = T,col.names = T)

DSCs.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'DSCs'))
DSCs.seu.data <- GetAssayData(DSCs.seu, assay = 'RNA',slot = 'data') %>% as.data.frame()
DSCs.seu.data <- DSCs.seu.data[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()
write.table(DSCs.seu.data, file = '/data/02.project/00.other.lab/01.hualun/02.mus.HP.placenta/03.analysis/add.analysis/rm.lpts/24.20.newplot/15.meta.compass/DSCs/DSCs.seu.data.txt',sep = '\t',row.names = T,col.names = T)

MSCs.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'MSCs'))
MSCs.seu.data <- GetAssayData(MSCs.seu, assay = 'RNA',slot = 'data') %>% as.data.frame()
MSCs.seu.data <- MSCs.seu.data[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()
write.table(MSCs.seu.data, file = '/data/02.project/00.other.lab/01.hualun/02.mus.HP.placenta/03.analysis/add.analysis/rm.lpts/24.20.newplot/15.meta.compass/MSCs/MSCs.seu.data.txt',sep = '\t',row.names = T,col.names = T)

SpTs.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'SpTs'))
SpTs.seu.data <- GetAssayData(SpTs.seu, assay = 'RNA',slot = 'data') %>% as.data.frame()
SpTs.seu.data <- SpTs.seu.data[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()
write.table(SpTs.seu.data, file = '/data/02.project/00.other.lab/01.hualun/02.mus.HP.placenta/03.analysis/add.analysis/rm.lpts/24.20.newplot/15.meta.compass/SpTs/SpTs.seu.data.txt',sep = '\t',row.names = T,col.names = T)

S.TGCs.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'S-TGCs'))
S.TGCs.seu.data <- GetAssayData(S.TGCs.seu, assay = 'RNA',slot = 'data') %>% as.data.frame()
S.TGCs.seu.data <- S.TGCs.seu.data[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()
write.table(S.TGCs.seu.data, file = '/data/02.project/00.other.lab/01.hualun/02.mus.HP.placenta/03.analysis/add.analysis/rm.lpts/24.20.newplot/15.meta.compass/S.TGCs/S.TGCs.seu.data.txt',sep = '\t',row.names = T,col.names = T)

macrop.seu.data <- GetAssayData(macrop.seu, assay = 'RNA',slot = 'data') %>% as.data.frame()
macrop.seu.data <- macrop.seu.data[mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol,] %>% na.omit()
write.table(macrop.seu.data, file = '/data/02.project/00.other.lab/01.hualun/02.mus.HP.placenta/03.analysis/add.analysis/rm.lpts/24.20.newplot/15.meta.compass/macrop/macrop.seu.data.txt',sep = '\t',row.names = T,col.names = T)

# 16.cellchat -----
library(CellChat)

# LP
split.seu <- SplitObject(Clean_sct.inte.rm.edc,split.by = 'condition')
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 

lp.cc <- createCellChat(object = split.seu$LP, meta = split.seu$LP@meta.data, group.by = "cell.type.percise.new")
lp.cc@DB <- CellChatDB.use
lp.cc <- subsetData(lp.cc) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

### HP
hp.cc <- createCellChat(object = split.seu$HP, meta = split.seu$HP@meta.data, group.by = "cell.type.percise.new")
hp.cc@DB <- CellChatDB.use
hp.cc <- subsetData(hp.cc) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

### merge
cc.list <- list(LP = lp.cc, HP = hp.cc)
cc.merge <- mergeCellChat(cc.list, add.names = names(cc.list))
gg1 <- compareInteractions(cc.merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cc.merge, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cc.merge, weight.scale = T,color.use = cc.cols)
netVisual_diffInteraction(cc.merge, weight.scale = T, measure = "weight",color.use = cc.cols)
netVisual_diffInteraction(cc.merge, weight.scale = T, measure = "weight",color.use = cc.cols,sources.use = 'Adipocytes')
netVisual_diffInteraction(cc.merge, weight.scale = T, measure = "weight",color.use = cc.cols,sources.use = 'Macrophages')


### plot
cc.cols <- c(use.cols,npg.cols)[1:13]
names(cc.cols) <- levels(Clean_sct.inte.rm.edc$cell.type.percise.new)
netVisual_heatmap(lp.cc,  color.heatmap = "Reds",color.use = cc.cols)
netVisual_heatmap(hp.cc,  color.heatmap = "Reds")

netVisual_heatmap(lp.cc,  color.heatmap = "Reds",measure = 'weight')
netVisual_heatmap(hp.cc,  color.heatmap = "Reds",measure = 'weight')

gg1 <- netVisual_heatmap(cc.merge,color.use = cc.cols,)
gg2 <- netVisual_heatmap(cc.merge, measure = "weight",color.use = cc.cols)
gg1 + gg2


# 17.DC genest ----
# DC geneset ----
DC.geneset <- rio::import_list('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/DC.geneset.xlsx')
genage.use <- list(c(str_to_title(str_to_lower(DC.geneset[["GenAge"]][["Symbol"]]))))
infm.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["Inflammatory response genes"]]))))
dna_da.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["DNA damage genes"]]))))
dna_re.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["DNA repair genes"]]))))
dna_jc.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["Cell junction genes"]]))))
ros.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["ROS genes"]]))))
nfkb.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["NF-κB pathway genes"]]))))
wnt.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["WNT pathway genes"]]))))
hed.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["Hedgehog pathway genes"]]))))
bmp.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["BMP pathway genes"]]))))
tgf.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["TGF-β pathway genes"]]))))


# DNARE
Clean_sct.inte.rm.edc <- AddModuleScore(Clean_sct.inte.rm.edc, features = dna_re.use, name = 'DNA_da',slot = 'data',assay = 'SCT')

colnames(Clean_sct.inte.rm.edc@meta.data)
DNA_da.score <- Clean_sct.inte.rm.edc@meta.data[,c('cell.type.percise.new','DNA_da1','condition')] %>% na.omit()
DNA_da.score$condition <- factor(DNA_da.score$condition, levels = c('LP','HP'))

ggplot(DNA_da.score, aes(x = DNA_da1, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  ggtitle("DNA DA") + xlab('AUC') + ylab("") + 
  xlim(-0.05,0.08) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

## AUC
library(AUCell)
expr_mat <- GetAssayData(Clean_sct.inte.rm.edc, assay="SCT", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('DNA_DA' = c(dna_da.use[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings)))
Clean_sct.inte.rm.edc$DNA_RE_AUC <- as.numeric(getAUC(cells_AUC)['DNA_DA',])

infm.auc <- Clean_sct.inte.rm.edc@meta.data[,c('DNA_RE_AUC','condition')] %>% na.omit()
table(infm.auc$condition)

ggplot(infm.auc, aes(x = DNA_RE_AUC, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  ggtitle("DNA damage") + xlab('AUC') + ylab("") + 
  # xlim(0.0,0.20) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 



# nfkb
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = nfkb.use, name = 'nfkb.use',slot = 'data')
summary(Age3.hd.8um$nfkb.use1)
SpatialFeaturePlot(Age3.hd.8um,features = 'nfkb.use1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.5)

colnames(Age3.hd.8um@meta.data)
plot.data <- Age3.hd.8um@meta.data[,c(15,28)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="nfkb.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NFKB score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')


Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = nfkb.use, name = 'nfkb.use',slot = 'data')
summary(Young1.hd.8um$genage.score1)
SpatialFeaturePlot(Young1.hd.8um,features = 'nfkb.use1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.5)
colnames(Young1.hd.8um@meta.data)
plot.data <- Young1.hd.8um@meta.data[,c(15,22)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')


KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = nfkb.use, name = 'nfkb.use',slot = 'data')
summary(KO2.hd.8um$nfkb.use1)
SpatialFeaturePlot(KO2.hd.8um,features = 'nfkb.use1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.5)
colnames(KO2.hd.8um@meta.data)
plot.data <- KO2.hd.8um@meta.data[,c(15,22)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')



## wnt
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = wnt.use, name = 'wnt.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = wnt.use, name = 'wnt.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = wnt.use, name = 'wnt.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,29)], Young1.hd.8um@meta.data[,c(19,28)], KO2.hd.8um@meta.data[,c(19,28)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-wnt
DB.wnt <- plot.data[plot.data$sub_region == 'DB',]
DB.wnt %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(wnt.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.wnt, x="sample", y="wnt.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('wnt score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-wnt
JZ.wnt <- plot.data[plot.data$sub_region == 'JZ',]
JZ.wnt %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(wnt.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.wnt, x="sample", y="wnt.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('wnt score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-wnt
LZ.wnt <- plot.data[plot.data$sub_region == 'LZ',]
LZ.wnt %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(wnt.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.wnt, x="sample", y="wnt.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('wnt score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## hed
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = hed.use, name = 'hed.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = hed.use, name = 'hed.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = hed.use, name = 'hed.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,30)], Young1.hd.8um@meta.data[,c(19,29)], KO2.hd.8um@meta.data[,c(19,29)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-hed
DB.hed <- plot.data[plot.data$sub_region == 'DB',]
DB.hed %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(hed.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.hed, x="sample", y="hed.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('hed score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-hed
JZ.hed <- plot.data[plot.data$sub_region == 'JZ',]
JZ.hed %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(hed.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.hed, x="sample", y="hed.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('hed score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-hed
LZ.hed <- plot.data[plot.data$sub_region == 'LZ',]
LZ.hed %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(hed.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.hed, x="sample", y="hed.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('hed score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")


## bmp
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = bmp.use, name = 'bmp.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = bmp.use, name = 'bmp.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = bmp.use, name = 'bmp.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,31)], Young1.hd.8um@meta.data[,c(19,30)], KO2.hd.8um@meta.data[,c(19,30)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-bmp
DB.bmp <- plot.data[plot.data$sub_region == 'DB',]
DB.bmp %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(bmp.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.bmp, x="sample", y="bmp.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('bmp score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-bmp
JZ.bmp <- plot.data[plot.data$sub_region == 'JZ',]
JZ.bmp %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(bmp.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.bmp, x="sample", y="bmp.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('bmp score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-bmp
LZ.bmp <- plot.data[plot.data$sub_region == 'LZ',]
LZ.bmp %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(bmp.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.bmp, x="sample", y="bmp.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('bmp score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## tgf
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = tgf.use, name = 'tgf.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = tgf.use, name = 'tgf.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = tgf.use, name = 'tgf.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,32)], Young1.hd.8um@meta.data[,c(19,31)], KO2.hd.8um@meta.data[,c(19,31)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-tgf
DB.tgf <- plot.data[plot.data$sub_region == 'DB',]
DB.tgf %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(tgf.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.tgf, x="sample", y="tgf.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('tgf score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-tgf
JZ.tgf <- plot.data[plot.data$sub_region == 'JZ',]
JZ.tgf %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(tgf.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.tgf, x="sample", y="tgf.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('tgf score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-tgf
LZ.tgf <- plot.data[plot.data$sub_region == 'LZ',]
LZ.tgf %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(tgf.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.tgf, x="sample", y="tgf.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('tgf score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

# 18.compass ----
adipo.mat <- GetAssayData(adipo.seu,assay = 'SCT',slot = 'data')
write.table(adipo.mat, file = '25.6.add/adipo.mat.txt',quote = F,row.names = T,col.names = T,sep = '\t')

macrop.mat <- GetAssayData(macrop.seu,assay = 'SCT',slot = 'data')
write.table(macrop.mat, file = '25.6.add/macrop.mat.txt',quote = F,row.names = T,col.names = T,sep = '\t')

# 19.GeneSwitches DISCARD-----
library(GeneSwitches)

plot_cells(cds,
           reduction_method="UMAP", 
           color_cells_by="SCT_snn_res.1",
           label_cell_groups = T,
           show_trajectory_graph = F,
           cell_stroke = .2,
           group_label_size = 5,
           cell_size = .3) +
  ggtitle(plot.ct) +
  # scale_color_manual(values = scales::alpha(colour = c('#48cae4',"#e9c46a"),alpha = 0.7)) + 
  NoAxes()

macrop.seu@meta.data$pseudotime <- cds@principal_graph_aux$UMAP$pseudotime
Seurat2GSsce <- function(Seuratrds, clusters){
  Seuratrds <- subset(Seuratrds, seurat_clusters %in% clusters)
  sce <- SingleCellExperiment(assays = List(expdata = Seuratrds@assays$SCT@data))
  colData(sce)$Pseudotime <- Seuratrds@meta.data$pseudotime
  UMAP <- Seuratrds@reductions$umap@cell.embeddings
  reducedDims(sce) <- SimpleList(UMAP = UMAP)
  return (sce)
}

# 根据monocle3的轨迹判断有两条发育路线
sce_p1 <- Seurat2GSsce(rds, c(11,4,13,10,2,1))
sce_p2 <- Seurat2GSsce(rds, c(11,4,13,12,15,6,5))

sce_p1 <- binarize_exp(sce_p1, ncores = 2)
sce_p2 <- binarize_exp(sce_p2, ncores = 2)
# 后面可以比较 sce_p1和sce_p2

# 20.GeneTrajectory ----
library(GeneTrajectory)
macrop.seu <- macrop.seu %>% RunUMAP(reduction = 'pca',dims = 1:30)
p1 <- DimPlot(macrop.seu,group.by = 'condition',cols = c('#0B996F',"#D6570D")) + ggtitle("Macrophages") + NoAxes() 
p2 <- FeaturePlot(macrop.seu,features = 'Cd38',cols = c('gray90','red3'),order = T) + ggtitle("Macrophages") + NoAxes() 
p1 + p2

macrop.seu.sub500 <- subset(macrop.seu, cells = sample(colnames(macrop.seu),1000))
DimPlot(macrop.seu.sub500,group.by = 'condition',cols = c('#0B996F',"#D6570D")) 

DefaultAssay(macrop.seu.sub500) <- 'RNA'
macrop.seu.sub500 <- NormalizeData(macrop.seu.sub500)
macrop.seu.sub500 <- FindVariableFeatures(macrop.seu.sub500, nfeatures = 2500)
all_genes <- VariableFeatures(macrop.seu.sub500)
match('Cd38',all_genes)
expr_percent <- apply(as.matrix(macrop.seu.sub500@assays$RNA$data[all_genes, ]) > 0, 1, sum)/ncol(macrop.seu.sub500)
genes <- all_genes[which(expr_percent > 0.01 & expr_percent < 0.5)]
length(genes)
match('Cd38',genes)

macrop.seu.sub500 <- GeneTrajectory::RunDM(macrop.seu.sub500)
cell.graph.dist <- GetGraphDistance(macrop.seu.sub500, K = 10)
cg_output <- CoarseGrain(macrop.seu.sub500, cell.graph.dist, genes, N = 100)

##
library(plot3D)
library(reticulate)
conda_list()
use_condaenv("/data/00.software/00.sysoftware/anaconda3/envs/rnaseq", required = TRUE)
cal_ot_mat_from_numpy <- reticulate::import('gene_trajectory.compute_gene_distance_cmd')$cal_ot_mat_from_numpy
gene.dist.mat <- cal_ot_mat_from_numpy(ot_cost = cg_output[["graph.dist"]], gene_expr = cg_output[["gene.expression"]], num_iter_max = 50000, show_progress_bar = TRUE)
rownames(gene.dist.mat) <- cg_output[["features"]]
colnames(gene.dist.mat) <- cg_output[["features"]]
dim(gene.dist.mat)

gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))

gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(4,7,7), K = 5)
table(gene_trajectory$selected)

# Extract the ordered list of genes along each gene trajectory
gene_list <- list()
for (i in 1:3){
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  gene_list[[i]] <- genes
}

saveRDS(macrop.seu.sub500,file = 'macrop.seu.sub500.rds',compress = F)
macrop.seu.sub500.back <- macrop.seu.sub500

library(SeuratWrappers)
macrop.seu.sub500 <- RunALRA(macrop.seu.sub500)
macrop.seu.sub500.win <- readRDS('TC_macrop.seu.sub500.rds')
macrop.seu.sub500.win <- AddGeneBinScore(macrop.seu.sub500.win, gene_trajectory, N.bin = 5,assay = 'alra', trajectories = 1:2, reverse = c(F, F))

# Visualize gene bin plots for each gene trajectory
FeaturePlot(macrop.seu.sub500.win, pt.size = 0.05, features = paste0("Trajectory",2,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

gene_trajectory['Cd38',]

gene_list <- list()
for (i in 1:2){
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  gene_list[[i]] <- genes
}

gene_labels <- paste("-----", rownames(gene_embedding))
names(gene_labels) <- rownames(gene_embedding)
genes_selected <- c("Cd38") 
gene_labels[names(gene_labels) %notin% genes_selected] <- ""
par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))
text3D(gene_embedding[,1],
       gene_embedding[,2],
       gene_embedding[,3],  labels = gene_labels,
       add = T, colkey = FALSE, cex = 2)

scatter3D(gene_embedding[,1], gene_embedding[,2], gene_embedding[,3],
          colvar = as.integer(as.factor(gene_trajectory$selected)) - 1,
          col = ramp.col(c(hue_pal()(3))),
          pch=19, cex=1, main="trajectory", theta=45, phi=0)

# 21.USS ----

## uss by GSVA ----
library(GSVA)
library(GSEABase)

uss.geneset <- rio::import('/data/02.project/00.other.lab/01.hualun/06.uterus/mus.aging.geneset.xlsx')
uss.geneset$CA <- uss.geneset$CA %>% tolower() %>% stringr::str_to_title() 
uss.geneset$SE <- uss.geneset$SE %>% tolower() %>% stringr::str_to_title() 

gene_sets <- GeneSetCollection(list(
  GeneSet(na.omit(uss.geneset$SM), setName="SM"),
  GeneSet(na.omit(uss.geneset$CA), setName="CA"),
  GeneSet(na.omit(uss.geneset$GA), setName="GA"),
  GeneSet(na.omit(uss.geneset$SE), setName="SE")
))

expr_mat <- GetAssayData(Clean_sct.inte.rm.edc, slot="data",assay = 'SCT')
ssgsea_res <- gsva(as.matrix(expr_mat),
                   gene_sets,
                   method="ssgsea",
                   ssgsea.norm=TRUE,
                   parallel.sz=40)

Clean_sct.inte.rm.edc[["ssgsea_scores"]] <- CreateAssayObject(ssgsea_res)
meta_add <- t(ssgsea_res) %>% as.data.frame()
Clean_sct.inte.rm.edc <- AddMetaData(Clean_sct.inte.rm.edc, metadata = meta_add)

## add SASP & INFM
gene_sets <- GeneSetCollection(list(
  GeneSet(na.omit(sasp.geneuse[[1]]), setName="SASP"),
  GeneSet(na.omit(infm.geneuse[[1]]), setName="INFM")
))
ssgsea_res <- gsva(as.matrix(expr_mat),
                   gene_sets,
                   method="ssgsea",
                   ssgsea.norm=TRUE,
                   parallel.sz=40)

rownames(Clean_sct.inte.rm.edc@meta.data) <- Clean_sct.inte.rm.edc$cells
meta_add <- t(ssgsea_res) %>% as.data.frame()
Clean_sct.inte.rm.edc <- AddMetaData(Clean_sct.inte.rm.edc, metadata = meta_add)

# ratio
md <- Clean_sct.inte.rm.edc@meta.data
type_meds <- md %>%
  dplyr::group_by(cell.type.percise.new) %>%
  dplyr::summarise(
    med_SM = median(SM, na.rm = TRUE),
    # med_CA = median(CA, na.rm = TRUE),
    # med_GA = median(GA, na.rm = TRUE),
    # med_SE = median(SE, na.rm = TRUE),
    med_SASP = median(SASP, na.rm = TRUE),
    med_INFM = median(INFM, na.rm = TRUE),
    .groups = "drop"
  )

md2 <- md %>%
  left_join(type_meds, by = "cell.type.percise.new") %>%
  dplyr::mutate(
    is_senescent = 
      (SM > med_SM) &
      # (CA > med_CA) &
      # (GA > med_GA) &
      # (SE > med_SE) 
      (SASP > med_SASP) &
      (INFM > med_INFM)
  )

prop_tbl <- md2 %>%
  group_by(cell.type.percise.new, condition) %>%
  dplyr::summarise(
    total_cells     = n(),
    sen_cells       = sum(is_senescent, na.rm = TRUE),
    prop_senescent  = sen_cells / total_cells,
    .groups = "drop"
  )

ggplot(prop_tbl,
       aes(x = condition,
           y = prop_senescent,
           fill = condition)) +                  
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Condition",
    y = "Proportion of Senescent Cells",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )
rio::export(prop_tbl, file = '25.8.add/uss.3geneset.cell.ratio.table.xlsx')

## deg of USS cells (by GSVA) ----
library(DESeq2)
Clean_sct.inte.rm.edc$is_senescent <- md2$is_senescent
agg <- AggregateExpression(
  object     = Clean_sct.inte.rm.edc,
  assay      = "SCT",
  slot       = "counts", 
  group.by   = c("sample", "cell.type.percise.new", "is_senescent")
)
pb_counts <- agg$SCT %>% as.data.frame() %>% round(0)

coldata <- data.frame(
  sample_celltype = colnames(pb_counts),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(
    col = sample_celltype,
    into = c("sample", "cell.type.percise.new", "is_senescent"),
    sep  = "_",
    remove = FALSE
  )
pb_list <- split(seq_len(ncol(pb_counts)), coldata$cell.type)
pseudo_bulk_by_ct <- lapply(pb_list, function(idxs){
  list(
    counts  = pb_counts[, idxs, drop = FALSE],
    coldata = coldata[idxs, ]
  )
})

results_cross <- map(pseudo_bulk_by_ct, function(x){
  mat <- x$counts
  cd  <- x$coldata
  
  # keep <- cd$is_senescent %in% c("TURE", "FALSE")
  # mat  <- mat[, keep]
  # cd   <- cd[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = cd,
    design    = ~ is_senescent
  )
  
  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]
  
  dds <- DESeq(dds, quiet = TRUE)
  # TURE VS FALSE
  res1 <- results(dds,
                  contrast = c("is_senescent","TRUE", "FALSE"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  up1   <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down1 <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  
  list(
    up_intersect   = up1,
    down_intersect = down1
  )
})

summary_df <- imap_dfr(results_cross, ~ tibble(
  cell.type      = .y,
  direction      = c("up","down"),
  n.genes        = c(length(.x$up_intersect),
                     length(.x$down_intersect)),
  genes          = c(paste(.x$up_intersect, collapse = ","),
                     paste(.x$down_intersect, collapse = ","))
))
rio::export(summary_df, file = './25.8.add/uss.3geneset.sn-DEGs.pseudobulk.xlsx')

## uss cell type umap ----
genes <- c("SM","SASP","INFM",'CA','GA','SE')

plot.list1 <- list()
plot.list2 <- list()
c(levels(Clean_sct.inte.rm.edc$cell.type.percise.new))
DefaultAssay(Clean_sct.inte.rm.edc) <- 'RNA'
for (ct in 'Adipocytes') {
  print(ct)
  Idents(Clean_sct.inte.rm.edc) <- Clean_sct.inte.rm.edc$cell.type.percise.new
  sub.seu.1 <- subset(Clean_sct.inte.rm.edc, idents = ct)
  
  # sub.seu.1.split <- SplitObject(sub.seu.1, split.by = 'condition')
  # 
  # sub.seu.1.merge <-  merge(sub.seu.1.split$LP, y = c(sub.seu.1.split$HP))
  # sub.seu.1.merge <-  sub.seu.1.merge %>%
  #   JoinLayers() %>%
  #   NormalizeData() %>% 
  #   FindVariableFeatures() %>%
  #   ScaleData() %>%
  #   RunPCA() %>% 
  #   RunUMAP( reduction = "pca", 
  #            dims = 1:12)
  # sub.seu.1.merge$is_senescent <- factor(sub.seu.1.merge$is_senescent, levels = c('TRUE','FALSE'))
  # ct.plot <- DimPlot(sub.seu.1.merge, group.by = 'is_senescent',cols = c('red3','gray90'),pt.size = .3,order = T) + ggtitle(ct)
  # tmp <- FetchData(sub.seu.1, vars = c('umap_1','umap_2', 'is_senescent'))
  # tmp <- tmp[order(tmp$is_senescent,decreasing = T),]
  # colnames(tmp) <- c('umap_1','umap_2',ct)
  
  # ct.plot <- ggplot(tmp) +
  #   geom_point(aes(x = umap_1, y = umap_2, color = is_senescent), size = 0.2) +
  #   theme_classic() +
  #   labs(color = ct) +
  #   scale_color_manual(values = c('red3','gray90')) +
  #   theme(axis.title = element_blank(),
  #         axis.text = element_blank(),
  #         axis.line = element_blank(),
  #         axis.ticks = element_blank(),
  #         legend.title = element_text(size = 12),
  #         legend.position = 'top',
  #         legend.key.width = unit(.3, "cm"))
  # 
  # ct.name <- paste0(ct, '_plot')
  # plot.list1[[ct.name]] <- ct.plot
  # 
  for (gene in genes) {
    tmp <- FetchData(sub.seu.1, vars = c('umap_1','umap_2', gene))
    tmp <- tmp[order(tmp[[gene]]), ]
    
    if (!is.numeric(tmp[[gene]])) {
      warning(paste0("基因列 '", gene, "' 不是数值型，尝试转换..."))
      tmp[[gene]] <- as.numeric(as.character(tmp[[gene]]))
    }
    
    gene.plot <-
      ggplot(tmp) +
      geom_point(aes(x = umap_1, y = umap_2, color = !!sym(gene)), size = 0.2) +
      scale_color_distiller(
        palette = 'RdBu',
        direction = -1,
        limits = c(NA, max(tmp[[gene]])*0.85),  # 下限自动取数据最小值，上限固定为2
        oob = scales::squish
      ) +
      ggtitle(paste0(ct,'_',gene)) +
      theme_classic() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(size = 12, face = "bold.italic"),
            legend.key.width = unit(.3, "cm"))
    erich.name <- paste0(ct,gene, '_plot')
    plot.list2[[erich.name]] <- gene.plot
  }
  
}
patchwork::wrap_plots(c(plot.list1),nrow = 3)
patchwork::wrap_plots(c(plot.list2),ncol = 3)
rm(list = ls(pattern = 'plot.list'))

# boxplot
plot.data <- Clean_sct.inte.rm.edc@meta.data[,c(22,36:41)]
plot.data.long <- reshape2::melt(plot.data)
ggplot(plot.data.long, aes(x=condition, y=value, fill=condition)) + 
  geom_boxplot(alpha = .7,outlier.size = .5,width = .7) + 
  theme_bw()+
  labs(x = '', y = 'Adipocytes score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  facet_wrap(~variable, scale="free",nrow = 2) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")


# cytotrace ----
library(CytoTRACE2)

adipo.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'Adipocytes'))
dim(adipo.seu)
DimPlot(adipo.seu)

adipo.seu.cyto <- cytotrace2(adipo.seu,is_seurat = T,ncores = 50,species = 'mouse')
ggboxplot(adipo.seu.cyto@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
          color = "black",
          fill="condition",
          xlab = F,
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=.1,
          outlier.shape=NA,
          legend = "right",
          alpha = 0.8) + 
  ylab('Potency score')  + ggtitle('Adipocytes') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  scale_fill_manual(values = con_colors)

library(CytoTRACE)
fibro.top2a.seu <- subset(Clean_sct.inte.rm.edc, cells = rownames(fibro.top2a.seu@meta.data))
dim(fibro.top2a.seu)

fibro.top2a.cyto <- cytotrace2(fibro.top2a.seu,is_seurat = T,ncores = 50,species = 'mouse')
ggboxplot(fibro.top2a.cyto@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
          color = "black",
          fill="condition",
          xlab = F,
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=.1,
          outlier.shape=NA,
          legend = "right",
          alpha = 0.8) + 
  ylab('Potency score')  + ggtitle('Top2a+ Fibroblast') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  scale_fill_manual(values = con_colors)

# 19.young cell type score ----

# SCT
LP.seu <- subset(Clean_sct.inte.rm.edc, condition == 'LP')
table(LP.seu$cell.type.percise.new)
DefaultAssay(LP.seu) <- 'SCT'
Idents(LP.seu) <- LP.seu$cell.type.percise.new
LP.ct.hvgs <- FindAllMarkers(LP.seu, only.pos = T,min.pct = 0.1,logfc.threshold = 0.75,recorrect_umi = F) %>% dplyr::filter(p_val < 0.05)
LP.ct.hvgs <- LP.ct.hvgs[LP.ct.hvgs$avg_log2FC > 1,]

for (ct in as.character(levels(LP.ct.hvgs$cluster))) {
  cat("working on", ct, "\n")
  Clean_sct.inte.rm.edc <- AddModuleScore(Clean_sct.inte.rm.edc, assay = 'SCT',features = list(c(LP.ct.hvgs[LP.ct.hvgs$cluster == ct,]$gene)),name = paste0('LP_',ct,'_score'))
}

LP.score.ct <- data.frame(condition = 0,cell.type = 0, score = 0)
for (ct in as.character(levels(LP.ct.hvgs$cluster))) {
  cat("working on", ct, "\n")
  tmp <- Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$cell.type.percise.new == ct, c('condition','cell.type.percise.new',paste0('LP_',ct,'_score1'))]
  colnames(tmp) <- c('condition','cell.type','score')
  LP.score.ct <- rbind(LP.score.ct,tmp)
}
LP.score.ct <- LP.score.ct[-1,]
LP.score.ct$condition <- factor(LP.score.ct$condition, levels = c('LP','HP'))

ggplot(LP.score.ct, aes(x=condition, y=score, fill=condition)) + 
  geom_boxplot(alpha = .7,outlier.size = .5,width = .7) + 
  theme_bw()+
  labs(x = '', y = 'LP score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = con_colors) +
  facet_wrap(~cell.type, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")
