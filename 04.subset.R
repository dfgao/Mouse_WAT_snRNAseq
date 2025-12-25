# 01. adipoblast subtypes ----

### umap ----
adipo.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'Adipocytes'))
adipo.seu <- adipo.seu %>% RunUMAP( reduction = "integrated.dr", dims = 1:30)
adipo.seu <- adipo.seu %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1),graph.name = 'SCT_snn')
DimPlot(adipo.seu, reduction = "umap", label = F,group.by = 'SCT_snn_res.0.2') + scale_color_d3() & NoAxes() & ggtitle('')
DimPlot(adipo.seu, reduction = "umap", label = F,group.by = 'SCT_snn_res.0.2',split.by = 'condition') + scale_color_d3() & NoAxes() & ggtitle('')
d3.cols <- pal_d3("category10")(10)
table(adipo.seu$condition)
# adipo.seu.split <- SplitObject(adipo.seu, split.by = 'condition')
# adipo.seu.merge <-  merge(adipo.seu.split$Young, y = c(adipo.seu.split$Aged_9,adipo.seu.split$Aged_13))
# adipo.seu.merge <-  adipo.seu.merge %>%
#   JoinLayers() %>% 
#   FindVariableFeatures() %>% 
#   ScaleData() %>% 
#   RunPCA() %>% 
#   RunUMAP( reduction = "pca", 
#            dims = 1:20)
# adipo.seu.merge$condition <- factor(adipo.seu.merge$condition, levels = c('Young','Aged_9','Aged_13'))
# DimPlot(adipo.seu.merge, reduction = "umap", label = F,group.by = 'RNA_snn_res.0.2') + scale_color_d3() & NoAxes() & ggtitle('')
# DimPlot(adipo.seu.merge, reduction = "umap", label = F,group.by = 'RNA_snn_res.0.2',split.by = 'condition') + scale_color_d3() & NoAxes() & ggtitle('')

### dotplot ----
Idents(adipo.seu) <- adipo.seu$SCT_snn_res.0.2
DefaultAssay(adipo.seu) <- 'SCT'
adipo.deg<- FindAllMarkers(adipo.seu, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
adipo.deg.top10 <- adipo.deg %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
com.dot.new(adipo.seu, feature = c('Acsl1','B3galt2','Cdo1',
                                   'Dock2','Rftn1','Diaph3',
                                   'Nnat','Npr3','Kirrel',
                                   'Igf1','Acaca','Cish',
                                   'Cfd','Hp','Apoe')
            ,groups = "condition",strip.color = d3.cols[1:5])
rio::export(adipo.deg, file = '25.8.add/adipo.cell.type.degs.xlsx')

### heatmap ----
# ht.data <- prepareDataFromscRNA(adipo.seu,
#                                 diffData = adipo.deg,
#                                 showAverage = T,
#                                 assays = 'SCT',slot = 'data',
#                                 group.by = 'cell.type.percise.new',keep.uniqGene = F,
#                                 scale.data = T)
# enrich.go <- enrichCluster(object = ht.data,
#                            OrgDb = org.Mm.eg.db,
#                            type = "BP",
#                            organism = "mmu",
#                            pvalueCutoff = 0.01,
#                            topn = 10,
#                            seed = 1234)
# enrich.go$ratio <- -log(enrich.go$pvalue)
# 
# pdf('./2025.7.add/adipo.cell.type.degs.heatmap.pdf',height = 12,width = 14,onefile = F)
# visCluster(object = ht.data,
#            ht.col.list = list(col_range = c(-4, 0, 4)),
#            plot.type = "both",
#            column_names_rot = 45,
#            show_row_dend = F,
#            # markGenes = markGenes,
#            markGenes.side = "left",
#            annoTerm.data = enrich.go,
#            line.side = "left",
#            cluster.order = c(1:4),
#            go.col = rep(d3.cols[1:4],each = 10),
#            sample.col = d3.cols[1:4],
#            ctAnno.col = d3.cols[1:4],
#            go.size = 8,bar.width = 5,
#            add.bar = T)
# dev.off()
# 
# enrich.go <- enrichCluster(object = ht.data,
#                            OrgDb = org.Mm.eg.db,
#                            type = "BP",
#                            organism = "mmu",
#                            pvalueCutoff = 0.01,
#                            topn = 50,
#                            seed = 1234)
# enrich.go[enrich.go$group == 'C1',]$group <- 'Nalf1'
# enrich.go[enrich.go$group == 'C2',]$group <- 'Celc3b'
# enrich.go[enrich.go$group == 'C3',]$group <- 'Ccl2'
# enrich.go[enrich.go$group == 'C4',]$group <- 'Top2a'
# rio::export(enrich.go, file = '2025.7.add/adipo.cell.type.degs.GOBP.top50.xlsx')

### LTSR ----
res <- run_celltype_composition(
  seu                   = adipo.seu,
  sample_col            = "sample",
  celltype_col          = "SCT_snn_res.0.2",
  treatment_levels      = c("HP",'LP'),
  reps                  = c("1","2"),
  ltsr_vars             = c("Treatment","Rep"),
  FC = 1.2
)
print(res$combined)

Clean_sct.inte.rm.edc$cell.type.minor.add.adipo <- as.character(Clean_sct.inte.rm.edc$cell.type.percise.new)
Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$cell.type.percise.new == 'Adipocytes',]$cell.type.minor.add.adipo <- as.character(adipo.seu$SCT_snn_res.0.2)

res <- run_celltype_composition(
  seu                   = Clean_sct.inte.rm.edc,
  sample_col            = "sample",
  celltype_col          = "cell.type.minor.add.adipo",
  treatment_levels      = c("LP",'HP'),
  reps                  = c("1","2"),
  ltsr_vars             = c("Treatment","Rep"),
  FC = 1.5
)
print(res$combined)

adipo.seu$SCT_snn_res.0.2 <- as.character(adipo.seu$SCT_snn_res.0.2)
plot_cell_fraction_sankey(obj = adipo.seu,
                          condition = 'condition',
                          cell.type = 'SCT_snn_res.0.2',
                          cols = d3.cols)

### ca_degs using findmarkers ----
comlist <- t(combn(c("Young","Aged_9","Aged_13"), 3))
deg_corss_condition_age_up <- list()
Idents(adipo.seu) <- 'cell.type.percise.new'

#### ca_DEG-up
for (ct in unique(adipo.seu$cell.type.percise.new)) {
  cat("==> Processing cell type:", ct, "\n")
  
  sub.seu.1 <- subset(adipo.seu, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    Young <- comlist[gp,1]
    Aged_9 <- comlist[gp,2]
    Aged_13 <- comlist[gp,3]
    DEGs_young_age9 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_9,
                                   ident.2 = Young,
                                   logfc.threshold = 0.1,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) 
    
    DEGs_age9_age13 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_13,
                                   ident.2 = Aged_9,
                                   logfc.threshold = 0.2,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    ca_DEGs_age_up <- intersect(rownames(DEGs_young_age9), rownames(DEGs_age9_age13))
    deg_corss_condition_age_up[[ct]] <- ca_DEGs_age_up
  }
}
DEGs_young_age9['Nme1',]
DEGs_age9_age13['Nme1',]
FeaturePlot(adipo.seu, features = 'Nme1',split.by = 'condition',order = T,min.cutoff = 'q50')

library(purrr)
max_len <- max(map_int(deg_corss_condition_age_up, length))
deg_corss_condition_age_up.df <- map_dfr(deg_corss_condition_age_up, ~ { length(.x) <- max_len; .x })
deg_corss_condition_age_up.df <- list2DF(lapply(deg_corss_condition_age_up, `length<-`, max(lengths(deg_corss_condition_age_up))))
colnames(deg_corss_condition_age_up.df)
rio::export(deg_corss_condition_age_up.df, file = '2025.7.add/adipo.ca_DGEs_age_up_by_findmarkers.xlsx')

## down
deg_corss_condition_age_down <- list()
for (ct in unique(adipo.seu$cell.type.percise.new)) {
  cat("==> Processing cell type:", ct, "\n")
  
  sub.seu.1 <- subset(adipo.seu, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    Young <- comlist[gp,1]
    Aged_9 <- comlist[gp,2]
    Aged_13 <- comlist[gp,3]
    DEGs_young_age9 <- FindMarkers(sub.seu.1,
                                   ident.1 = Young,
                                   ident.2 = Aged_9,
                                   logfc.threshold = 0.1,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) 
    
    DEGs_age9_age13 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_9,
                                   ident.2 = Aged_13,
                                   logfc.threshold = 0.2,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    ca_DEGs_age_down <- intersect(rownames(DEGs_young_age9), rownames(DEGs_age9_age13))
    deg_corss_condition_age_down[[ct]] <- ca_DEGs_age_down
  }
}

max_len <- max(map_int(deg_corss_condition_age_down, length))
deg_corss_condition_age_down.df <- map_dfr(deg_corss_condition_age_down, ~ { length(.x) <- max_len; .x })
deg_corss_condition_age_down.df <- list2DF(lapply(deg_corss_condition_age_down, `length<-`, max(lengths(deg_corss_condition_age_down))))
colnames(deg_corss_condition_age_down.df)
rio::export(deg_corss_condition_age_down.df, file = '2025.7.add/adipo.ca_DGEs_age_down_by_findmarkers.xlsx')

### ca-DEGs from pseudobulk ----
library(stringr)

agg <- AggregateExpression(
  object     = adipo.seu,
  assay      = "RNA",
  slot       = "counts", 
  group.by   = c("sample", "cell.type.percise.new", "condition")
)
pb_counts <- agg$RNA %>% as.data.frame() %>% round(0)

coldata <- data.frame(
  sample_celltype = colnames(pb_counts),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(
    col = sample_celltype,
    into = c("sample", "cell.type.percise.new", "condition"),
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

library(DESeq2)
results_cross <- map(pseudo_bulk_by_ct, function(x){
  mat <- x$counts
  cd  <- x$coldata
  
  keep <- cd$condition %in% c("Young", "Aged-9", "Aged-13")
  mat  <- mat[, keep]
  cd   <- cd[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = cd %>% mutate(condition = factor(condition, 
                                                 levels = c("Young","Aged-9","Aged-13"))),
    design    = ~ condition
  )
  
  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]
  
  dds <- DESeq(dds, quiet = TRUE)
  # Aged_9 vs Young
  res1 <- results(dds,
                  contrast = c("condition","Aged-9","Young"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  # Aged_13 vs Aged_9
  res2 <- results(dds,
                  contrast = c("condition","Aged-13","Aged-9"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  up1   <- dplyr::filter(res1, log2FoldChange > 0.5) %>% pull(gene)
  down1 <- dplyr::filter(res1, log2FoldChange < -0.5) %>% pull(gene)
  up2   <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down2 <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  
  list(
    up_intersect   = intersect(up1,   up2),
    down_intersect = intersect(down1, down2)
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
rio::export(summary_df, file = '2025.7.add//adipo.ca_DGEs_by_pseudobulk.xlsx')
FeaturePlot(adipo.seu,features = 'Blnk',order = T,split.by = 'condition')

### ca_DEG young or aged_9 vs aged13 ----
results_cross <- map(pseudo_bulk_by_ct, function(x){
  mat <- x$counts
  cd  <- x$coldata
  
  keep <- cd$condition %in% c("Young", "Aged-9", "Aged-13")
  mat  <- mat[, keep]
  cd   <- cd[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = cd %>% mutate(condition = factor(condition, 
                                                 levels = c("Young","Aged-9","Aged-13"))),
    design    = ~ condition
  )
  
  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]
  
  dds <- DESeq(dds, quiet = TRUE)
  # Aged_9 vs Young
  res1 <- results(dds,
                  contrast = c("condition","Aged-9","Young"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  # Aged_13 vs Aged_9
  res2 <- results(dds,
                  contrast = c("condition","Aged-13","Young"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  up1   <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down1 <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  up2   <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down2 <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  
  list(
    up_intersect   = intersect(up1,   up2),
    down_intersect = intersect(down1, down2)
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
rio::export(summary_df, file = '2025.7.add//adipo.ca_DGEs_by_pseudobulk_age13-young & age13-age9 jiaoji.xlsx')
FeaturePlot(adipo.seu,features = 'Blnk',order = T,split.by = 'condition')

### cytotrace ----
library(CytoTRACE2)
adipo.top2a.seu <- subset(adipo.seu, subset = (RNA_snn_res.0.2 == '3'))
dim(adipo.top2a.seu)

adipo.top2a.cyto <- cytotrace2(adipo.top2a.seu,is_seurat = T,ncores = 50,species = 'mouse')
ggboxplot(adipo.top2a.cyto@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
          color = "black",
          fill="condition",
          xlab = F,
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=.1,
          outlier.shape=NA,
          legend = "right",
          alpha = 0.8) + 
  ylab('Potency score')  + ggtitle('Top2a+ adipoblast') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  scale_fill_manual(values = con_colors)

adipo.top2a.seu <- subset(Clean_sct.inte.rm.pt, cells = rownames(adipo.top2a.seu@meta.data))
dim(adipo.top2a.seu)

adipo.top2a.cyto <- cytotrace2(adipo.top2a.seu,is_seurat = T,ncores = 50,species = 'mouse')
ggboxplot(adipo.top2a.cyto@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
          color = "black",
          fill="condition",
          xlab = F,
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=.1,
          outlier.shape=NA,
          legend = "right",
          alpha = 0.8) + 
  ylab('Potency score')  + ggtitle('Top2a+ adipoblast') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  scale_fill_manual(values = con_colors)
### pathway ----
library(ggridges)
# library(KEGGREST)
# estrogen.pathway <- keggGet('mmu04915')
# estrogen.pathway <- sapply(seq(2,268,2), function(x){
#   gns <- unlist(strsplit(estrogen.pathway[[1]][['GENE']][x],";"))[1]
# })
# estrogen.geneuse <- list(c(estrogen.pathway))

#### addmodulescore
# NAD

adipo.seu <- AddModuleScore(adipo.seu, features = nad.geneuse,name = 'NAD')
# adipo.seu$condition <- factor(adipo.seu$condition, levels = c('HP','LP'))
adipo.seu@meta.data %>% 
  ggplot( aes(x = NAD1,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "NAD score",
    y = "CAondition",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

# milo ----
library(miloR)
library(scales)
library(SingleCellExperiment)
library(ggbeeswarm)

table(Clean_sct.inte.rm.edc$condition)
milo.seu.add.adipo <- subset(Clean_sct.inte.rm.edc, cells = c(
  rownames(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'HP',] ),
  sample(rownames(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$condition == 'LP',]), 18436)
))

milo.seu.add.adipo <- as.SingleCellExperiment(milo.seu.add.adipo,assay = 'SCT') %>%
  Milo() %>%
  miloR::buildGraph(k = 30, d = 50) %>% 
  makeNhoods(
    prop = 0.2,                 
    k = 30,         
    d=50,refinement_scheme="graph",               
    refined = T)
milo.seu.add.adipo <- countCells(milo.seu.add.adipo, 
                                  meta.data = data.frame(colData(milo.seu.add.adipo)),                     
                                  sample="sample")

milo.traj_design <- data.frame(colData(milo.seu.add.adipo))[,c("sample", "condition")]#分别是重复样本ID和分组
milo.traj_design$sample <- as.factor(milo.traj_design$sample)
milo.traj_design <- distinct(milo.traj_design)
rownames(milo.traj_design) <- milo.traj_design$sample
milo.traj_design$condition <- factor(milo.traj_design$condition, levels = c('LP','HP'))

milo.da_results <- testNhoods(milo.seu.add.adipo,                 
                              design = ~ condition,  
                              fdr.weighting="graph-overlap",          
                              design.df = milo.traj_design)
milo.seu.add.adipo <- buildNhoodGraph(milo.seu.add.adipo)

miloR::plotNhoodGraphDA(milo.seu.add.adipo, milo.da_results, alpha=0.1,) +  
  scale_fill_gradient2(low="#4cc9f0",#修改颜色            
                       mid="#F2F2F2",
                       high="#FF5E5B",                  
                       name="log2FC",                   
                       limits=c(-5,5),                  
                       oob=squish)

milo.da_results <- annotateNhoods(milo.seu.add.adipo, milo.da_results, coldata_col = "cell.type.minor.add.adipo")
milo.da_results$cell.type.minor.add.adipo <- factor(milo.da_results$cell.type.minor.add.adipo, levels = rev(c(seq(0,4), levels(Clean_sct.inte.rm.edc$cell.type.percise.new)[-1])))

plotbee (milo.da_results, alpha = 0.2,
         group.by = "cell.type.minor.add.adipo") +  
  scale_color_gradient2(low="#4cc9f0",           
                        mid="#F2F2F2",                  
                        high="#FF5E5B",                  
                        # limits=c(-5,5),                  
                        oob=squish) + 
  labs(x="", y="Log2 Fold Change") +  
  theme_bw(base_size=10)+  
  theme(axis.text = element_text(colour = 'black')) + 
  scale_size_continuous(range = .1) + ggtitle('HP vs LP')

### mito ----
mitocarta <- rio::import_list('../../Mouse.MitoCarta3.0.xls')
mitocarta.gmt <- GeneSetCollection(
  apply(mitocarta$`C MitoPathways`, 1, function(row) {
    genes <- trimws(strsplit(row["Genes"], ",")[[1]])
    GeneSet(geneIds       = genes,
            setName       = row["MitoPathway"],
            shortDescription = row["MitoPathway.Hierarchy"])
  })
)

mitocarta$`C MitoPathways`$Genes <- as.character(mitocarta$`C MitoPathways`$Genes)
mitocarta.list <- setNames(
  strsplit(mitocarta$`C MitoPathways`$Genes , ",\\s*"),    
  mitocarta$`C MitoPathways`$MitoPathway              
)

# GSVA
gene_sets <- GeneSetCollection(list(
  GeneSet(na.omit(mitocarta.list$OXPHOS), setName="OXPHOS.GSVA"),
  GeneSet(na.omit(mitocarta.list$`Mitochondrial central dogma`), setName="Central.GSVA"),
  GeneSet(na.omit(mitocarta.list$Translation), setName="Translation.GSVA"),
  GeneSet(na.omit(mitocarta.list$Metabolism), setName="Mito.meta.GSVA"),
  GeneSet(na.omit(mitocarta.list$`NAD biosynthesis and metabolism`), setName="Mito.NAD.GSVA"),
  GeneSet(na.omit(mitocarta.list$`Mitochondrial dynamics and surveillance`), setName="Mito.dyna.GSVA")
))

expr_mat <- GetAssayData(adipo.seu, slot="data",assay = 'RNA')
ssgsea_res <- gsva(as.matrix(expr_mat),
                   gene_sets,
                   method="ssgsea",
                   ssgsea.norm=TRUE,
                   parallel.sz=40)

meta_add <- t(ssgsea_res) %>% as.data.frame()
adipo.seu <- AddMetaData(adipo.seu, metadata = meta_add)

## OXPHOS
adipo.seu <- AddModuleScore(adipo.seu, features = list(c(mitocarta.list$OXPHOS)),name = 'OXPHOS')
p1 <- adipo.seu@meta.data %>% 
  ggplot( aes(x = OXPHOS.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  labs(
    x = "Condition",
    y = "OXPHOS GSVA",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

## central
adipo.seu <- AddModuleScore(adipo.seu, features = list(c(mitocarta.list$`Mitochondrial central dogma`)),name = 'Mito_central')
p2 <- adipo.seu@meta.data %>% 
  ggplot( aes(x = Central.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  labs(
    x = "Condition",
    y = "Mitochondrial central dogma GSVA",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

## Translation
adipo.seu <- AddModuleScore(adipo.seu, features = list(c(mitocarta.list$Translation)),name = 'Mito_translation')
p3 <- adipo.seu@meta.data %>% 
  ggplot( aes(x = Translation.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Mitochondrial translation GSVA",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

## Translation
adipo.seu <- AddModuleScore(adipo.seu, features = list(c(mitocarta.list$Metabolism)),name = 'Mito_metabolism')
p4 <- adipo.seu@meta.data %>% 
  ggplot( aes(x = Mito.meta.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Mitochondrial metabolism GSVA",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

## NAD
adipo.seu <- AddModuleScore(adipo.seu, features = list(c(mitocarta.list$`NAD biosynthesis and metabolism`)),name = 'Mito_NAD')
p5 <- adipo.seu@meta.data %>% 
  ggplot( aes(x = Mito.NAD.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "NAD biosynthesis and metabolism GSVA",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

## dynamic
adipo.seu <- AddModuleScore(adipo.seu, features = list(c(mitocarta.list$`Mitochondrial dynamics and surveillance`)),name = 'Mito_dyna')
p6 <- adipo.seu@meta.data %>% 
  ggplot( aes(x = Mito.dyna.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Mitochondrial dynamics and surveillance GSVA",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )
(p1 + p2 + p3) / (p4 + p5 + p6)

adipo.seu@meta.data %>% 
  ggplot( aes(x = Mito.meta.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  labs(
    x = "Condition",
    y = "Mitochondrial central dogma GSVA",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

# mt genes 
mt.genes <- mus.gene.info %>% filter(V3 == "protein_coding") %>% pull(V2) %>% str_subset("^mt-")
adipo.seu <- AddModuleScore(adipo.seu, features = list(c(mt.genes)),name = 'Mito_mt_genes')
adipo.seu@meta.data %>% 
  ggplot( aes(x = Mito_mt_genes1,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Mitochondrial dynamics and surveillance score",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

#### ADD SCORE
plot.list <- list()
for (genelist in c('OXPHOS1','Mito_central1','Mito_translation1','Mito_metabolism1','Mito_NAD1','Mito_dyna1')) {
  plot.score <- adipo.seu@meta.data[,colnames(adipo.seu@meta.data) %in% c('condition', genelist) ]
  colnames(plot.score) <- c('condition','genelist')
  plot.obj <- ggplot(plot.score, aes(x = genelist,
                                     y = condition,
                                     fill = condition)) +                  
    geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
    scale_fill_manual(values = con_colors) + 
    labs(
      x = "Condition",
      y = paste0(genelist,'_addscore'),
      title = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
      panel.grid.major.x = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
    )
  plot.name <- paste0(genelist,'_plot')
  plot.list[[plot.name]] <- plot.obj
}
patchwork::wrap_plots(c(plot.list),nrow = 2)

### uss ----
table(Clean_sct.inte.rm.pt@meta.data[Clean_sct.inte.rm.pt$cell.type.percise.new == 'Stromal cells', ]$cells == adipo.seu$cells)
adipo.seu$is_senescent <- Clean_sct.inte.rm.pt@meta.data[Clean_sct.inte.rm.pt$cell.type.percise.new == 'Stromal cells', ]$is_senescent

prop_tbl <- adipo.seu@meta.data %>%
  as_tibble() %>%                  
  group_by(RNA_snn_res.0.2, condition) %>%
  dplyr::summarise(
    total_cells    = n(),                           
    sen_cells      = sum(is_senescent, na.rm = TRUE),
    prop_senescent = sen_cells / total_cells,
    .groups        = "drop"
  )
ggplot(prop_tbl,
       aes(x = condition,
           y = prop_senescent,
           fill = condition)) +                  
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ RNA_snn_res.0.2, nrow = 3) +
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
rio::export(prop_tbl, file = '2025.8.5.add/adipo.uss.ratio.xlsx')

# deg of USS cells
adipo.seu.sct <- subset(Clean_sct.inte.rm.pt, cells = colnames(adipo.seu))
adipo.seu.sct@meta.data <- adipo.seu@meta.data

agg <- AggregateExpression(
  object     = adipo.seu.sct,
  assay      = "SCT",
  slot       = "counts", 
  group.by   = c("sample", "RNA_snn_res.0.2", "is_senescent")
)
pb_counts <- agg$SCT %>% as.data.frame() %>% round(0)

coldata <- data.frame(
  sample_celltype = colnames(pb_counts),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(
    col = sample_celltype,
    into = c("sample", "RNA_snn_res.0.2", "is_senescent"),
    sep  = "_",
    remove = FALSE
  )
pb_list <- split(seq_len(ncol(pb_counts)), coldata$RNA_snn_res.0.2)
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
rio::export(summary_df, file = './2025.8.5.add/adipo.sn-DEGs.pseudobulk.xlsx')

### compass ----
table(adipo.seu$condition)
adipo.banlance.sub <- subset(adipo.seu, cells = c(
  rownames(adipo.seu@meta.data[adipo.seu$condition == 'HP',] ),
  sample(rownames(adipo.seu@meta.data[adipo.seu$condition == 'LP',]), 1376)
))
table(adipo.banlance.sub$condition)

adipo.mat <- GetAssayData(adipo.banlance.sub,assay = 'SCT',slot = 'data')
# dir.create('25.8.add/compass/adipo')
write.table(adipo.mat, file = '25.8.add/compass/adipo/adipo.mat.txt',quote = F,row.names = T,col.names = T,sep = '\t')

adipo.meta <- adipo.banlance.sub@meta.data[,c(21,22)]
write.table(adipo.meta, file = '25.8.add/compass/adipo/adipo.meta.txt',quote = F,row.names = T,col.names = T,sep = '\t')

# pseudobulk ----
pseudo.adipo <- 


# 02. macrophages ----

macrop.seu <- subset(Clean_sct.inte.rm.edc, subset = (cell.type.percise.new == 'Macrophages'))
macrop.seu <- macrop.seu %>% RunUMAP( reduction = "pca", dims = 1:30)
DimPlot(macrop.seu,split.by = 'condition')
macrop.seu.split <- SplitObject(macrop.seu,split.by = 'condition')
macrop.seu.merge <- merge(x = macrop.seu.split$Young, y = c(macrop.seu.split$Aged_9, 
                                                            macrop.seu.split$Aged_13))
macrop.seu.merge <- macrop.seu.merge %>% 
  JoinLayers() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP( reduction = "pca", 
           dims = 1:20)
macrop.seu.merge$condition <- factor(macrop.seu.merge$condition, levels = c('Young','Aged_9','Aged_13'))
DimPlot(macrop.seu.merge, group.by = 'condition', label = F,cols = con_colors) + ggtitle("") + NoAxes()
VlnPlot(macrop.seu,group.by = 'condition',features = 'Cd38',cols = con_colors,pt.size = 0) + theme(axis.title.x = element_blank())
FeaturePlot(macrop.seu,features = 'Cd38',cols = c('gray90','red3'),order = T,split.by = 'condition') & NoAxes()
FeaturePlot(macrop.seu.merge,features = 'Cd38',cols = c('gray90','red3'),order = T,split.by = 'condition') & NoAxes()

macrop.seu$cd38 <- FetchData(macrop.seu,vars = 'Cd38')
summary(macrop.seu@meta.data[macrop.seu$condition == 'Young',]$cd38)
cd38.fraction <- data.frame(table(macrop.seu@meta.data[macrop.seu$cd38 > 0,]$condition) / table(macrop.seu@meta.data$condition))

ggplot(cd38.fraction,
       aes(x = Var1,
           y = Freq,
           fill = Var1)) +                   
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = con_colors) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Condition",
    y = "Proportion of Cd38+ cells",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

macrop.seu$is_cd38 <- ifelse(macrop.seu$cd38 > 0,'Cd38 +','Cd38-')
DimPlot(macrop.seu, group.by = 'is_cd38',split.by = 'condition') + scale_color_simpsons() + NoAxes()


# 03.Epithelium ----
Epi.sub <- subset(Clean_sct.inte.rm.edc,
                  cells = rownames(Clean_sct.inte.rm.edc@meta.data[Clean_sct.inte.rm.edc$cell.type.percise.new == 'Epithelial cells',]))
Epi.sub <- Epi.sub %>%
  FindNeighbors( reduction = "harmony", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "harmony", 
           dims = 1:30)

DimPlot(Epi.sub, reduction = "umap", label = T,group.by = 'RNA_snn_res.0.2')
DimPlot(Epi.sub, reduction = "umap", label = T,group.by = 'condition')

library(CytoTRACE2)
Epi.sub.cyto <- cytotrace2(Epi.sub,is_seurat = T,ncores = 50,species = 'mouse')
ggboxplot(Epi.sub.cyto@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
          color = "black",
          fill="condition",
          xlab = F,
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=.1,
          outlier.shape=NA,
          legend = "right",
          alpha = 0.8) + 
  ylab('Potency score')  + ggtitle('Epithelium') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  scale_fill_manual(values = con_colors)

### pathway ----
library(ggridges)
# library(KEGGREST)
# estrogen.pathway <- keggGet('mmu04915')
# estrogen.pathway <- sapply(seq(2,268,2), function(x){
#   gns <- unlist(strsplit(estrogen.pathway[[1]][['GENE']][x],";"))[1]
# })
# estrogen.geneuse <- list(c(estrogen.pathway))

#### addmodulescore
# estrogen
est.recp.gene <- getGO("GO:0030520")
Epi.sub <- AddModuleScore(Epi.sub, features = list(c(est.recp.gene$`intracellular estrogen receptor signaling pathway`)),name = 'ESTR')
Epi.sub$condition <- factor(Epi.sub$condition, levels = rev(levels(Epi.sub$condition)))
p1 <- Epi.sub@meta.data %>% 
  ggplot( aes(x = ESTR1,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Estrogen receptor pathway score",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

# progesterone 
prog.recp.gene <- getGO("GO:0050847")
Epi.sub <- AddModuleScore(Epi.sub, features = list(c(prog.recp.gene$`progesterone receptor signaling pathway`)),name = 'PROG')

p2 <- Epi.sub@meta.data %>% 
  ggplot( aes(x = PROG1,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Progesterone receptor pathway score",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )
p1 + p2

#### GSVA
gene_sets <- GeneSetCollection(list(
  GeneSet(na.omit(est.recp.gene$`intracellular estrogen receptor signaling pathway`), setName="ESTR.GSVA"),
  GeneSet(na.omit(prog.recp.gene$`progesterone receptor signaling pathway`), setName="PROG.GSVA")
))

expr_mat <- GetAssayData(Epi.sub, slot="data",assay = 'RNA')
ssgsea_res <- gsva(as.matrix(expr_mat),
                   gene_sets,
                   method="ssgsea",
                   ssgsea.norm=TRUE,
                   parallel.sz=40)

meta_add <- t(ssgsea_res) %>% as.data.frame()
Epi.sub <- AddMetaData(Epi.sub, metadata = meta_add)

# estrogen
p1 <- Epi.sub@meta.data %>% 
  ggplot( aes(x = ESTR.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Estrogen receptor pathway score",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )

# progesterone 
p2 <- Epi.sub@meta.data %>% 
  ggplot( aes(x = PROG.GSVA,
              y = condition,
              fill = condition)) +                  
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = con_colors) + 
  facet_wrap(~ cell.type.percise.new, nrow = 3) +
  labs(
    x = "Condition",
    y = "Progesterone receptor pathway score",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "#EEEEEE", colour = NA),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),legend.position = ''
  )
p1 + p2
