# 1.in tissue integrate -----
options(future.globals.onReference = "error")

# integrated----

plan(multisession, workers=1)
Clean_qc.merge.filtered <- Clean_qc.merge.filtered %>%
  NormalizeData()

plan(multisession, workers=10)
Clean_qc.merge.filtered <- Clean_qc.merge.filtered %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('mitoRatio',"S.Score", "G2M.Score")) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(resolution = .5, cluster.name = 'unintegrated_clusters') %>%
  RunUMAP(verbose = F,dims = 1:30, reduction = "pca")

Clean_qc.merge.filtered <- RunUMAP(Clean_qc.merge.filtered, umap.method = 'uwot-learn', dims = 1:30 ,metric = "correlation")

DimPlot(Clean_qc.merge.filtered, reduction = "umap", group.by = 'unintegrated_clusters')
DimPlot(Clean_qc.merge.filtered, reduction = "umap", group.by = 'sample') # batch effect is slight

library(reticulate)
use_condaenv(condaenv = "base", conda = '/data/00.software/00.sysoftware/anaconda3/bin/conda')
Sys.setenv(RETICULATE_PYTHON="/data/00.software/00.sysoftware/anaconda3/bin/python")

Clean_qc.merge.filtered <- RunUMAP(Clean_qc.merge.filtered, umap.method = 'umap-learn', dims = 1:30, reduction.name ='py.umap')


Clean_inte <- Clean_qc.merge.filtered %>% 
  IntegrateLayers(method = RPCAIntegration,orig.reduction = "pca", new.reduction = "inte.rpca", k.anchor = 10) %>% 
  FindNeighbors( reduction = "inte.rpca", dims = 1:30) %>% 
  FindClusters(resolution = .5, cluster.name = 'rpca.cluster') %>% 
  RunUMAP( reduction = "inte.rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(Clean_inte, reduction = "umap.rpca", group.by = 'rpca.cluster',label = T)
DimPlot(Clean_inte, reduction = "umap.rpca", group.by = 'sample',label = F)

Clean_inte.merge <- JoinLayers(Clean_inte)
Clean_inte.markers <- FindAllMarkers(Clean_inte.merge, only.pos = T, min.pct = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
Clean_inte.markers.top50 <- Clean_inte.markers %>% group_by(cluster) %>% top_n(n = 50,wt = avg_log2FC)
mad(Clean_inte.merge@meta.data[Clean_inte.merge$rpca.cluster == '0',]$mitoRatio, constant=1)

FeaturePlot(Clean_inte.merge, features = c('Cd79a','Ms4a1','Sdc1','Cd19','Mki67'), reduction = 'umap.rpca',
            cols = c('gray90','red3'),
            order = T,
            min.cutoff = 0,raster = T) 

FeaturePlot(Clean_inte.merge, features = c('Pparg','Ppara', # 09
                                           'Ghr','Tshr', # 09
                                           'Pdgfra','Fbn1', # 6
                                           'Pecam1','Ptprb', # 13 
                                           'Reln','Abcc9', # 13
                                           'Kcnq5','Notch3', # 
                                           'Myh8','Myh11',
                                           'Nkain2','Grid2',
                                           'Hbb-bs','Hba-a1',
                                           'Adgre1','Apoe',  # 1 11 8
                                           'Dock2','Pax5', # B 7 12 
                                           'Krt17','Krt14', # Keratinocytes 14
                                           'Skap1','Ms4a4b', # T memory 4 
                                           'Upk3b','Nkain4','Msln', # Mesothelial  2
                                           'C1qc', 'Mrc1', # marcop 1
                                           'H2-Eb1','Ciita', # DCs 5
                                           'Sirpb1c', # mono 8
                                           'S100a9', # neut 11
                                           'Fscn1', # DC pro 15
                                           ))

Idents(Clean_inte) <- Clean_inte$rpca.cluster
Clean_inte <- RenameIdents(Clean_inte, 
                           '0' = 'Adipocytes',
                           '9' = 'Adipocytes',
                           '6' = 'FAPs',
                           '13' = 'Endothelial cells',
                           '' = 'Fibroblasts',
                           '19' = 'Fibroblasts',
                           '1' = 'Immune cells',
                           '3' = 'Immune cells',
                           '8' = 'Immune cells',
                           '13' = 'Immune cells',
                           '17' = 'Immune cells',
                           '18' = 'Immune cells',
                           # '18' = 'Immune cells',
                           # '16' = 'Immune cells',
                           '4' = 'FAPs',
                           '7' = 'FAPs',
                           '11' = 'FAPs',
                           '12' = 'FAPs',
                           '5' = 'Adipocytes',
                           '10' = 'Adipocytes',
                           '15' = 'SMCs',
                           '14' = 'Low-quanlity cells'
)

Clean_inte$large.ct <- Idents(Clean_inte)
DimPlot(Clean_inte, reduction = "umap",label = T,group.by = 'large.ct',repel = T) + scale_color_simpsons()

Clean_inte.merge <- JoinLayers(Clean_inte)
Clean_inte.markers <- FindAllMarkers(Clean_inte.merge, only.pos = T, min.pct = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
Clean_inte.markers.top20 <- Clean_inte.markers %>% group_by(cluster) %>% top_n(n = 20,wt = avg_log2FC)

lq.tmp <- subset(Clean_inte, idents = 'Low-quanlity cells')
# FeaturePlot(lq.tmp, features = c('nCount_RNA'),max.cutoff = 750,min.cutoff = 500)

VlnPlot(Clean_inte, features = c('nFeature_RNA','nCount_RNA','mitoRatio'),pt.size =  0)

# table(Clean_inte@meta.data[Clean_inte$large.ct == 'Neurons',]$sample)

Clean_inte.merge.remove.lq <- subset(Clean_inte.merge, cells = colnames(lq.tmp), invert = T)
DimPlot(Clean_inte.merge.remove.lq, reduction = "umap",label = T,group.by = 'large.ct') + scale_color_simpsons()

