library(Seurat)
library(SeuratWrappers)
library(DropletUtils)
library(ggplot2)
library(RColorBrewer)
library(knitr)
library(dplyr)
library(harmony)
library(glmGamPoi)
library(patchwork)
library(cowplot)

#standard workflow steps QC and filtering and normalisation 

# i will use my data that i already have. 
set.seed(02023)



ifnb <- merged_seurat

# QC and filtering
# explore QC

# filter
ifnb
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200)

# standard workflow steps
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'orig.ident')
before 


# run Harmony -----------
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = TRUE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:20,1:20]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = ifnb.harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = ifnb.harmony, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)

# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>% 
  identity()


# visualize 
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'Day')

before|after

#this doesn't seem to give me the best results. will try a different approach i found here https://github.com/immunogenomics/harmony/issues/41



