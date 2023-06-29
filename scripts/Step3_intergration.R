# script to integrate scRNA-Seq datasets to correct for batch effects
# setwd("~/Desktop/demo/single_cell_integrate")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(scuttle)

# get data location
C1_STC <- read.csv('data/C1_STC_gene_counts.csv',header=T, row.names = 1)
C4_STC <- read.csv('data/C4_STC_gene_count.csv',header=T, row.names = 1)
C2_Day25 <- read.csv('data/C2_Day25_gene_counts.csv',header=T, row.names = 1)
C5_Day25 <- read.csv('data/C5_Day25_gene_count.csv',header=T, row.names = 1)
C2_Day55 <- read.csv('data/C2_Day55_gene_count.csv',header=T, row.names = 1)
C3_Day55 <- read.csv('data/C3_Day55_gene_counts.csv',header=T, row.names = 1)
C3_Day80 <- read.csv('data/C3_Day80_gene_count.csv',header=T, row.names = 1)
C5_Day80 <- read.csv('data/C5_Day80_gene_counts.csv',header=T, row.names = 1)

#Make serat objects not filterd
C1_STC_umap_object <- CreateSeuratObject(counts = C1_STC, project = "C1_STC", min.features = 1000,max.feature = 3000)
C4_STC_umap_object <- CreateSeuratObject(counts = C4_STC, project = "C4_STC",  min.features = 1000,max.feature = 3000)
C2_Day25_umap_object <- CreateSeuratObject(counts = C2_Day25, project = "C2_Day25", min.features = 1000,max.feature = 3000)
C2_Day55_umap_object <- CreateSeuratObject(counts = C2_Day55, project = "C2_Day55", min.features = 700,max.feature = 3000)
C5_Day25_umap_object <- CreateSeuratObject(counts = C5_Day25, project = "C5_Day25", min.features = 1000,max.feature = 3000)
C3_Day55_umap_object <- CreateSeuratObject(counts = C3_Day55, project = "C3_Day55", min.features = 500,max.feature = 3000)
C3_Day80_umap_object <- CreateSeuratObject(counts = C3_Day80, project = "C3_Day80", min.features = 500,max.feature = 3000)
C5_Day80_umap_object <- CreateSeuratObject(counts = C5_Day80, project = "C5_Day80", min.features = 400,max.feature = 3000)

#start from here if doublets have been removed 

# merge datasets
merged_seurat <- merge(C1_STC_umap_object, y = c(C4_STC_umap_object, C2_Day25_umap_object,C5_Day25_umap_object, C2_Day55_umap_object, C3_Day55_umap_object, C5_Day80_umap_object, C3_Day80_umap_object),
                       add.cell.ids = c("C1_STC", "C4_STC", "C2_Day25","C5_Day25", "C2_Day55", "C3_Day55", "C3_Day80", "C5_Day80"),
                       project = 'kolf2.1_gene')

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column to makw a batch col
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('batch', 'Day', 'Barcode'), 
                                    sep = '_')

#Make some extra meta data cols 
View(merged_seurat@meta.data)

#Would do mt content and cell cycle content here if a could
#Something to mention in the discussion and limitations of the project


#QC filtering together
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Calculate outliers and filter merged data
ouliers_Counts <- isOutlier(merged_seurat@meta.data[["nCount_RNA"]], type = "both", log = TRUE, batch=merged_seurat@meta.data[["batch"]])
ouliers_features <- isOutlier(merged_seurat@meta.data[["nFeature_RNA"]], type = "both", log = TRUE, batch=merged_seurat@meta.data[["batch"]])


#havn't worked out exactly how to filter for each batch
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800)
#QC filtering together
VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


######
#SCT code would go here if using 
###

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered) # if using SCT dont run this
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered) # if using SCT dont run this
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered) # if using SCT dont run this
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.1)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.3)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.5)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.7)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.9)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


#plots
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'orig.ident')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'batch')
p3 <- FeaturePlot(merged_seurat_filtered, reduction = 'umap', features = 'nCount_RNA')
p4 <- FeaturePlot(merged_seurat_filtered, reduction = "umap", features = 'nFeature_RNA')

pdf(file = "merged_objects.pdf", width = 6, height = 6) 
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off()


# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'orig.ident')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], nfeatures = 2000)
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list, )

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)



p5 <- DimPlot(seurat.integrated, reduction = 'umap', split.by = 'Day', group.by = "orig.ident")
p6 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'seurat_clusters')
p7 <- FeaturePlot(seurat.integrated, reduction = "umap", features = 'nCount_RNA')
p8 <- FeaturePlot(seurat.integrated, reduction = "umap", features = 'nFeature_RNA')


#can plot based on day
p9 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Day')
p10 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'batch')

pdf(file = "intergrated_objects.pdf", width = 10, height = 10) 
grid.arrange(p5, p9, p7, p8, ncol = 2, nrow = 2)
dev.off()

grid.arrange(p9, p10, ncol = 2, nrow = 2)


p11 <-DimPlot(seurat.integrated, reduction = 'umap', group.by = "RNA_snn_res.0.7")


p5 | p11

VlnPlot(seurat.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

FeaturePlot(seurat.integrated, features = c('ENSG00000164649.20', "ENSG00000127152.18"))


###save intergrated file
saveRDS(seurat.integrated, file = "seurat.intergrated.rds")

###### SCT code 
#######
#here will try an normalize with SCT as i can remove batch effects and then i will see how it looks

merged_seurat_filtered.list <- SplitObject(merged_seurat, split.by = "orig.ident")


for (i in 1:length(merged_seurat_filtered.list)) {
  merged_seurat_filtered.list[[i]] <- SCTransform(merged_seurat_filtered.list[[i]], method = "glmGamPoi", verbose = T) # have removed vars.to.regress = c("batch"), didn't seem to work but might go back and trouble shoot
}

options(future.globals.maxSize = 8000 * 1024^2)
merged_seurat_filtered.features <- SelectIntegrationFeatures(object.list = merged_seurat_filtered.list, nfeatures = 2000)
merged_seurat_filtered.list <- PrepSCTIntegration(object.list = merged_seurat_filtered.list, anchor.features = merged_seurat_filtered.features, 
                                                  verbose = FALSE)

merged_seurat.anchors <- FindIntegrationAnchors(object.list = merged_seurat_filtered.list, normalization.method = "SCT", 
                                                anchor.features = merged_seurat_filtered.features, verbose = FALSE)
merged_seuratSCT.integrated <- IntegrateData(anchorset = merged_seurat.anchors, normalization.method = "SCT", 
                                             verbose = FALSE)

merged_seuratSCT.integrated <- RunPCA(merged_seuratSCT.integrated, verbose = FALSE)
merged_seuratSCT.integrated <- RunUMAP(merged_seuratSCT.integrated, dims = 1:30)
merged_seuratSCT.integrated <- FindNeighbors(object = merged_seuratSCT.integrated, dims = 1:20)
merged_seuratSCT.integrated <- FindClusters(object = merged_seuratSCT.integrated, resolution = 0.1)
merged_seuratSCT.integrated <- FindClusters(object = merged_seuratSCT.integrated, resolution = 0.3)
merged_seuratSCT.integrated <- FindClusters(object = merged_seuratSCT.integrated, resolution = 0.5)
merged_seuratSCT.integrated <- FindClusters(object = merged_seuratSCT.integrated, resolution = 0.7)
merged_seuratSCT.integrated <- FindClusters(object = merged_seuratSCT.integrated, resolution = 0.9)


#plots
p1 <- DimPlot(merged_seuratSCT.integrated, reduction = 'umap', group.by = 'orig.ident')
p2 <- DimPlot(merged_seuratSCT.integrated, reduction = 'umap', group.by = 'Day')
p3 <- FeaturePlot(merged_seuratSCT.integrated, reduction = 'umap', features = 'nCount_RNA')
p4 <- FeaturePlot(merged_seuratSCT.integrated, reduction = "umap", features = 'nFeature_RNA')
p5 <- DimPlot(merged_seuratSCT.integrated, reduction = 'umap', group.by = 'seurat_clusters')

grid.arrange(p1, p2, p3, p4, p5)
######

FeaturePlot(seurat.integrated, features = c("ENSG00000067225.18", "ENSG00000092841.18"),  min.cutoff = 'q10')
#PKM MYL6

