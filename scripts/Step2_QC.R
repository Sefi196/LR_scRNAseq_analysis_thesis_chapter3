library(scDblFinder)
library(DoubletFinder)
library(scater)
library(scran)
library(BiocSingular)
library(gridExtra)
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(stringr)
library(gridExtra)
library(grid)
library(cowplot)

###########QC and filtering #######

# get data location for filt CSV objects
C1_STC <- read.csv('geneSymbol_C1STC_gene_count.csv', header=T, row.names = 1)
C4_STC <- read.csv('geneSymbol_C4STC_gene_count.csv',header=T, row.names = 1)
C2_Day25 <- read.csv('geneSymbol_C2Day25_gene_count.csv',header=T, row.names = 1)
C5_Day25 <- read.csv('geneSymbol_C5Day25_gene_count.csv',header=T, row.names = 1)
C2_Day55 <- read.csv('geneSymbol_C2Day55_gene_count.csv',header=T, row.names = 1)
C3_Day55 <- read.csv('geneSymbol_C3Day55_gene_count.csv',header=T, row.names = 1)
C3_Day80 <- read.csv('geneSymbol_C3Day80_gene_count.csv',header=T, row.names = 1)
C5_Day80 <- read.csv('geneSymbol_C5Day80_gene_count.csv',header=T, row.names = 1)


### plot UMAPS and sumamry stats from filtered objects 
##UMAP plotting fucntion
#fucntion takes a single cell count matrix -> outputputs sumamry plots and UMAP object and removes doublets.
plot_umap<- function(count.matrix, min.features = 1000, max.features = 999999999, max.counts =10000, min.counts=10000, npc = 20, cluster_res = 0.7, fig_name = '', project = "", MT=10){
  rst_figures <- list()
  rst_table = data.frame()
  
  # init 
  seurat_object <- CreateSeuratObject(counts = count.matrix, project = project, min.cells = 5, min.features=1)
  rst_table <- rbind(rst_table, data.frame("Cells"=dim(seurat_object)[2],
                                           "median feature per Cell"=median(seurat_object$nFeature_RNA), "Meadian reads per feature"= median(seurat_object$nCount_RNA), row.names = paste0('No filter'),check.names = FALSE))
  
  # normalize data
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  plot(VlnPlot(seurat_object, features = c("nFeature_RNA","nCount_RNA","percent.mt")))
  x = median(seurat_object@meta.data[["percent.mt"]])
  
  #remove unwanted cells. below are default settings but you can modify these
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < MT & nCount_RNA < max.counts & nCount_RNA > min.counts) 
  rst_table <- rbind(rst_table, data.frame("Cells"=dim(seurat_object)[2],
                                           "median feature per Cell"=median(seurat_object$nFeature_RNA), "Meadian reads per feature"= median(seurat_object$nCount_RNA), row.names = paste0('filter(min-ft %MT)', min.features,paste0(' '),MT),check.names = FALSE))
  #plot_violin
  plot(VlnPlot(seurat_object, features = c("nFeature_RNA","nCount_RNA","percent.mt")))
  
  vln1 <- VlnPlot(seurat_object, features = c("nFeature_RNA"))
  vln2 <- VlnPlot(seurat_object, features = c("nCount_RNA"))
  vln3 <- VlnPlot(seurat_object, features = c("percent.mt"))
  
  #now you have removed unwanted cells, it is time to normalize the data. By default, Seurat employs a global-scaling normalization method "LogNormalize" that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  #you can alternatively input seurat_object <- NormalizeData(seurat_object) instead of above.
  #we now identify highly variable features in order to determine a subset of features that exhibit high cell-to-cell 
  
  #variation in the dataset.
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  
  #now we apply a linear transformation (scaling) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  #we can visualise both cells and features that define the PCA
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  #this can be visualised using an Elbow Plot
  #rst_figures <- append(rst_figures, ElbowPlot(seurat_object))
  
  #to cluster cells, I have used the number of sig. PCs that I observed in the above plots. The findneighbors function is for constuction of a KNN graph based on euclidean distance in PCA space and refines edge weights between any two cells based on the shared overlap in their local neighborhoods (jaccard similarity). It uses the input of previously defined dimensionality of the dataset.
  seurat_object <- FindNeighbors(seurat_object, dims = 1:npc)
  #now to actually cluster the cells, we apply modularity optimisation techniques (default is Louvain algorithm). The findclusters function contains a resolution parameter which sets the granularity of downstream clustering. Settings are recommended between 0.4-1.2, but may need to be increased for larger datasets.
  seurat_object <- FindClusters(seurat_object, resolution = cluster_res)
  
  #run non-linear dimensional reduction (UMAP/tSNE)
  seurat_object <- RunUMAP(seurat_object, dims = 1:npc)
  
  ### filter out doublets (remember  to modify doublet rate if samples have variable target cells)
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_pbmc <- paramSweep_v3(seurat_object, PCs = 1:20, sct = FALSE)
  sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
  bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
  
  pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- seurat_object@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.039*nrow(seurat_object@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  
  # run doubletFinder 
  seurat_object <- doubletFinder_v3(seurat_object, 
                                    PCs = 1:20, 
                                    pN = 0.25, 
                                    pK = pK, 
                                    nExp = nExp_poi.adj,
                                    reuse.pANN = FALSE, sct = FALSE)
  
  #remove string afet _
  colnames(seurat_object@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(seurat_object@meta.data))
  
  
  #summary doublets
  statsDoublets <- seurat_object@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(median_nCount_RNA = median(nCount_RNA),
              median_nFeature_RNA = median(nFeature_RNA),
              count = n())
  
  # visualize doublets
  doublets <- DimPlot(seurat_object, reduction = 'umap', group.by = "DF.classifications")
  
  seurat_object <- subset(seurat_object, subset = DF.classifications == 'Singlet')
  
  
  #### appending in figure 
  rst_figures <- append(rst_figures, list(DimPlot(seurat_object, reduction = "umap") + labs(color = "cluster \n(from PCA)", title = '') + theme(text = element_text(size = 10))  )) 
  
  
  rst_figures <- append(rst_figures, list(
    FeaturePlot(seurat_object, reduction = "umap", features = 'nCount_RNA')+labs(color = "UMI count",title = '')+ theme(text = element_text(size = 10)),
    FeaturePlot(seurat_object, reduction = "umap", features = 'nFeature_RNA')+labs(color = str_wrap("Feature count (isoform/gene)",15),title = '') + theme(text = element_text(size = 10))
  ))
  
  plot_pc <- ElbowPlot(seurat_object)+labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10))
  plot_umap <- grid.arrange( plot_pc, tableGrob(rst_table), rst_figures[[1]], rst_figures[[2]], rst_figures[[3]], vln1, vln2, vln3,
                             ncol=2, top=textGrob(fig_name))
  
  plot(doublets)
  tbl_sts1 <- tableGrob(statsDoublets)
  grid.newpage()
  grid.draw(tbl_sts1)
  
  #### write out on new page final sumamry stats after all filetreing is complete 
  stats_sumary <- rbind("sanple ID" = project,
                        "Cells"=dim(seurat_object)[2],
                        "median feature per Cell"=median(seurat_object$nFeature_RNA),
                        "Meadian reads per gene/isoform"= median(seurat_object$nCount_RNA),
                        "Max features" = max.features,
                        "Min features" = min.features,
                        "Min Counts" = min.counts,
                        "Max Counts" = max.counts,
                        "MT percentage" = MT,
                        "NPCs" = npc,
                        "median_percentMT_before_filter" = x,
                        "median_percentMT_after_filter" = median(seurat_object@meta.data[["percent.mt"]]))
  
  tbl_sts2<- tableGrob(stats_sumary)
  
  write.table(stats_sumary, file = paste0(project,"_stats.csv")) 
  
  grid.newpage()
  grid.draw(tbl_sts2)
  
  list(plot_umap, seurat_object, 
       statsDoublets, stats_sumary)
  
}

#### Generate QC plots and output umap object

###C1STC
pdf(file = "C1STC_QC.pdf", width = 10, height = 10) 
plotsC1_STC <- plot_umap(C1_STC, min.features = 1000, max.features = 10000, max.counts =100000, min.counts=800, npc = 20, cluster_res = 0.7,fig_name = 'C1_STC (gene counts, Kolf2.1)', project = "C1_STC", MT=11)
dev.off()

C1_STC_umap_object <- plotsC1_STC[[2]]

###C4_STC
pdf(file = "C4_STC_QC.pdf", width = 10, height = 10) 
plotsC4_STC <- plot_umap(C4_STC, min.features = 1000, max.features = 10000, max.counts =100000, min.counts=800,  npc = 20, cluster_res = 0.7,fig_name = 'C4_STC (gene counts, Kolf2.1)', project = "C4_STC", MT=11)
dev.off()

C4_STC_umap_object <- plotsC4_STC[[2]]

###C2Day25
pdf(file = "C2Day25_QC.pdf", width = 10, height = 10) 
plotsC2Day25 <- plot_umap(C2_Day25, min.features = 1000, max.features = 10000, max.counts =100000, min.counts=800, npc = 20, cluster_res = 0.7,fig_name = 'C2_Day25 (gene counts, Kolf2.1)', project = "C2_Day25", MT=11)
dev.off()

C2_Day25_umap_object <- plotsC2Day25[[2]]

###C2Day55
pdf(file = "C2Day55_QC.pdf", width = 10, height = 10) 
plotsC2Day55 <- plot_umap(C2_Day55, min.features = 1000, max.features = 10000, max.counts =100000, min.counts=800, npc = 20, cluster_res = 0.7,fig_name = 'C2_Day55 (gene counts, Kolf2.1)', project = "C2_Day55", MT=11)
dev.off()

C2_Day55_umap_object <- plotsC2Day55[[2]]

###C3Day55
pdf(file = "C3Day55_QC.pdf", width = 10, height = 10) 
plotsC3Day55 <- plot_umap(C3_Day55, min.features = 1000, max.features = 10000, max.counts =100000, min.counts=800, npc = 20, cluster_res = 0.7,fig_name = 'C3_Day55 (gene counts, Kolf2.1)', project = "C3_Day55", MT=11)
dev.off()

C3_Day55_umap_object <- plotsC3Day55[[2]]

###C3Day80
pdf(file = "C3Day80_QC.pdf", width = 10, height = 10) 
plotsC3Day80 <- plot_umap(C3_Day80, min.features = 1000, max.features = 10000, max.counts =100000, min.counts=800, npc = 20, cluster_res = 0.7,fig_name = 'C3_Day80 (gene counts, Kolf2.1)', project = "C3_Day80", MT=11)
dev.off()

C3_Day80_umap_object <- plotsC3Day80[[2]]


###C5Day25
pdf(file = "C5Day25_QC.pdf", width = 10, height = 10) 
plotsC5Day25 <- plot_umap(C5_Day25, min.features = 1500, max.features = 10000, max.counts =100000, min.counts=800, npc = 20, cluster_res = 0.7,fig_name = 'C5_Day25 (gene counts, Kolf2.1)', project = "C5_Day25", MT=11)
dev.off()

C5_Day25_umap_object <- plotsC5Day25[[2]]


###C5Day80
pdf(file = "C5Day80_QC.pdf", width = 10, height = 10) 
plotsC5Day80 <- plot_umap(C5_Day80, min.features = 1000, max.features= 10000, max.counts =100000, min.counts=800, npc = 20, cluster_res = 0.7,fig_name = 'C5Day80 (gene counts, Kolf2.1)', project = "C5_Day80", MT=11)
dev.off()

C5_Day80_umap_object <- plotsC5Day80[[2]]

####### now onto intergration. 