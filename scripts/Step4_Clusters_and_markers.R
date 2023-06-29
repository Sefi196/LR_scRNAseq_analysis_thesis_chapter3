#Marker genes and clusters 

library(Seurat)
library(tidyverse)
library(dplyr)

# visualize data
clusters <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'seurat_clusters')

condition|clusters

# findAll markers -----------------

p6 <- DimPlot(seurat.integrated, reduction = 'umap', split.by = 'seurat_clusters')

all.markers <- FindAllMarkers(seurat.integrated,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)

Markers_C12 <- all.markers %>%
  dplyr::filter(cluster == 12 & p_val_adj < 0.05)

FeaturePlot(seurat.integrated, features = c('ENSG00000147862.17', 'ENSG00000121897.14', 'ENSG00000115526.11', 'ENSG00000079931.15'),  min.cutoff = 'q10') 

# 'NFIB' 'LIAS' 'CHST10 -> nevous system devlopmet', 'MOXD1'   


#not the correct way of using this function but still intersting
markers_cluster10 <- FindConservedMarkers(seurat.integrated,
                                         ident.1 = 10, 
                                         grouping.var = 'orig.ident')


# so we will use several marker genes to find out what we have. using this website for refercece 
#http://xteam.xbio.top/CellMarker/search.jsp?species=Human&tissue=Brain&cellname=Stem%20cell

#StemCells 

FeaturePlot(seurat.integrated, features = c('ENSG00000181449.4'),  min.cutoff = 'q10') 
#Sox2, #Nestin 

#progenitors
FeaturePlot(seurat.integrated, features = c('ENSG00000105281.12'))

#SLC1A5 (ATB5)



#Found some de genes so thought i could plot them 
FeaturePlot(seurat.integrated, features = c('ENSG00000170035.16')) | p9
UBE2E3

#Neurotransmitter secretion	LIN7B RAB3A
FeaturePlot(seurat.integrated, features = c('ENSG00000104863.12', 'ENSG00000105649.9'),min.cutoff = 'q10')

