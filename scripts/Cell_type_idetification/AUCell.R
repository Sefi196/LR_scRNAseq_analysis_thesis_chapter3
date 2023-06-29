if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AUCell")

library(AUCell)
library(Seurat)
library(GSEABase)


markers <- read.csv("AUcell/PanglaoDB_markers_27_Mar_2020.tsv", sep = "\t")
markers_gluta <- markers[markers$cell.type == "Interneurons",]
genes <- markers_gluta$official.gene.symbol

counts.intergrated <- GetAssayData(object = seurat.intergrated.test, slot = "scale.data")
cell_rankings <- AUCell_buildRankings(counts.intergrated)
cells_AUC <- AUCell_calcAUC(sig, cell_rankings)


cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign=TRUE)


cells_assignment$geneSet$assignment

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
assignmentTable <- assignmentTable[!duplicated(assignmentTable$cell), ]
row.names(assignmentTable) <- assignmentTable$cell
assignmentTable$cell <- NULL
colnames(assignmentTable)[2] <- "geneSet"

head(assignmentTable)

seurat.intergrated.test <- AddMetaData(
  object = seurat.intergrated.test,
  metadata = assignmentTable,
  col.name = 'AUC')


new_cells <- names(which(getAUC(cells_AUC)["geneSet",]>0.065))

length(new_cells)

seurat.intergrated.test$AUC_test <- ifelse(colnames(seurat.intergrated.test) %in% new_cells, "neurons", "otherC")


DimPlot(object = seurat.intergrated.test, group.by = "AUC", label = F)
