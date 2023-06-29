library("SingleR")
library("celldex")
library("scuttle")
library("Seurat")
library(ExperimentHub)
library(scRNAseq)
devtools::install_github('dviraran/SingleR')


ref <- celldex::HumanPrimaryCellAtlasData()

results <- SingleR(test = as.SingleCellExperiment(allen_m1c_2019_ssv4), ref = ref, labels = ref$label.main)

allen_m1c_2019_ssv4$singlr_labels <- results$labels

DimPlot(allen_m1c_2019_ssv4, reduction = 'umap', group.by = 'singlr_labels', label = TRUE)

ExperimentHub::listResources()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TabulaMurisData")

library(ExperimentHub)

######

results <- SingleR(test = as.SingleCellExperiment(seurat.intergrated), ref = datasets, labels = listData[["clust_all_neurons"]])

###
counts

#####
library(SingleR)
SingleR::   
singler = CreateSinglerSeuratObject(counts, seurat.integrated.test, project.name,
                                    min.genes = 500, technology, species = "Human", seurat.integrated.test$RNA_snn_res.0.3,
                                    normalize.gene.length = T, min.cells = 2, npca = 20,
                                    regress.out = "nUMI", reduce.seurat.object = T)


seurat.integrated.test$RNA_snn_res.0.3
