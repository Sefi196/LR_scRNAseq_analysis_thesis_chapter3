library("Seurat")
library('SCINA')

#import counts
counts <- read.csv("seurat.intergrated.norm.counts.csv", row.names =1)

#import signitures
sig <- read.csv("SCINA/example_signatures.csv")

#sig[sig == ""] <- NA
#sig <- sig[, 8:10]
sig <- sig[-c(90:221), ]
sig$Stem.Cells <- NULL
sig$X <- NULL

results = SCINA(exp = counts, signatures = sig, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.99999, sensitivity_cutoff = 1.00, rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')


#View(results$cell_labels)
#View(results$probabilities)

seurat.intergrated.test <- seurat.intergrated
Idents(seurat.intergrated.test) <- results$cell_labels

DimPlot(seurat.intergrated.test, group.by = "ident", label=F) | DimPlot(seurat.intergrated.test, group.by = "RNA_snn_res.0.3") | DimPlot(seurat.intergrated.test, group.by = "customclassif") 

#### lets s
# Define a vector of cell identities that you want to dimplot
cell_subset <- c("radial.glial.cells", "Neural.stem/precursor.cells", "unknown")

# Dimplot only the subset of cell identities
DimPlot(object = seurat.integrated.test, reduction = "umap", group.by = "ident")


subset.test <- subset(x = seurat.intergrated.test, idents = c("Neural.stem.precursor.cells", "GABAergic.neurons", "Glutaminergic.neurons"), invert = FALSE)
DimPlot(object = subset.test, reduction = "umap", group.by = "ident")

seurat.intergrated.test -> seurat.intergrated





