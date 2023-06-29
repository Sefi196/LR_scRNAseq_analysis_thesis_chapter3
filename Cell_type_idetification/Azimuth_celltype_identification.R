library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

datarefs <- AvailableData()

pbmcsca <- RunAzimuth(
  merged_seurat,
  reference="bonemarrowref",
)


DimPlot(pbmcsca, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3)
DimPlot(pbmcsca, group.by = "customclassif", label = TRUE, label.size = 3)

# Extract counts data from Seurat object
counts_data <- as.data.frame(seurat.intergrated@assays$RNA@counts)

# Save counts data to CSV file
write.csv(counts_data, file = "seurat.intergrated.raw.counts.csv", row.names = TRUE)

saveRDS(merged_seurat, "merged_seurat.rds")

######
#Will try and grab both the Day 80 sample together

# merge datasets
merged_seurat_day80 <- merge(C5_Day80_umap_object, y = c( C3_Day80_umap_object),
                       add.cell.ids = c("C3_Day80", "C5_Day80"),
                       project = 'kolf2.1_gene')

saveRDS(merged_seurat_day80, "merged_seurat_Day80.rds")
