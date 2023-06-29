library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)

# ...1 Convert to cell_data_set object ------------------------

cds <- as.cell_data_set(seurat.integrated)
cds


# to get cell metadata
colData(cds)
cds <- estimate_size_factors(cds)

#this doesn't work for me as i don't have any meta data
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)

# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- seurat.integrated@meta.data[["orig.ident"]]
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- seurat.integrated@reductions$umap@cell.embeddings



# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'orig.ident',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "Day",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

plot_cells(cds,
           color_cells_by = 'orig.ident',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = NULL) #need to specify the root cluster

p100 <- plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)


p101 <- plot_cells(cds,
           color_cells_by = 'Day',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE) +
  theme(legend.position = "right")


pdf(file = "trajectory_objects.pdf", width = 8, height = 4) 
p100 | p101
dev.off()

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(cds$monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = Day)) +
  geom_boxplot()




# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells_for_GO <- deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK')

write.csv(deg_bcells, "genes_change_with_psedotime.csv")

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

pdf(file = "genes_that_change_with_pseudo.pdf", width = 10, height = 10)
FeaturePlot(seurat.integrated, features = c('FP671120.4' , 'TUBA1A', 'MALAT1',
                                            'ACTB', 'HAPLN1','MT-ND4'), min.cutoff = 'q10')
dev.off()

p102 <- FeaturePlot(seurat.integrated, features = c("BCL11B"), min.cutoff = 'q10')

pdf(file = "BCL11B.pdf", width = 4, height = 4)
p100
dev.off()

# visualizing pseudotime in seurat

seurat.integrated$pseudotime <- pseudotime(cds)
FeaturePlot(seurat.integrated, features = "pseudotime", label = T)

p100 | p101 | p102


########
pdf(file = "Some_interesting_markers.pdf", width = 10, height = 10) 
FeaturePlot(seurat.integrated, features =  c("VIM", "INA", "TUB3", "CDC20", "PAX6", "HES1", "FGFBS",
                                             "TCF7L2", "BCL11B", "NRG1"),  min.cutoff = 'q10')
dev.off()

FeaturePlot(seurat.integrated, features = c('ENSG00000157168.20')) #NRG1

#VIM INA TUB3, CDC20( Cell polifiaration), pax6 HES1 FGFBS, TCF7L2, BCL11B 


