library(dplyr)
library(Seurat)

# get data location
C1_STC <- read.csv('../C5Day80_transcript_count.csv',header=T, row.names = 1)
C1_STC <- subset(C1_STC, select = -c(gene_id))

sum(C1_STC)

###############
###Im concerend that the objects now are huge in comaprsion to running flames indivadually. 
genes.percent.expression <- rowMeans(C1_STC>0 )*100   #percent of cells expressing each gene
df <- as.data.frame(genes.percent.expression)
boxplot(df)

C1_STC_seurat <- CreateSeuratObject(C1_STC)
VlnPlot(C1_STC_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#However, you can not filter out certain genes unless you create a new Seurat object, like this.
#Filter
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])  #select isoforms expressed in at least 1% of cells
counts.sub <- C1_STC[genes.filter,]

sum(counts.sub)

####create another seurat object to see if the distribution is ok
C1_STC_seurat.filt <- CreateSeuratObject(counts.sub)
VlnPlot(C1_STC_seurat.filt, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#### write out filtered matrix
write.csv(counts.sub, "filt_C5Day80_transcript_count.csv")

#plot filetrind count matrix
filt.genes.percent.expression <- rowMeans(counts.sub>0 )*100   
df.2 <- as.data.frame(filt.genes.percent.expression)

#plot boxplot of filtered and raw. 
boxplot(df) 
boxplot(df.2)



#not loving this even with minsupcount of 20. 
#lets check how the C1STC sample ran on its own 

#C1_STC_old <- read.csv('../C1_STC/flames_outs/C1_STC_transcript_count.csv',header=T, row.names = 1)
#C1_STC_old <- subset(C1_STC_old, select = -c(gene_id))

#sum(C1_STC_old)

#genes.percent.expression <- rowMeans(C1_STC_old>0 )*100   #percent of cells expressing each gene
#df <- as.data.frame(genes.percent.expression)
#boxplot(df)

########
#lets see what the data should look like for just known isoforms
#Known <- C1_STC[grep("ENST", row.names(C1_STC)), ]
#novel <- C1_STC[grep("ENSG", row.names(C1_STC)), ]

#genes.percent.expression <- rowMeans(novel>0 )*100   #percent of cells expressing each gene
#df <- as.data.frame(genes.percent.expression)
#boxplot(df) 
#summary(df)


