library(Seurat)
library(DESeq2)
library(tidyverse)
library(d3heatmap)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(data.table)
library(dplyr)

#can start with integrated or non integrated data. pseubulk using raw counts/ 
#will use seurat intergrated object. 

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

DefaultAssay(seurat.integrated)

cts <- AggregateExpression(seurat.integrated, 
                           group.by = c('orig.ident'),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA
countdata <- as.matrix(cts)
head(countdata)
#So tutorial continues down bellow but i think a tthis stage all im interested in is DE between groupd and not DE between clusters as many clusters won't be in all groups'
#So i think we can go straint from the cts to DE analysis 

# Normal condition 
Time <- c('early', 'early', 'late', 'late', 'late', 'early', 'early', 'late')
Day <- c('STC', 'Day25', 'Day55', 'Day55', 'Day80', 'STC', 'Day25', 'Day80')
batch <- as.factor(c(2, 2, 2, 3, 3, 3, 5, 5))
diff <- as.factor(c(3, 4, 3, 4, 3, 5, 4, 5))

deseq1<-data.frame(Time,Day,batch,diff)

# Make DESeq dataset
(coldata <- data.frame(row.names=colnames(countdata), deseq1))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~Time + batch + diff)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)


# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- as.data.frame(results(dds, contrast = c("Time","early", "late")))
res

res_filt<- res %>% dplyr::filter(padj < 0.05) %>%
  dplyr::filter(log2FoldChange < -2)

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot", genecol = "black", fitcol = "cyan", finalcol = "blue", legend = TRUE)

# RLD for viewing
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Plot residual p-values
hist(res$pvalue, breaks=50, col="grey")

#Set colours for plotting
# for full set of all 15 samples
mycols <- brewer.pal(8, "Accent")[1:length(unique(Day))]

# Heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
pdf("heatmap-samples.pdf", 18, 18, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[Day], RowSideColors=mycols[batch],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# PCA
rld_pca <- function (rld, intgroup = "Day", ntop = 500, colors=NULL, main="Principal Component Analysis", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1) # change the [X,this] to get other PC components 
  pc1lab <- paste0("PC1: ",as.character(pc1var),"% variance")
  pc2lab <- paste0("PC2: ",as.character(pc2var),"% variance")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx, offset=1)) #makes sample labells
  #legend(legendpos, legend=levels(fac), col=colors, pch=6)
  legend("topright", legend=levels(fac), col=colors, pch = 6,
         xpd=TRUE, horiz=FALSE, bty="n"
  )
}
pdf("pca-genes 5y_1-2.pdf", 6, 6, pointsize=13)
rld_pca(rld, colors=mycols, intgroup="Day", xlim=c(-150, 150), ylim=c(-150, 150))
dev.off()

# MA Plot
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="blue", pch=20, cex=1))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
pdf("diffexpr-maplot-0.05.pdf", 18, 18, pointsize=20)
maplot(res, main="MA Plot")
dev.off()

# Volcano Plot
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, xlab="log2(Fold Change)", legendpos="topright", labelsig=FALSE, textcx=1.5, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlab=xlab, cex.axis=1.8, cex.lab=1.5, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p-adj<",sigthresh,sep=""), paste("|log2(FC)|>",lfcthresh,sep=""), "both"), cex=1.5, pch=20, col=c("blue","orange","green"))
}
pdf("diffexpr-volcanoplot-hi-res.pdf", 18, 18, pointsize=20)
volcanoplot(res, lfcthresh=2, sigthresh=0.05, legendpos="bottomright")
dev.off()





#########


# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)
cts.t[1:8,1:20]

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

cts.split$`5`[1:15,1:15]

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

cts.split.modified$Day25[1:10,1:10]

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')

#Let's run DE analysis with B cells
# 1. Get counts matrix
counts_bcell <- cts.split.modified$`B cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_bcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('stim', samples), 'Stimulated', 'Control')) %>%
  column_to_rownames(var = 'samples')

# get more information from metadata



