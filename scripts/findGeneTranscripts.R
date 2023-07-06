library(Matrix)
library(dplyr)

# Function to find gene ID and transcripts for a given gene symbol
findGeneTranscripts <- function(geneSymbol, mappingFilePath, countMatrix) {
  # Read the mapping file to get gene ID for the provided gene symbol
  mapping <- read.csv(mappingFilePath, header = TRUE, stringsAsFactors = FALSE)
  geneID <- mapping %>% 
    filter(GeneSymbol == geneSymbol) %>% 
    pull(gene_id)
  
  # Read the count matrix from the provided file path
  countMatrix <- read.csv(countMatrixPath, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract the transcript IDs for the given gene ID
  transcripts <- unique(countMatrix[countMatrix[, 2] == geneID, 1])
  
  # Return the gene ID and transcript IDs
  return(list(GeneID = geneID, Transcripts = transcripts))
}

###### Running the fucntion 

# Mapping file path (GeneSymbol to GeneID)
mappingFilePath <- "/data/gpfs/projects/punim1441/genomes/GTF/gencode.v31_gene_ID_Gene_symbol.csv"


# Mapping file path (GeneSymbol to GeneID)
countMatrixPath <- "/data/gpfs/projects/punim1441/Project_Kolf_SCLR/rebase/flames_all_together_supcnt_20/gene_isoform_seurat_objects/geneid_ENSID_raw.all.transcript.counts.csv"


# Find gene ID and transcript IDs for the given gene symbol
result <- findGeneTranscripts("SNAP25", mappingFilePath, countMatrixPath)

#if you want to use them for ploting with feature plot as it stands we need to change _ to -
result$Transcripts <- gsub("_", "-", result$Transcripts)


###can print out the 
# Print the gene ID and transcript IDs
cat("Gene ID for", geneSymbol, "is:", result$GeneID, "\n")
cat("Transcript IDs for", geneSymbol, "are:", result$Transcripts, "\n")

