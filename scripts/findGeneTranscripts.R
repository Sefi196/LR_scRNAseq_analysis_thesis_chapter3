# Function to find gene ID and transcripts for a given gene symbol
findGeneTranscripts <- function(geneSymbol, mappingFilePath, countMatrixPath) {
  # Read the mapping file to get gene ID for the provided gene symbol
  mapping <- read.csv(mappingFilePath, header = TRUE, stringsAsFactors = FALSE)
  geneID <- mapping %>% 
    filter(GeneSymbol == geneSymbol) %>% 
    pull(gene_id)
  
  # Read the count matrix from the provided file path
  countMatrix <- read.csv(countMatrixPath, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract the transcript IDs for the given gene ID
  transcripts <- countMatrix[countMatrix[, 2] == geneID, 1]
  
  # Return the gene ID and transcript IDs
  return(list(GeneID = geneID, Transcripts = transcripts))
}

# Find gene ID and transcript IDs for the given gene symbol
result <- findGeneTranscripts(geneSymbol, mappingFilePath, countMatrixPath)

# Print the gene ID and transcript IDs
cat("Gene ID for", geneSymbol, "is:", result$GeneID, "\n")
cat("Transcript IDs for", geneSymbol, "are:", result$Transcripts, "\n")