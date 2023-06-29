library(dplyr)

####change ENG1D to gene symbol so i can use packages for idetifaction of cell types


# Load the reference table
ref_table <- read.csv("/data/gpfs/projects/punim1441/genomes/GTF/gencode.v31_gene_ID_Gene_symbol.csv", header = TRUE)

# Load the data object with ENSG row names
data_obj <- read.csv('../C5Day80_gene_count.csv', header=T)

merged.file = merge(data_obj, ref_table, by.x='gene_id', by.y='gene_id')

#remove first col
merged.file <- merged.file[, -1]

#see how many duplicated items
sum(duplicated(merged.file$GeneSymbol))

#### remove duplicated rows 
duplicates <- duplicated(merged.file$GeneSymbol)

# Subset the data frame to exclude the duplicated rows
df_unique <- subset(merged.file, !duplicates)
#df_duplicates <- subset(merged.file, duplicates)

#remove first col
#df_duplicates <- df_duplicates[, -1]

# Set the last column as row names
row.names(df_unique) <- df_unique[, ncol(df_unique)]
df_unique <- df_unique[, -ncol(df_unique)]

#write out new data frame
write.csv(df_unique, "geneSymbol_C5Day80_gene_count.csv")

###if interested in duplicated rows use this to investifate the data frame 
### i foind the majority of duplicates are on the Y chromosme and constiute about 2% f the total data genes in my data.
### could check the number of counts -> looks like its less than .1% of counts so im not going to worry about it.
#sum(df_duplicates)
#sum(df_unique)

#### i removed them as i think this is a small proportion but it could be possible to change the ids but then they won't match across the samples 
#dup_rows <- merged.file[duplicated(merged.file[["GeneSymbol"]]), ]
#dup_rows <- dup_rows %>% select(c('gene_id',"GeneSymbol"))
#########