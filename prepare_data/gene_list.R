# Combine the datasets to create a list of differentially expressed genes

library(data.table)

PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"

# Load the gene list datasets and only the required columns
edge_r <- data.frame(fread(
  file.path(PATH, "EdgeR.csv"),
  select=c("external_gene_name", "ensembl_gene_id_version")
))

final_gene_list <- read.csv(
  file.path(PATH, "FinalGeneList_1.5_fold_NoTrkA.txt"),
  col.names=c("ensembl_gene_id_version", "external_gene_name"),
  sep="\t"
)

de_seq2 <- data.frame(fread(
  file.path(PATH, "DESeq2.csv"),
  select=c("external_gene_name", "ensembl_gene_id_version")
))

limma_voom <- data.frame(fread(
  file.path(PATH, "limma-voom.csv"),
  select=c("external_gene_name", "ensembl_gene_id_version")
))

# Add genes from https://www.sciencedirect.com/science/article/pii/S2405580821001758
external_gene_name = c("ACE2", "CD147", "PPIA", "PPIB")
ensembl_gene_id_version = c(
  "ENSG00000130234",
  "ENSG00000172270",
  "ENSG00000196262",
  "ENSG00000166794"
)

# Create a dataframe using the two vectors from above
covid_genes = data.frame(external_gene_name, ensembl_gene_id_version)

# Concatenate the dataframes
gene_list <- rbind(edge_r, final_gene_list, de_seq2, limma_voom, covid_genes)

# Remove the period and numbers following it
# e.g. "ENSG00000123.45" -> "ENSG00000123"
gene_list$ensembl_gene_id_version <- gsub("\\..*","", gene_list$ensembl_gene_id_version)

# Remove the duplicated entries
gene_list <- gene_list[!duplicated(gene_list), ]

# Rename column
names(gene_list)[names(gene_list) == "ensembl_gene_id_version"] <- "ensembl_gene_id"

# Save file
write.csv(gene_list, file.path(PATH, "gene_list.csv"), row.names=FALSE)
