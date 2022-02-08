library(data.table)
library(GEOquery)


path <- "~/Imperial/neuroblastoma_gene_signature/data/"

# Load differentially expressed genes
gene_list <- read.csv(file.path(path, "gene_list.csv"))

# Load the patients' data
geo_id <- "GSE49711"
patients <- load_and_prepare_patient_data(geo_id)

gse <- getGEO(geo_id)[[1]]

expression_data <- read.csv("~/Downloads/GSE49711_SEQC_NB_TUC_G_log2.txt", sep='\t')

expression_data$X00gene_id

expression_data2 <- expression_data[expression_data$X00gene_id %in% gene_list$external_gene_name,]
expression_data2 <- expression_data2[,2:length(expression_data)]

expression_data2 <- expression_data2[!duplicated(expression_data2$X00gene_id),]

# Select the genes for future assignment
genes <- expression_data2$X00gene_id

# Filter to only the expression features
expression_data2 <- expression_data2[,grepl("SEQC", colnames(expression_data2))]

# Undo log2 transformation and Z-transform the expression values
expression_data2 <- data.frame(scale(2^expression_data2))

# Transpose the expression data
expression_data2_t <- transpose(expression_data2)

# Label the columns with the gene names
colnames(expression_data2_t) <- genes

# Add sequence_id
expression_data2_t$sequence_id <- colnames(expression_data2)

patients <- pData(gse)

merge(
  patients,
  expression_data2_t,
  on="sequence_id",
)




load_and_prepare_patient_data <- function(geo_id) {
  
  # Load GEO data
  gse <- getGEO(geo_id)[[1]]
  
  # Get the patient data
  patients <- pData(gse)
  
  # Clear the row names
  rownames(patients) <- c()
  
  # Filter to relevant columns
  patients <- patients[c(
    "title",
    "geo_accession",
    "age at diagnosis:ch1",
    "Sex:ch1",
    "d_fav_all:ch1",
    "efs bin:ch1",
    "efs day:ch1",
    "high risk:ch1"
  )]
  
  patients <- rename_columns(patients)
  patients <- correct_data_types(patients)
  
  return(patients)
}


