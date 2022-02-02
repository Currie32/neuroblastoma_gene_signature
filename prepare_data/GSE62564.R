library("biomaRt")
library(data.table)
library(GEOquery)


# Load the patients' data
patients <- load_and_prepare_patient_data()

# Load the patients' expression data
expression_data <- read.table(
  "../data/GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt",
  header=TRUE,
  sep="\t"
)

# Load differentially expressed genes
gene_list <- read.csv("../data/gene_list.csv")

# Create join table for merging the gene_list with the
# gene expression data
join_table <- create_join_table(gene_list)

# Add refseq_mrna values to gene_list
gene_list <- merge(gene_list, join_table, on="ensembl_gene_id_version")

# Add external_gene_name to expression data
expression_data <- merge(
  gene_list[c("refseq_mrna", "external_gene_name")],
  expression_data,
  by.x="refseq_mrna",
  by.y="RefSeqID"
)

expression_data <- prepare_expression_data(expression_data)

patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

write.csv(
  patients,
  "../data/GSE62564.csv",
  row.names=FALSE
)


load_and_prepare_patient_data <- function() {
  
  # Load GEO data
  geo_id <- "GSE62564"
  gse <- getGEO(geo_id)[[1]]
  
  # Get the patient data
  patients <- pData(gse)
  
  # Clear the row names
  rownames(patients) <- c()
  
  # Filter to relevant columns
  patients <- patients[c(
    "title",
    "geo_accession",
    "age:ch1",
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


rename_columns <- function(df) {
  names(df)[names(df) == "title"] <- "sequence_id"
  names(df)[names(df) == "geo_accession"] <- "geo_id"
  names(df)[names(df) == "Sex:ch1"] <- "male"
  names(df)[names(df) == "age:ch1"] <- "age_days"
  names(df)[names(df) == "high risk:ch1"] <- "high_risk"
  names(df)[names(df) == "d_fav_all:ch1"] <- "favourable"
  names(df)[names(df) == "efs bin:ch1"] <- "event_free_survival"
  names(df)[names(df) == "efs day:ch1"] <- "event_free_survival_days"
  
  return(df)
}


correct_data_types <- function(df) {
  
  # Trim unnecessary characters from string
  df$sequence_id <- substr(df$sequence_id, 1, nchar(df$sequence_id) - 4)
  
  # Convert string values to integers
  df$high_risk[df$high_risk == "HR"] <- 1
  df$high_risk[is.na(df$high_risk)] <- 0
  
  # Identify the male patients
  df$male[df$male == "M"] <- 1
  df$male[df$male == "F"] <- 0
  
  # Convert null values to -1, so the feature can be an integer
  df$favourable[is.na(df$favourable)] <- -1
  
  df$age_days <- as.integer(df$age_days)
  df$male <- as.integer(df$male)
  df$high_risk <- as.integer(df$high_risk)
  df$favourable <- as.integer(df$favourable)
  df$event_free_survival <- as.integer(df$event_free_survival)
  df$event_free_survival_days <- as.integer(df$event_free_survival_days)
  
  return(df)
}


create_join_table <- function(gene_list) {
  
  # Connect to ensembl database
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  
  # Create join table for merging the gene_list with the
  # gene expression data
  join_table <- getBM(
    attributes=c("refseq_mrna", "ensembl_gene_id"),
    filters="ensembl_gene_id",
    values=gene_list$ensembl_gene_id,
    mart=ensembl
  )
  
  # Filter out rows missing a refseq_mrna value
  join_table <- join_table[join_table$refseq_mrna != "", ]
  
  return(join_table)
}


prepare_expression_data <- function(expression_data){
  
  # Remove duplicate genes
  expression_data <- expression_data[!duplicated(expression_data$external_gene_name),]
  
  # Select the genes for future assignment
  genes <- expression_data$external_gene_name
  
  # Filter to only the expression features
  expression_data <- expression_data[,grepl("SEQC", colnames(expression_data))]
  
  # Transpose the expression data
  expression_data_t <- transpose(expression_data)
  
  # Label the columns with the gene names
  colnames(expression_data_t) <- genes
  
  # Add sequence_id
  expression_data_t$sequence_id <- colnames(expression_data)
  
  return(expression_data_t)
}
