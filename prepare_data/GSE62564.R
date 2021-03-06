# Prepare GSE62564 for analysis

library("biomaRt")
library(GEOquery)


# Update path to data as needed
PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


create_join_table <- function(gene_list) {
  #' Create a join table between the list of differentially expressed genes
  #' and the gene expression data.
  #' 
  #' gene_list List(str): a list of ensembl gene IDs
  #' 
  #' return data.frame: the join table
  
  # Connect to ensembl database
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  
  # Create join table for merging the gene_list with the
  # gene expression data
  join_table <- getBM(
    attributes=c("refseq_mrna", "ensembl_gene_id"),
    filters="ensembl_gene_id",
    values=gene_list,
    mart=ensembl
  )
  
  # Filter out rows missing a refseq_mrna value
  join_table <- join_table[join_table$refseq_mrna != "", ]
  
  return(join_table)
}


load_and_prepare_expression_data <- function(gene_list){
  #' Load and prepare the GEO expression data
  #' 
  #' gene_list data.frame: contains information about the list of
  #'                       differentially expressed genes
  #'                       
  #' returns data.frame: the prepared expression data
  
  # Load the patients' expression data
  # Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62564
  # Needs to be downloaded as it is too large to store on github
  expression_data <- read.table(
    file.path(PATH, "GSE62564/GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt"),
    header=TRUE,
    sep="\t"
  )
  
  # Convert corrupt data to NA, then filter out rows with NA values
  expression_data <- expression_data[rowSums(is.na(expression_data)) == 0,]
  
  # Add external_gene_name to expression data
  expression_data <- merge(
    gene_list[c("refseq_mrna", "external_gene_name")],
    expression_data,
    by.x="refseq_mrna",
    by.y="RefSeqID"
  )
  
  # Average expression data from duplicated genes
  ## Identify which genes have multiple entries
  genes_duplicated <- unique(
    expression_data$external_gene_name[duplicated(expression_data$external_gene_name)]
  )
  
  for (gene in genes_duplicated) {
    
    # Get the mean expression values for each gene
    mean_values <- data.frame(t(colMeans(
      expression_data[expression_data$external_gene_name == gene, 3:length(expression_data)]
    )))
    
    # The information to identify the genes and patients
    row_info <- expression_data[expression_data$external_gene_name == gene, 1:2][1,]
    
    # Combine the expression data and identification info
    row_data <- data.frame(row_info, mean_values)
    
    # Remove the original expression data from the dataframe
    expression_data <- expression_data[expression_data$external_gene_name != gene,]
    
    # Add the averaged expression data to the dataframe
    expression_data <- rbind(expression_data, row_data)
  }
  
  # Get the gene names for future assignment
  genes <- expression_data$external_gene_name
  
  # Filter to only the expression features
  expression_data <- expression_data[,grepl("SEQC", colnames(expression_data))]
  
  # Undo log2 transformation and Z-transform the expression values
  expression_data <- data.frame(scale(2^expression_data))
  
  # Transpose the expression data
  expression_data_t <- data.frame(t(expression_data))
  
  # Label the columns with the gene names
  colnames(expression_data_t) <- genes
  
  # Add sequence_id
  expression_data_t$sequence_id <- colnames(expression_data)
  
  return(expression_data_t)
}


load_and_prepare_patient_data <- function() {
  #' Load and prepare the patient data for analysis
  #' 
  #' return data.frame: the patient data
  
  # Load GEO data
  gse <- getGEO("GSE62564")[[1]]
  
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
  #' Rename the features of the patient data
  #' 
  #' df data.frame: the patient data to be renamed
  #' 
  #' returns data.frame: the renamed patient data
  
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
  #' Format the patient data so that it is ready for analysis
  #' 
  #' df data.frame: the patient data
  #' 
  #' return data.frame: the formatted patient data
  
  # Trim unnecessary characters from the sequence_id
  df$sequence_id <- substr(df$sequence_id, 1, nchar(df$sequence_id) - 4)
  
  # Convert string values to integers
  df$high_risk[df$high_risk == "HR"] <- 1
  df$high_risk[is.na(df$high_risk)] <- 0
  
  # Identify the male patients
  df$male[df$male == "M"] <- 1
  df$male[df$male == "F"] <- 0
  
  # Convert null values to -1, so the feature can be an integer
  df$favourable[is.na(df$favourable)] <- -1
  
  # Convert features to integers
  df$age_days <- as.integer(df$age_days)
  df$male <- as.integer(df$male)
  df$high_risk <- as.integer(df$high_risk)
  df$favourable <- as.integer(df$favourable)
  df$event_free_survival <- as.integer(df$event_free_survival)
  df$event_free_survival_days <- as.integer(df$event_free_survival_days)
  
  return(df)
}


# Load differentially expressed genes
gene_list <- read.csv(file.path(PATH, "gene_list/gene_list.csv"))

# Create join table for merging the gene_list with the
# gene expression data
join_table <- create_join_table(gene_list$ensembl_gene_id)

# Add refseq_mrna values to gene_list
gene_list <- merge(gene_list, join_table, on="ensembl_gene_id_version")

# Get expression data
expression_data <- load_and_prepare_expression_data(gene_list)

# Get patient data
patients <- load_and_prepare_patient_data()

# Merge the patient data with their expression data
patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

# Save the patient data
write.csv(patients, file.path(PATH, "processed/GSE62564.csv"), row.names=FALSE)
