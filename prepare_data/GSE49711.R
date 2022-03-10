# Prepare GSE49711 for analysis

library(GEOquery)
library(fastDummies)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


load_and_prepare_patient_data <- function() {
  #' Load and prepare the patient data for analysis
  #' 
  #' return data.frame: the patient data
  
  # Load GEO data
  gse <- getGEO("GSE49711")[[1]]
  
  # Get the patient data
  patients <- pData(gse)
  
  # Clear the row names
  rownames(patients) <- c()
  
  # Filter to relevant columns
  patients <- patients[c(
    "title",
    "geo_accession",
    "Sex:ch1",
    "age at diagnosis:ch1",
    "high risk:ch1",
    "death from disease:ch1",
    "mycn status:ch1",
    "inss stage:ch1"
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
  names(df)[names(df) == "age at diagnosis:ch1"] <- "age_at_diagnosis_days"
  names(df)[names(df) == "high risk:ch1"] <- "high_risk"
  names(df)[names(df) == "death from disease:ch1"] <- "death_from_disease"
  names(df)[names(df) == "mycn status:ch1"] <- "mycn_amplification"
  names(df)[names(df) == "inss stage:ch1"] <- "inss_stage"
  
  return(df)
}


correct_data_types <- function(df) {
  #' Format the patient data so that it is ready for analysis
  #' 
  #' df data.frame: the patient data
  #' 
  #' return data.frame: the formatted patient data
  
  # Identify the male patients
  df$male[df$male == "M"] <- 1
  df$male[df$male == "F"] <- 0
  
  # Convert null values to -1, so the feature can be an integer
  df <- df[df$mycn_amplification != "N/A", ]

  # Convert features to integers
  df$male <- as.integer(df$male)
  df$age_at_diagnosis_days <- as.integer(df$age_at_diagnosis_days)
  df$high_risk <- as.integer(df$high_risk)
  df$death_from_disease <- as.integer(df$death_from_disease)
  df$mycn_amplification <- as.integer(df$mycn_amplification)

  # Create dummy columns from inss_stage
  df <- dummy_cols(
    df,
    select_columns="inss_stage",
    remove_first_dummy=TRUE,
    remove_selected_columns=TRUE
  )
  
  return(df)
}


load_and_prepare_expression_data <- function(gene_list) {
  #' Load and prepare the GEO expression data
  #' 
  #' gene_list List(str): the differentially expressed gene names to filter the
  #'                      the expression data to.
  #'                       
  #' returns data.frame: the prepared expression data
  
  # Load expression data
  # Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711
  # Needs to be downloaded as it is too large to store on github
  df <- read.csv(file.path(PATH, "GSE49711/GSE49711_SEQC_NB_TUC_G_log2.txt"), sep='\t')
  
  # Filter to genes in the gene list
  df <- df[df$X00gene_id %in% gene_list,]
  
  # Convert corrupt data to NA, then filter out rows will NA values
  df <- df %>% dplyr::na_if(0)
  df <- df[rowSums(is.na(df)) == 0,]
  
  # Select the genes for future assignment
  genes <- df$X00gene_id
  
  # Filter to only the expression features
  df <- df[,grepl("SEQC", colnames(df))]
  
  # Transpose the expression data
  df_t <- t(df)
  
  # Label the columns with the gene names
  colnames(df_t) <- genes
  
  # Undo log2 transformation and Z-transform the expression values
  df_t <- data.frame(scale(2^df_t))
  
  # Add sequence_id
  df_t$sequence_id <- colnames(df)
  df_t$sequence_id <- substr(df_t$sequence_id, 1, 10)
  
  return(df_t)
}


# Load differentially expressed genes
gene_list <- read.csv(file.path(PATH, "gene_list/gene_list.csv"))

# Get the patient data
patients <- load_and_prepare_patient_data()

# Get the expression data
expression_data <- load_and_prepare_expression_data(gene_list$external_gene_name)

# Merge patient data with their expression data
patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

# Save the patient data
write.csv(patients, file.path(PATH, "processed/GSE49711.csv"), row.names=FALSE)
