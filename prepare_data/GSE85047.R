library(data.table)
library(GEOquery)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


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
    "age_at_diagnosis(days):ch1",
    "event_progression_free:ch1",
    "profression_free_survival_time(days):ch1",
    "mycn_amplification:ch1",
    "inss:ch1"
  )]
  
  patients <- rename_columns(patients)
  # patients <- correct_data_types(patients)
  
  return(patients)
}


rename_columns <- function(df) {
  
  names(df)[names(df) == "title"] <- "sequence_id"
  names(df)[names(df) == "geo_accession"] <- "geo_id"
  names(df)[names(df) == "age_at_diagnosis(days):ch1"] <- "age_at_diagnosis_days"
  names(df)[names(df) == "event_progression_free:ch1"] <- "event_free_survival"
  names(df)[names(df) == "profression_free_survival_time(days):ch1"] <- "event_free_survival_days"
  names(df)[names(df) == "mycn_amplification:ch1"] <- "mycn_status"
  names(df)[names(df) == "inss:ch1"] <- "inss_stage"
  
  return(df)
}


correct_data_types <- function(df) {
  
  # Convert null values to -1, so the feature can be an integer
  df$mycn_status[df$mycn_status == "N/A"] <- -1
  
  df$male <- as.integer(df$male)
  df$age_at_diagnosis_days <- as.integer(df$age_at_diagnosis_days)
  df$high_risk <- as.integer(df$high_risk)
  df$death_from_disease <- as.integer(df$death_from_disease)
  df$mycn_status <- as.integer(df$mycn_status)
  
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
  
  # Load expression data
  # Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711
  df <- read.csv(file.path(PATH, "GSE49711_SEQC_NB_TUC_G_log2.txt"), sep='\t')
  
  # Filter to genes in the gene list
  df <- df[df$X00gene_id %in% gene_list,]
  
  # Remove duplicate genes
  df <- df[!duplicated(df$X00gene_id),]
  
  # Select the genes for future assignment
  genes <- df$X00gene_id
  
  # Filter to only the expression features
  df <- df[,grepl("SEQC", colnames(df))]
  
  # Transpose the expression data
  df_t <- transpose(df)
  
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
gene_list <- read.csv(file.path(PATH, "gene_list.csv"))

# Load the patients' data
geo_id <- "GSE85047"
patients <- load_and_prepare_patient_data(geo_id)

head(patients)

colSums(is.na(patients))


expression_data <- load_and_prepare_expression_data(gene_list$external_gene_name)

# Merge patient data with their expression data
patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

# Save patient data
write.csv(patients, file.path(PATH, "GSE49711.csv"), row.names=FALSE)
