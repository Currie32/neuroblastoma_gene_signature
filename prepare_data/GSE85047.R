library(GEOquery)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


prepare_patient_data <- function(gse) {
  
  # Extract patient data
  patients <- gse[["GSE85047_series_matrix.txt.gz"]]@phenoData@data
  
  # Filter to relevant columns
  patients <- patients[c(
    "title",
    "geo_accession",
    "age_at_diagnosis(days):ch1",
    "inss:ch1",
    "mycn_amplification:ch1",
    "event_progression_free:ch1",
    "profression_free_survival_time(days):ch1"
  )]
  
  patients <- rename_columns(patients)
  patients <- correct_data_types(patients)
  
  return(patients)
}


rename_columns <- function(df) {
  
  names(df)[names(df) == "title"] <- "sequence_id"
  names(df)[names(df) == "geo_accession"] <- "geo_id"
  names(df)[names(df) == "age_at_diagnosis(days):ch1"] <- "age_at_diagnosis_days"
  names(df)[names(df) == "inss:ch1"] <- "inss_stage"
  names(df)[names(df) == "mycn_amplification:ch1"] <- "mycn_amplification"
  names(df)[names(df) == "event_progression_free:ch1"] <- "event_free_survival"
  names(df)[names(df) == "profression_free_survival_time(days):ch1"] <- "event_free_survival_days"
  
  return(df)
}

correct_data_types <- function(df) {
  
  # Remove rows with null values
  df <- df[df$age_at_diagnosis_days != "NA", ]
  df <- df[df$event_free_survival_days != "NA", ]
  df <- df[df$mycn_amplification != "NA", ]
  df <- df[df$event_free_survival != "NA", ]
  df <- df[df$inss_stage != "NA", ]
  
  # Convert columns to integers
  df$age_at_diagnosis_days <- as.integer(df$age_at_diagnosis_days)
  df$event_free_survival_days <- as.integer(df$event_free_survival_days)
  
  # Update values to 1 or 0
  df$mycn_amplification <- ifelse(df$mycn_amplification == "yes", 1, 0)
  df$event_free_survival <- ifelse(df$event_free_survival == "yes", 1, 0)
  
  # Update values of inss_stage
  df$inss_stage[df$inss_stage == "st1"] <- "1"
  df$inss_stage[df$inss_stage == "st2"] <- "2"
  df$inss_stage[df$inss_stage == "st3"] <- "3"
  df$inss_stage[df$inss_stage == "st4"] <- "4"
  df$inss_stage[df$inss_stage == "st4s"] <- "4S"
  
  # Create dummy columns from inss_stage
  df <- dummy_cols(
    df,
    select_columns="inss_stage",
    remove_first_dummy=TRUE,
    remove_selected_columns=TRUE
  )
  
  return(df)
}


prepare_expression_data <- function(gse, gene_list) {
  
  # Get gene names for expression data
  genes <- gene_names(gse)
  
  # Extract the expression data
  expression_data <- gse[["GSE85047_series_matrix.txt.gz"]]@assayData[["exprs"]]
  expression_data <- data.frame(expression_data)
  
  # Remove rows that are missing a gene name
  expression_data <- expression_data[genes != "",]
  genes <- genes[genes != ""]
  
  # Remove rows with duplicate genes
  expression_data <- expression_data[!duplicated(genes),]
  genes <- genes[!duplicated(genes)]
  
  # Transpose the data to have the genes as columns
  expression_data_t <- data.frame(t(expression_data))
  
  # Label columns with gene names
  colnames(expression_data_t) <- genes
  
  # Filter to genes in gene_list
  expression_data_t <- expression_data_t[colnames(expression_data_t) %in% gene_list]
  
  # Add geo_id for merging with patients data
  expression_data_t$geo_id <- colnames(expression_data)
  
  return(expression_data_t)
}


gene_names <- function(df) {
  
  feature_data <- gse[["GSE85047_series_matrix.txt.gz"]]@featureData@data
  genes <- str_split_fixed(feature_data$gene_assignment, " // ", n=3)[,2]
  
  return(genes)
}


# Load differentially expressed genes
gene_list <- read.csv(file.path(PATH, "gene_list.csv"))

# Load GEO data
gse <- getGEO("GSE85047")

# Prepare datasets for merging
patients <- prepare_patient_data(gse)
expression_data <- prepare_expression_data(gse, gene_list$external_gene_name)

# Merge patient data with their expression data
patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

# Save patient data
write.csv(patients, file.path(PATH, "GSE85047.csv"), row.names=FALSE)
