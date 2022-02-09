library(data.table)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


load_and_prepare_patient_data <- function() {
  
  # Load header row
  headers <- read.csv(
    file.path(PATH, "data_clinical_patient.txt"),
    nrows = 1,
    sep="\t"
  )
  # Load data, but skip first five rows since they do not contain the data
  patients <- read.csv(
    file.path(PATH, "data_clinical_patient.txt"),
    skip=5,
    header=F,
    sep="\t"
  )
  # Rename the columns using the headers
  colnames(patients) <- headers
  
  # Filter to relevant columns
  patients <- patients[c(
    "#Identifier to uniquely specify a patient.",
    "Sex",
    "Age at which a condition or disease was first diagnosed.",
    "Staging according to the International Neuroblastoma Staging System",
    "Risk group.",
    "Tumor Sample Histology",
    "Event Free Survival Censored",
    "EFS time."
  )]
  
  patients <- rename_columns(patients)
  patients <- correct_data_types(patients)
  
  return(patients)
}


rename_columns <- function(df) {
  
  names(df)[names(df) == "#Identifier to uniquely specify a patient."] <- "sequence_id"
  names(df)[names(df) == "Sex"] <- "male"
  names(df)[names(df) == "Age at which a condition or disease was first diagnosed."] <- "age_at_diagnosis_days"
  names(df)[names(df) == "Staging according to the International Neuroblastoma Staging System"] <- "inss_stage"
  names(df)[names(df) == "Risk group."] <- "high_risk"
  names(df)[names(df) == "Tumor Sample Histology"] <- "favourable"
  names(df)[names(df) == "Event Free Survival Censored"] <- "event_free_survival"
  names(df)[names(df) == "EFS time."] <- "event_free_survival_days"
  
  return(df)
}


correct_data_types <- function(df) {
  
  # Identify the male patients
  df$male[df$male == "Male"] <- 1
  df$male[df$male == "Female"] <- 0
  df$male <- as.integer(df$male)
  
  # Update values of inss_stage
  df$inss_stage[df$inss_stage == "Stage 1"] <- "1"
  df$inss_stage[df$inss_stage == "Stage 2a"] <- "2"
  df$inss_stage[df$inss_stage == "Stage 2b"] <- "2"
  df$inss_stage[df$inss_stage == "Stage 3"] <- "3"
  df$inss_stage[df$inss_stage == "Stage 4"] <- "4"
  df$inss_stage[df$inss_stage == "Stage 4s"] <- "4S"
  
  # Drop rows with an unknown inss_stage value
  df <- df[df$inss_stage != "Unknown",]
  
  # Create dummy columns from inss_stage
  df <- dummy_cols(
    df,
    select_columns="inss_stage",
    remove_first_dummy=TRUE,
    remove_selected_columns=TRUE
  )
  
  # Set high risk values to 1, all else is 0
  df$high_risk <- ifelse(df$high_risk == "High Risk", 1, 0)
  
  # Set favourable values to 1, all else is 0
  df$favourable <- ifelse(df$favourable == "Favourable", 1, 0)
  
  # Drop rows with a null event_free_survival value
  df <- df[!is.na(df$event_free_survival),]
  
  return(df)
}


load_and_prepare_expression_data <- function(gene_list) {
  
  # Load expression data
  # Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711
  df <- read.csv(
    file.path(PATH, "data_mrna_seq_rpkm_zscores_ref_all_samples.txt"),
    sep="\t"
  )
  
  # Filter to genes in the gene list
  df <- df[df$Hugo_Symbol %in% gene_list,]
  
  # Remove duplicate genes
  df <- df[!duplicated(df$Hugo_Symbol),]
  
  # Select the genes for future assignment
  genes <- df$Hugo_Symbol
  
  # Filter to only the expression features
  df <- df[,grepl("TARGET", colnames(df))]
  
  # Transpose the expression data
  df_t <- transpose(df)
  
  # Label the columns with the gene names
  colnames(df_t) <- genes
  
  # Undo log2 transformation and Z-transform the expression values
  df_t <- data.frame(scale(2^df_t))
  
  # Add sequence_id
  df_t$sequence_id <- colnames(df)
  df_t$sequence_id <- gsub("\\.", "-", df_t$sequence_id)
  df_t$sequence_id <- substr(df_t$sequence_id, 1, 16)
  
  return(df_t)
}

patients <- load_and_prepare_patient_data()

# Load differentially expressed genes
gene_list <- read.csv(file.path(PATH, "gene_list.csv"))

expression_data <- load_and_prepare_expression_data(gene_list$external_gene_name)

# Merge patient data with their expression data
patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

# Save patient data
write.csv(patients, file.path(PATH, "target_2018.csv"), row.names=FALSE)


