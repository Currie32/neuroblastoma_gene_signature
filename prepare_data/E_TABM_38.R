PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


load_and_prepare_patient_data <- function() {
  
  # Data is downloaded from: https://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-38/files/
  patients <- read.csv("~/Downloads/E-TABM-38.sdrf.txt", sep="\t")
  
  # Filter to relevant columns
  patients <- patients[c(
    "Hybridization.Name",
    "Characteristics..Sex.",
    "Characteristics..Age.",
    "Characteristics..DiseaseStaging.",
    "Factor.Value..ClinicalInformation.",
    "Factor.Value..EventFreeSurvival."
  )]
  
  patients <- rename_columns(patients)
  patients <- correct_data_types(patients)
  
  return(patients)
}


rename_columns <- function(df) {
  
  names(df)[names(df) == "Hybridization.Name"] <- "sequence_id"
  names(df)[names(df) == "Characteristics..Sex."] <- "male"
  names(df)[names(df) == "Characteristics..Age."] <- "age_at_diagnosis_days"
  names(df)[names(df) == "Characteristics..DiseaseStaging."] <- "inss_stage"
  names(df)[names(df) == "Factor.Value..ClinicalInformation."] <- "event_free_survival"
  names(df)[names(df) == "Factor.Value..EventFreeSurvival."] <- "event_free_survival_days"
  
  return(df)
}


correct_data_types <- function(df) {
  
  df <- df[!is.na(df$age_at_diagnosis_days),]
  
  # Identify the male patients
  df$male[df$male == "male"] <- 1
  df$male[df$male == "female"] <- 0
  
  # Update values of inss_stage
  df$inss_stage[df$inss_stage == "Stage 1"] <- "1"
  df$inss_stage[df$inss_stage == "Stage 2"] <- "2"
  df$inss_stage[df$inss_stage == "Stage 2a"] <- "2"
  df$inss_stage[df$inss_stage == "Stage 2b"] <- "2"
  df$inss_stage[df$inss_stage == "Stage 3"] <- "3"
  df$inss_stage[df$inss_stage == "Stage 4"] <- "4"
  df$inss_stage[df$inss_stage == "Stage 4s"] <- "4S"
  
  # Create dummy columns from inss_stage
  df <- dummy_cols(
    df,
    select_columns="inss_stage",
    remove_first_dummy=TRUE,
    remove_selected_columns=TRUE
  )
  
  # Update values of event_free_survival
  df$event_free_survival[df$event_free_survival == "death from disease"] <- 1
  df$event_free_survival[df$event_free_survival == "relapse of disease"] <- 1
  df$event_free_survival[df$event_free_survival == "no event"] <- 0
  
  return(df)
}


load_and_prepare_expression_data <- function() {
  
  expression_data <- read.csv(file.path(PATH, "Warnat19092005_E-NBDE-1_FGEDM.txt"), sep="\t")
  
  # Remove unncessary row
  expression_data <- expression_data[-c(1),]
  
  # Convert expression values to numerics and z-transform them
  expression_data <- data.frame(sapply(expression_data, as.numeric))
  expression_data <- data.frame(sapply(expression_data, scale))
  
  # Add gene names to the expression data
  gene_names <- get_gene_names()
  expression_data$gene <- gene_names
  
  # Filter to the relevant genes
  gene_list <- read.csv(file.path(PATH, "gene_list.csv"))
  expression_data <- expression_data[expression_data$gene %in% gene_list$external_gene_name,]
  
  # Transpose the dataframe
  expression_data_t <- data.frame(t(expression_data))
  
  # Name the columns after the gene names
  colnames(expression_data_t) <- expression_data$gene
  
  # Add and format sequence_id
  expression_data_t$sequence_id <- colnames(expression_data)
  expression_data_t$sequence_id <- str_replace(substr(expression_data_t$sequence_id, 1, 10), "\\.", " ")
  
  return(expression_data_t)
}


get_gene_names <- function() {
  
  # Data is downloaded from: https://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-38/files/
  gene_names <- read.csv(file.path(PATH, "A-MEXP-255_comments.txt"), sep="\t")
  controls_and_experiments <- read.csv(file.path(PATH, "A-MEXP-255.adf.txt"), sep="\t", skip=23)
  
  # Identify the controls rows for which we do not have expression data
  n_control_rows <- sum(controls_and_experiments$Reporter.Group.role. == "CONTROL") + 1
  
  # Filter out the control rows from the list of gene names
  gene_names <- gene_names[n_control_rows:length(gene_names$Comment.AEReporterName.),]
  
  return(gene_names)
}


# Load the prepare the data
patients <- load_and_prepare_patient_data()
expression_data <- load_and_prepare_expression_data()

patients <- merge(
  patients,
  expression_data,
  on="sequence_id"
)

write.csv(patients, file.path(PATH, "E_TABM_38.csv"), row.names=FALSE)
