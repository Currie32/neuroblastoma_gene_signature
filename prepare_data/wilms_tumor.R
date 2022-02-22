PATH <- "~/Downloads/wt_target_2018_pub/"


load_and_prepare_patient_data <- function() {

  headers <- read.csv(
    file.path(PATH, "data_clinical_patient.txt"),
    nrows = 1,
    header=F,
    sep="\t"
  )
  # Load data, but skip first five rows since they do not contain the data
  patients <- read.csv(
    file.path(PATH, "data_clinical_patient.txt"),
    skip=5,
    header=F,
    sep="\t"
  )
  
  colnames(patients) <- headers
  
  patients <- patients[c(
    "#Patient Identifier",
    "Sex",
    "Diagnosis Age (days)",
    "Neoplasm American Joint Committee on Cancer Clinical Group Stage",
    "Event Type",
    "Time To Event (days)"
  )]
  
  patients <- rename_columns(patients)
  patients <- correct_data_types(patients)
  
  return(patients)
}


rename_columns <- function(df) {
  
  names(df)[names(df) == "#Patient Identifier"] <- "sequence_id"
  names(df)[names(df) == "Sex"] <- "male"
  names(df)[names(df) == "Diagnosis Age (days)"] <- "age_at_diagnosis_days"
  names(df)[names(df) == "Neoplasm American Joint Committee on Cancer Clinical Group Stage"] <- "inss_stage"
  names(df)[names(df) == "Event Type"] <- "event_free_survival"
  names(df)[names(df) == "Time To Event (days)"] <- "event_free_survival_days"
  
  return(df)
}


correct_data_types <- function(df) {
  
  # Identify the male patients
  df$male[df$male == "Male"] <- 1
  df$male[df$male == "Female"] <- 0
  df$male <- as.integer(df$male)
  
  # Update values of inss_stage
  df$inss_stage[df$inss_stage == "I"] <- "1"
  df$inss_stage[df$inss_stage == "II"] <- "2"
  df$inss_stage[df$inss_stage == "III"] <- "3"
  df$inss_stage[df$inss_stage == "IV"] <- "4"

  # Drop rows with an unknown inss_stage value
  df <- df[!df$inss_stage %in% c("", "II/V", "III/V", "IIIB", "IIIB/V", "IV/V", "U", "V"),]
  
  # Create dummy columns from inss_stage
  df <- dummy_cols(
    df,
    select_columns="inss_stage",
    remove_first_dummy=TRUE,
    remove_selected_columns=TRUE
  )
  
  df$event_free_survival[df$event_free_survival == "None"] <- 1
  df$event_free_survival[df$event_free_survival %in% c("Progression", "Relapse")] <- 0
  
  return(df)
}

df <- read.csv(
  file.path(PATH, "data_mrna_seq_rpkm_zscores_ref_all_samples.txt"),
  sep="\t"
)
df <- df[df$Hugo_Symbol %in% gene_list$external_gene_name,]
df <- df[!duplicated(df$Hugo_Symbol),]
genes <- df$Hugo_Symbol
df <- df[,grepl("TARGET", colnames(df))]
df_t <- t(df)
colnames(df_t) <- genes
df_t <- data.frame(df_t)
df_t$sequence_id <- colnames(df)
df_t$sequence_id <- gsub("\\.", "-", df_t$sequence_id)
df_t$sequence_id <- substr(df_t$sequence_id, 1, 16)

df


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
  df_t <- t(df)
  
  # Label the columns with the gene names
  colnames(df_t) <- genes
  
  # Undo log2 transformation and Z-transform the expression values
  df_t <- data.frame(df_t)
  
  # Add sequence_id
  df_t$sequence_id <- colnames(df)
  df_t$sequence_id <- gsub("\\.", "-", df_t$sequence_id)
  df_t$sequence_id <- substr(df_t$sequence_id, 1, 16)
  
  return(df_t)
}

# Load differentially expressed genes
gene_list <- read.csv("~/Imperial/neuroblastoma_gene_signature/data/gene_list.csv")

patients <- load_and_prepare_patient_data()

expression_data <- load_and_prepare_expression_data(gene_list$external_gene_name)

# Merge patient data with their expression data
patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

write.csv(patients, "~/Imperial/neuroblastoma_gene_signature/data/wilms_tumor.csv", row.names=FALSE)
