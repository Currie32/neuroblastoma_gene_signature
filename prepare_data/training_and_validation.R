library(dplyr)
library(tibble)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


extract_patient_traits <- function(GSE49711, GSE62564) {
  
  patient_traits <- GSE49711[c(
    "sequence_id",
    "geo_id",
    "male",
    "age_at_diagnosis_days",
    "high_risk",
    "mycn_amplification",
    "inss_stage_2",
    "inss_stage_3",
    "inss_stage_4",
    "inss_stage_4S"
  )]
  
  patient_traits <- merge(
    patient_traits,
    GSE62564[c(
      "sequence_id",
      "event_free_survival",
      "event_free_survival_days"
    )],
    on="sequence_id"
  )
  
  return(patient_traits)
}


format_GSE_patient_data <- function(gene_list, patients_GSE, patients_traits) {
  
  patients_GSE <- patients_GSE[
    colnames(patients_GSE) %in% append("sequence_id", gene_list)
  ]
  
  patients_GSE <- merge(
    patients_traits,
    patients_GSE,
    on="sequence_id"
  )
  
  return(patients_GSE)
}


combine_datasets <- function(GSE1, GSE2) {
  
  patients <- bind_rows(GSE1, GSE2)
  patients[is.na(patients)] <- 0
  
  return(patients)
}


add_missing_features <- function(train, df) {
  
  missing_features <- setdiff(names(train), names(df))
  df[missing_features] <- 0
  
  return (df)
}


correlated_genes <- function(genes, df) {
  
  correlations <- cor(df[colnames(df) %in% genes & colSums(df == 0) != length(df)])
  correlations <- abs(data.frame(correlations))
  
  corr_threshold <- 0.7
  max_corr <- corr_threshold
  max_feature <- ""
  max_feature_paired <- ""
  
  for (feature in colnames(correlations)) {
    corr <- max(correlations[feature][correlations[feature] < 1])
    
    if (!is.na(corr)) {
      if (corr > max_corr) {
        max_corr <- corr
        max_feature <- feature
      }
      else if (corr == max_corr) {
        max_feature_paired <- feature
      }
    }
  }
  
  if (max_corr > corr_threshold) {
    
    joined_feature <- paste(c(max_feature, max_feature_paired), collapse='_')
    
    genes <- append(genes, joined_feature)
    genes <- genes[!genes %in% c(max_feature, max_feature_paired)]
    
    df[joined_feature] <- rowMeans(df[c(max_feature, max_feature_paired)])
    df <- df[!colnames(df) %in% c(max_feature, max_feature_paired)]
    
    correlated_genes(genes, df)
  }
  else {
    return(df)
  }
}


average_correlated_genes <- function(genes, df, df_corr) {
  
  genes_to_average <- colnames(df_corr[!colnames(df_corr) %in% colnames(df)])
  
  for (genes_correlated in genes_to_average) {
    genes_correlated_split <- strsplit(genes_correlated, '_')[[1]]
    df[genes_correlated] <- scale(rowMeans(df[genes_correlated_split]))
  }
  
  return(df)
}


GSE49711 <- read.csv(file.path(PATH, "GSE49711.csv"))
GSE62564 <- read.csv(file.path(PATH, "GSE62564.csv"))
GSE85047 <- read.csv(file.path(PATH, "GSE85047.csv"))
target_2018 <- read.csv(file.path(PATH, "target_2018.csv"))

# Load differentially expressed genes
gene_list <- read.csv(file.path(PATH, "gene_list.csv"))
genes <- gene_list$external_gene_name

patients_traits <- extract_patient_traits(GSE49711, GSE62564)

GSE49711 <- format_GSE_patient_data(
  genes, 
  GSE49711,
  patients_traits
)

GSE62564 <- format_GSE_patient_data(
  genes, 
  GSE62564,
  patients_traits
)

# Create training dataset
train <- combine_datasets(GSE49711, GSE62564)
train_corr <- correlated_genes(genes, train)
train_corr <- average_correlated_genes(genes, train, train_corr)

# Add missing features to testing datasets
GSE85047 <- add_missing_features(train_corr, GSE85047)
target_2018 <- add_missing_features(train_corr, target_2018)

GSE85047_corr <- average_correlated_genes(genes, GSE85047, GSE85047_corr)
target_2018_corr <- average_correlated_genes(genes, target_2018, target_2018_corr)

# Save datasets
write.csv(train_corr, file.path(PATH, "train.csv"), row.names=FALSE)
write.csv(GSE85047_corr, file.path(PATH, "test_GSE85047.csv"), row.names=FALSE)
write.csv(target_2018_corr, file.path(PATH, "test_target.csv"), row.names=FALSE)
