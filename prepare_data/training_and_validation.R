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


add_correlated_genes <- function(genes, df) {
  
  correlations <- cor(df[colnames(df) %in% genes & colSums(df == 0) != length(df)])
  correlations <- abs(data.frame(correlations))
  
  corr_threshold <- 0.85
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
    
    add_correlated_genes(genes, df)
  }
  else {
    return(df)
  }
}


average_correlated_genes_train <- function(genes, df, train) {
  
  for (genes_correlated in genes) {
    genes_correlated_split <- strsplit(genes_correlated, '_')[[1]]
    df[genes_correlated] <- scale(rowMeans(train[genes_correlated_split]))
  }
  
  return(df)
}


average_correlated_genes <- function(genes, df) {

  for (genes_correlated in genes) {
    genes_correlated_split <- strsplit(genes_correlated, '_')[[1]]
    
    if (all(df[genes_correlated_split] == 0)) {
      df[genes_correlated] <- 0
    }
    else {
      df[genes_correlated] <- scale(rowMeans(df[genes_correlated_split]))
    }
  }
  
  return(df)
}


GSE49711 <- read.csv(file.path(PATH, "GSE49711.csv"))
GSE62564 <- read.csv(file.path(PATH, "GSE62564.csv"))
GSE85047 <- read.csv(file.path(PATH, "GSE85047.csv"))
target_2018 <- read.csv(file.path(PATH, "target_2018.csv"))
wilms_tumor <- read.csv(file.path(PATH, "wilms_tumor.csv"))

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
train_corr <- add_correlated_genes(genes, train)

correlated_genes <- colnames(train_corr)[!colnames(train_corr) %in% colnames(train)]

train_corr <- average_correlated_genes_train(correlated_genes, train_corr, train)


# Add missing features to testing datasets
GSE85047_corr <- add_missing_features(train, GSE85047)
target_2018_corr <- add_missing_features(train, target_2018)
wilms_tumor_corr <- add_missing_features(train, wilms_tumor)

GSE85047_corr <- average_correlated_genes(correlated_genes, GSE85047_corr)
target_2018_corr <- average_correlated_genes(correlated_genes, target_2018_corr)
wilms_tumor_corr <- average_correlated_genes(correlated_genes, wilms_tumor_corr)

# Save datasets
write.csv(train_corr, file.path(PATH, "train.csv"), row.names=FALSE)
write.csv(GSE85047_corr, file.path(PATH, "test_GSE85047.csv"), row.names=FALSE)
write.csv(target_2018_corr, file.path(PATH, "test_target.csv"), row.names=FALSE)
write.csv(wilms_tumor_corr, file.path(PATH, "test_wilms_tumor.csv"), row.names=FALSE)
