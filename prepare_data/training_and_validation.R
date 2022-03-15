# Create the training and validation datasets

library(dplyr)
library(tibble)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


extract_patient_traits <- function(GSE49711, GSE62564) {
  #' Get the patient traits from the two datasets that will be used
  #' as the training data
  #' 
  #' GSE49711 data.frame: one of the training datasets
  #' GSE62564 data.frame: the other training dataset
  #' 
  #' return data.frame: phenotype data of the patients used to train the models
  
  # The patient traits from GSE49711
  patient_traits_1 <- GSE49711[c(
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
  
  # The patient traits from GSE62564
  patient_traits_2 <- GSE62564[c(
    "sequence_id",
    "event_free_survival",
    "event_free_survival_days"
  )]
  
  # Merge both sets of patient traits
  patient_traits <- merge(
    patient_traits_1,
    patient_traits_2,
    on="sequence_id"
  )
  
  return(patient_traits)
}


merge_expression_and_patient_data <- function(genes, gse_data, patients_traits) {
  #' Get the expression data, then merge it with the patient data.
  #'
  #' genes List(str): names of the genes to extract from gse_data
  #' gse_data data.frame: contains the expression data for the patients
  #' patients_traits data.frame: contains the patients' phenotype data
  #' 
  #' return data.frame: The combined expression and phenotype data for a
  #'                    set of patients
  
  # Get the expression data
  expression_data <- gse_data[
    colnames(gse_data) %in% append("sequence_id", genes)
  ]
  
  # Merge the phenotype and expression data
  patients_data <- merge(
    patients_traits,
    expression_data,
    on="sequence_id"
  )
  
  return(patients_data)
}


training_data <- function(df1, df2) {
  #' Create the training dataset by merging two datasets
  #' 
  #' df1 data.frame: one of the training datasets
  #' df2 data.frame: the other training dataset
  #' 
  #' return data.frame: the training data
  
  # Concatenate the dataframes
  patients <- bind_rows(df1, df2)
  
  # Set null values to 0
  # This happens since not all genes are present in both datasets
  patients[is.na(patients)] <- 0
  
  return(patients)
}


find_and_average_correlated_genes <- function(genes, df) {
  #' Average the expression values of correlated genes
  #' 
  #' genes List(str): names of differential expressed genes
  #' df data.frame: contains the expression values
  #' 
  #' return data.frame: the patient and expression data, but with averaged
  #'                    expression values for correlated genes
  
  # Measure the correlation between differentially expressed genes
  # that do not contain only 0 values
  correlations <- cor(df[colnames(df) %in% genes & colSums(df == 0) != length(df)])
  
  # Get the absolute values of the correlations
  correlations <- abs(data.frame(correlations))
  
  # Preset values for finding correlated genes
  corr_threshold <- 0.9
  max_corr <- corr_threshold
  max_feature <- ""
  max_feature_paired <- ""
  
  # Loop through each gene
  for (gene in colnames(correlations)) {
    
    # Find the most correlated gene
    corr <- max(correlations[gene][correlations[gene] < 1])
    
    if (!is.na(corr)) {
      if (corr > max_corr) {
        max_corr <- corr
        max_gene <- gene
      }
      else if (corr == max_corr) {
        max_gene_paired <- gene
      }
    }
  }
  
  # If the max correlation is greater than the correlation threshold,
  # average expression values for the two genes and
  # remove the original values from the data
  if (max_corr > corr_threshold) {
    
    joined_genes <- paste(c(max_gene, max_gene_paired), collapse='_')
    
    genes <- append(genes, joined_genes)
    genes <- genes[!genes %in% c(max_gene, max_gene_paired)]
    
    df[joined_genes] <- rowMeans(df[c(max_gene, max_gene_paired)])
    df <- df[!colnames(df) %in% c(max_gene, max_gene_paired)]
    
    # Repeat the workflow to try to find another pair of highly correlated genes
    find_and_average_correlated_genes(genes, df)
  }
  else {
    return(df)
  }
}


average_correlated_genes_train <- function(genes, df_correlated, train) {
  #' Average the expression values of the correlated genes in the training dataset
  #' 
  #' genes List(str): names of the correlated genes
  #' df_correlated data.frame: contains the correlated gene names
  #' train data.frame: contains the original expression values to be averaged
  #' 
  #' return data.frame: contains the correlated genes with averaged expression values
  
  for (genes_correlated in genes) {
    
    # Get the names of the correlated genes to be averaged
    genes_correlated_split <- strsplit(genes_correlated, '_')[[1]]
    
    # Average the expression values of the correlated genes
    df_correlated[genes_correlated] <- scale(rowMeans(train[genes_correlated_split]))
  }
  
  return(df_correlated)
}


add_missing_features <- function(train, df) {
  #' Add missing features to a testing dataframe
  #' 
  #' train data.frame: contains all of the required features
  #' df data.frame: a testing dataframe that might be missing features
  #' 
  #' return data.frame: a testing dataframe that contains all of the required features
  
  # Find the missing features and set their values to 0
  missing_features <- setdiff(names(train), names(df))
  df[missing_features] <- 0
  
  return (df)
}


average_correlated_genes <- function(genes, df) {
  #' Average the expression values of the correlated genes in a testing dataset
  #' 
  #' genes List(str): names of the correlated genes
  #' df data.frame: a testing dataset containing expression data
  #' 
  #' return data.frame: a testing dataset with averaged expression values for
  #'                    correlated genes
  
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


is_under_18_months <- function(df) {
  #' Identify which patients are under 18 months old
  #' 
  #' df data.frame: a dataframe missing this feature
  #' 
  #' return List(int): identifies which patients are under 18 months old
  
  under_18_months <- df$age_at_diagnosis_days < 548 # 548 = 365 * 1.5
  under_18_months[under_18_months == TRUE] <- 1
  under_18_months[under_18_months == FALSE] <- 0
  
  return(under_18_months)
}


# Load the datasets used in the analysis
GSE49711 <- read.csv(file.path(PATH, "processed/GSE49711.csv"))
GSE62564 <- read.csv(file.path(PATH, "processed/GSE62564.csv"))
GSE85047 <- read.csv(file.path(PATH, "processed/GSE85047.csv"))
target_2018 <- read.csv(file.path(PATH, "processed/target_2018.csv"))
E_TABM_38 <- read.csv(file.path(PATH, "processed/E_TABM_38.csv"))
wilms_tumor <- read.csv(file.path(PATH, "processed/wilms_tumor.csv"))

paste(colnames(GSE49711)[
  (colnames(GSE49711) %in% colnames(GSE62564))
  & (colnames(GSE49711) %in% colnames(GSE85047))
  & (colnames(GSE49711) %in% colnames(target_2018))
  & (colnames(GSE49711) %in% colnames(E_TABM_38))
], collapse=', ')

# Load the names of the differentially expressed genes
genes <- read.csv(file.path(PATH, "gene_list/gene_list.csv"))$external_gene_name

# Get the patient traits from the training datasets
patients_traits <- extract_patient_traits(GSE49711, GSE62564)

# Combine the expression and patient data
GSE49711 <- merge_expression_and_patient_data(genes, GSE49711, patients_traits)
GSE62564 <- merge_expression_and_patient_data(genes, GSE62564, patients_traits)

# Create training dataset
train <- training_data(GSE49711, GSE62564)

length(colnames(train)[colnames(train) %in% genes])

# Average the expression values of correlated genes
train_corr <- find_and_average_correlated_genes(genes, train)

# Identify the names of the correlated genes
correlated_genes <- colnames(train_corr)[!colnames(train_corr) %in% colnames(train)]

# Average the expression values of the correlated genes
train_corr <- average_correlated_genes_train(correlated_genes, train_corr, train)

# Add missing features to testing datasets
GSE85047_corr <- add_missing_features(train, GSE85047)
target_2018_corr <- add_missing_features(train, target_2018)
E_TABM_38_corr <- add_missing_features(train, E_TABM_38)
wilms_tumor_corr <- add_missing_features(train, wilms_tumor)

# Average the expression values of the correlated genes
GSE85047_corr <- average_correlated_genes(correlated_genes, GSE85047_corr)
target_2018_corr <- average_correlated_genes(correlated_genes, target_2018_corr)
E_TABM_38_corr <- average_correlated_genes(correlated_genes, E_TABM_38_corr)
wilms_tumor_corr <- average_correlated_genes(correlated_genes, wilms_tumor_corr)

# Add under_18_months feature
train_corr$under_18_months <- is_under_18_months(train_corr)
GSE85047_corr$under_18_months <- is_under_18_months(GSE85047_corr)
target_2018_corr$under_18_months <- is_under_18_months(target_2018_corr)
E_TABM_38_corr$under_18_months <- is_under_18_months(E_TABM_38_corr)
wilms_tumor_corr$under_18_months <- is_under_18_months(wilms_tumor_corr)

# Save datasets
write.csv(train_corr, file.path(PATH, "modelling/train.csv"), row.names=FALSE)
write.csv(GSE85047_corr, file.path(PATH, "modelling/test_GSE85047.csv"), row.names=FALSE)
write.csv(target_2018_corr, file.path(PATH, "modelling/test_target.csv"), row.names=FALSE)
write.csv(E_TABM_38_corr, file.path(PATH, "modelling/test_E_TABM_38.csv"), row.names=FALSE)
write.csv(wilms_tumor_corr, file.path(PATH, "modelling/test_wilms_tumor.csv"), row.names=FALSE)
