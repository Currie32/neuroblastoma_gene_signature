# Prepare GSE85047 for analysis

library(GEOquery)
library(fastDummies)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"


prepare_patient_data <- function(gse) {
  #' Prepare the patient data for analysis
  #' 
  #' gse List: contains information about the study
  #' 
  #' return data.frame: the patient data
  
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
  #' Rename the features of the patient data
  #' 
  #' df data.frame: the patient data to be renamed
  #' 
  #' returns data.frame: the renamed patient data
  
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
  #' Format the patient data so that it is ready for analysis
  #' 
  #' df data.frame: the patient data
  #' 
  #' return data.frame: the formatted patient data
  
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
  #' Prepare the Target expression data
  #' 
  #' gse List: contains information about the study
  #' gene_list List(str): the differentially expressed gene names to filter the
  #'                      the expression data to.
  #'                       
  #' returns data.frame: the prepared expression data
  
  # Get gene names for expression data
  genes <- gene_names(gse)
  
  # Extract the expression data
  expression_data <- gse[["GSE85047_series_matrix.txt.gz"]]@assayData[["exprs"]]
  expression_data <- data.frame(expression_data)
  
  # Remove rows that are missing a gene name
  expression_data <- expression_data[genes != "",]
  genes <- genes[genes != ""]
  
  # Average expression data from duplicated gens
  ## Identify which genes have multiple entries
  genes_duplicated <- unique(genes[duplicated(genes)])
  
  # Stored the averaged expression data
  expression_data_mean <- data.frame()

  for (duplicated_gene in genes_duplicated) {
    
    # Identify the gene's rows
    indicies <- duplicated_gene == genes
    
    # Get the mean expression values for the gene
    mean_values <- data.frame(t(colMeans(expression_data[indicies, ])))
    
    # Add the averaged expression data to the dataframe
    expression_data_mean <- rbind(expression_data_mean, mean_values)
  }
  
  # Identify the indicies of the genes with multiple readings
  duplicated_indicies <- genes %in% genes_duplicated
  
  # Remove the original expression data from the data
  expression_data <- expression_data[!duplicated_indicies,]
  
  # Remove the duplicated gene entries from the list of gene names
  genes <- genes[!duplicated_indicies]
  
  # Add the averaged expression data to the dataframe
  expression_data <- rbind(expression_data, expression_data_mean)
  
  # Add the names of the duplicated genes to the list of gene names
  genes <- c(genes, genes_duplicated)
  
  # Transpose the data to have the genes as columns
  expression_data_t <- data.frame(t(expression_data))
  
  # Label columns with gene names
  colnames(expression_data_t) <- genes
  
  # Filter to genes in gene_list
  expression_data_t <- expression_data_t[colnames(expression_data_t) %in% gene_list]
  
  # Z-transform the expression data
  expression_data_t <- data.frame(scale(expression_data_t))

  # Add geo_id for merging with patients data
  expression_data_t$geo_id <- colnames(expression_data)
  
  return(expression_data_t)
}


gene_names <- function(gse) {
  #' Get the gene names from the GEO data
  #' 
  #' gse List: contains information about the study
  #' 
  #' return List(str): Names of the genes used in the study
  
  feature_data <- gse[["GSE85047_series_matrix.txt.gz"]]@featureData@data
  genes <- str_split_fixed(feature_data$gene_assignment, " // ", n=3)[,2]
  
  return(genes)
}


# Load differentially expressed genes
gene_list <- read.csv(file.path(PATH, "gene_list/gene_list.csv"))

# Load GEO data
gse <- getGEO("GSE85047")

# Prepare datasets for merging
patients <- prepare_patient_data(gse)
expression_data <- prepare_expression_data(gse, gene_list$external_gene_name)

# Merge the patient data with their expression data
patients <- merge(
  patients,
  expression_data,
  on="sequence_id",
)

# Save the patient data
write.csv(patients, file.path(PATH, "processed/GSE85047.csv"), row.names=FALSE)
