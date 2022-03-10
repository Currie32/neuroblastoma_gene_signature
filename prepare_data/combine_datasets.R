# Define Functions
extract_patient_traits <- function(GSE49711, GSE62564) {
  
  patient_traits <- GSE49711[c(
    "sequence_id",
    "geo_id",
    "male",
    "age_at_diagnosis_days",
    "high_risk",
    "mycn_status",
    "inss_stage"
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


format_GSE_patient_data <- functions(gene_list, patients_GSE, patients_traits) {
  
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


combine_GSE_datasets <- function(GSE1, GSE2) {
  
  patients <- bind_rows(GSE1, GSE2)
  patients[is.na(patients_all)] <- 0
  
  return(patients)
}


path <- "~/Imperial/neuroblastoma_gene_signature/data/"

GSE49711 <- read.csv(file.path(path, "GSE49711.csv"))
GSE62564 <- read.csv(file.path(path, "GSE62564.csv"))

# Load differentially expressed genes
gene_list <- read.csv(file.path(path, "gene_list.csv"))

GSE49711 <- format_GSE_patient_data(
  gene_list$external_gene_name, 
  GSE49711,
  patients_traits
)

GSE62564 <- format_GSE_patient_data(
  gene_list$external_gene_name, 
  GSE62564,
  patients_traits
)

training <- combine_datasets(GSE49711, GSE62564)

write.csv(training, file.path(path, "training.csv"), row.names=FALSE)
