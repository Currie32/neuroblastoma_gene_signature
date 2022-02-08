PATH = "~/Downloads/nbl_target_2018_pub/"


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


patients <- read.csv(file.path(PATH, "data_clinical_patient.txt"), sep="\t")

data <- read.csv(file.path(PATH, "data_mrna_seq_rpkm_zscores_ref_all_samples.txt"), sep="\t")
data2 <- read.csv(file.path(PATH, "data_mrna_agilent_microarray_zscores_ref_all_samples.txt"), sep="\t")

# Load differentially expressed genes
gene_list <- read.csv(file.path(PATH, "gene_list.csv"))

headers <- read.csv(
  file.path(PATH, "data_clinical_patient.txt"),
  nrows = 1,
  sep="\t"
)
patients <- read.csv(
  file.path(PATH, "data_clinical_patient.txt"),
  skip=5,
  header=F,
  sep="\t"
)
colnames(patients) <- headers

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

pa



head(patients$`Tumor Sample Histology`, 10)


"high_risk",
"mycn_status",
"inss_stage_2",
"inss_stage_3",
"inss_stage_4",
"inss_stage_4S"

colnames(patients)

head(patients)

sum(!is.na(patients$`Event Free Survival Censored`))


data
