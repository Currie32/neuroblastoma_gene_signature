PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"

NAMES_DATA <- c("train", "GSE85047", "target", "E_TABM_38", "wilms_tumor")

load_data <- function() {
  #' Load the datasets used in the analysis.
  #' 
  #' return List(data.frame): List of the dataframes used in the analysis
  
  # Load the datasets
  train <- read.csv(file.path(PATH, "train.csv"))
  GSE85047 <- read.csv(file.path(PATH, "test_GSE85047.csv"))
  target <- read.csv(file.path(PATH, "test_target.csv"))
  E_TABM_38 <- read.csv(file.path(PATH, "test_E_TABM_38.csv"))
  wilms_tumor <- read.csv(file.path(PATH, "test_wilms_tumor.csv"))
  
  # Create a list to store all of the data
  data <- list(train, GSE85047, target, E_TABM_38, wilms_tumor)
  names(data) <- NAMES_DATA
  
  return(data)
}

# Load the data used in the analysis
data <- load_data()

# Vectors for the stats to keep
n_patients <- c()
is_male <- c()
under_18_months <- c()
event_free_survival <- c()
event_free_survival_days <- c()
event_free_survival_days_0 <- c()
event_free_survival_days_1 <- c()
age_at_diagnosis_days <- c()
age_at_diagnosis_days_0 <- c()
age_at_diagnosis_days_1 <- c()
inss_stage_1 <- c()
inss_stage_2 <- c()
inss_stage_3 <- c()
inss_stage_4 <- c()
inss_stage_4S <- c()


# Loop through each dataset and calculate the summary stats
for (i in seq(length(data))) {
  
  n_patients[i] <- nrow(data[[i]])
  
  is_male[i] <- round(mean(data[[i]]$male), 2)
  
  under_18_months[i] <- round(mean(data[[i]]$under_18_months), 2)
  
  event_free_survival[i] <- round(mean(data[[i]]$event_free_survival), 2)
  
  event_free_survival_days[i] <- round(mean(data[[i]]$event_free_survival_days), 0)
  event_free_survival_days_0[i] <- round(mean(data[[i]]$event_free_survival_days[data[[i]]$event_free_survival == 0]), 0)
  event_free_survival_days_1[i] <- round(mean(data[[i]]$event_free_survival_days[data[[i]]$event_free_survival == 1]), 0)
  
  age_at_diagnosis_days[i] <- round(mean(data[[i]]$age_at_diagnosis_days), 0)
  age_at_diagnosis_days_0[i] <- round(mean(data[[i]]$age_at_diagnosis_days[data[[i]]$event_free_survival == 0]), 0)
  age_at_diagnosis_days_1[i] <- round(mean(data[[i]]$age_at_diagnosis_days[data[[i]]$event_free_survival == 1]), 0)
  
  inss_stage_1[i] <- round(
    1
    - mean(data[[i]]$inss_stage_2)
    - mean(data[[i]]$inss_stage_3)
    - mean(data[[i]]$inss_stage_4)
    - mean(data[[i]]$inss_stage_4S)
  , 2)
  inss_stage_2[i] <- round(mean(data[[i]]$inss_stage_2), 2)
  inss_stage_3[i] <- round(mean(data[[i]]$inss_stage_3), 2)
  inss_stage_4[i] <- round(mean(data[[i]]$inss_stage_4), 2)
  inss_stage_4S[i] <- round(mean(data[[i]]$inss_stage_4S), 2)
}

# Put all the stats in a dataframe
stats <- data.frame(
  n_patients,
  is_male,
  under_18_months,
  event_free_survival,
  event_free_survival_days,
  event_free_survival_days_0,
  event_free_survival_days_1,
  age_at_diagnosis_days,
  age_at_diagnosis_days_0,
  age_at_diagnosis_days_1,
  inss_stage_1,
  inss_stage_2,
  inss_stage_3,
  inss_stage_4,
  inss_stage_4S
)

# Add the names of the dataframe
rownames(stats) <- NAMES_DATA

# Replace 0 with N/A since this metric is missing from the dataset
stats$is_male[2] <- "N/A"

# View the summary data
View(t(stats))
