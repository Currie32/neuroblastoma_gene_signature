library(ggplot2)
library(MuMIn)
library(pheatmap)
library(precrec)
library(rms)
library(stringr)
library(survival)
library(survminer)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"
P_VALUE = 0.01

# The day to split the training and testing data so that
# the proportional hazards are consistent over time. More info: 
# https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
CUT_DAY <- c(420)


load_training_data <- function() {
  
  train <- read.csv(file.path(PATH, "train.csv"))
  
}

load_validation_data <- function() {
  
  # Load the datasets
  test_GSE85047 <- read.csv(file.path(PATH, "test_GSE85047.csv"))
  test_target <- read.csv(file.path(PATH, "test_target.csv"))
  test_E_TABM_38 <- read.csv(file.path(PATH, "test_E_TABM_38.csv"))
  test_wilms_tumor <- read.csv(file.path(PATH, "test_wilms_tumor.csv"))
  
  # Create a list to store all of the data
  data <- list(test_GSE85047, test_target, test_E_TABM_38, test_wilms_tumor)
  
  return(data)
}


survival_data_train <- function(data) {
  
  survival <- Surv(data$event_free_survival_days, data$event_free_survival)
  return(survival)
  
}


survival_data_validation <- function(data) {
  
  # For storing the survival objects
  survival <- list()
  
  # Loop through the datasets and add the survival object
  for (i in seq(length(data))) {
    
    df <- data[[i]]
    df_survival <- Surv(df$event_free_survival_days, df$event_free_survival)
    
    survival[[i]] <- df_survival
  }
  
  return(survival)
}


split_data_train <- function(train, train_survival) {
  
  train_split <- survSplit(
    train_survival ~ ., 
    data=train,
    cut=CUT_DAY,
    zero=-1,
    episode="time_group",
  )
  
  return(train_split)
}


split_data_validation <- function(data, survival) {
  
  # Store the newly splitted dataframes
  splits <- c()
  
  # Split the data and add it to the list
  for (i in seq(length(data))) {
    df <- data[[i]]
    df_survival <- survival[[i]]
    
    df_split <- survSplit(
      df_survival ~ ., 
      data=df,
      cut=CUT_DAY,
      zero=-1,
      episode="time_group",
    )
    
    splits[[i]] <- df_split
  }
  
  return(splits)
}


modelling_genes <- function(genes, train, survival, interaction) {
  
  # Find the combination of features that are most predictive of survival
  genes_candidate <- univariate_cox_regression(genes, train, survival)
  genes_significant <- multivariate_cox_regression(genes_candidate, survival, train, FALSE, FALSE)
  
  if (interaction) {
    genes_significant <- multivariate_cox_regression(genes_significant, survival, train, FALSE, TRUE)
  }
  
  genes_best <- aic_feature_selection(genes_significant, survival, train, '')
  
  return(genes_best)
}


univariate_cox_regression <- function(genes, data, survival) {
  #' Find the genes that are statistically significant in a univariate cox regression.
  #'
  #' genes (List(str)): Names of genes
  #' df (data.frame): Contains survival and gene expression data
  #'
  #' return List(str): Names of statistically significant genes
  
  # Store the p-values
  candidate_genes <- c()
  
  for (gene in genes) {
    
    # Perform cox regression
    result <- coxph(
      as.formula(paste("survival ~ ", gene)),
      data=data,
      na.action=na.pass
    )
    
    # Only add a gene to the list of candidates if its p-value is
    # under the threshold
    p_value <- summary(result)$coefficients[5]
    if (p_value <= P_VALUE) {
      candidate_genes <- append(candidate_genes, gene)
    }
  }
  
  return(candidate_genes)
}


multivariate_cox_regression <- function(
  features, survival, data, use_phenotype_features, use_interaction
) {
  #' Find which features (genes and patient traits) are statistically
  #' significant in combination
  #' 
  #' features (List(str)): genes that were significant during the
  #'                              univariate cox regression
  #' survival (double): survival information for each patient
  #' df (data.frame): contains the expression and survival data for each patient
  #' 
  #' returns (List(str)): the statistically significant features
  
  if (use_phenotype_features) {
    
    # List of phenotype features to use
    features_all = c(
      "under_18_months",
      "mycn_amplification",
      "inss_stage_2",
      "inss_stage_3",
      "inss_stage_4",
      "inss_stage_4S",
      features
    )
    features_joined <- paste(features_all, collapse=" + ")
  
    # Create a model with the phenotype features
    model <- coxph(
      as.formula(paste("survival ~ ", features_joined)),
      data=data,
      na.action=na.pass
    )
  }
  
  else if (use_interaction) {
    
    # Format the candidate genes the cox regression
    candidate_genes_joined <- paste(features, collapse=" + ")
    
    # Create an interaction model
    model <- coxph(
      as.formula(paste("survival ~ (", candidate_genes_joined, ")^2")),
      data=data,
      na.action=na.pass
    )
  }
  
  else {
    
    # Format the candidate genes the cox regression
    candidate_genes_joined <- paste(features, collapse=" + ")

    # Create a multivariate cox-regression model  
    model <- coxph(
      as.formula(paste("survival ~ ", candidate_genes_joined)),
      data=data,
      na.action=na.pass
    )
  }
  
  # Find the significant features
  result <- summary(model)
  p_values <- result$coefficients[,5]
  significant_features <- names(p_values[p_values < P_VALUE])
 
  return(significant_features)
}


aic_feature_selection <- function(features, survival, data, delimiter) {
  #' Use AIC to find the optimal combination of features.
  #' 
  #' features (List(str)): the statistically significant features from the
  #'                    multivariate cox regression
  #' survival (double): survival information for each patient
  #' df (data.frame): contains the expression and survival data for each patient
  #' 
  #' returns (List(str)): the best combination of features according to AIC 
  
  # Format the features for the cox regression
  features_joined <- paste(features, collapse = " + ")
  
  # Fit top model
  top_model <- coxph(
    as.formula(paste("survival ~ ", features_joined)),
    data=data,
    na.action=na.pass
  )
  
  # "Dredge" all subsets of the top model
  models <- dredge(top_model)
  dredge_results <- summary(model.avg(models, fit=TRUE))
  
  # Filter to the best performing features
  features <- rownames(dredge_results$coefmat.full)
  best_model_features_integers <- rownames(dredge_results$msTable)[1]
  best_model_features_integers <- as.integer(strsplit(best_model_features_integers, delimiter)[[1]])
  best_model_features <- features[best_model_features_integers]
  
  return(best_model_features)
}


baseline_model <- function(train_split, train_survival) {
  
  # The baseline features for the model, i.e. all features, but the
  # gene expression ones
  features <- c(
    "under_18_months",
    "mycn_amplification",
    "inss_stage_2",
    "inss_stage_3",
    "inss_stage_4",
    "inss_stage_4S"
  )
  features_joined <- paste(features, collapse = " + ")

  top_model_baseline <- coxph(
    as.formula(paste("train_survival ~", features_joined)),
    data=train_split,
    na.action=na.pass
  )
  
  # "Dredge" all subsets of the top model
  models_baseline <- dredge(top_model_baseline)
  
  # Get the dredge results
  dredge_results_baseline <- summary(model.avg(models_baseline, fit=TRUE))
  
  # Identify the features of the dredged models
  features_baseline <- rownames(dredge_results_baseline$coefmat.full)
  
  # Identify the best combination of features
  best_model_baseline_features_integers <- rownames(dredge_results_baseline$msTable)[1]
  best_model_baseline_features_integers <- as.integer(strsplit(best_model_baseline_features_integers, '/')[[1]])
  best_model_baseline_features <- features_baseline[as.integer(strsplit(as.character(best_model_baseline_features_integers), split="")[[1]])]
  
  # Train a model using the best combination of features
  features_joined <- paste(paste(best_model_baseline_features, ":strata(time_group)"), collapse = " + ")

  model_baseline <- coxph(
    train_survival ~ mycn_amplification + under_18_months:strata(time_group) + inss_stage_2:strata(time_group) + inss_stage_3:strata(time_group) + inss_stage_4:strata(time_group) + inss_stage_4S:strata(time_group),
    data=train_split,
    na.action=na.pass
  )
  
  return(model_baseline)
}


no_interaction_model <- function(train_split, train_survival) {
  
  model <- coxph(
    train_survival ~ CD147:strata(time_group) + CD9:strata(time_group) + DLG2 + HEBP2:strata(time_group) + HSD17B12 + NXT2:strata(time_group) + RACK1 + TXNDC5:strata(time_group),
    data=train_split,
    na.action=na.pass
  )
  
  return(model)
}


interaction_model <- function(train_split, train_survival) {
  
  model <- coxph(
    train_survival ~ CD9:strata(time_group) + DLG2 + HEBP2:strata(time_group) + HSD17B12 + NXT2:strata(time_group) + TXNDC5:strata(time_group) + CD9:NXT2:strata(time_group) + NXT2:TXNDC5,
    data=train_split,
    na.action=na.pass
  )
  
  return(model)
}


predict_risk_score <- function(model, data, risk_score) {
  
  # Create the risk score for the training dataframe
  data[[1]][risk_score] <- predict(model, data_split[[1]])
  data[[1]][paste(risk_score, "low", sep="_")] <- data[[1]][risk_score] < median(data[[1]][risk_score][[1]])
  
  # Calculate the median train risk score to be used as the high/low threshold
  median_train_risk_score <- median(data[[1]][risk_score][[1]])
  
  # Add risk scores to the remaining dataframes
  for (i in seq(2, length(data))) {
    data[[i]][risk_score] <- predict(model, data[[i]])
    data[[i]][paste(risk_score, "low", sep="_")] <- data[[i]][risk_score] < median_train_risk_score
  }
  
  return(data)
}


prognostic_model <- function(data, survival, risk_score) {
  
  features_significant <- multivariate_cox_regression(c(risk_score), survival, data, TRUE, FALSE)
  features_best <- aic_feature_selection(features_significant, survival, data, '')
  
  if (risk_score == "risk_score_interaction") {
    data$risk_score <- data$risk_score_interaction
  }
  else if (risk_score == "risk_score_no_interaction") {
    data$risk_score <- data$risk_score_no_interaction
  }
  
  model_prognostic <- coxph(
    survival ~ inss_stage_2 + inss_stage_3 + inss_stage_4:strata(time_group) + inss_stage_4S + under_18_months:strata(time_group) + risk_score:strata(time_group),
    data=data,
  )
  
  return(model_prognostic)
}


compare_risk_scores <- function(data) {
  
  # Names of the predicted risk scores
  risk_scores <- c(
    "risk_score_baseline",
    "risk_score_no_interaction",
    "risk_score_interaction",
    "prognostic_score_no_interaction",
    "prognostic_score_interaction"
  )
  # Names of the datasets
  data_names <-  c(
    "train", "GSE85047", "target", "E_TABM-38", "wilms_tumor"
  )
  
  results <- list()
  
  for (i in seq(length(risk_scores))) {
    
    auc_scores <- c()
    
    for (ii in seq(length(data))) {
      precrec_obj <- evalmod(
        scores=data[[ii]][risk_scores[i]],
        labels=data[[ii]]$event_free_survival,
      )
      auc_score <- round(auc(precrec_obj)$auc[1], 4)
      auc_scores <- append(auc_scores, auc_score)
    }
    results[[i]] <- auc_scores
  }
  names(results) <- risk_scores
  results <- data.frame(results)
  rownames(results) <- data_names
  
  pheatmap(
    results,
    cluster_rows=F,
    cluster_cols=F,
    display_numbers=T
  )
  
  return(results)
}


kaplan_meier_plot <- function() {
  
  # Need to create this object because of bug in R.
  survival_object <- survival_split[[1]]
  
  # Kaplan Meier Train
  fit <- survfit(
    survival_object ~ prognostic_score_interaction_low,
    data=data_split[[1]],
    na.action=na.pass
  )
  ggsurvplot(
    fit=fit,
    conf.int=TRUE,
    pval=TRUE,
    risk.table="abs_pct",
    xlab="Years", 
    ylab="Overall survival probability",
    xscale=365,
    break.x.by=365,
    xlim=c(0, 1826)
  )
  
  # Other possible validation and diagnostic tools
  # summary(fit_train)$table
  # ggsurvplot(fit_train, fun="cloglog")
  # ggcoxdiagnostics(model, type="schoenfeld")
}


plot_roc_curve <- function(data, risk_score) {
  
  precrec_obj <- evalmod(
    scores=data[[risk_score]],
    labels=data$event_free_survival,
  )
  autoplot(precrec_obj)
}


compare_auc_values <- function(data, samples) {
  
  auc_scores <- vector()
  auc_scores_baseline <- vector()
  
  for (i in 1:samples) {
    df_sample <- data[sample(nrow(data), length(data), replace=TRUE), ]
    
    precrec_obj <- evalmod(
      scores=df_sample$prognostic_score,
      labels=df_sample$event_free_survival,
    )
    precrec_obj_baseline <- evalmod(
      scores=df_sample$risk_score_baseline,
      labels=df_sample$event_free_survival,
    )
    
    auc_score <- auc(precrec_obj)$aucs[1]
    auc_scores <- append(auc_scores, auc_score)
    
    auc_score_baseline <- auc(precrec_obj_baseline)$aucs[1]
    auc_scores_baseline <- append(auc_scores_baseline, auc_score_baseline)
  }
  
  t.test(auc_scores, auc_scores_baseline, alternative="greater")
}


# Prepare the data for training
train <- load_training_data()
train_survival <- survival_data_train(train)
train_split <- split_data_train(train, train_survival)

validation <- load_validation_data()
validation_survival <- survival_data_validation(validation)
validation_split <- split_data_validation(validation, validation_survival)

# Get the genes from the training dataset
genes <- colnames(train)[14:length(colnames(train)) - 1]

# Identify the genes for modelling
genes_modelling_no_interaction <- modelling_genes(genes, train, train_survival, FALSE)
genes_modelling_interaction <- modelling_genes(genes, train, train_survival, TRUE)

# Create models
model_baseline <- baseline_model(train_split, train_survival)
model_no_interaction <- no_interaction_model(train_split, train_survival)
model_interaction <- interaction_model(train_split, train_survival)

# Check that all models meet the assumptions of the Cox proportional hazards model
zph_model_baseline <- cox.zph(model_baseline)
zph_model_no_interaction <- cox.zph(model_no_interaction)
zph_model_interaction <- cox.zph(model_interaction)

# All values should be below 0.05
print(zph_model_baseline)
print(zph_model_no_interaction)
print(zph_model_interaction)

# Plot the diagnostic plots to confirm there is no trend
ggcoxzph(zph_model_baseline)
ggcoxzph(zph_model_no_interaction)
ggcoxzph(zph_model_interaction)

# Add risk scores
data_split <- predict_risk_score(model_baseline, data_split, "risk_score_baseline")
data_split <- predict_risk_score(model_no_interaction, data_split, "risk_score_no_interaction")
data_split <- predict_risk_score(model_interaction, data_split, "risk_score_interaction")

# Update data_split_train with the risk scores
data_split_train <- data_split[[1]]

# Train the prognostic models
model_prognostic_no_interaction <- prognostic_model(
  data_split_train, survival_train_split, "risk_score_no_interaction"
)
model_prognostic_interaction <- prognostic_model(
  data_split_train, survival_train_split, "risk_score_interaction"
)

# Check that all models meet the assumptions of the Cox proportional hazards model
zph_model_prognostic_no_interaction <- cox.zph(model_prognostic_no_interaction)
zph_model_prognostic_interaction <- cox.zph(model_prognostic_interaction)

# All values should be below 0.05
print(zph_model_prognostic_no_interaction)
print(zph_model_prognostic_interaction)

# Plot the diagnostic plots to confirm there is no trend
ggcoxzph(zph_model_prognostic_no_interaction)
ggcoxzph(zph_model_prognostic_interaction)

# Add prognostic scores
data_split <- predict_risk_score(model_prognostic_no_interaction, data_split, "prognostic_score_no_interaction")
data_split <- predict_risk_score(model_prognostic_interaction, data_split, "prognostic_score_interaction")

# Create a heatmap comparing the risk and prognostic scores
results_risk_scores <- compare_risk_scores(data_split)

# Can't use input parameters because of limitations in R
kaplan_meier_plot()

















# Check if difference in AUC scores is significant
# (between baseline and prognostic models)
compare_auc_values(test_GSE85047_split, 500)
compare_auc_values(test_target_split, 500)
compare_auc_values(test_E_TABM_38_split, 500)
compare_auc_values(test_all_split, 500)


train_split$prognostic_score_low <- as.numeric(train_split$prognostic_score_low)
train_split$survival_train_2 <- as.numeric(as.factor(train_split$survival_train))

# Kaplan Meier Train
fit_train <- survfit(
  Surv(train_split$event_free_survival_days, train_split$event_free_survival) ~ prognostic_score_low,
  data=train_split,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_train,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)
ggsurvplot(fit_train, conf.int=TRUE, fun="cloglog")

# Kaplan Meier test_GSE85047
fit_test_GSE85047 <- survfit(
  Surv(test_GSE85047_split$event_free_survival_days, test_GSE85047_split$event_free_survival) ~ prognostic_score_low,
  data=test_GSE85047_split,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_test_GSE85047,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

# Kaplan Meier test_target
fit_test_target <- survfit(
  Surv(test_target_split$event_free_survival_days, test_target_split$event_free_survival) ~ prognostic_score_low,
  data=test_target_split,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_test_target,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

fit_test_all <- survfit(
  Surv(test_all_split$event_free_survival_days, test_all_split$event_free_survival) ~ prognostic_score_low,
  data=test_all_split,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_test_all,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

fit_test_wilms_tumor <- survfit(
  Surv(test_wilms_tumor_split$event_free_survival_days, test_wilms_tumor_split$event_free_survival) ~ prognostic_score_low,
  data=test_wilms_tumor_split,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_test_wilms_tumor,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

# Plot hazard scores
summary(model_prognostic)

train_split$`inss_stage_4_time_group_1` <- ifelse(train_split$time_group == 1, train_split$inss_stage_4, 0)
train_split$`inss_stage_4_time_group_2` <- ifelse(train_split$time_group == 2, train_split$inss_stage_4, 0)
train_split$`risk_score_time_group_1` <- ifelse(train_split$time_group == 1, train_split$risk_score, 0)
train_split$`risk_score_time_group_2` <- ifelse(train_split$time_group == 2, train_split$risk_score, 0)
train_split$`under_18_months_time_group_1` <- ifelse(train_split$time_group == 1, train_split$under_18_months, 0)
train_split$`under_18_months_time_group_2` <- ifelse(train_split$time_group == 2, train_split$under_18_months, 0)

names(model_prognostic)

names(model_prognostic$coefficients) %in% colnames(train_split)

trace(ggforest, edit = TRUE)

library(broom)

tidy(model_prognostic)
gmodel <- glance(model)

eval(model_prognostic$call$data)
attr(model_prognostic$terms, "dataClasses")[-1] <- c(
  "inss_stage_2" = "numeric",
  "inss_stage_3" = "numeric",
  "inss_stage_4S" = "numeric",
  "inss_stage_4_time_group_1" = "numeric",
  "inss_stage_4_time_group_2" = "numeric",
  "risk_score_time_group_1" = "numeric",
  "risk_score_time_group_2" = "numeric",
  "under_18_months_time_group_1" = "numeric",
  "under_18_months_time_group_2" = "numeric"
)
  



terms1 <- terms[1]

names(terms)

as.data.frame(tidy(model_prognostic, conf.int = TRUE))

names(model_prognostic$coefficients) <- c(
  "inss_stage_2",
  "inss_stage_3",
  "inss_stage_4S",
  "inss_stage_4_time_group_1",
  "inss_stage_4_time_group_2",
  "risk_score_time_group_1",
  "risk_score_time_group_2",
  "under_18_months_time_group_1",
  "under_18_months_time_group_2"
)
train_split$under_18_months_time_group_1

ggforest(model, data=train)
ggforest(model_prognostic, data=train_split)
ggforest(model_baseline, data=train)




# Plot gene expression heatmap
pheatmap(
  train_split[order(train_split$risk_score),][colnames(train_split) %in% genes_best],
  show_rownames=FALSE,
  cluster_rows=F,
  cluster_cols=F,
  breaks=seq(-2, 3, by = 0.06),
)




# Check number of high and low prognostic risk scores
table(train_split$prognostic_score_low)
table(test_GSE85047_split$prognostic_score_low)
table(test_target_split$prognostic_score_low)


# Compare event_free_survival_days between risk groups
ggplot(train_split, aes(x=prognostic_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

ggplot(test_GSE85047_split, aes(x=prognostic_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

ggplot(test_target_split, aes(x=prognostic_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

ggplot(train_split, aes(x=prognostic_score_low, y=CD9_NXT2)) + 
  geom_boxplot(aes(ymin=-1, ymax=3))

t.test(train_split$CD9_NXT2[train_split$risk_score_low == 0], train_split$CD9_NXT2[train_split$risk_score_low == 1])

train_split$CD9_NXT2[train_split$risk_score_low == 0]


shapiro.test(train_split$CD9_NXT2[train_split$risk_score_low == 0])
shapiro.test(train_split$CD9_NXT2[train_split$risk_score_low == 1])

train_split$CD9_NXT2 <- train_split$CD9 * train_split$NXT2

train_split$risk_score_low_factor <- as.factor(train_split$risk_score_low)

wilcox.test(CD9_NXT2 ~ risk_score_low, data=train_split)

asd <- data.frame(t(train[colnames(train) %in% c(genes_best)][1:20,]))
asd$gene_name <- row.names(asd)

ed_target <- read.csv(file.path(PATH, "expression_data_target.csv"))

ed_target$KLF4

test_target$GRHL1 <- ed_target$GRHL1
test_target$HDAC5 <- ed_target$HDAC5

genes_found <- c("GRHL1", "KLF4")

genes_best

coxph(as.formula(paste(
  "Surv(event_free_survival_days, event_free_survival) ~ ", "HDAC5")),
  data = test_target)

univariate_cox_regression(c("GRHL1"), test_target)

cor(test_target$event_free_survival_days, test_target$HDAC5)

par(mfrow = c(1, 1))
plot(test_target$event_free_survival_days, test_target$HDAC5)



library(gtools)
combinations(4, 4)


