library(fastDummies)
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


is_under_18_months <- function(df) {
  under_18_months <- df$age_at_diagnosis_days < 548
  under_18_months[under_18_months == TRUE] <- 1
  under_18_months[under_18_months == FALSE] <- 0
  
  return(under_18_months)
}


univariate_cox_regression <- function(genes, df) {
  #' Find the genes that are statistically significant in a univariate cox regression.
  #'
  #' genes (List(str)): Names of genes
  #' df (data.frame): Contains survival and gene expression data
  #'
  #' return List(str): Names of statistically significant genes
  
  # Store the p-values
  p_values <- c()
  
  for (gene in genes) {
    
    # Perform cox regression
    result <- coxph(as.formula(paste(
      "Surv(event_free_survival_days, event_free_survival) ~ ", gene)),
      data = df)
    
    # Get a store p_value
    p_value <- summary(result)$coefficients[5]
    p_values[gene] <- p_value
  }
  
  # Filter to p-values that are statistically significant
  p_values <- data.frame(genes, p_values)
  candidate_genes <- p_values[p_values$p_values < P_VALUE,]$genes
  length(candidate_genes)
  
  return(candidate_genes)
}


multivariate_cox_regression <- function(candidate_genes, survival, df, phenotype_features, interaction) {
  #' Find which features (genes and patient traits) are statistically
  #' significant in combination
  #' 
  #' candidate_genes (List(str)): genes that were significant during the
  #'                              univariate cox regression
  #' survival (double): survival information for each patient
  #' df (data.frame): contains the expression and survival data for each patient
  #' 
  #' returns (List(str)): the statistically significant features
  
  # Format the candidate genes the cox regression
  canadidate_genes_joined <- paste(candidate_genes, collapse = " + ")

  # Multivariate cox regression
  if (phenotype_features) {
    model <- coxph(
      as.formula(paste(
        "survival ~ under_18_months + mycn_amplification + inss_stage_2 + inss_stage_3 + inss_stage_4 + inss_stage_4S +",
        canadidate_genes_joined)),
      data=df,
      na.action=na.pass
    )
  }
  else if (interaction) {
    model <- coxph(
      as.formula(paste("survival_train ~ (", canadidate_genes_joined, ")^2")),
      data=train,
      na.action=na.pass
    )
  }
  else {
    model <- coxph(
      as.formula(paste("survival ~ ", canadidate_genes_joined)),
      data=df,
      na.action=na.pass
    )
  }
  
  # Find the significant features
  result <- summary(model)
  p_values <- result$coefficients[,5]
  significant_features <- names(p_values[p_values < P_VALUE])
  
  return(significant_features)
}




aic_feature_selection <- function(features, survival, df, delimiter) {
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
    data=df,
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


train_model <- function(features, df) {
  
  features_joined <- paste(features, collapse = " + ")
  
  # Train model
  model <- coxph(
    as.formula(paste("survival_train ~ ", features_joined)),
    data=df,
    na.action=na.pass
  )
  
  return(model)
}


baseline_model <- function(df, survival) {
  
  # Fit the top model
  top_model_baseline <- coxph(
    survival_train ~  under_18_months + mycn_amplification + inss_stage_2 + inss_stage_3 + inss_stage_4 + inss_stage_4S,
    data=df,
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
  model_baseline <- train_model(best_model_baseline_features, train)
  
  return(model_baseline)
}


plot_roc_curve <- function(df, risk_score) {
  
  precrec_obj <- evalmod(
    scores=df[[risk_score]],
    labels=df$event_free_survival,
  )
  autoplot(precrec_obj)
  auc(precrec_obj)
}


compare_auc_values <- function(df, samples) {
  
  auc_scores <- vector()
  auc_scores_baseline <- vector()
  
  for (i in 1:samples) {
    df_sample <- df[sample(nrow(df), length(df), replace=TRUE), ]
    
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


# Load training data
train <- read.csv(file.path(PATH, "train.csv"))
test_GSE85047 <- read.csv(file.path(PATH, "test_GSE85047.csv"))
test_target <- read.csv(file.path(PATH, "test_target.csv"))

# Add under_18_months feature
train$under_18_months <- is_under_18_months(train)
test_GSE85047$under_18_months <- is_under_18_months(test_GSE85047)
test_target$under_18_months <- is_under_18_months(test_target)

# Get the genes from the training dataset
genes <- colnames(train)[14:length(colnames(train)) - 1]

# Create survival for cox regressions
survival_train <- Surv(train$event_free_survival_days, train$event_free_survival)
survival_test_GSE85047 <- Surv(test_GSE85047$event_free_survival_days, test_GSE85047$event_free_survival)
survival_test_target <- Surv(test_target$event_free_survival_days, test_target$event_free_survival)

# Find the combination of features that are most predictive of survival
genes_candidate <- univariate_cox_regression(genes, train)
genes_significant <- multivariate_cox_regression(genes_candidate, survival_train, train[genes], FALSE, FALSE)
genes_interaction <- multivariate_cox_regression(genes_significant, survival_train, train[genes], FALSE, TRUE)
genes_best <- aic_feature_selection(genes_interaction, survival_train, train, '')

# Train the model using the best gene set
model <- train_model(genes_best, train)

# Create risk scores
train$risk_score <- predict(model, train)
train$risk_score_low <- train$risk_score < median(train$risk_score)

test_GSE85047$risk_score <- predict(model, test_GSE85047)
test_GSE85047$risk_score_low <- test_GSE85047$risk_score < median(train$risk_score)

test_target$risk_score <- predict(model, test_target)
test_target$risk_score_low <- test_target$risk_score < median(train$risk_score)

# Create baseline model risk scores
model_baseline <- baseline_model(train, survival_train)
train$risk_score_baseline <- predict(model_baseline, train, type="lp")
train$risk_score_low_baseline <- train$risk_score_baseline < median(train$risk_score_baseline)

test_GSE85047$risk_score_baseline <- predict(model_baseline, test_GSE85047, type="lp")
test_GSE85047$risk_score_low_baseline <- test_GSE85047$risk_score_baseline < median(train$risk_score_baseline)

test_target$risk_score_baseline <- predict(model_baseline, test_target, type="lp")
test_target$risk_score_low_baseline <- test_target$risk_score_baseline < median(train$risk_score_baseline)


plot_roc_curve(train, "risk_score")
plot_roc_curve(test_GSE85047, "risk_score")
plot_roc_curve(test_target, "risk_score")

plot_roc_curve(train, "risk_score_baseline")
plot_roc_curve(test_GSE85047, "risk_score_baseline")
plot_roc_curve(test_target, "risk_score_baseline")


# Kaplan Meier Train
fit_train <- survfit(
  survival_train ~ risk_score_low,
  data=train,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_train,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

# Kaplan Meier test_GSE85047
fit_test_GSE85047 <- survfit(
  survival_test_GSE85047 ~ risk_score_low,
  data=test_GSE85047,
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
  survival_test_target ~ risk_score_low,
  data=test_target,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_test_target,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

# Train model with risk score and pheontype features
features_significant <- multivariate_cox_regression(c('risk_score'), survival_train, train, TRUE, FALSE)
features_best <- aic_feature_selection(features_significant, survival_train, train, '')
model <- train_model(features_best, train)

# Create prognostic score (i.e. the more advanced risk score)
train$prognostic_score <- predict(model, train)
train$prognostic_score_low <- train$prognostic_score < median(train$prognostic_score)

test_GSE85047$prognostic_score <- predict(model, test_GSE85047)
test_GSE85047$prognostic_score_low <- test_GSE85047$prognostic_score < median(train$prognostic_score)

test_target$prognostic_score <- predict(model, test_target)
test_target$prognostic_score_low <- test_target$prognostic_score < median(train$prognostic_score)

# Plot results
plot_roc_curve(train, "prognostic_score")
plot_roc_curve(test_GSE85047, "prognostic_score")
plot_roc_curve(test_target, "prognostic_score")

# Kaplan Meier Train
fit_train <- survfit(
  survival_train ~ prognostic_score_low,
  data=train,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_train,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

# Kaplan Meier test_GSE85047
fit_test_GSE85047 <- survfit(
  survival_test_GSE85047 ~ prognostic_score_low,
  data=test_GSE85047,
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
  survival_test_target ~ risk_score_low_baseline,
  data=test_target,
  na.action=na.pass
)
ggsurvplot(
  fit = fit_test_target,
  conf.int=TRUE,
  pval = TRUE,
  xlab = "Days", 
  ylab = "Overall survival probability"
)

# Plot hazard scores
ggforest(model, data=train)
ggforest(model_baseline, data=train)

# Check if difference in AUC scores is significant
# (between baseline and prognostic models)
compare_auc_values(test_GSE85047, 500)
compare_auc_values(test_target, 500)

# Plot gene expression heatmap
pheatmap(
  train[order(train$risk_score),][colnames(train) %in% genes_best],
  show_rownames=FALSE,
  cluster_rows=F,
  cluster_cols=F,
  breaks=seq(-3, 3, by = 0.1),
)


# Check number of high and low prognostic risk scores
table(train$prognostic_score_low)
table(test_GSE85047$prognostic_score_low)
table(test_target$prognostic_score_low)


# Compare event_free_survival_days between risk groups
ggplot(train, aes(x=prognostic_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

ggplot(test_GSE85047, aes(x=prognostic_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

ggplot(test_target, aes(x=prognostic_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

