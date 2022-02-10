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
  candidate_genes <- p_values[p_values$p_values < 0.01,]$genes
  length(candidate_genes)
  
  return(candidate_genes)
}


multivariate_cox_regression <- function(candidate_genes, survival, df) {
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
  # canadidate_genes_joined <- paste("(", canadidate_genes_joined, ")^2")
  
  # Multivariate cox regression
  model <- coxph(
    as.formula(paste(
      "survival ~ under_18_months + male + high_risk + mycn_status + inss_stage_2 + inss_stage_3 + inss_stage_4 + inss_stage_4S +",
      canadidate_genes_joined)),
    data=df,
    na.action=na.pass
  )
  
  # Find the significant features
  result <- summary(model)
  p_values <- result$coefficients[,5]
  significant_features <- names(p_values[p_values < 0.01])
  
  return(significant_features)
}


aic_feature_selection <- function(features, survival, df) {
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
  dredge_results <- summary(model.avg(models, subset=delta < 4, fit=TRUE))
  
  # Filter to the best performing features
  best_model_features_integers <- rownames(dredge_results$msTable)[1]
  best_model_features_integers <- as.integer(strsplit(best_model_features_integers, '/')[[1]])
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


plot_roc_curve <- function(df) {
  
  precrec_obj <- evalmod(
    scores=df$risk_score,
    labels=df$event_free_survival,
  )
  autoplot(precrec_obj)
  auc(precrec_obj)
}


# Load training data
train <- read.csv(file.path(PATH, "train.csv"))
test_GSE85047 <- read.csv(file.path(PATH, "test_GSE85047.csv"))
test_target <- read.csv(file.path(PATH, "test_target.csv"))

# Add under_18_months feature
train$under_18_months <- is_under_18_months(train)
test_GSE85047$under_18_months <- is_under_18_months(test_GSE85047)
test_target$under_18_months <- is_under_18_months(test_target)

# Load gene list
gene_list <- read.csv(file.path(PATH, "gene_list.csv"))
genes <- colnames(train[, (colnames(train) %in% gene_list$external_gene_name)])

# Create survival for cox regressions
survival_train <- Surv(train$event_free_survival_days, train$event_free_survival)
survival_test_GSE85047 <- Surv(test_GSE85047$event_free_survival_days, test_GSE85047$event_free_survival)
survival_test_target <- Surv(test_target$event_free_survival_days, test_target$event_free_survival)

# Find the combination of features that are most predictive of survival
candidate_genes <- univariate_cox_regression(genes, train)
features_significant <- multivariate_cox_regression(candidate_genes, survival_train, train)
features_best <- aic_feature_selection(features_significant, survival_train, train)
model <- train_model(features_best, train)

# Create risk scores
train$risk_score <- predict(model, type="lp")
train$risk_score_low <- train$risk_score < median(train$risk_score)

test_GSE85047$risk_score <- predict(model, test_GSE85047, type="lp")
test_GSE85047$risk_score_low <- test_GSE85047$risk_score < median(train$risk_score)

test_target$risk_score <- predict(model, test_target, type="lp")
test_target$risk_score_low <- test_target$risk_score < median(train$risk_score)

plot_roc_curve(train)
plot_roc_curve(test_GSE85047)
plot_roc_curve(test_target)

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

# Plot hazard scores
ggforest(model, data=train)

# Plot gene expression heatmap
pheatmap(
  train[order(train$risk_score),][c(features_best)][,2:length(features_best)],
  show_rownames=FALSE,
  cluster_rows=F,
  cluster_cols=F
)

# Compare event_free_survival_days between risk groups
ggplot(train, aes(x=risk_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

ggplot(test_GSE85047, aes(x=risk_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

ggplot(test_target, aes(x=risk_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

