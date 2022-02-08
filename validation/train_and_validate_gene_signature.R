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
    as.formula(paste("survival ~ ", features_joined)),
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


plot_kaplan_meier <- function(training, survival) {
  
  fit1 <- survfit(survival ~ risk_score_low, data=training, na.action=na.pass)
  ggsurvplot(
    fit = fit1,
    conf.int=TRUE,
    pval = TRUE,
    xlab = "Days", 
    ylab = "Overall survival probability"
  )
}


# Load training data
training <- read.csv(file.path(PATH, "training.csv"))
training$under_18_months <- is_under_18_months(training)

# Load gene list
gene_list <- read.csv(file.path(PATH, "gene_list.csv"))
genes <- colnames(training[, (colnames(training) %in% gene_list$external_gene_name)])

# Create survival for cox regressions
survival <- Surv(training$event_free_survival_days, training$event_free_survival)

# Find the combination of features that are most predictive of survival
candidate_genes <- univariate_cox_regression(genes, training)
features_significant <- multivariate_cox_regression(candidate_genes, survival, training)
features_best <- aic_feature_selection(features_significant, survival, training)
model <- train_model(features_best, training)

# Create risk scores
training$risk_score <- predict(model, type="lp")
training$risk_score_low <- training$risk_score < median(training$risk_score)

plot_roc_curve(training)
plot_kaplan_meier(training, survival)

# Plot hazard scores
ggforest(model, data = training)

# Plot gene expression heatmap
pheatmap(
  training[order(training$risk_score),][c(features_best)][,2:length(features_best)],
  show_rownames=FALSE,
  cluster_rows=F,
  cluster_cols=F
)

# Compare event_free_survival_days between risk groups
ggplot(training, aes(x=risk_score_low, y=event_free_survival_days)) + 
  geom_boxplot()

