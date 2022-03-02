source("~/Imperial/neuroblastoma_gene_signature/validation/plot_ggforest.R")

library(ggplot2)
library(MuMIn)
library(pheatmap)
library(pROC)
library(rms)
library(stringr)
library(survival)
library(survivalROC)
library(survminer)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/"
P_VALUE = 0.01

# The day to split the training and testing data so that
# the proportional hazards are consistent over time. More info: 
# https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
CUT_DAY <- c(410)

# Names of the predicted risk scores
RISK_SCORES <- c(
  "risk_score_baseline",
  "risk_score_no_interaction",
  "risk_score_interaction",
  "prognostic_score_no_interaction",
  "prognostic_score_interaction"
)

# Names of the datasets in each combination
VALIDATION_COMBOS_NAMES <- c(
  "GSE85047", "target", "E-TABM-38",
  "GSE85047_target", "GSE85047_E-TABM-38", "target_E-TABM-38",
  "all"
)


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


univariate_cox_regression <- function(genes, train, survival) {
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
      data=train,
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
  features, train_survival, train, use_phenotype_features, use_interaction
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
      as.formula(paste("train_survival ~ ", features_joined)),
      data=train,
      na.action=na.pass
    )
  }
  
  else if (use_interaction) {
    
    # Format the candidate genes the cox regression
    candidate_genes_joined <- paste(features, collapse=" + ")
    
    # Create an interaction model
    model <- coxph(
      as.formula(paste("train_survival ~ (", candidate_genes_joined, ")^2")),
      data=train,
      na.action=na.pass
    )
  }
  
  else {
    
    # Format the candidate genes the cox regression
    candidate_genes_joined <- paste(features, collapse=" + ")

    # Create a multivariate cox-regression model  
    model <- coxph(
      as.formula(paste("train_survival ~ ", candidate_genes_joined)),
      data=train,
      na.action=na.pass
    )
  }
  
  # Find the significant features
  result <- summary(model)
  p_values <- result$coefficients[,5]
  significant_features <- names(p_values[p_values < P_VALUE])
 
  return(significant_features)
}


aic_feature_selection <- function(features, train_survival, train, delimiter) {
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
    as.formula(paste("train_survival ~ ", features_joined)),
    data=train,
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


baseline_model <- function(train, train_survival) {
  
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
    data=train,
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
    data=train,
    na.action=na.pass
  )
  
  return(model_baseline)
}


no_interaction_model <- function(train, train_survival) {
  
  model <- coxph(
    train_survival ~ CD147:strata(time_group) + CD9:strata(time_group) + DLG2 + HEBP2:strata(time_group) + HSD17B12 + NXT2:strata(time_group) + RACK1 + TXNDC5:strata(time_group),
    data=train,
    na.action=na.pass
  )
  
  return(model)
}


interaction_model <- function(train, train_survival) {
  
  model <- coxph(
    train_survival ~ CD9:strata(time_group) + DLG2 + HEBP2:strata(time_group) + HSD17B12 + NXT2:strata(time_group) + TXNDC5:strata(time_group) + CD9:NXT2:strata(time_group) + NXT2:TXNDC5,
    data=train,
    na.action=na.pass
  )
  
  return(model)
}

predict_risk_score <- function(model, train, validation, risk_score) {
  
  # Create the risk score for the training dataframe
  train[risk_score] <- predict(model, train)
  train[paste(risk_score, "low", sep="_")] <- train[risk_score] < median(train[risk_score][[1]])
  
  # Calculate the median train risk score to be used as the high/low threshold
  median_train_risk_score <- median(train[risk_score][[1]])
  
  # Add risk scores to the remaining dataframes
  for (i in seq(length(validation))) {
    validation[[i]][risk_score] <- predict(model, validation[[i]])
    validation[[i]][paste(risk_score, "low", sep="_")] <- validation[[i]][risk_score] < median_train_risk_score
  }
  
  results <- list()
  results$train <- train
  results$validation <- validation
  
  return(results)
}


prognostic_model <- function(train, train_survival, risk_score) {
  
  features_significant <- multivariate_cox_regression(risk_score, train_survival, train, TRUE, FALSE)
  features_best <- aic_feature_selection(features_significant, train_survival, train, '')
  
  if (risk_score == "risk_score_no_interaction") {

    model <- coxph(
      train_survival ~ inss_stage_2 + inss_stage_3 + inss_stage_4:strata(time_group) + inss_stage_4S + under_18_months:strata(time_group) + risk_score_no_interaction:strata(time_group),
      data=train,
      na.action=na.pass
    )
    
  }
  else if (risk_score == "risk_score_interaction") {
    
    model <- coxph(
      train_survival ~ inss_stage_2 + inss_stage_3 + inss_stage_4:strata(time_group) + inss_stage_4S + under_18_months:strata(time_group) + risk_score_interaction:strata(time_group),
      data=train,
      na.action=na.pass
    )
  }
  return(model)
}


add_five_year_cutoff <- function(results) {
  
  results$train$event_free_survival_5y <- results$train$event_free_survival
  results$train$event_free_survival_days_5y <- results$train$event_free_survival_days
  
  results$train$event_free_survival_5y[
    (results$train$event_free_survival_days_5y > 365*5) & (results$train$event_free_survival_5y == 1)
  ] <- 0
  results$train$event_free_survival_days_5y[results$train$event_free_survival_days_5y > 365*5] <- 365*5
  
  for (i in seq(length(results$validation))) {
    results$validation[[i]]$event_free_survival_5y <- results$validation[[i]]$event_free_survival
    results$validation[[i]]$event_free_survival_days_5y <- results$validation[[i]]$event_free_survival_days
    
    results$validation[[i]]$event_free_survival_5y[
      (results$validation[[i]]$event_free_survival_days_5y > 365*5) & (results$validation[[i]]$event_free_survival_5y == 1)
    ] <- 0
    results$validation[[i]]$event_free_survival_days_5y[results$validation[[i]]$event_free_survival_days_5y > 365*5] <- 365*5
  }
  return(results)
}


compare_risk_scores <- function(data, data_names) {
  
  results <- list()
  
  for (risk_score in RISK_SCORES) {
    
    print(risk_score)
    
    auc_scores <- c()
    
    for (ii in seq(length(data))) {

      precrec_obj <- evalmod(
        scores=data[[ii]][risk_score],
        labels=data[[ii]]$event_free_survival_5y,
      )
      auc_score <- round(precrec::auc(precrec_obj)$auc[1], 4)
      auc_scores <- append(auc_scores, auc_score)
    }
    results[[risk_score]] <- auc_scores
  }
  
  results <- data.frame(results)
  rownames(results) <- data_names
  
  pheatmap(
    results,
    cluster_rows=F,
    cluster_cols=F,
    display_numbers=T
  )
}


kaplan_meier_plot <- function() {
  
  results$validation[[1]]$event_free_survival <- as.numeric(results$validation[[1]]$event_free_survival)
  survival_object <- Surv(
    results$validation[[1]]$event_free_survival_days, results$validation[[1]]$event_free_survival
  )

  # Kaplan Meier Train
  fit <- survfit(
    survival_object ~ prognostic_score_interaction_low,
    data=results$validation[[1]],
    na.action=na.pass
  )
  ggsurvplot(
    fit=fit,
    conf.int=TRUE,
    pval=TRUE,
    risk.table="abs_pct",
    xlab="Years", 
    ylab="Overall survival probability",
    xscale=365, # Divide x-axis into years
    break.x.by=365, # Assume each year is only 365 days
    xlim=c(0, 1826) # Limit to five years
  )
  
  # Other possible validation and diagnostic tools
  # summary(fit_train)$table
  # ggsurvplot(fit_train, fun="cloglog")
  # ggcoxdiagnostics(model, type="schoenfeld")
}


compare_auc_values <- function(data, samples) {
  
  auc_scores <- vector()
  auc_scores_baseline <- vector()
  
  for (i in 1:samples) {
    df_sample <- data[sample(nrow(data), length(data), replace=TRUE), ]

    precrec_obj <- evalmod(
      scores=df_sample$prognostic_score_interaction,
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


evaluate_validation_datasets_combinations <- function(results) {
  
  # Identify the matching columns between datasets so that
  # they can be concatenated
  required_columns <- c(RISK_SCORES, "event_free_survival_5y")
  
  # Add the different combinations of validation datasets to a list
  validation_combos <- list()
  validation_combos[[1]] <- results$validation[[1]][required_columns]
  validation_combos[[2]] <- results$validation[[2]][required_columns]
  validation_combos[[3]] <- results$validation[[3]][required_columns]
  validation_combos[[4]] <- rbind(validation_combos[[1]], validation_combos[[2]])
  validation_combos[[5]] <- rbind(validation_combos[[1]], validation_combos[[3]])
  validation_combos[[6]] <- rbind(validation_combos[[2]], validation_combos[[3]])
  validation_combos[[7]] <- rbind(validation_combos[[1]], validation_combos[[2]], validation_combos[[3]])
  
  # Create a heatmap to compare the AUC scores across the
  # combination of datasets and risk scores
  compare_risk_scores(validation_combos, VALIDATION_COMBOS_NAMES)
  
  return(validation_combos)
}


plot_time_dependent_roc_curves <- function(data, risk_score) {
  
  roc_object_1 <- survivalROC(data$event_free_survival_days, 
                              data$event_free_survival, 
                              data[[risk_score]],
                              predict.time = 365,
                              method = 'KM')
  roc_object_3 <- survivalROC(data$event_free_survival_days, 
                              data$event_free_survival, 
                              data[[risk_score]],
                              predict.time = 365*3,
                              method = 'KM')
  roc_object_5 <- survivalROC(data$event_free_survival_days, 
                              data$event_free_survival, 
                              data[[risk_score]],
                              predict.time = 365*5,
                              method = 'KM')
  
  df1 <- data.frame(fp=roc_object_1$FP, tp=roc_object_1$TP)
  df3 <- data.frame(fp=roc_object_3$FP, tp=roc_object_3$TP)
  df5 <- data.frame(fp=roc_object_5$FP, tp=roc_object_5$TP)
  
  auc_label1 <- sprintf('AUC (1Y) = %.4f', roc_object_1$AUC)
  auc_label3 <- sprintf('AUC (3Y) = %.4f', roc_object_3$AUC)
  auc_label5 <- sprintf('AUC (5Y) = %.4f', roc_object_5$AUC)
  
  (ggplot()
    + geom_line(data = df1, aes(fp, tp, color="AUC (1Y)"))
    + geom_line(data=df3, aes(fp, tp, color="AUC (3Y)"))
    + geom_line(data=df5, aes(fp, tp, color="AUC (5Y)"))
    + geom_abline(intercept = 0, slope = 1)
    + theme_light()
    + annotate('text', x = 0.7, y = 0.4, label = auc_label1, size = 4)
    + annotate('text', x = 0.7, y = 0.3, label = auc_label3, size = 4)
    + annotate('text', x = 0.7, y = 0.2, label = auc_label5, size = 4)
    + scale_colour_manual("Legend", 
                          breaks = c("AUC (1Y)", "AUC (3Y)", "AUC (5Y)"),
                          values = c("AUC (1Y)"="dark green", "AUC (3Y)"="red", "AUC (5Y)"="blue"))
    + ylab('Sensitivity')
    + xlab('1-Specificity')
    + ggtitle(sprintf('Event-free Survival at years 1, 3, 5')))
}


measure_auc_p_values <- function(data) {
  
  prognostic_scores <- c(
    "prognostic_score_no_interaction",
    "prognostic_score_interaction"
  )
  
  results <- list()
  
  for (score in prognostic_scores) {
    
    p_values <- c()
    
    for (i in seq(length(data))) {
      roc_model <- roc(data[[i]]$event_free_survival_5y, data[[i]][[score]])
      roc_baseline <- roc(data[[i]]$event_free_survival_5y, data[[i]]$risk_score_baseline)
      
      result <- roc.test(roc_model, roc_baseline)
      p_value <- round(result$p.value, 4)
      p_values <- append(p_values, p_value)
    }
    results[[score]] <- p_values
  }
  
  results <- data.frame(results)
  rownames(results) <- VALIDATION_COMBOS_NAMES
  
  pheatmap(
    results,
    cluster_rows=F,
    cluster_cols=F,
    display_numbers=T,
    legend=F,
    main="P-values of AUC scores between prognostic models and baseline model",
    breaks=c(
      seq(0, 0.05, length.out=ceiling(10)),
      seq(0.1, 1, length.out=ceiling(10))
    ),
    color=colorRampPalette(c("light blue", "white", "orange"))(20)
  )
}


# Prepare the data for training
train <- load_training_data()
train_survival <- survival_data_train(train)

validation <- load_validation_data()
validation_survival <- survival_data_validation(validation)

# Get the genes from the training dataset
genes <- colnames(train)[14:length(colnames(train)) - 1]

# Identify the genes for modelling
genes_modelling_no_interaction <- modelling_genes(genes, train, train_survival, FALSE)
genes_modelling_interaction <- modelling_genes(genes, train, train_survival, TRUE)

# Split the data to help keep the proportional hazards consistent over time
train <- split_data_train(train, train_survival)
validation <- split_data_validation(validation, validation_survival)

# Create models
model_baseline <- baseline_model(train, train_survival)
model_no_interaction <- no_interaction_model(train, train_survival)
model_interaction <- interaction_model(train, train_survival)

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
results <- predict_risk_score(model_baseline, train, validation, "risk_score_baseline")
results <- predict_risk_score(model_no_interaction, results$train, results$validation, "risk_score_no_interaction")
results <- predict_risk_score(model_interaction, results$train, results$validation, "risk_score_interaction")

# Update data_split_train with the risk scores
train <- results$train
validation <- results$validation

# Train the prognostic models
model_prognostic_no_interaction <- prognostic_model(train, train_survival, "risk_score_no_interaction")
model_prognostic_interaction <- prognostic_model(train, train_survival, "risk_score_interaction")

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
results <- predict_risk_score(model_prognostic_no_interaction, train, validation, "prognostic_score_no_interaction")
results <- predict_risk_score(model_prognostic_interaction, results$train, results$validation, "prognostic_score_interaction")

# Truncate the survival data to five years for measuring the performance
# of the models
results <- add_five_year_cutoff(results)

# Combine train and validation into one list
data <- list()
data[[1]] <- results$train
data <- append(data, results$validation)

# Create a heatmap comparing the risk and prognostic scores
data_names <-  c("train", "GSE85047", "target", "E_TABM-38", "wilms_tumor")
compare_risk_scores(data, data_names)

# Can't use input parameters because of limitations in survfit function
kaplan_meier_plot()

# Plot gene expression heatmap
pheatmap(
  train[order(train$risk_score_interaction),][colnames(train) %in% genes_modelling_interaction],
  show_rownames=FALSE,
  cluster_rows=F,
  cluster_cols=F,
  breaks=seq(-2, 3, by = 0.06),
)

# Plot hazard scores from models
plot_ggforest(model_prognostic_interaction, data=train)

# Check number of high and low prognostic risk scores
table(results$train$prognostic_score_interaction_low)
table(results$validation[[1]]$prognostic_score_interaction_low)
table(results$validation[[2]]$prognostic_score_interaction_low)
table(results$validation[[3]]$prognostic_score_interaction_low)

# Compare event_free_survival_days between risk groups
ggplot(results$train, aes(x=prognostic_score_interaction_low, y=event_free_survival_days)) + geom_boxplot()
ggplot(results$validation[[1]], aes(x=prognostic_score_interaction_low, y=event_free_survival_days)) + geom_boxplot()
ggplot(results$validation[[2]], aes(x=prognostic_score_interaction_low, y=event_free_survival_days)) + geom_boxplot()
ggplot(results$validation[[3]], aes(x=prognostic_score_interaction_low, y=event_free_survival_days)) + geom_boxplot()
ggplot(results$validation[[4]], aes(x=prognostic_score_interaction_low, y=event_free_survival_days)) + geom_boxplot()

# Create a heatmap comparing the AUC scores using different combinations
# of validation datasets and risk scores
validation_combos <- evaluate_validation_datasets_combinations(results)

# Plot the 1, 3, and 5 year ROC curves for a dataset and risk score
plot_time_dependent_roc_curves(results$validation[[1]], "prognostic_score_interaction")

# Measure the p-values for the ROC curves of the prognostic and baseline models
# across the combinations of validation datasets
measure_auc_p_values(validation_combos)
