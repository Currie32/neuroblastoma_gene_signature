source("~/Imperial/neuroblastoma_gene_signature/validation/plot_ggforest.R")

library(grid)
library(ggplot2)
library(MuMIn)
library(pheatmap)
library(precrec)
library(pROC)
library(randomForestSRC)
library(rms)
library(stringr)
library(survival)
library(survivalROC)
library(survminer)
library(VennDiagram)


PATH <- "~/Imperial/neuroblastoma_gene_signature/data/modelling"
P_VALUE = 0.005

# The day to split the training and testing data so that
# the proportional hazards are consistent over time. More info: 
# https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
CUT_DAY <- c(410)

# Names of the predicted risk scores
RISK_SCORES <- c(
  "risk_score_baseline",
  "risk_score_no_interaction",
  "risk_score_interaction",
  "risk_score_rf",
  "prognostic_score_no_interaction",
  "prognostic_score_interaction",
  "prognostic_score_rf"
)

# Names of the datasets in each combination
VALIDATION_COMBOS_NAMES <- c(
  "GSE85047", "target", "E-TABM-38", "wilms_tumor",
  "GSE85047_target", "GSE85047_E-TABM-38", "target_E-TABM-38",
  "GSE85047_target_E-TABM-38"
)


load_training_data <- function() {
  
  train <- read.csv(file.path(PATH, "train.csv"))
  
  return(train)
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


add_five_year_cutoff <- function(train, validation) {
  #' Add the five-year survival data to the training and validation datasets.
  #' 
  #' train data.frame: the training dataset
  #' validation data.frame: the validation datasets
  #' 
  #' return List(data.frame): the training and validation datasets with
  #'                          five-year survival results
  
  # Get the survival features
  train$event_free_survival_5y <- train$event_free_survival
  train$event_free_survival_days_5y <- train$event_free_survival_days
  
  # Update their values at the five-year mark
  train$event_free_survival_5y[
    (train$event_free_survival_days_5y > 365*5) & (train$event_free_survival_5y == 1)
  ] <- 0
  train$event_free_survival_days_5y[
    train$event_free_survival_days_5y > 365*5
  ] <- 365*5
  
  # Do the same thing as above, but on the validation datasets
  for (i in seq(length(validation))) {
    validation[[i]]$event_free_survival_5y <- validation[[i]]$event_free_survival
    validation[[i]]$event_free_survival_days_5y <- validation[[i]]$event_free_survival_days
    
    validation[[i]]$event_free_survival_5y[
      (validation[[i]]$event_free_survival_days_5y > 365*5) & (validation[[i]]$event_free_survival_5y == 1)
    ] <- 0
    validation[[i]]$event_free_survival_days_5y[validation[[i]]$event_free_survival_days_5y > 365*5] <- 365*5
  }
  
  # Add the training and validation datasets to a list of results
  results <- list()
  results$train <- train
  results$validation <- validation
  
  return(results)
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


modelling_genes <- function(genes, train, train_survival, interaction) {
  #' Perform multiple cox regressions and AIC feature selection to determine
  #' the optimal set of genes for predicting the risk of neuroblastoma
  #' 
  #' genes List(str): names of the genes to analyse
  #' train data.frame: contains the phenotype and expression data
  #' train_survival double: the survival object for the training dataset
  #' interaction bool: TRUE if the interaction step should be used
  #' 
  #' return List(str): names of the most predictive genes
  
  genes_candidate <- univariate_cox_regression(genes, train, train_survival)
  genes_significant <- multivariate_cox_regression(genes_candidate, train_survival, train, FALSE, FALSE)
  
  if (interaction) {
    genes_significant <- multivariate_cox_regression(genes_significant, train_survival, train, FALSE, TRUE)
  }
  
  genes_best <- aic_feature_selection(genes_significant, train_survival, train, '')
  
  return(genes_best)
}


univariate_cox_regression <- function(genes, train, train_survival) {
  #' Find the genes that are statistically significant in a univariate cox regression.
  #'
  #' genes List(str): Names of genes to analyse
  #' train data.frame: Contains the gene expression data
  #' train_survival double: The survival object for the training data
  #'
  #' return List(str): Names of statistically significant genes
  
  # Store the significant genes
  candidate_genes <- c()
  
  for (gene in genes) {
    
    # Perform cox regression
    result <- coxph(
      as.formula(paste("train_survival ~ ", gene)),
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
  #' significant using a multivariate cox regression
  #' 
  #' features List(str): genes that were significant during the
  #'                     univariate cox regression
  #' train_survival double: The survival object for the training data
  #' train data.frame: contains the expression and survival data for each patient
  #' use_phenotype_features bool: TRUE if the specificed phenotype features
  #'                              should be used
  #' use_interaction bool: TRUE if an interaction term for each pair of features
  #'                       should be included
  #' 
  #' returns List(str): the statistically significant features
  
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
  #' features List(str): the statistically significant features from a
  #'                     multivariate cox regression
  #' train_survival double: survival information for each patient
  #' train data.frame: contains the expression and survival data for each patient
  #' delimiter str: the character to split the index of the best features
  #' 
  #' returns List(str): the best combination of features according to AIC 
  
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
  #' Build a baseline model to measure the predictive performance of only
  #' phenotype features
  #' 
  #' train data.frame: contains the phenotype features to analysis
  #' train_survival double: the survival object for the training data
  #' 
  #' return List: the trained baseline model
  
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

  # # Fit top model
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

  # Train the baseline model
  # Below are the features_joined and have been stratified across the time groups
  model_baseline <- coxph(
    train_survival ~ mycn_amplification + under_18_months:strata(time_group) + inss_stage_2:strata(time_group) + inss_stage_3:strata(time_group) + inss_stage_4:strata(time_group) + inss_stage_4S:strata(time_group),
    data=train,
    na.action=na.pass
  )
  
  return(model_baseline)
}


no_interaction_model <- function(train, train_survival) {
  #' Train a cox proportional-hazard model, which does not include any
  #' interacting features
  #' 
  #' train data.frame: contains the expression features to model
  #' train_survival double: the survival object for the training data
  #' 
  #' return List: the trained model
  
  # The features below are genes_modelling_no_interaction, with the
  # necessary stratification
  model <- coxph(
    train_survival ~ CD147:strata(time_group) + HEBP2:strata(time_group) + HSD17B12 + NXT2:strata(time_group) + RACK1 + TXNDC5:strata(time_group),
    data=train,
    na.action=na.pass
  )
  
  return(model)
}


interaction_model <- function(train, train_survival) {
  #' Train a cox proportional-hazard model, using interacting features
  #' 
  #' train data.frame: contains the expression features to model
  #' train_survival double: the survival object for the training data
  #' 
  #' return List: the trained model
  
  # The features below are genes_modelling_interaction, with the
  # necessary stratification
  model <- coxph(
    train_survival ~ HEBP2:strata(time_group) + HSD17B12 + NXT2:strata(time_group) + TXNDC5:strata(time_group) + NXT2:TXNDC5,
    data=train,
    na.action=na.pass
  )
  
  return(model)
}


rf_model <- function(train, genes_rf){
  #' Train randfom forest survival model
  #' 
  #' train data.frame: contains the expression features to model
  #' genes_rf List(str): names of the genes to train the random forest model
  #' 
  #' return List: the trained model
  
  # Subset to GSE49711
  train_rf <- train[1:498, ]
  set.seed(1)

  # Join the random forest genes into a single string
  genes_rf_joined <- paste(genes_rf, collapse=" + ")
  
  # Train the model
  model <- rfsrc(
    as.formula(paste("Surv(event_free_survival_days_5y, event_free_survival_5y) ~", genes_rf_joined)),
    data=train_rf, 
    ntree=200, 
    mtry=2,
    nodesize=3,
    nsplit=1,
    importance=TRUE
  )

  return(model)
}


differentiation_model <- function(train, train_survival) {
  #' Train a cox proportional-hazard model using the genes that affect the
  #' differentiation pathway
  #' 
  #' train data.frame: contains the expression features to model
  #' train_survival double: the survival object for the training data
  #' 
  #' return List: the trained model
  
  model <- coxph(
    train_survival ~ B3GAT1 + CD147 + FNBP1:strata(time_group) + IGSF10 + RACK1 + TXNDC5:strata(time_group),
    data=train,
    na.action=na.pass
  )
  return(model)
}


tumourigenesis_model <- function(train, train_survival) {
  #' Train a cox proportional-hazard model using the genes that affect the
  #' tumourigenesis pathway
  #' 
  #' train data.frame: contains the expression features to model
  #' train_survival double: the survival object for the training data
  #' 
  #' return List: the trained model
  
  model <- coxph(
    train_survival ~ HEBP2:strata(time_group) + IGSF10 + IQCE + RACK1 + TXNDC5:strata(time_group),
    data=train,
    na.action=na.pass
  )
  return(model)
}


metastasis_model <- function(train, train_survival) {
  #' Train a cox proportional-hazard model using the genes that affect the
  #' metastasis pathway
  #' 
  #' train data.frame: contains the expression features to model
  #' train_survival double: the survival object for the training data
  #' 
  #' return List: the trained model
  
  model <- coxph(
    train_survival ~ HAPLN4:strata(time_group) + CCDC125 + CD147 + CNKSR3 +  FNBP1:strata(time_group) + RACK1,
    data=train,
    na.action=na.pass
  )
  return(model)
}


predict_risk_score_lm <- function(model, train, validation, risk_score) {
  #' Predict a risk score, using a linear model, on the training and validation datasets
  #' 
  #' model List: the model to use to make the predictions
  #' train data.frame: the training dataset
  #' validation List(data.frame): the validation datasets
  #' risk_score str: name of the risk score to predict
  #' 
  #' return List(data.frame): the training and validation datasets with their
  #'                          newly predicted risk scores
  
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
  
  # Add the training and validation datasets to a list of results
  results <- list()
  results$train <- train
  results$validation <- validation
  
  return(results)
}


predict_risk_score_rf <- function(model_rf, train, validation, risk_score, genes_rf) {
  #' Predict a risk score, using a random forest model, on the training
  #' and validation datasets
  #' 
  #' model_rf List: the model to use to make the predictions
  #' train data.frame: the training dataset
  #' validation List(data.frame): the validation datasets
  #' risk_score str: name of the risk score to predict
  #' genes_rf List(str): names of the genes used to train the random forest model
  #' 
  #' return List(data.frame): the training and validation datasets with their
  #'                          newly predicted risk scores
  
  # Predict the risk scores for the training dataframe
  train[risk_score] <- randomForestSRC::predict.rfsrc(model_rf, train[genes_rf])$predicted
  train[paste(risk_score, "low", sep="_")] <- train[risk_score] < median(train[risk_score][[1]])
  
  # Calculate the median train risk score to be used as the high/low threshold
  median_train_risk_score <- median(train[risk_score][[1]])
  
  # Add risk scores to the remaining dataframes
  for (i in seq(length(results$validation))) {
    validation[[i]][risk_score] <- randomForestSRC::predict.rfsrc(model_rf, validation[[i]][genes_rf])$predicted
    validation[[i]][paste(risk_score, "low", sep="_")] <- validation[[i]][risk_score] < median_train_risk_score
  }
  
  # Add the training and validation datasets to a list of results
  results <- list()
  results$train <- train
  results$validation <- validation
  
  return(results)
}

scale_risk_scores <- function(train) {
  #' Z-transform the risk scores that will be used to train another model
  #' 
  #' train data.frame: the training dataframe that is used to train models
  #' 
  #' return data.frame: the training dataframe with scaled risk scores
  
  train$risk_score_differentiation <- c(scale(train$risk_score_differentiation))
  train$risk_score_immune_disregulation <- c(scale(train$risk_score_immune_disregulation))
  train$risk_score_metabolism <- c(scale(train$risk_score_metabolism))
  train$risk_score_metastasis <- c(scale(train$risk_score_metastasis))
  train$risk_score_nuclear_transport <- c(scale(train$risk_score_nuclear_transport))
  train$risk_score_tumourigenesis <- c(scale(train$risk_score_tumourigenesis))
  train$risk_score_rf <- c(scale(train$risk_score_rf))
  
  return(train)
}


combined_functional_model <- function(train, train_survival) {
  #' Train a cox proportional-hazard model using the risk scores from the
  #' functional models to create a combined functional model
  #' 
  #' train data.frame: contains the risk scores to model
  #' train_survival double: the survival object for the training data
  #' 
  #' return List: the trained model
  
  model_combine <- coxph(
    train_survival ~ risk_score_differentiation + risk_score_immune_disregulation + risk_score_metabolism + risk_score_metastasis + risk_score_nuclear_transport + risk_score_tumourigenesis,
    data=train,
    na.action=na.pass
  )
  
  return(model_combine)
}


prognostic_model <- function(train, train_survival, risk_score) {
  #' Train a prognostic model
  #' 
  #' train data.frame: the data to train the prognostic model
  #' train_survival double: the survival object for the training data
  #' risk_score str: the name of the risk score to use to train the prognostic model
  #' 
  #' return List: the trained prognostic model
  
  # Identify the significant features for training the model
  features_significant <- multivariate_cox_regression(risk_score, train_survival, train, TRUE, FALSE)
  features_best <- aic_feature_selection(features_significant, train_survival, train, '')
  
  if (risk_score == "risk_score_no_interaction") {

    model <- coxph(
      train_survival ~ mycn_amplification + inss_stage_2 + inss_stage_3 + inss_stage_4:strata(time_group) + inss_stage_4S + under_18_months:strata(time_group) + risk_score_no_interaction:strata(time_group),
      data=train,
      na.action=na.pass
    )
  }
  else if (risk_score == "risk_score_interaction") {
    
    model <- coxph(
      train_survival ~ mycn_amplification + inss_stage_2 + inss_stage_3 + inss_stage_4:strata(time_group) + inss_stage_4S + under_18_months:strata(time_group) + risk_score_interaction:strata(time_group),
      data=train,
      na.action=na.pass
    )
  }
  else if (risk_score == "risk_score_rf") {
    
    model <- coxph(
      train_survival ~ mycn_amplification + inss_stage_2 + inss_stage_3 + inss_stage_4:strata(time_group) + inss_stage_4S + under_18_months:strata(time_group) + risk_score_rf:strata(time_group),
      data=train,
      na.action=na.pass
    )
  }
  return(model)
}


compare_risk_scores <- function(data, data_names, title, x_label, y_label) {
  #' Compare the performance of a set of risk scores on a list of datasets
  #' using the AUC metric
  #' 
  #' data List(data.frame): datasets on which to measure the performance
  #' data_names List(str): names of the datasets
  #' title str: title for the heatmap
  #' x_label str: x-axis label
  #' y_label str: y-axis label
  
  results <- list()
  
  for (risk_score in RISK_SCORES) {

    auc_scores <- c()
    
    for (i in seq(length(data))) {
      
      # Set random forest AUC values to 0 when predicting on the E-TABM-38
      # dataset since it is missing required genes
      if (
        risk_score %in% c("risk_score_rf", "prognostic_score_rf")
        & (grepl("E-TABM-38", data_names[i], fixed = TRUE))
      ) {
        auc_scores <- append(auc_scores, 0)
      }
      # The RF risk score is only trained on one of the training datasets
      else if ((risk_score == "risk_score_rf") & (i == 1)) {
        roc_model <- roc(data[[i]][1:986, ]$event_free_survival_5y, data[[i]][1:986, ][[risk_score]])
        auc_score <- round(roc_model$auc[1], 4)
        auc_scores <- append(auc_scores, auc_score)
      }
      else {
        roc_model <- roc(data[[i]]$event_free_survival_5y, data[[i]][[risk_score]])
        auc_score <- round(roc_model$auc[1], 4)
        auc_scores <- append(auc_scores, auc_score)
      }
    }
    results[[risk_score]] <- auc_scores
  }
  
  results <- data.frame(results)
  rownames(results) <- data_names
  
  setHook(
    "grid.newpage",
    function() pushViewport(viewport(
      x=0.92, y=0.9, width=0.9, height=0.8, name="vp", just=c("right","top"))
    ),
    action="prepend"
  )
  pheatmap(
    results,
    cluster_rows=F,
    cluster_cols=F,
    display_numbers=T,
    angle_col=315,
    legend=FALSE,
    fontsize=14
  )
  setHook("grid.newpage", NULL, "replace")
  grid.text(title, y=1.05, gp=gpar(fontsize=18))
  grid.text(x_label, y=-0.02, gp=gpar(fontsize=16))
  grid.text(y_label, x=1.04, y=0.7, rot=90, gp=gpar(fontsize=16))
}


kaplan_meier_plot <- function() {
  #' Create a Kaplan Meier plot to measure the performance of a risk score
  #' on a dataset
  
  # Kaplan Meier Train
  fit <- survfit(
    Surv(
      event_free_survival_days_5y, event_free_survival_5y
    ) ~ prognostic_score_interaction_low,
    data=results$validation[[4]],
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


plot_gene_expression <- function(data, risk_score, genes, title, x_label, y_label) {
  #' Plot the Z-transformed gene expression data ranked by a risk score.
  #' 
  #' data data.frame: contains the expression data and risk score
  #' risk_score str: name of the risk score by which to rank
  #' genes List(str): names of the genes to plot
  #' title str: title of the heatmap
  #' x_label str: x-axis label of the heatmap
  #' y_label str: y-axis label of the heatmap
  
  # Create a "grid space" for the heatmap to be plotted within
  setHook(
    "grid.newpage",
    function() pushViewport(viewport(
      x=0.98, y=0.9, width=0.9, height=0.8, name="vp", just=c("right","top"))
    ),
    action="prepend"
  )
  # Plot the heatmap
  pheatmap(
    data[order(data[risk_score]), genes],
    show_rownames=FALSE,
    cluster_rows=F,
    cluster_cols=F,
    breaks=seq(-2.5, 2.5, by = 0.07),
    fontsize=12
  )
  # Add the titles and labels
  setHook("grid.newpage", NULL, "replace")
  grid.text(title, y=1.05, gp=gpar(fontsize=18))
  grid.text(x_label, y=-0.04, gp=gpar(fontsize=16))
  grid.text(y_label, x=-0.04, rot=90, gp=gpar(fontsize=16))
}


evaluate_validation_datasets_combinations <- function(results, title, x_label, y_label) {
  #' Create a heatmap comparing the AUC scores using different combinations
  #' of validation datasets and risk scores
  #' 
  #' results List(data.frame): contains the validation datasets
  #' #' title str: title for the heatmap
  #' x_label str: x-axis label
  #' y_label str: y-axis label
  #' 
  #' return List(data.frame): different combinations of the validation datasets
  
  # Identify the matching columns between datasets so that
  # they can be concatenated
  required_columns <- c(RISK_SCORES, "event_free_survival_5y", "event_free_survival_days_5y")
  
  # Add the different combinations of validation datasets to a list
  validation_combos <- list()
  validation_combos[[1]] <- results$validation[[1]][required_columns]
  validation_combos[[2]] <- results$validation[[2]][required_columns]
  validation_combos[[3]] <- results$validation[[3]][required_columns]
  validation_combos[[4]] <- results$validation[[4]][required_columns]
  validation_combos[[5]] <- rbind(validation_combos[[1]], validation_combos[[2]])
  validation_combos[[6]] <- rbind(validation_combos[[1]], validation_combos[[3]])
  validation_combos[[7]] <- rbind(validation_combos[[2]], validation_combos[[3]])
  validation_combos[[8]] <- rbind(validation_combos[[1]], validation_combos[[2]], validation_combos[[3]])
  
  # Create a heatmap to compare the AUC scores across the
  # combination of datasets and risk scores
  compare_risk_scores(validation_combos, VALIDATION_COMBOS_NAMES, title, x_label, y_label)
  
  return(validation_combos)
}


plot_time_dependent_roc_curves <- function(data, risk_score) {
  #' Plot the 1, 3, and 5 year ROC curves for a dataset and risk score
  #' 
  #' data data.frame: contains the survival data and risk score
  #' risk_score str: the risk score to be measured
  
  # Create the 1, 3, and 5 year ROC objects
  roc_object_1 <- survivalROC(data$event_free_survival_days_5y, 
                              data$event_free_survival_5y, 
                              data[[risk_score]],
                              predict.time = 365,
                              method = 'KM')
  roc_object_3 <- survivalROC(data$event_free_survival_days_5y, 
                              data$event_free_survival_5y, 
                              data[[risk_score]],
                              predict.time = 365*3,
                              method = 'KM')
  roc_object_5 <- survivalROC(data$event_free_survival_days_5y, 
                              data$event_free_survival_5y, 
                              data[[risk_score]],
                              predict.time = 365*5,
                              method = 'KM')
  
  # Move the false positive and true positive data into dataframes
  df1 <- data.frame(fp=roc_object_1$FP, tp=roc_object_1$TP)
  df3 <- data.frame(fp=roc_object_3$FP, tp=roc_object_3$TP)
  df5 <- data.frame(fp=roc_object_5$FP, tp=roc_object_5$TP)
  
  # Create annotations for the plot
  auc_label1 <- sprintf('AUC (1Y) = %.4f', roc_object_1$AUC)
  auc_label3 <- sprintf('AUC (3Y) = %.4f', roc_object_3$AUC)
  auc_label5 <- sprintf('AUC (5Y) = %.4f', roc_object_5$AUC)
  
  # Create the plot
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
  #' Test if the difference in AUC values between the models' prognostic scores
  #' and the baseline model's risk score is statistically significant
  #' 
  #' data List(data.frame): the datasets on which to measure the difference in AUC
  
  # The prognostic scores to compare against the baseline risk score
  prognostic_scores <- c(
    "prognostic_score_no_interaction",
    "prognostic_score_interaction",
    "prognostic_score_rf"
  )
  
  results <- list()
  
  for (score in prognostic_scores) {
    
    p_values <- c()
    
    for (i in seq(length(data))) {
      
      if ((score == "prognostic_score_rf") & (i %in% c(3, 6, 7, 8))) {
        p_values <- append(p_values, 1)
      }
      else {
        roc_model <- roc(data[[i]]$event_free_survival_5y, data[[i]][[score]])
        roc_baseline <- roc(data[[i]]$event_free_survival_5y, data[[i]]$risk_score_baseline)
        
        result <- roc.test(roc_model, roc_baseline, alternative="greater")
        p_value <- round(result$p.value, 4)
        p_values <- append(p_values, p_value)
      }
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
    fontsize=15,
    angle_col=315,
    main="P-values of AUC scores between prognostic models and baseline model",
    breaks=c(
      seq(0, 0.05, length.out=ceiling(10)),
      seq(0.1, 1, length.out=ceiling(10))
    ),
    color=colorRampPalette(c("light blue", "white", "orange"))(20)
  )
}


plot_overlapping_genes <- function() {
  #' Plot a venn diagram to show the overlapping genes between
  #' the three models. Use global variable since this function is only used once.
  
  # Create Venn Diagram 
  venn_diagram <- venn.diagram(
    x = list(
      genes_modelling_no_interaction,
      genes_modelling_interaction,
      genes_rf
    ),
    col=c("#440154ff", '#21908dff', '#fde725ff'),
    category.names=c("No interaction genes" , "Interaction genes", "Random Forest genes"),
    fill=c(alpha("#440154ff", 0.3), alpha('#21908dff', 0.3), alpha('#fde725ff', 0.3)),
    main="The overlapping genes between the three models",
    filename = NULL,
    cex=1.2,
    cat.cex=1.1,
    main.cex=1.3
  )
  
  # Plot it
  grid.newpage()
  grid.draw(venn_diagram)
}


# Load the training and validation data
train <- load_training_data()
validation <- load_validation_data()

# Get the training and validation survival data
train_survival <- survival_data_train(train)
validation_survival <- survival_data_validation(validation)

# Truncate the survival data to five years for measuring the performance
# of the models
results <- add_five_year_cutoff(train, validation)
train <- results$train
validation <- results$validation

# Get the genes from the training dataset
genes <- colnames(train)[16:length(colnames(train)) - 3]

# Identify the genes for modelling
genes_modelling_no_interaction <- modelling_genes(genes, train, train_survival, FALSE)
genes_modelling_interaction <- modelling_genes(genes, train, train_survival, TRUE)
genes_rf <- c(
"B3GAT1", "CCDC125", "FNBP1", "HAPLN4", "HEBP2", "HSD17B12", "IGSF10", "IQCE", "KCNQ3", "TOX2"
)

# Train the random forest model with the split data
model_rf <- rf_model(train, genes_rf)

# Split the data to help keep the proportional hazards consistent over time
train <- split_data_train(train, train_survival)
validation <- split_data_validation(validation, validation_survival)

# Create models
model_baseline <- baseline_model(train, train_survival)
model_no_interaction <- no_interaction_model(train, train_survival)
model_interaction <- interaction_model(train, train_survival)

## Functional models
model_differentiation <- differentiation_model(train, train_survival)
model_immune_disregulation <- coxph(train_survival ~ TOX2, data=train, na.action=na.pass)
model_metabolism <- coxph(train_survival ~ HSD17B12, data=train, na.action=na.pass)
model_metastasis <- metastasis_model(train, train_survival)
model_nuclear_transport <- coxph(train_survival ~ NXT2:strata(time_group), data=train, na.action=na.pass)
model_tumourigenesis <- tumourigenesis_model(train, train_survival)

# Check that all models meet the assumptions of the Cox proportional hazards model
zph_model_baseline <- cox.zph(model_baseline)
zph_model_no_interaction <- cox.zph(model_no_interaction)
zph_model_interaction <- cox.zph(model_interaction)
zph_model_differentiation <- cox.zph(model_differentiation)
zph_model_immune_disregulation <- cox.zph(model_immune_disregulation)
zph_model_metabolism <- cox.zph(model_metabolism)
zph_model_metastasis <- cox.zph(model_metastasis)
zph_model_nuclear_transport <- cox.zph(model_nuclear_transport)
zph_model_tumourigenesis <- cox.zph(model_tumourigenesis)

# All values should be below 0.05
print(zph_model_baseline)
print(zph_model_no_interaction)
print(zph_model_interaction)
print(zph_model_differentiation)
print(zph_model_immune_disregulation)
print(zph_model_metabolism)
print(zph_model_metastasis)
print(zph_model_nuclear_transport)
print(zph_model_tumourigenesis)

# Plot the diagnostic plots to confirm there is no trend
ggcoxzph(zph_model_baseline)
ggcoxzph(zph_model_no_interaction)
ggcoxzph(zph_model_interaction)
ggcoxzph(zph_model_tumourigenesis)
ggcoxzph(zph_model_metastasis)
ggcoxzph(zph_model_differentiation)
ggcoxzph(zph_model_immune_disregulation)
ggcoxzph(zph_model_metabolism)
ggcoxzph(zph_model_nuclear_transport)

# Add risk scores
results <- predict_risk_score_lm(model_baseline, train, validation, "risk_score_baseline")
results <- predict_risk_score_lm(model_no_interaction, results$train, results$validation, "risk_score_no_interaction")
results <- predict_risk_score_lm(model_interaction, results$train, results$validation, "risk_score_interaction")
results <- predict_risk_score_lm(model_differentiation, results$train, results$validation, "risk_score_differentiation")
results <- predict_risk_score_lm(model_immune_disregulation, results$train, results$validation, "risk_score_immune_disregulation")
results <- predict_risk_score_lm(model_metabolism, results$train, results$validation, "risk_score_metabolism")
results <- predict_risk_score_lm(model_metastasis, results$train, results$validation, "risk_score_metastasis")
results <- predict_risk_score_lm(model_nuclear_transport, results$train, results$validation, "risk_score_nuclear_transport")
results <- predict_risk_score_lm(model_tumourigenesis, results$train, results$validation, "risk_score_tumourigenesis")
results <- predict_risk_score_rf(model_rf, results$train, results$validation, "risk_score_rf", genes_rf)

# Update data_split_train with the risk scores
train <- results$train
validation <- results$validation

# Scale the risk scores
# train <- scale_risk_scores(train)

# Train the combined functional model
model_combine <- combined_functional_model(train, train_survival)

# Train the prognostic models
model_prognostic_no_interaction <- prognostic_model(train, train_survival, "risk_score_no_interaction")
model_prognostic_interaction <- prognostic_model(train, train_survival, "risk_score_interaction")
model_prognostic_rf <- prognostic_model(train, train_survival, "risk_score_rf")

# Check that all models meet the assumptions of the Cox proportional hazards model
zph_model_combine <- cox.zph(model_combine)
zph_model_prognostic_no_interaction <- cox.zph(model_prognostic_no_interaction)
zph_model_prognostic_interaction <- cox.zph(model_prognostic_interaction)
zph_model_prognostic_rf <- cox.zph(model_prognostic_rf)

# All values should be below 0.05
print(zph_model_combine)
print(zph_model_prognostic_no_interaction)
print(zph_model_prognostic_interaction)
print(zph_model_prognostic_rf)

# Plot the diagnostic plots to confirm there is no trend
ggcoxzph(zph_model_combine)
ggcoxzph(zph_model_prognostic_no_interaction)
ggcoxzph(zph_model_prognostic_interaction)
ggcoxzph(zph_model_prognostic_rf)

# Add prognostic scores
results <- predict_risk_score_lm(model_prognostic_no_interaction, train, validation, "prognostic_score_no_interaction")
results <- predict_risk_score_lm(model_prognostic_interaction, results$train, results$validation, "prognostic_score_interaction")
results <- predict_risk_score_lm(model_prognostic_rf, results$train, results$validation, "prognostic_score_rf")

# Combine train and validation into one list
data <- list()
data[[1]] <- results$train
data <- append(data, results$validation)

# Create a heatmap comparing the risk and prognostic scores
data_names <-  c("train", "GSE85047", "target", "E-TABM-38", "wilms_tumor")
compare_risk_scores(data, data_names, 'AUC value of risk scores on datasets', 'Risk score', 'Dataset')

# Can't use input parameters because of limitations in survfit function
kaplan_meier_plot()

# Plot gene expression heatmap
## Low scores at the top
genes_modelling_no_interaction_sorted <- c("CD147", "HEBP2", "NXT2", "HSD17B12", "TXNDC5", "RACK1")
genes_modelling_interaction_sorted <- c("HEBP2", "NXT2", "HSD17B12", "TXNDC5")
genes_modelling_rf_sorted <- c("FNBP1", "TOX2", "KCNQ3", "IQCE", "B3GAT1", "IGSF10", "HEBP2", "CCDC125", "HAPLN4", "HSD17B12")

plot_gene_expression(
  train,
  "risk_score_rf",
  genes_modelling_rf_sorted,
  "Gene expression ranked by random forest risk score",
  "Gene",
  "Expression level (Z-transformed)"
)

# Plot hazard scores from models
plot_ggforest(model_no_interaction, data=train, main="Hazard ratios of prognostic model, no interaction")
plot_ggforest(model_interaction, data=train, main="Hazard ratios of prognostic model, interaction")
plot_ggforest(model_rf, data=train, main="Hazard ratios of prognostic model, random forest")

plot_ggforest(model_prognostic_no_interaction, data=train, main="Hazard ratios of prognostic model, no interaction")
plot_ggforest(model_prognostic_interaction, data=train, main="Hazard ratios of prognostic model, interaction")
plot_ggforest(model_prognostic_rf, data=train, main="Hazard ratios of prognostic model, random forest")

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
validation_combos <- evaluate_validation_datasets_combinations(
  results,
  "AUC values of risk scores on validation datasets",
  "Risk scores from models",
  "Validation datasets"
)

# Plot the 1, 3, and 5 year ROC curves for a dataset and risk score
plot_time_dependent_roc_curves(validation_combos[[6]], "prognostic_score_no_interaction")
plot_time_dependent_roc_curves(validation_combos[[6]], "prognostic_score_interaction")
plot_time_dependent_roc_curves(train[1:986, ], "risk_score_rf")

# Measure the p-values for the ROC curves of the prognostic and baseline models
# across the combinations of validation datasets
measure_auc_p_values(validation_combos)

# Identify all of the genes used in the prognostic models
genes_all <- c(genes_modelling_interaction, genes_modelling_no_interaction, genes_rf)
genes_all <- genes_all[!duplicated(genes_all)]
genes_all <- Filter(function(x) !any(grepl(":", x)), genes_all)

# A venn diagram to show overlapping genes between models
plot_overlapping_genes()

# Inspect the coefficients of models' features
summary(model_no_interaction)
summary(model_interaction)
plot(model_rf)
summary(model_combine)
summary(model_prognostic_no_interaction)
summary(model_prognostic_interaction)
summary(model_prognostic_rf)

summary(model_no_interaction)

# To be removed
i <- 3

qwe <- as.data.frame(apply(results[[2]][[i]], 2, rev))
qwe <- qwe[!duplicated(qwe$sequence_id), ]

asd <- data.frame(
  qwe$sequence_id,
  as.numeric(qwe$risk_score_interaction),
  qwe$risk_score_interaction_low,
  as.numeric(qwe$risk_score_no_interaction),
  qwe$risk_score_no_interaction_low,
  as.numeric(qwe$risk_score_rf),
  qwe$risk_score_rf_low,
  as.numeric(qwe$prognostic_score_no_interaction),
  qwe$prognostic_score_no_interaction_low,
  as.numeric(qwe$prognostic_score_interaction),
  qwe$prognostic_score_interaction_low,
  as.numeric(qwe$prognostic_score_rf),
  qwe$prognostic_score_rf_low
)

colnames(asd) <- c(
  "patient_id",
  "risk_score_interaction", "risk_score_interaction_low",
  "risk_score_no_interaction", "risk_score_no_interaction_low",
  "risk_score_rf", "risk_score_rf_low",
  "prognostic_score_no_interaction", "prognostic_score_no_interaction_low",
  "prognostic_score_interaction", "prognostic_score_interaction_low",
  "prognostic_score_rf", "prognostic_score_rf_low"
)
data_names
write.csv(asd, "~/Imperial/neuroblastoma_gene_signature/data/temp/scores_E-TABM-38.csv", row.names=FALSE)


# To make a kaplan meier plot
fit <- survfit(
  Surv(
    event_free_survival_days, event_free_survival
  ) ~ risk_score_interaction_low,
  data=results$train,
  na.action=na.pass
)
ggsurvplot(
  fit=fit,
  conf.int=TRUE,
  pval=TRUE,
  #risk.table="abs_pct",
  xlab="Years", 
  ylab="Overall survival probability",
  xscale=365, # Divide x-axis into years
  break.x.by=365, # Assume each year is only 365 days
  xlim=c(0, 1826) # Limit to five years
)

summary(model_prognostic_rf)
