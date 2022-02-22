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
  candidate_genes <- c()
  
  for (gene in genes) {
    
    # gene_strata <- paste(gene, ":strata(tgroup)", sep="")
    
    # Perform cox regression
    result <- coxph(
      as.formula(paste("survival_train ~ ", gene)),
      data=df,
      na.action=na.pass
    )
    
    # Get a store p_value
    # result_p_values <- summary(result)$coefficients[,5]
    # significant_genes <- names(result_p_values[result_p_values < P_VALUE])
    # if (length(significant_genes) > 0) {
    #   candidate_genes <- append(candidate_genes, gene)
    # }
    
    p_value <- summary(result)$coefficients[5]
    if (p_value <= P_VALUE) {
      candidate_genes <- append(candidate_genes, gene)
    }
  }
  
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
  
  

  # Multivariate cox regression
  if (phenotype_features) {
    
    # Format the candidate genes the cox regression
    features = c(
      "under_18_months",
      "mycn_amplification",
      "inss_stage_2",
      "inss_stage_3",
      "inss_stage_4",
      "inss_stage_4S",
      candidate_genes
    )
    # features_joined <- paste(
    #   paste(features, ":strata(tgroup)"), collapse = " + "
    # )
    features_joined <- paste(features, collapse=" + ")
    
    model <- coxph(
      as.formula(paste("survival_train ~ ", features_joined)),
      data=df,
      na.action=na.pass
    )
  }
  else if (interaction) {
    
    # Format the candidate genes the cox regression
    # features_joined <- paste(
    #   paste(features, ":strata(tgroup)"), collapse = " + "
    # )
    candidate_genes_joined <- paste(candidate_genes, collapse=" + ")
    
    model <- coxph(
      as.formula(paste("survival_train ~ (", candidate_genes_joined, ")^2")),
      data=df,
      na.action=na.pass
    )
  }
  else {
    
    # canadidate_genes_joined <- paste(
    #   paste(candidate_genes, ":strata(tgroup)"), collapse = " + "
    # )
    candidate_genes_joined <- paste(candidate_genes, collapse=" + ")
    
    model <- coxph(
      as.formula(paste("survival_train ~ ", candidate_genes_joined)),
      data=df,
      na.action=na.pass
    )
  }
  
  # Find the significant features
  result <- summary(model)
  p_values <- result$coefficients[,5]
  significant_features <- names(p_values[p_values < P_VALUE])
  
  # significant_features <- str_remove(significant_features, "tgroup=1")
  # significant_features <- str_remove(significant_features, "tgroup=2")
  # significant_features <- str_remove(significant_features, "tgroup=3")
  # significant_features <- str_remove(significant_features, ":strata\\(tgroup\\)")
  # significant_features <- str_remove(significant_features, "strata\\(tgroup\\):")
  # significant_features <- significant_features[!duplicated(significant_features)]
  
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
    as.formula(paste("survival_train ~ ", features_joined)),
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
  
  # best_model_features <- str_remove(best_model_features, "tgroup=1")
  # best_model_features <- str_remove(best_model_features, "tgroup=2")
  # best_model_features <- str_remove(best_model_features, "tgroup=3")
  # best_model_features <- str_remove(best_model_features, ":strata\\(tgroup\\)")
  # best_model_features <- str_remove(best_model_features, "strata\\(tgroup\\):")
  # best_model_features <- best_model_features[!duplicated(best_model_features)]
  
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
    as.formula(paste("survival_train ~ ", features_joined)),
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
  
  # best_model_baseline_features <- str_remove(best_model_baseline_features, "tgroup=1")
  # best_model_baseline_features <- str_remove(best_model_baseline_features, "tgroup=2")
  # best_model_baseline_features <- str_remove(best_model_baseline_features, "tgroup=3")
  # best_model_baseline_features <- str_remove(best_model_baseline_features, ":strata\\(tgroup\\)")
  # best_model_baseline_features <- best_model_baseline_features[!duplicated(best_model_baseline_features)]
  
  # Train a model using the best combination of features
  features_joined <- paste(paste(best_model_baseline_features, ":strata(tgroup)"), collapse = " + ")
  
  model_baseline <- coxph(
    survival_train ~ under_18_months:strata(tgroup) + mycn_amplification + inss_stage_2:strata(tgroup) + inss_stage_3:strata(tgroup) + inss_stage_4:strata(tgroup) + inss_stage_4S:strata(tgroup),
    data=train_split,
    na.action=na.pass
  )
  
  temp <- cox.zph(model_baseline)
  print(temp)
  plot(temp)

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
test_wilms_tumor <- read.csv(file.path(PATH, "test_wilms_tumor.csv"))

# Add under_18_months feature
train$under_18_months <- is_under_18_months(train)
test_GSE85047$under_18_months <- is_under_18_months(test_GSE85047)
test_target$under_18_months <- is_under_18_months(test_target)
test_wilms_tumor$under_18_months <- is_under_18_months(test_wilms_tumor)

# Get the genes from the training dataset
genes <- colnames(train)[14:length(colnames(train)) - 1]


splits <- c(400, 600)

train_split <- survSplit(survival_train ~ ., 
                         data=train,
                         cut = splits,
                         zero=-1,
                         episode= "tgroup"
)

test_GSE85047_split <- survSplit(survival_test_GSE85047 ~ ., 
                                 data=test_GSE85047,
                                 cut = splits,
                                 zero=-1,
                                 episode= "tgroup"
)

test_target_split <- survSplit(survival_test_target ~ ., 
                               data=test_target,
                               cut = splits,
                               zero=-1,
                               episode= "tgroup"
)

test_wilms_tumor_split <- survSplit(survival_wilms_tumor ~ ., 
                                    data=test_wilms_tumor,
                                    cut = splits,
                                    zero=-1,
                                    episode= "tgroup"
)


# Create survival for cox regressions
survival_train <- Surv(train$event_free_survival_days, train$event_free_survival)
survival_test_GSE85047 <- Surv(test_GSE85047$event_free_survival_days, test_GSE85047$event_free_survival)
survival_test_target <- Surv(test_target$event_free_survival_days, test_target$event_free_survival)
survival_wilms_tumor <- Surv(wilms_tumor$event_free_survival_days, wilms_tumor$event_free_survival)

survival_asd <- Surv(asd$event_free_survival_days, asd$event_free_survival)

# Find the combination of features that are most predictive of survival

genes <- genes[!genes %in% c("TXNDC5")]

genes_candidate <- univariate_cox_regression(genes, train)
genes_significant <- multivariate_cox_regression(genes_candidate, survival_train, train, FALSE, FALSE)
genes_interaction <- multivariate_cox_regression(genes_significant, survival_train, train, FALSE, TRUE)
genes_best <- aic_feature_selection(genes_significant, survival_train, train_split, '/')
genes_best <- genes_best[!genes_best %in% c(
  "HEBP2", "TXNDC5", "NXT2:TXNDC5", "HEBP2:RACK1", "RACK1:TXNDC5"
  )]
genes_best <- append(genes_best, c("VAT1L", "NUDT10_NUDT11"))

model <- train_model(genes_best, train)




# Train model
model <- coxph(
  survival_train ~ CD9:strata(tgroup) + DLG2:strata(tgroup) + HSD17B12 + NUDT10_NUDT11 + NXT2:strata(tgroup) + PPIA + RACK1,
  data=train_split,
  na.action=na.pass
)

# Train the model using the best gene set
#model <- train_model(genes_best, train)
model_baseline <- baseline_model(train_split, survival_train)

temp <- cox.zph(model_baseline)
temp <- cox.zph(model)
print(temp)
plot(temp)

asd <- read.csv("~/Downloads/GSE108474_processed.csv")

asd[colnames(asd) %in% genes_best]
asd$NXT2 <- 0
asd$TXNDC5 <- 0
asd$CD147 <- 0
asd$RACK1 <- 0
asd$`CD9:NXT2` <- 0

asd$risk_score <- predict(model, asd)
asd$risk_score_low <- asd$risk_score < median(train$risk_score)

asd$event_free_survival <- asd$EVENT_OS
asd$event_free_survival_days <- asd$OVERALL_SURVIVAL_MONTHS * 30


# Create risk scores
# train$risk_score <- predict(model, train)
# train$risk_score_low <- train$risk_score < median(train$risk_score)
# 
# test_GSE85047$risk_score <- predict(model, test_GSE85047)
# test_GSE85047$risk_score_low <- test_GSE85047$risk_score < median(train$risk_score)
# 
# test_target$risk_score <- predict(model, test_target)
# test_target$risk_score_low <- test_target$risk_score < median(train$risk_score)
# 
# test_wilms_tumor$risk_score <- predict(model, test_wilms_tumor)
# test_wilms_tumor$risk_score_low <- test_wilms_tumor$risk_score < median(train$risk_score)

train_split$risk_score <- predict(model, train_split)
train_split$risk_score_low <- train_split$risk_score < median(train_split$risk_score)

test_GSE85047_split$risk_score <- predict(model, test_GSE85047_split)
test_GSE85047_split$risk_score_low <- test_GSE85047_split$risk_score < median(train_split$risk_score)

test_target_split$risk_score <- predict(model, test_target_split)
test_target_split$risk_score_low <- test_target_split$risk_score < median(train_split$risk_score)

test_wilms_tumor_split$risk_score <- predict(model, test_wilms_tumor_split)
test_wilms_tumor_split$risk_score_low <- test_wilms_tumor_split$risk_score < median(train_split$risk_score)

# Create baseline model risk scores
model_baseline <- baseline_model(train_split, survival_train)
train_split$risk_score_baseline <- predict(model_baseline, train_split)
train_split$risk_score_low_baseline <- train_split$risk_score_baseline < median(train_split$risk_score_baseline)

test_GSE85047_split$risk_score_baseline <- predict(model_baseline, test_GSE85047_split)
test_GSE85047_split$risk_score_low_baseline <- test_GSE85047_split$risk_score_baseline < median(train_split$risk_score_baseline)

test_target_split$risk_score_baseline <- predict(model_baseline, test_target_split)
test_target_split$risk_score_low_baseline <- test_target_split$risk_score_baseline < median(train_split$risk_score_baseline)

test_wilms_tumor_split$risk_score_baseline <- predict(model_baseline, test_wilms_tumor_split)
test_wilms_tumor_split$risk_score_low_baseline <- test_wilms_tumor_split$risk_score_baseline < median(train_split$risk_score_baseline)


# Plot ROC curves
plot_roc_curve(train_split, "risk_score")
plot_roc_curve(test_GSE85047_split, "risk_score")
plot_roc_curve(test_target_split, "risk_score")
plot_roc_curve(test_wilms_tumor_split, "risk_score")
plot_roc_curve(asd, "risk_score")

plot_roc_curve(train_split, "risk_score_baseline")
plot_roc_curve(test_GSE85047_split, "risk_score_baseline")
plot_roc_curve(test_target_split, "risk_score_baseline")
plot_roc_curve(test_wilms_tumor_split, "risk_score_baseline")


# Kaplan Meier Train
fit_train <- survfit(
  survival_train ~ risk_score_low,
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
features_significant <- multivariate_cox_regression(c('risk_score'), survival_train, train_split, TRUE, FALSE)
features_best <- aic_feature_selection(features_significant, survival_train, train_split, '')
model_prognostic <- train_model(features_best, train_split)

features_significant <- c("inss_stage_2", "inss_stage_3", "inss_stage_4", "risk_score")

model_prognostic <- coxph(
  survival_train ~ inss_stage_2 + inss_stage_3 + inss_stage_4:strata(tgroup) + inss_stage_4S +risk_score,
  # as.formula(paste("survival_train ~ ", paste(paste(features_best, ":strata(tgroup)"), collapse = " + "))),
  data=train_split,
)


temp <- cox.zph(model_prognostic)
print(temp)
plot(temp)

# Create prognostic score (i.e. the more advanced risk score)
train_split$prognostic_score <- predict(model_prognostic, train_split)
train_split$prognostic_score_low <- train_split$prognostic_score < median(train_split$prognostic_score)

test_GSE85047_split$prognostic_score <- predict(model_prognostic, test_GSE85047_split)
test_GSE85047_split$prognostic_score_low <- test_GSE85047_split$prognostic_score < median(train_split$prognostic_score)

test_target_split$prognostic_score <- predict(model_prognostic, test_target_split)
test_target_split$prognostic_score_low <- test_target_split$prognostic_score < median(train_split$prognostic_score)

# Plot results
plot_roc_curve(train_split, "prognostic_score")
plot_roc_curve(test_GSE85047_split, "prognostic_score")
plot_roc_curve(test_target_split, "prognostic_score")

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
  survival_test_target ~ prognostic_score_low,
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
compare_auc_values(test_GSE85047_split, 500)
compare_auc_values(test_target_split, 500)

genes_best

# Plot gene expression heatmap
pheatmap(
  train[order(train$risk_score),][colnames(train) %in% genes_best],
  show_rownames=FALSE,
  cluster_rows=F,
  cluster_cols=F,
  breaks=seq(-2, 3, by = 0.06),
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



model




