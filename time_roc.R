library(survival)
library(survivalROC)
library(ggplot2)

train <- read.csv('training_data.csv')
target <- read.csv('target_risk_df.csv')


train_roc_s <-  survivalROC(train$event_free_survival_days, 
                            train$event_free_survival, 
                            train$risk_score, predict.time = 365,
                            method = 'KM')

target_roc_s <- survivalROC(target$event_free_survival_days, 
                            target$event_free_survival, 
                            target$risk_score, predict.time = 365, 
                            method = 'KM')

target_roc_s3 <- survivalROC(target$event_free_survival_days, 
                            target$event_free_survival, 
                            target$risk_score, predict.time = 3*365, 
                            method = 'KM')

target_roc_s5 <- survivalROC(target$event_free_survival_days, 
                             target$event_free_survival, 
                             target$risk_score, predict.time = 5*365, 
                             method = 'KM')


# function takes survivalROC output and time at which prediction was made

tROC_plot <- function(survivalroc_obj, predtime){
  df <- data.frame(survivalroc_obj$FP, survivalroc_obj$TP)
  auc_label <- sprintf('AUC = %.4f', survivalroc_obj$AUC)
  p <- ggplot(data = df, aes(x = df[,1], y = df[,2])) + geom_line(size = 1.5) +
    geom_abline(intercept = 0, slope = 1) + theme_light() +
    annotate('text', x = 0.1, y = 0.9, label = auc_label, size = 5) +
    ylab('Sensitivity') + xlab('1-Specificity') + 
    ggtitle(sprintf('Survival at year %.0f', predtime))
  return(p)
}






