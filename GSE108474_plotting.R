# create Kaplan Meier curves using the dataset of patients with other brain cancers
# stratify the dataset by trkA activation or not


library(ggplot2)
library(survival)
library(survminer)

# Get data
PATH = "/Users/isabellezhang/R_projects/MSc_neuroblastoma"
data <- read.csv(file.path(PATH, "datasets/GSE108474_33_genes.csv"))

# Make column for trkA = 1 if activated
data$trkA_on <- ifelse(data$NTRK1 > 0 & data$TP53 < 0 
                       & data$PTPN6 < 0, 1, 0)

trkA_on <- data[c(which(data$trkA_on==1)),]

# Group into different brain cancers
astrocyt <- data[c(which(trkA_on$DISEASE_TYPE=="ASTROCYTOMA")),]
gbm <- data[c(which(trkA_on$DISEASE_TYPE=="GBM")),]
oligoden <- data[c(which(trkA_on$DISEASE_TYPE=="OLIGODENDROGLIOMA")),]




cat_exp <- function(gene){
  
  #input a col of z-scored gene expression
  #returns a categorical vector where:
  # H = high, z-score > 0.5
  # M = moderate, z-score >=-0.5 & <= 0.5
  # L = low, z-score < -0.5
  
  gene <- as.numeric(gene)
  gene_cat <- c()
  for (i in 1:length(gene)){
    
    if (gene[i] > 0.5){
      gene_cat <- append(gene_cat, "H")
    }
    else if (gene[i] < -0.5){
      gene_cat <- append(gene_cat, "L")
    }
    else{
      gene_cat <- append(gene_cat, "M")
    }
  }
  
  return(gene_cat)
}

plot_list <- list()


for (i in 2:(length(trkA_on)-4)){
  
  exp_status <- cat_exp(trkA_on[,i])
  d <- trkA_on
  d$exp_status <- exp_status
  
  par(mfrow=c(6,6))
  
  
  km_model <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, EVENT_OS)
                      ~ exp_status, 
                      data = trkA_on, na.action=na.pass)
  
  km_plot <- ggsurvplot(fit = km_model,
             conf.int=FALSE,
             pval = TRUE,
             xlab = "efs_MONTHS", 
             ylab = "Overall survival probability", 
             risk.table = TRUE,
             title = colnames(trkA_on[i])
  )
  
  plot_list[[i-1]] <- km_plot 
  
}




trkA_status <- data$trkA_on

km_model <- survfit(Surv(OVERALL_SURVIVAL_MONTHS, EVENT_OS)
                      ~ trkA_status, 
                      data = data, na.action=na.pass)
  
km_plot <- ggsurvplot(fit = km_model,
                      data = data,
                        conf.int=FALSE,
                        pval = TRUE,
                        xlab = "efs_MONTHS", 
                        ylab = "Overall survival probability", 
                        risk.table = TRUE)
  
  



