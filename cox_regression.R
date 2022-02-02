library(dplyr)
library(survminer)

gset <- getGEO('GSE49711')
gset <- gset[[1]]

gse2 <- getGEO('GSE62564')
gse2 <- gse2[[1]]


# load microarray data
data <- read.csv('~/Desktop/NB project/GSE49711_SEQC_NB_TUC_G_log2.txt', sep = '\t')

# select rows of the 167 DE genes
# genelist is list of 167 DE genes
test <- data[data$X00gene_id %in% genelist,]

# remove gene name column and z-transform
test2 <- test[,2:499]
test2 <- 2^test2
test2 <- t(scale(t(test2)))

# select for columns of pData: here select for age, MYCN status, death
idx <- which(colnames(pData(gset)) %in% c('characteristics_ch1.3', 
                                          'characteristics_ch1.4', 
                                          'characteristics_ch1.9'))


# create metadata table
metadata <- data.frame(pData(gset)[,idx], row.names = pData(gset)$title)

# add EFS data obtained from GSE62564
metadata$EFS <- pData(gse2)[,13]

# ensure colnames of expression = rownames of metadata (sample number)
colnames(test2) <- rownames(metadata)

# create table with cox data
coxdata <- data.frame(metadata, t(test2))
colnames(coxdata)[1:3] <- c('Age', 'MYCN', 'Death')

# prepare phenotype columns for cox
coxdata$Age <- as.numeric(gsub("age at diagnosis: ", "", coxdata$Age))
coxdata$MYCN <- as.numeric(gsub("mycn status: ", "", coxdata$MYCN))
coxdata$Death <- as.numeric(gsub("death from disease: ", "", coxdata$Death))
coxdata$EFS <- as.numeric(gsub("efs day: ", "", coxdata$EFS))


res <- RegParallel(
  data = coxdata,
  formula = 'Surv(EFS, Death) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  blocksize = 50,
  variables = colnames(coxdata)[5:ncol(coxdata)])
  
res$Variable <- test$X00gene_id
res <- res[order(res$LogRank, decreasing = FALSE),]
final <- subset(res, LogRank < 0.05)


# kaplan-meier plot

# create table
survplotdata <- coxdata[,c('EFS', 'Death', 'X18580', 'X14580')]
colnames(survplotdata) <- c('EFS', 'Death', 'SLC22A4', 'NXT2')

# set high/low expression boundaries as described in arthur's paper
highExpr <- 0.5
lowExpr <- -0.5

# categorise into low/mid/high
survplotdata$SLC22A4 <- ifelse(survplotdata$SLC22A4 >= highExpr, 'High', 
                               ifelse(survplotdata$SLC22A4 <= lowExpr, 'Low', 'Mid'))
survplotdata$NXT2 <- ifelse(survplotdata$NXT2 >= highExpr, 'High', 
                            ifelse(survplotdata$NXT2 <= lowExpr, 'Low', 'Mid'))

# set mid to reference level
survplotdata$SLC22A4 <- factor(survplotdata$SLC22A4, 
                               levels = c('Mid', 'Low', 'High'))
survplotdata$NXT2 <- factor(survplotdata$NXT2, 
                            levels = c('Mid', 'Low', 'High'))

SLC22A4_plot <- ggsurvplot(survfit(Surv(EFS, Death) ~ SLC22A4,
                   data = survplotdata),
           data = survplotdata,
           risk.table = F,
           pval = TRUE,
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)


NXT2_plot <- ggsurvplot(survfit(Surv(EFS, Death) ~ NXT2,
                   data = survplotdata),
           data = survplotdata,
           risk.table = F,
           pval = TRUE,
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)
