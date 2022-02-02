library(dplyr)

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

