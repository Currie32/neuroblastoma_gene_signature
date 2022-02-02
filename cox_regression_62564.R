library(dplyr)
library(RegParallel)
library(survival)
library(GEOquery)

data <- read.csv('geo64_clean_nonz.csv') #file already filtered for 167

# remove gene name column and z-transform
test2 <- data[,2:499]
test2 <- t(scale(t(test2)))

#gse64 pData does not contain mycn or death so geo11 is used
gset <- getGEO('GSE49711')
gset <- gset[[1]]
#gse11 does not include EFS data, so gse64 is used
gse2 <- getGEO('GSE62564')
gse2 <- gse2[[1]]


# select for columns of pData: here select for age, MYCN status, death
idx <- which(colnames(pData(gset)) %in% c('characteristics_ch1.3', 
                                          'characteristics_ch1.4', 
                                          'characteristics_ch1.9'))


# create metadata table
metadata <- data.frame(pData(gset)[,idx], row.names = pData(gset)$title)

# add EFS data obtained from GSE62564
metadata$EFS <- pData(gse2)[,13]


#check that sample names match exactly between pdata and Z-scores 
#all((colnames(test2) == rownames(metadata)) == TRUE)  
# TRUE

# create a merged pdata and Z-scores object
coxdata <- data.frame(metadata, t(test2))
colnames(coxdata)[1:3] <- c('Age', 'MYCN', 'Death')

# prepare phenotype columns for cox
coxdata$Age <- as.numeric(gsub("age at diagnosis: ", "", coxdata$Age))
coxdata$MYCN <- as.numeric(gsub("mycn status: ", "", coxdata$MYCN)) #NA for 5 rows, not an issue at this stage
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

res$Variable <- data$external_gene_name
res <- res[order(res$LogRank, decreasing = FALSE),]
final <- subset(res, LogRank < 0.05)

write.csv(final, "cox_regression_results_62564.csv")
