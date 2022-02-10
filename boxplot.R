library(dplyr)
library(ggplot2)

# create stratification table 

pheno <- getGEO('GSE49711')
pheno <- pheno[[1]]

# select for age and mycn columns
idx <- which(colnames(pData(pheno)) %in% c('characteristics_ch1.3', 
                                          'characteristics_ch1.4'))

# create table
base <- data.frame(pData(pheno)[idx], rownames = pData(pheno)$title)

# reformat into 0 and 1 only
base$mycn <- as.numeric(gsub("mycn status: ", "", base$mycn))
base$over_18m <- as.numeric(gsub("age at diagnosis: ", "", base$over_18m))
base$over_18m <- ifelse(base$over_18m > 548, 1, 0)

# load microarray expression and risk scores data
expr <- read.csv('~/Desktop/GSE49711_SEQC_NB_MAV_G_log2.20121127.txt', sep = '\t')
risk <- read.csv('~/Desktop/training_data.csv')

# z-transform data
expr_t <- expr[,10:ncol(expr)]
expr_t <- scale(t(2^expr_t))

expr_t <- data.frame(expr_t)
colnames(expr_t) <- expr$X.Gene

# subset for TrkA pathway genes
cols <- which(colnames(expr_t) %in% c('TP53', 'PTPN6', 'NTRK1'))
trka_data <- data.frame(expr_t[cols])

# TrkA expression present if high NTRK & low TP53 & low PTPN6
base$trkA <- ifelse(trka_data$NTRK1 > 0 & trka_data$TP53 < 0 
                    & trka_data$PTPN6 < 0, 1, 0)

# load expression data with peaks removed, subsetted for target genes
expr_target <- read.csv('~/Desktop/gse49711_peaks_removed_ztransformed.csv')

# reformat data
rownames(expr_target) <- expr_target$external_gene_name
expr_target <- expr_target[,3:ncol(expr_target)]

# create dataframe of stratification levels and expression
data <- data.frame(base, t(expr_target))
data$risk_score <- risk$risk_score[499:nrow(risk)]

# create interaction factors to subset into stratification + survival status
inter_age <- interaction(data$over_18m, data$death)

# rename factors
levels(inter_age) <- c('Under.Alive', 'Over.Alive', 'Under.Dead','Over.Dead')

# reorder factor levels for graphing
inter_age <- factor(inter_age, levels = c('Under.Alive', 'Under.Dead', 'Over.Alive',
                                  'Over.Dead'))

inter_MYCN <- interaction(data$mycn, data$death)
levels(inter_MYCN) <- c('No.Alive', 'Yes.Alive', 'No.Dead','Yes.Dead')
inter_MYCN <- factor(inter_MYCN, levels = c('No.Alive', 'No.Dead', 'Yes.Alive',
                                          'Yes.Dead'))

inter_TrkA <- interaction(data$trkA, data$death)
levels(inter_TrkA) <- c('No.Alive', 'Yes.Alive', 'No.Dead','Yes.Dead')
inter_TrkA <- factor(inter_TrkA, levels = c('No.Alive', 'No.Dead', 'Yes.Alive',
                                            'Yes.Dead'))

# create labeller text
status_age <- c(`0` = 'Under 18 months', `1` = 'Above 18 months')
status_MYCN <- c(`0` = 'No MYCN amplification', `1` = 'MYCN amplification')
status_TrkA <- c(`0` = 'No TrkA pathway', `1` = 'TrkA pathway present')

p1 <- ggplot(aes(y = risk_score, x = inter_age, fill = factor(death)), data = data) + 
  geom_boxplot() + theme_light() + ggtitle('Stratified by age') +
  facet_grid(~factor(over_18m), scales = 'free_x',
             labeller = as_labeller(status_age)) +
  scale_fill_discrete(labels = c("Alive", 'Deceased')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab('Risk score') +
  labs(fill = 'Status')

p2 <- ggplot(aes(y = risk_score, x = inter_MYCN[-c(112, 189, 255, 312, 361)], fill = factor(death)), 
             data = data[-c(112, 189, 255, 312, 361),]) + 
  geom_boxplot() + theme_light() + ggtitle('Stratified by MYCN') +
  facet_grid(~factor(mycn), scales = 'free_x',
             labeller = as_labeller(status_MYCN)) +
  scale_fill_discrete(labels = c("Alive", 'Deceased')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab('Risk score') +
  labs(fill = 'Status')

p3 <- ggplot(aes(y = risk_score, x = inter_TrkA, fill = factor(death)), data = data) + 
  geom_boxplot() + theme_light() + ggtitle('Stratified by TrkA') +
  facet_grid(~factor(data$trkA), scales = 'free_x',
             labeller = as_labeller(status_TrkA)) +
  scale_fill_discrete(labels = c("Alive", 'Deceased')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab('Risk score') +
  labs(fill = 'Status')

grid.arrange(p1, p2, p3, nrow=3)



######

# t-testing

shapiro.test(data$risk_score)

wilcox.test(data[which(inter_age == 'Under.Alive'),]$risk_score, 
            data[which(inter_age == 'Under.Dead'),]$risk_score)

wilcox.test(data[which(inter_age == 'Over.Alive'),]$risk_score, 
            data[which(inter_age == 'Over.Dead'),]$risk_score)

wilcox.test(data[which(inter_MYCN == 'No.Alive'),]$risk_score, 
            data[which(inter_MYCN == 'No.Dead'),]$risk_score)

wilcox.test(data[which(inter_MYCN == 'Yes.Alive'),]$risk_score, 
            data[which(inter_MYCN == 'Yes.Dead'),]$risk_score)

wilcox.test(data[which(inter_TrkA == 'No.Alive'),]$risk_score, 
            data[which(inter_TrkA == 'No.Dead'),]$risk_score)

wilcox.test(data[which(inter_TrkA == 'Yes.Alive'),]$risk_score, 
            data[which(inter_TrkA == 'Yes.Dead'),]$risk_score)

