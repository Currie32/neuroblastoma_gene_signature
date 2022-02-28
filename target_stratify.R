phenot <- read.csv('~/Desktop/NB project/nbl_target_2018_pub_clinical_data.tsv',
                   sep = '\t')
exprt <- read.csv('~/Desktop/data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt',
                  sep = '\t')


samples <- colnames(exprt)[3:ncol(exprt)]
phenot <- phenot[phenot$Sample.ID %in% samples,]

cols <- which(colnames(phenot) %in% c('Diagnosis.Age..days.', 
                                      'Tumor.Sample.Histology', 
                                      'MYCN', 'Overall.Survival.Status'))

base <- data.frame(phenot[cols], row.names = phenot$Sample.ID)
colnames(base) <- c('over_18m', 'mycn', 'death', 'histology')

base$mycn <- ifelse(base$mycn == 'Not Amplified', 0, 1)
base$over_18m <- ifelse(base$over_18m > 548, 1, 0)
base$histology <- ifelse(base$histology == 'Unfavorable', 1,
                        ifelse(base$histology == 'Favorable', 0, NA))


base$death <- ifelse(base$death == '1:DECEASED', 1,
                     ifelse(base$death == '0:LIVING', 0, NA))

genelist <- append(genelist, c('PTPN6', 'TP53', 'NTRK1'))


exprt <- exprt[exprt$Hugo_Symbol %in% genelist,]
exprt <- exprt[!duplicated(exprt$Hugo_Symbol),]
trka_cols <- which(exprt$Hugo_Symbol %in% c('TP53', 'PTPN6', 'NTRK1'))

trka_data <- data.frame(exprt[trka_cols,])
rownames(trka_data) <- c('NTRK1', 'PTPN6', 'TP53')
trka_data <- trka_data[,3:ncol(trka_data)]

trka_data <- t(trka_data)
trka_data <- data.frame(trka_data)

base$trkA <- ifelse(trka_data$NTRK1 > 0 & trka_data$TP53 < 0 
                    & trka_data$PTPN6 < 0, 1, 0)


rownames(exprt) <- exprt$Hugo_Symbol
exprt <- exprt[,3:ncol(exprt)]
exprt <- data.frame(t(exprt))

expr_stratify <- data.frame(base, exprt)

              