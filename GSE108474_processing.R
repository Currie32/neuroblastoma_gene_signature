# process the GSE108474 dataset 
# return a dataset of the patients with expression data of the 33 genes + trkA 

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108474

library(GEOquery)
library(biomaRt)
library(dplyr)
library(batman)

#PATH = "C:/Users/tomri/OneDrive - Imperial College London/neuroblastoma/scripts/data"

# load series and platform data from GEO
Sys.setenv(VROOM_CONNECTION_SIZE = 100000000)
gset <- getGEO("GSE108474", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)


######################


#To convert RefSeqID into external_gene_name
ensembl <- useMart("ensembl") #database
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl) # dataset

array_to_gname = getBM(attributes = c("affy_hg_u133_plus_2", "external_gene_name"),
                       filter= c("affy_hg_u133_plus_2"),
                       values=rownames(ex), mart=ensembl) #query


#table to relevant genes

genes <- read.csv(file.path("datasets/common_genes.csv"), row.names = 1)
genes <- rbind(genes, "NTRK1") 
genes <- rbind(genes, "PTPN6") 
genes <- rbind(genes, "TP53") 

array_to_gname <- array_to_gname[array_to_gname$external_gene_name %in% genes$x, ]
length(unique(array_to_gname$external_gene_name))#23 unique genes out of the 33 univariate cox 

#subset ex to only probes in array_to_gname
ex <- ex[rownames(ex) %in% array_to_gname$affy_hg_u133_plus_2, ]
ex <- data.frame(ex)

#add column for gene names
g <- c()
for (i in rownames(ex)){
  g <- append(g, array_to_gname$external_gene_name[which(
    array_to_gname$affy_hg_u133_plus_2== i)])
}
ex$name <- g

#get average gene exp.
dummy_frame <- data.frame()
dup <- unique(g[duplicated(g)])

for (i in dup){
  r <- colMeans(ex[which(ex$name == i),1:length(ex)-1 ])
  dummy_frame <- rbind(dummy_frame, r)
}
dummy_frame$name <- dup
colnames(dummy_frame) <- colnames(ex)

dummy_frame <- rbind(dummy_frame, ex[which(!ex$name %in% dup), ])
ex <- dummy_frame
rownames(ex) <- ex$name
ex <- ex[, 1:length(ex)-1]

ex <- scale(t(ex))

#


info <- read.csv(file.path("datasets/GSE108474_REMBRANDT_clinical.data.txt"), sep = "\t")
features <- c("SUBJECT_ID", "EVENT_OS", "OVERALL_SURVIVAL_MONTHS","DISEASE_TYPE")
#filter for patients with expression data
have_exp <- to_logical(info$GENE_EXP, custom_true = c("YES"))
info <- info[which(have_exp == TRUE), features]
#remove rows with incomplete survival info
info <- na.omit(info)

#equalise patient ID's
pid <- read.csv(file.path("datasets/GSE108474_patient_id.csv"))
pid <- pid[, 1:2]

sample <- c()
for (i in info$SUBJECT_ID){
  sample <- append(sample, pid$Accession[which(pid$Title == i)])
}
#remove old ID, append new
info <- info[, -c(1)]
info$sample <- sample

#remove the last 8 patients which have N/A expression 
ex <- ex[1:(nrow(ex)-8), ]

#join survival features with expression
sample <- unlist(rownames(ex))
ex <- cbind(ex, sample)
ex <- merge(ex, info, by="sample")
rownames(ex) <- ex$sample
ex <- ex[, -c(1)]



#write.csv(ex, "GSE108474_processed.csv")

