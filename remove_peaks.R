library(dplyr)

gse <- read.csv('~/Desktop/geo64_clean_nonz.csv')

# par(mfrow = c(12,12))
# par(mar = c(1,1,1,1))

# for(i in 1:nrow(gse)){
  #hist(as.numeric(gse[i,2:ncol(gse)]), main = "", label=F, plot = T)
#}
# peaks at -6 visible

# remove presumed sequencing error at -6.665
gse_removed <- gse %>% dplyr::na_if(-6.665)

#for(i in 1:nrow(gse_removed)){
#hist(as.numeric(gse_removed[i,2:ncol(gse)]), main = "", label = F ,plot = T)
#}
# peaks at -6 removed

# transform data by undoing log2, z-transform
expr_62 <- gse_removed[,2:ncol(gse_removed)]
expr_62 <- t(scale(t(2^expr_62)))

expr_62 <- data.frame(gse_removed$external_gene_name, expr_62)
colnames(expr_62)[1] <- 'external_gene_name'
  

#########

gse2 <- read.csv('~/Desktop/gse11_DEgenes_nonz.csv')

# remove sequencing errors, this time with value 0
gse2_removed <- gse2 %>% dplyr::na_if(0)

# transform
expr_49 <- gse2_removed[,2:ncol(gse2_removed)]
expr_49 <- t(scale(t(2^expr_49)))

expr_49 <- data.frame(gse2_removed$X.Gene, expr_49)
colnames(expr_49)[1] <- 'external_gene_name'

