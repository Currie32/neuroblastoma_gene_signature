library(dplyr)
library(ggpubr)
library(sm)
library(RVAideMemoire)


PATH = "/Users/isabellezhang/R_projects/MSc_neuroblastoma"

gene_exp <- read.csv(file.path(PATH, "datasets/stratification_table_gse.csv"))
gene_exp <- data.frame(gene_exp)

# non-parametric ancova
x <- gene_exp$CD9
y <- gene_exp$NXT2
g <- gene_exp$death
sm.ancova(y, x, g, h = 1, model = "none",
          xlim = c(-2,4),
          xlab = "CD9 expression",
          ylim = c(-2,4),
          ylab = "NXT2 expression")
legend(2.5, 3.5, legend=c("Survival", "Death"),
       col=c("pink", "green"), lty=1:2, cex=0.8)




# stratify by age
over18 <- gene_exp[c(which(gene_exp$over_18m==1)),]
under18 <- gene_exp[c(which(gene_exp$over_18m==0)),]

# spearman's rank test 
# this is a separate analysis from the ancova below
over18_r <- cor.test(over18$CD9, over18$NXT2, method = 'spearman', exact = FALSE)
under18_r <- cor.test(under18$CD9, under18$NXT2, method = 'spearman', exact = FALSE)

spearman.ci(over18$CD9, over18$NXT2, nrep = 1000, conf.level = 0.95)


# non-parametric ancova
x <- gene_exp$CD9
y <- gene_exp$NXT2
g <- gene_exp$over_18
sm.ancova(y, x, g, h = 1, model = "none",
          xlim = c(-2,4),
          xlab = "CD9 expression",
          ylim = c(-2,4),
          ylab = "NXT2 expression")
legend(2.5, 3.5, legend=c("Under 18 m", "Over 18 m"),
       col=c("pink", "green"), lty=1:2, cex=0.8)
