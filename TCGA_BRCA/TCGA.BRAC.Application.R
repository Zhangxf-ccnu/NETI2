# Run NETI2 on TCGA breast cancer datasets

rm(list=ls())
library(NETI2)

data("TCGA.BRCA")

# run NETI2 using NETI2
TCGA.BRCA.NETI2 = NETI2(TCGA.BRCA$BRCA.data,TCGA.BRCA$BRCA.purity, lambda = 1.6, tau = 0.4,delta = 0.2)