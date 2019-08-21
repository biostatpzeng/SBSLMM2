library(data.table)
rm(list=ls())
##neale result
neale <- data.frame(fread("/net/mulan/Biobank/rawdata/EGAD00010001226/ukbbSummaryData/21001.assoc.tsv.gz?dl=0%20-O%2021001.assoc.tsv.gz"))

##my result
pheno <- data.frame(fread("/net/mulan/disk2/yasheng/output/summary_pheno4.assoc.txt"))

##intersect
inter <- intersect(neale[, 2], pheno[, 2])
neale_inter <- neale[match(inter, neale[, 2]), ]
pheno_inter <- pheno[match(inter, pheno[, 2]), ]
idx_pheno <- which(as.numeric(pheno_inter[, 11]) == 0)

plot(-log10(neale_inter$pval[-idx_pheno]), -log10(as.numeric(pheno_inter[-idx_pheno, 11])))