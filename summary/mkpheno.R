library(data.table)
library(parallel)
## total sample for GWAS
pheno_tot <- data.frame(fread("/net/mulan/disk2/yasheng/phenotype_file/pheno_tot_v2.txt"))
fam <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/genotype/chr1.fam"))
fam_pheno <- cbind(fam[, c(1:5)], pheno_tot)
fam_output <- mclapply(c(1:22), function(chr){
  write.table(fam_pheno, file = paste0("/net/mulan/disk2/yasheng/tot_plink_file/chr", chr, ".fam"), 
              col.names = F, row.names = F, quote = F)
  return(chr)
},  mc.cores=22)

## training sample for five flod cross validation
setwd("/net/mulan/home/yasheng/Biobank/data/pheno/pheno_alp_cross")
pheno_tot <- matrix(NA, nrow = 337198, ncol = 80)
for (pheno in c(1: 16)){
  b <- (pheno-1)*5+1
  e <- (pheno-1)*5+5
  pheno_tot[, c(b: e)] <- as.matrix(read.table(paste0("phenoRsd_", pheno, "_train.txt")))
}
fam <- data.frame(fread("/net/mulan/disk2/yasheng/tot_plink_file/chr1.fam"))
fam_pheno <- cbind(fam[, c(1:5)], pheno_tot)
fam_output <- mclapply(c(1:22), function(chr){
  write.table(fam_pheno, file = paste0("/net/mulan/disk2/yasheng/tot_plink_file/chr", chr, ".fam"), 
              col.names = F, row.names = F, quote = F)
  return(chr)
 },  mc.cores=22)
