library(parallel)
## output to fam file for cross validation 
pheno_num <- 16
pheno <- matrix(NA, nrow = 337198, ncol = pheno_num*5)
for (p in 1: pheno_num){
  begin <- (p-1)*5+1
  end <- (p-1)*5+5
  pheno_str <- paste0("/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_", p, "_train.txt")
  pheno[, c(begin:end)] <- as.matrix(fread(pheno_str))
}
fam <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/genotype/chr1.fam"))
fam_pheno <- cbind(fam[, c(1:5)], pheno)
fam_output <- mclapply(c(1:22), function(chr){
  write.table(fam_pheno, file = paste0("/net/mulan/disk2/yasheng/plink_file/genotype/chr", chr, ".fam"), 
              col.names = F, row.names = F, quote = F)
  return(chr)
},  mc.cores=22)

## output to fam file for GWAS analysis
fam <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/genotype/chr1.fam"))
pheno <- data.frame(fread("/net/mulan/disk2/yasheng/phenotype_file/pheno_tot_v2.txt"))
fam_pheno <- cbind(fam[, c(1:5)], pheno)
fam_output <- mclapply(c(1:22), function(chr){
  write.table(fam_pheno, file = paste0("/net/mulan/disk2/yasheng/plink_file/genotype/chr", chr, ".fam"), 
              col.names = F, row.names = F, quote = F)
  return(chr)
},  mc.cores=22)