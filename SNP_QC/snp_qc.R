library(data.table)
setwd("/net/mulan/Biobank/rawdata/EGAD00010001225/001")

## result
tot_list <- vector("numeric", 22)
duplicated_list <- vector("numeric", 22)
maf_list <- vector("numeric", 22)
INFO_list <- vector("numeric", 22)

for (i in 1:22){
  snp <- data.frame(fread(paste0("ukb_mfi_chr", i, "_v2.txt")))
  tot_list[i] <- nrow(snp)
  
  ## duplicated
  freq <- data.frame(table(snp[, 1]))
  snp_dup <- as.character(freq[freq[, 2] > 1, 1])
  duplicated_list[i] <- length(snp_dup)

  ## maf < 0.01
  snp_maf <- snp[snp[, 5] < 0.01, 1]
  maf_list[i] <- length(snp_maf)

  ## INFO < 0.8
  snp_INFO <- snp[snp[, 6] < 0.8, 1]
  INFO_list[i] <- length(snp_INFO)

  ## delete snp list
  snp_del <- unique(c(snp_dup, snp_maf, snp_INFO))
  write.table(snp_del, file = paste0("/net/mulan/disk2/yasheng/plink_file/removal_snp_list/chr", i, ".txt"), 
              row.names = F, col.names = F, quote = F)
}

## summary
tot_sum <- sum(tot_list)
duplicated_sum <- sum(duplicated_list)
maf_sum <- sum(maf_list)
INFO_sum <- sum(INFO_list)

cat ("total SNP: ", tot_sum, "\n")
cat ("deplicated SNP: ", duplicated_sum, "\n")
cat ("maf<0.01 SNP: ", maf_sum, "\n")
cat ("INFO<0.8 SNP: ", INFO_sum, "\n")
