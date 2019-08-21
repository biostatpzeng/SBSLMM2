#! /usr/bin/env Rscript
library(plyr)
library(dplyr)
library(data.table)
library(optparse)

tot <- 0

for (chr in 1: 21){
 
  bim_str <- paste0("/net/mulan/home/zhongshy/UKBiobank/datamanagement/rssnp/chr_", chr, ".bim")
  mfi_str <- paste0("/net/mulan/Biobank/rawdata/EGAD00010001225/001/ukb_mfi_chr", chr, "_v2.txt")

  bim <- data.frame(fread(bim_str))
  mfi <- data.frame(fread(mfi_str))
  bim_freq <- data.frame(table(bim[, 2]))
  bim_uni <- bim_freq[bim_freq[, 2] == 1, 1]

  mfi_high <- mfi[mfi[, 6] >= 0.8, 1]

  snp_inter <- intersect(bim_uni, mfi_high)
  write.table(snp_inter, file = paste0("/net/mulan/home/yasheng/Biobank/plink_file/chr", chr, ".txt"), 
              row.names = F, col.names = F, quote = F)
  
}


 fam <- read.table("/net/mulan/home/yasheng/Biobank/plink_file/chr5.fam")
 fam[, 6] <- 1
for (i in 1: 22){
  
 write.table(fam, file = paste0("/net/mulan/home/yasheng/Biobank/plink_file/chr", i, ".fam"), 
             row.names = F, col.names = F, quote = F)
  
}


