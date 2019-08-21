library(data.table)
library(plyr)
block <- read.table("/net/mulan/home/yasheng/Biobank/data/LDblock/chr1.txt")
chr <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/genotype/chr1.bim"))

block_num <- 1
getIdx <- function(block_num){
  start <- block[block_num, 2]
  end <- block[block_num, 3]
  in_idx <- ifelse(chr[, 4] - start > 0 & end - chr[, 4] > 0, 1, 0)
  return(in_idx)
}

x <- aaply(c(1: nrow(block)), 1, getIdx)

num_snp_block <- rowSums(as.matrix(x))
max(num_snp_block)
min(num_snp_block)