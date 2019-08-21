#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)
library(lassosum)

## Parameter setting
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary data prefix", metavar="character"), 
  make_option("--ref", type="character", default=NULL,
              help="INPUT: reference file prefix", metavar="character"), 
  make_option("--test", type="character", default=NULL,
              help="INPUT: max number prefix", metavar="character"),
  make_option("--beta", type="character", default=NULL,
              help="OUTPUT: output of coefficient", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/net/mulan/disk2/yasheng/simulation2/SLM/output/summary_block5_her0.5_cross8_dist1_ps0.001.assoc.txt",
#             ref = "/net/mulan/disk2/yasheng/sample500/frq/chr1",
#             test = "/net/mulan/disk2/yasheng/simulation2/SLM/data/test_block5_her0.5_cross8_dist1_ps0.001",
#             beta = "/net/mulan/disk2/yasheng/simulation2/SLM/lassosum/esteff_block5_her0.5_cross8_dist1_ps0.001.profile")

summ <- data.frame(fread(opt$summ))
if (length(unique(summ[, 3])) != length(summ[, 3])){
  freq <- data.frame(table(summ[, 3]))
  pos_nouni <- as.character(freq[freq[, 2] > 1, 1])
  summ <- summ[-which(summ[, 3] %in% pos_nouni), ]
}

#p to cor
n <- summ[1, 4] + summ[1, 5]
setwd(system.file("data", package="lassosum"))
LDblocks <- "EUR.hg19"
cor <- p2cor(p = summ[, 11], n = n, sign = summ[, 9])
idx_cor <- which(is.na(cor))
if (length(idx_cor) != 0){
  cor <- cor[-idx_cor]
  summ <- summ[-idx_cor, ]
}
out <- lassosum.pipeline(cor = cor, chr = summ[, 1], pos = summ[, 3], 
                         A1 = summ[, 6], A2 = summ[, 7],
                         ref.bfile = opt$ref, test.bfile = opt$test, 
                         LDblocks = LDblocks)
beta <- validate(out, plot = F)$best.beta
idx1 <- which(abs(beta - 0.0) < 1e-20)
inter_summ <- summ[summ[, 3] %in% out$sumstats$pos, ]
idx2 <- which(inter_summ[, 8] == 0)
idx <- unique(c(idx1, idx2))

beta_dat <- data.frame(inter_summ[, 2], inter_summ[, 6], beta)[-idx, ]
write.table(beta_dat, file = opt$beta, col.names = F, row.names = F, quote = F)