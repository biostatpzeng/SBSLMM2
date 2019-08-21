#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--pred", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--pheno", type="character", default=NULL,
              help="INPUT: phenotype number", metavar="character"), 
  make_option("--r2", type="character", default=NULL,
              help="OUTPUT: r2 file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

pred_pheno <- data.frame(fread(opt$pred, header = T))[, 6]
pheno <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/subsample/pheno_sub_v2.txt"))[, as.numeric(opt$pheno)]
r2 <- cor(pred_pheno, pheno)^2
write.table(r2, file = opt$r2, col.names = F, row.names = F, quote = F)
