#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--pheno", type="character", default=NULL,
              help="INPUT: phenotype file from plink", metavar="character"),
  make_option("--r2", type="character", default=NULL,
              help="OUTPUT: r2 file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

pred_pheno <- data.frame(fread(opt$pheno, header = T))
r2 <- cor(pred_pheno[, 3], pred_pheno[, 6])^2
write.table(r2, file = opt$r2, col.names = F, row.names = F, quote = F)