#! /usr/bin/env Rscript
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--gemma", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"),
  make_option("--plink", type="character", default=NULL,
              help="OUTPUT: plink file", metavar="character") 
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## gemma to plink file
gemma <- data.frame(fread(opt$gemma, header = F))
p <- ifelse(gemma[, 11] == 0, 1e-20, gemma[, 11])
plink <- data.frame(CHR = gemma[, 1], SNP = gemma[, 2], BP = gemma[, 3],
                    NMISS = gemma[, 4], BETA = gemma[, 9], SE = gemma[, 10],
                    R2 = 0, T = gemma[, 9] / gemma[, 10], P = gemma[, 11], BETAS=gemma[, 9])
write.table(plink, file = opt$plink, row.names = F, quote = F)