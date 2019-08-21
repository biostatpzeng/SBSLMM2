#! /usr/bin/env Rscript
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--gemma", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"),
  make_option("--ldpred", type="character", default=NULL,
              help="OUTPUT: ldpred file", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

gemma <- data.frame(fread(opt$gemma))
pval <- ifelse(gemma[, 11] == 0, 2e-50, gemma[, 11])
ldpred <- data.frame(chr = paste0("chr", gemma[, 1]), pos = gemma[, 3],
                     ref = gemma[, 7], alt = gemma[, 6], reffrq = gemma[, 8],
                     info = 0.9, rs = gemma[, 2], pval = pval,
                     effalt = gemma[, 9])
write.table(ldpred, file = opt$ldpred, row.names = F, quote = F)