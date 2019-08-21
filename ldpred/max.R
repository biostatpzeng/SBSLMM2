#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--r2", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--pbest", type="character", default=NULL,
              help="INPUT: plink file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## p threshold
pth <- c("1.0000e+00", "3.0000e-01", "1.0000e-01", "3.0000e-02", "1.0000e-02",
         "3.0000e-03", "1.0000e-03", "3.0000e-04", "1.0000e-04")
# pth <- c("3.0000e-01", "1.0000e-01")
r2_vec <- vector("numeric", length(pth))
for (i in 1: length(pth)){
  r2_vec[i] <- read.table(paste0(opt$r2, pth[i], ".txt"))[1, 1]
}
pbest_val <- pth[which(r2_vec == max(r2_vec))]
write.table(pbest_val, file = opt$pbest, 
            row.names = F, col.names = F, quote = F)
