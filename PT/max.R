#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--r2", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--pbest1", type="character", default=NULL,
              help="INPUT: plink file", metavar="character"),
  make_option("--pbest2", type="character", default=NULL,
              help="INPUT: plink file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## p threshold
pth1 <- c("5e-8", "1e-6", "1e-4", "1e-3", "1e-2", "5e-2", "1e-1", "2e-1", "5e-1", "1.0")
r2_vec1 <- vector("numeric", length(pth1))
for (i in 1: length(pth1)){
  r2_tmp <- try(read.table(paste0(opt$r2, pth1[i], ".txt"))[1, 1], silent = T)
  if (inherits(r2_tmp, "try-error") == F){
    r2_vec1[i] <- r2_tmp
  }
}
pth2 <- c("1e-8", "5e-8", "1e-7", "5e-7", "1e-6")
r2_vec2 <- vector("numeric", length(pth2))
for (i in 1: length(pth2)){
  r2_tmp <- try(read.table(paste0(opt$r2, pth2[i], ".txt"))[1, 1], silent = T)
  if (inherits(r2_tmp, "try-error") == F){
    r2_vec2[i] <- r2_tmp
  }
}
write.table(pth1[which.max(r2_vec1)], file = opt$pbest1, 
            row.names = F, col.names = F, quote = F)
write.table(pth2[which.max(r2_vec2)], file = opt$pbest2, 
            row.names = F, col.names = F, quote = F)