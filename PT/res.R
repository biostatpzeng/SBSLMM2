#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--gemma", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--plink", type="character", default=NULL,
              help="INPUT: plink file", metavar="character"), 
  make_option("--ref", type="character", default=NULL,
              help="INPUT: ref file", metavar="character"),
  make_option("--pth", type="character", default=NULL,
              help="INPUT: p threshold", metavar="character"),
  make_option("--clump", type="character", default=NULL,
              help="OUTPUT: clump file", metavar="character"), 
  make_option("--eff", type="character", default=NULL,
              help="OUTPUT: clump file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## clumping
system(paste0("plink-1.9 --bfile ", opt$ref, " --clump ", opt$plink, 
              " --clump-kb 1000 --clump-r2 0.1 --clump-p1 ", opt$pth, 
              " --clump-p2 ", opt$pth, " --out ", opt$clump))
## gemma
gemma <- data.frame(fread(opt$gemma, header = F))

clump_file <- data.frame(fread(paste0(opt$clump, ".clumped")))
sel_eff <- gemma[gemma[, 2] %in% clump_file[, 3], ]
write.table(sel_eff, file = opt$eff, 
            col.names = F, row.names = F, quote = F, sep = "\t")
