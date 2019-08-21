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
  make_option("--pred", type="character", default=NULL,
              help="OUTPUT: predict file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(gemma = "/net/mulan/disk2/yasheng/pheno1/summ/summary_cross1_chr22.assoc.txt", 
#             plink = "/net/mulan/disk2/yasheng/pheno1/summ/summary_cross1_chr22.assoc.plink", 
#             clump = "/net/mulan/disk2/yasheng/pheno1/PT_2000/summary_cross1_chr22", 
#             ref = "/net/mulan/disk2/yasheng/sample2000/frq/chr22", 
#             pth = "5e-8", 
#             pred = "/net/mulan/disk2/yasheng/pheno1/PT_2000/pheno_cross1_chr22")

## clumping
system(paste0("plink-1.9 --bfile ", opt$ref, " --clump ", opt$plink, 
              " --clump-kb 1000 --clump-r2 0.1 --clump-p1 ", opt$pth, " --clump-p2 ", opt$pth, 
              " --out ", opt$clump))
## gemma
gemma <- data.frame(fread(opt$gemma, header = F))

## prediction
clump_file <- data.frame(fread(paste0(opt$clump, ".clumped")))
sel_eff <- gemma[gemma[, 2] %in% clump_file[, 3], c(2, 6, 9)]
write.table(sel_eff, file = paste0(opt$clump, ".txt"), row.names = F, quote = F, col.names = F)
system(paste0("plink-1.9 --bfile ", opt$ref, " --score ", opt$clump, ".txt", 
              " 1 2 3 sum --out ", opt$pred))
