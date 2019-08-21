#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary data prefix", metavar="character"), 
  make_option("--ref", type="character", default=NULL,
              help="INPUT: reference file", metavar="character"), 
  make_option("--prop", type="character", default=NULL,
              help="INPUT: prop of fix effect", metavar="character"),
  make_option("--r2", type="character", default=NULL,
              help="INPUT: r2", metavar="character"),  
  make_option("--p", type="character", default=NULL,
              help="INPUT: p threshold", metavar="character"),  
  make_option("--effl", type="character", default=NULL,
              help="OUTPUT: large effect file", metavar="character"),
  make_option("--effs", type="character", default=NULL,
              help="OUTPUT: small effect file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/net/mulan/disk2/yasheng/pheno1/summ/summary_cross1_chr1",
#             ref = "/net/mulan/disk2/yasheng/sample500/frq/chr1",
#             prop = 0.0001,
#             effl = "/net/mulan/disk2/yasheng/pheno1/slmm_ukb/l_summary_cross1_chr1",
#             effs = "/net/mulan/disk2/yasheng/pheno1/slmm_ukb/s_summary_cross1_chr1")

## clumping significant SNPs
plink_str <- paste0(opt$summ, ".plink.txt")
system(paste0("plink-1.9 --bfile ", opt$ref, " --clump ", plink_str,
              " --clump-kb 1000 --clump-r2 ", opt$r2, " --clump-p1 ", opt$p, 
              " --clump-p2 ", opt$p, " --out ", opt$effl))
# ## independent SNPs 
# plink_file <- data.frame(fread(plink_str))
# snp_clump <- data.frame(fread(paste0(opt$effl, ".clumped")))[, 3]
# write.table(plink_file[plink_file[, 2] %in% snp_clump, ], 
#             file = paste0(opt$effl, "_clump.txt"), 
#             row.names = F, col.names = F, quote = F, sep = "\t")
# system(paste0("plink-1.9 --bfile ", opt$ref, " --clump ", plink_str,
#               " --indep-pairwise 2 2 0.5", 
#               " --out ", opt$effl))

## number of large effect SNP
summ_dat <- data.frame(fread(paste0(opt$summ, ".assoc.txt"), header = F))
clump_file <- try(data.frame(fread(paste0(opt$effl, ".clumped")))[, 3], silent = T)
if (inherits(clump_file, "try-error")){
  sig_idx <- which.min(summ_dat[, 11])
  write.table(summ_dat[sig_idx, ], file = paste0(opt$effl, ".txt"),
              row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(summ_dat[-sig_idx, ], file = paste0(opt$effs, ".txt"),
              row.names = F, col.names = F, quote = F, sep = "\t")
}else{
  snp_prop <- floor(as.numeric(opt$prop)*nrow(summ_dat))
  snp_inter <- length(clump_file)
  max_num <- ifelse(snp_prop > snp_inter, snp_inter, snp_prop)
  cat("snp_prop: ", snp_prop, ", snp_inter: ", snp_inter, "\n")
  snp_inter_max <- clump_file[c(1: max_num)]
  summ_indep <- summ_dat[summ_dat[, 2] %in% snp_inter_max, ]
  # summ_indep <- summ_indep[order(summ_indep[, 11]), ]
  # summ_indep <- summ_indep[c(1: max_num), ]
  summ_indep <- summ_indep[order(summ_indep[, 3]), ]
  snp_indep <- summ_indep[, 2]
  write.table(summ_indep, file = paste0(opt$effl, ".txt"),
              row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(summ_dat[!summ_dat[, 2] %in% snp_indep, ], file = paste0(opt$effs, ".txt"),
              row.names = F, col.names = F, quote = F, sep = "\t")
}


system(paste0("rm ", opt$effl, ".log"))
system(paste0("rm ", opt$effl, ".clumped"))
