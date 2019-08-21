#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)
library(lassosum)
library(parallel)

## Parameter setting
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary data prefix", metavar="character"), 
  make_option("--ref", type="character", default=NULL,
              help="INPUT: reference file prefix", metavar="character"), 
  make_option("--valid", type="character", default=NULL,
              help="INPUT: max number prefix", metavar="character"),
  make_option("--pheno", type="character", default=NULL,
              help="INPUT: phenotype number", metavar="character"),
  make_option("--thread", type="character", default=NULL,
              help="INPUT: threads", metavar="character"),
  make_option("--beta", type="character", default=NULL,
              help="OUTPUT: output of coefficient", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/net/mulan/disk2/yasheng/pheno5/summ/summary_cross1.assoc.txt",
#             ref = "/net/mulan/disk2/yasheng/sample500/frq/merge",
#             valid = "/net/mulan/disk2/yasheng/plink_file/subsample/frq/merge",
#             thread = "6",
#             pheno = "5",
#             beta = "/net/mulan/disk2/yasheng/pheno1/lassosum_ukb/esteff_cross1_chr")

summ <- data.frame(fread(opt$summ))

## sample size
n <- as.numeric(summ[1, 4]) + as.numeric(summ[1, 5])

## p to cor
pval <- ifelse(summ[, 11] == 0, 1e-100, summ[, 11]) 
t1 <- system.time(cor <- p2cor(p = pval, n = n, sign = summ[, 9]))

## block information
setwd(system.file("data", package="lassosum"))
LDblocks <- "EUR.hg19"

threads <- as.numeric(opt$thread)
cl <- makeCluster(threads, type="FORK")
t2 <- system.time(out <- lassosum.pipeline(cor = cor, chr = summ[, 1], pos = summ[, 3], 
                         A1 = summ[, 6], A2 = summ[, 7], destandardize = T,
                         ref.bfile = opt$ref, test.bfile = opt$valid,
                         LDblocks = LDblocks, cluster = cl))
ind_idx <- read.table("/net/mulan/disk2/yasheng/plink_file/subsample/idx_sub.txt")
pheno_tot <- read.table("/net/mulan/disk2/yasheng/plink_file/subsample/pheno_sub_v2.txt")
pheno_s <- data.frame(FID = ind_idx[, 1], IID = ind_idx[, 1], 
                      pheno_tot[, as.numeric(opt$pheno)]) 
t3 <- system.time(v <- validate(out, pheno = pheno_s, plot = F))
beta <- v$best.beta
idx_eff <- which(beta != 0) ##lasso zero
idx <- out$sumstats$order[idx_eff]
frq <- summ[idx, 8]
beta_dat <- data.frame(summ[idx, 1], summ[idx, 2], summ[idx, 6], 
                       beta[idx_eff])
cat("total time: ", t1[3]+t2[3]+t3[3], "\n")
beta_output <- mclapply(c(1:22), function(chr){
  beta_tmp <- beta_dat[beta_dat[, 1] == chr, ]
  write.table(beta_tmp[, -1], file = paste0(opt$beta, chr, ".txt"), 
              col.names = F, row.names = F, quote = F)
  return(chr)
},  mc.cores=threads)