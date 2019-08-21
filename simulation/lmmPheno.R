## LMM setting
library(data.table)
library(plyr)
library(rmutil)
library(optparse)

## parameter setting
args_list = list(
  make_option("--path", type="character", default=NULL,
              help="INPUT: pathway", metavar="numeric"),
  make_option("--block", type="character", default=NULL,
              help="INPUT: block", metavar="numeric"),
  make_option("--herit", type="character", default=NULL,
              help="INPUT: heritability", metavar="numeric"),
  make_option("--cross", type="character", default=NULL,
              help="INPUT: cross", metavar="numeric"), 
  make_option("--type", type="character", default=NULL,
              help="INPUT: selected blocks", metavar="numeric"), 
  make_option("--dist", type="character", default=NULL,
              help="INPUT: distribution assumption", metavar="numeric")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(path = "/net/mulan/disk2/yasheng/simulation2/LMM/data/",
#             block = 5,
#             herit = 0.5,
#             p = 0.002,
#             prop = 0.1,
#             cross = 1,
#             type = 0,
#             dist = 2)

## seed setting for replicate
seed_str <- c(1111, 529, 93, 725, 1988, 2017, 1932, 989, 20170529, 19890725)
snp_info_str <- "/net/mulan/disk2/yasheng/plink_file/simulation/sim2"

getblock <- function(end_loc){
  loc <- ifelse(snp_info[, 4] < end_loc, 1, 0)
  return(loc)
}

##
snp_str <- paste0(opt$path, "snp_block", opt$block, "_her", opt$herit, "_cross", opt$cross, "_dist", opt$dist, ".txt")
geno_str <- paste0(opt$path, "geno_block", opt$block, "_her", opt$herit, "_cross", opt$cross, "_dist", opt$dist)
eff_str <- paste0(opt$path, "eff_block", opt$block, "_her", opt$herit, "_cross", opt$cross, "_dist", opt$dist, ".txt")
pheno_str <- paste0(opt$path, "pheno_block", opt$block, "_her", opt$herit, "_cross", opt$cross, "_dist", opt$dist)

## STEP1: select blocks in chromosome 1
block_info <- as.matrix(read.table("/net/mulan/home/yasheng/Biobank/data/LDblock/chr1.txt")[, 3])
snp_info <- data.frame(fread(paste0(snp_info_str, ".bim")))
block <- t(aaply(c(1:nrow(block_info)), 1, function(a) getblock(block_info[a, 1]))) 
block <- nrow(block_info) - rowSums(block) + 1
if (opt$type == 0){
  set.seed(seed_str[as.numeric(opt$cross)])
  block_sel <- sample(c(1:nrow(block_info)), opt$block)
}
if (opt$type == 1){
  block_sel <- c(1: as.numeric(opt$block))
}
snp_sel <- snp_info[block%in%block_sel, 2]
write.table(snp_sel, file = snp_str, col.names = F, row.names = F, quote = F)

## STEP2: build genotype and phenotype data
sel_cmd <- paste0("plink-1.9 --bfile ", snp_info_str, " --extract ", snp_str, " --make-bed --out ", geno_str)
system(sel_cmd)
set.seed(seed_str[as.numeric(opt$cross)])
if(opt$dist == "1"){
  eff_var <- as.numeric(opt$herit) / length(snp_sel)
  eff <- rnorm(length(snp_sel), 0, sqrt(eff_var))
}
if(opt$dist == "2"){
  eff <- rt(length(snp_sel), 4) / sqrt(length(snp_sel)) * as.numeric(opt$herit)
}
if(opt$dist == "3"){
  eff <- rlaplace(length(snp_sel), 0, sqrt(as.numeric(opt$herit)/(2*length(snp_sel))))
}
eff <- cbind(snp_sel, eff)
write.table(eff, file = eff_str, col.names = F, row.names = F, quote = F)

gcta_str <- "/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64"
pheno_cmd <- paste0(gcta_str, " --bfile ", geno_str, " --simu-qt --simu-causal-loci ", eff_str,
                    " --simu-hsq ", opt$herit, " --simu-rep 1 --out ", pheno_str)
system(pheno_cmd)


##
pheno_sim <- data.frame(fread(paste0(pheno_str, ".phen")))[, 3]
fam <- data.frame(fread(paste0(geno_str, ".fam")))
fam[, 6] <- pheno_sim
write.table(fam, paste0(geno_str, ".fam"), col.names = F, row.names = F, quote = F)

##
system(paste0("rm ", geno_str, ".log"))
system(paste0("rm ", pheno_str, ".log"))
system(paste0("rm ", pheno_str, ".par"))
system(paste0("rm ", pheno_str, ".phen"))
system(paste0("rm ", snp_str))
system(paste0("rm ", eff_str))
