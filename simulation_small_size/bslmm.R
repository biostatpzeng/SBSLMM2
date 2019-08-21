## LMM setting
library(data.table)
library(plyr)
library(optparse)

## parameter setting
args_list = list(
  make_option("--bim", type="character", default=NULL,
              help="INPUT: bim file", metavar="numeric"),
  make_option("--eff", type="character", default=NULL,
              help="INPUT: effect of bslmm", metavar="numeric"),
  make_option("--effc", type="character", default=NULL,
              help="OUTPUT: result of bslmm", metavar="numeric")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## load data
bim_file  <- data.frame(fread(opt$bim))
eff_file <- data.frame(fread(opt$eff))

## intersect
inter_snp <- intersect(bim_file[, 2], eff_file[, 2])
bim_inter <- bim_file[match(inter_snp, bim_file[, 2]), c(2, 5)]
eff_inter <- eff_file[match(inter_snp, eff_file[, 2]), c(2, 5, 6, 7)]

## merge
effc_file <- data.frame(SNP = bim_inter[, 1], 
                        allele = bim_inter[, 2], 
                        effct = eff_inter[, 2] + eff_inter[, 3] * eff_inter[, 4])
write.table(effc_file, file = opt$effc, row.names = F, col.names = F, quote = F)