library(data.table)
library(optparse)

## Parameter setting
args_list <- list(
  make_option("--summary", type = "character", default = NULL,
              help = "INPUT: the summary statistics (gemma output format)", metavar = "character"),
  make_option("--tmpPath", type = "character", default = NULL,
              help = "INPUT: the temptorary path", metavar = "character"),
  make_option("--outPath", type="character", default=NULL,
              help="INPUT: the output path", metavar="character"),
  make_option("--plink", type = "character", default = NULL,
              help = "INPUT: the perfix of Plink software", metavar = "character"),
  make_option("--ref", type = "character", default = NULL,
              help = "INPUT: the perfix of reference panel", metavar = "character"),
  make_option("--r2", type = "character", default = NULL,
              help = "INPUT: the cutoff of SNPs clumping (default:0.1)", metavar = "character"),
  make_option("--pv", type = "character", default = NULL,
              help = "INPUT: the cutoff of SNPs pruning (default:1e-6)", metavar = "character"),
  make_option("--sbslmm", type = "character", default = NULL,
              help = "INPUT: the perfix of sbslmm software", metavar = "character"),
  make_option("--mafMax", type = "character", default = NULL,
              help = "INPUT: the maximium of the difference between reference panel and summary data", metavar = "character"),
  make_option("--nsnp", type = "character", default = NULL,
              help = "INPUT: the number of SNPs in whole genome", metavar = "character"),
  make_option("--block", type = "character", default = NULL,
              help = "INPUT: the block information (Berisa and Pickrell 2015)", metavar = "character"),
  make_option("--h2", type = "character", default = NULL,
              help = "INPUT: the heritability of trait", metavar = "character"),
  make_option("--thread", type = "character", default = NULL,
              help = "INPUT: the number of threads", metavar = "character")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- try(parse_args(opt_parser), silent = T)

## output the options
cat("Acccept Options: \n")
cat("--summary: ", opt$summary, "\n")
cat("--tmpPath: ", opt$tmpPath, "\n")
cat("--outPath: ", opt$outPath, "\n")
cat("--plink:   ", opt$plink, "\n")
r2 <- ifelse(is.null(opt$r2), 0.1, opt$r2)
cat("--r2:      ", r2, "\n")
pv <- ifelse(is.null(opt$pv), 1e-6, opt$pv)
cat("--p:       ", pv, "\n")
mafMax <- ifelse(is.null(opt$mafMax), 0.2, opt$mafMax)
cat("--mafMax:  ", mafMax, "\n")
cat("--sbslmm:  ", opt$sbslmm, "\n")
cat("--nsnp:    ", opt$nsnp, "\n")
cat("--block:   ", opt$block, "\n")
cat("--h2:      ", opt$h2, "\n")
cat("--thread:  ", opt$thread, "\n")

## check the options
if (!file.exists(opt$summary)){
  cat(paste0("ERROR: ", opt$summary, " does not exist!\n"))
  q()
}
if (!file.exists(opt$tmpPath)){
  cat(paste0("ERROR: ", opt$tmpPath, " does not exist!\n"))
  q()
}
if (!file.exists(opt$outPath)){
  cat(paste0("ERROR: ", opt$outPath, " does not exist!\n"))
  q()
}
if (!file.exists(opt$plink)){
  cat(paste0("ERROR: ", opt$plink, " does not exist!\n"))
  q()
}
if (!file.exists(opt$sbslmm)){
  cat(paste0("ERROR: ", opt$sbslmm, " does not exist!\n"))
  q()
}
if (!file.exists(opt$block)){
  cat(paste0("ERROR: ", opt$block, " does not exist!\n"))
  q()
}
## gemma format to plink format
summ_gemma <- data.frame(fread(opt$summary, header = T))
p <- ifelse(summ_gemma[, 11] == 0, 1e-100, summ_gemma[, 11])
summ_plink <- data.frame(CHR = summ_gemma[, 1], SNP = summ_gemma[, 2],
                         BP = summ_gemma[, 3], NMISS = summ_gemma[, 4],
                         BETA = summ_gemma[, 9], SE = summ_gemma[, 10],
                         R2 = 0, T = summ_gemma[, 9] / summ_gemma[, 10],
                         P = summ_gemma[, 11], BETAS=summ_gemma[, 9])
s <- strsplit(opt$summary, "\\/")[[1]]
s <- s[length(s)]
prefix_file <- strsplit(s, "\\.")[[1]][1]
write.table(summ_plink, file = paste0(opt$tmpPath, prefix_file, ".txt"), row.names = F, quote = F)

## Plink clumping
system(paste0(opt$plink, " --bfile ", opt$ref, " --clump ", opt$tmpPath, prefix_file, ".txt",
              " --clump-kb 1000 --clump-r2 ", r2, " --clump-p1 ", pv,
              " --clump-p2 ", pv, " --out ", opt$tmpPath, "l_", prefix_file))

## number of large effect SNP
snp_clump <- try(data.frame(fread(paste0(opt$tmpPath, "l_", prefix_file, ".clumped")))[, 3], silent = T)
if (inherits(snp_clump, "try-error")){
  cat ("WARNINGS: After clumping, we select no indepent large effect SNPs. We select top five SNPs.")
  idx_sig <- order(summ_gemma[, 11])[1:5]
}else{
  # maxmium of SNP number
  num_prop <- floor(0.01 * nrow(summ_gemma))
  num_clump <- length(snp_clump)
  num_max <- ifelse(num_prop > num_clump, num_clump, num_prop)
  snp_max <- snp_clump[c(1: num_max)]

  idx_sig <- which(summ_gemma[, 2] %in% snp_max)
  idx_sig <- sort(idx_sig)
}
write.table(summ_gemma[idx_sig, ], file = paste0(opt$tmpPath, "l_", prefix_file, ".txt"),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(summ_gemma[-idx_sig, ], file = paste0(opt$tmpPath, "s_", prefix_file, ".txt"),
            row.names = F, col.names = F, quote = F, sep = "\t")

## SBSLMM
sampleSize <- summ_gemma[1, 4] + summ_gemma[1, 5]
system(paste0(opt$sbslmm, " -s ", opt$tmpPath, "s_", prefix_file, ".txt",
                          " -l ", opt$tmpPath, "l_", prefix_file, ".txt",
                          " -r ", opt$ref,
                          " -n ", sampleSize,
                          " -mafMax ", mafMax,
                          " -nsnp ", opt$nsnp,
                          " -b ", opt$block,
                          " -h ", opt$h2,
                          " -t ", opt$thread,
                          " -eff ", opt$outPath, prefix_file, ".txt"))
system(paste0("rm ", opt$tmpPath))