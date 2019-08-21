###This program aims to processing dataset, including phenotye and QC. 
library(data.table)
library(plyr)
library(stringr)
library(parallel)

#########STEP1: Combine the sqc and fam
print("STEP1: Combine the sqc and fam...")
#File input
sqc <- data.frame(fread("/net/mulan/Biobank/rawdata/EGAD00010001225/001/ukb_sqc_v2.txt", 
                        stringsAsFactors = F))
header <- data.frame(read.table("/net/mulan/Biobank/rawdata/EGAD00010001225/001/ukb_sqc_v2_head.txt", 
                                stringsAsFactors = F))
fam <- data.frame(fread("/net/mulan/data/UKB/ukb30186_baf_chr9_v2_s488363.fam", 
                        stringsAsFactors = F))

#
if (dim(fam)[1] != dim(sqc)[1]){
  
  stop("ERROR: The number of rows is not the same!")
}else{
  
  sqc <- cbind.data.frame(fam[1], sqc)
  names(sqc) <- c("eid", header[, 1])
  cat("QC sample:", dim(sqc)[1], "\n")
}

#########STEP2: Select the samples

##in.Phasing.Input.chr1_22==1
##in.white.British.ancestry.subset==1
##used.in.pca.calculation==1
##excess.relatives==0
##putative.sex.chromosome.aneuploidy==0
##eid > 0
print("STEP2: Select the samples...")
##Process the QC results
#Genetyping
sqc.i <- sqc[which(sqc$in.Phasing.Input.chr1_22==1), ]
eid.i <- sqc.i[, 1]
cat("Genotyping success:", dim(sqc.i)[1], "\n")
#Others
cnd.i <- sqc.i$in.white.British.ancestry.subset == 1
cat("White British ancestry subset:", sum(cnd.i), "\n")
cnd.ii <- sqc.i$excess.relatives == 0 
cat("Excess relatives:", sum(!cnd.ii), "\n")
cnd.iii <- sqc.i$putative.sex.chromosome.aneuploidy == 0 
cat("Sex chromosome aneuploidy:", sum(!cnd.iii), "\n")
cnd.ix <- sqc.i$used.in.pca.calculation == 1
cat("Used in PCA calculation:", sum(cnd.ix), "\n")
cnd.x <- sqc.i$eid > 0
cat("Redacted:", sum(!cnd.x), "\n")

#Index of selected samples
cnd <- cnd.i & cnd.ii & cnd.iii & cnd.ix & cnd.x
cat("Samples Remaining:", sum(cnd), "\n")

## the phenotype data 
#########STEP3: phenotype data
print("STEP3: Process phenotype data...")
data <- data.frame(fread("/net/mulan/data/UKB/ukb10683.csv", header = T))
phenoList <- read.table("/net/mulan/home/yasheng/Biobank/data/pheno/nly_pheno.txt",
                        header = T, sep = "\t")
label <- paste0("X", phenoList$UKBB_code, ".0.0")
#Check the phenotype
if (all(label %in% colnames(data)) == F){

  stop("ERROR: The labels of phenotype are wrong!")
}
#Select and sort the phenotype
pheno <- data[which(data$eid %in% eid.i), c(1, match(label, colnames(data)))]
idx <- match(eid.i, pheno[, 1])
pheno <- pheno[idx, ]
pheno[which(is.na(pheno[, 1])), 1] <- c(-1:-14)

#########STEP4: regress between pheno and sqc
print("STEP4: Regress between pheno and sqc...")
#Check the size of sqc and pheno
if (all(pheno[, 1] == sqc.i[, 1]) == F){

  stop("ERROR: phenotype do not match sqc!")
}

#########
###phenotype from Alkes Price paper
#########
pheno_ap <- matrix(NA, ncol = 16, nrow = nrow(pheno))
pheno_ap[, 1] <- pheno[, 2]
pheno_ap[, 2] <- pheno[, 3]
pheno_ap[, 3] <- pheno[, 4]
pheno_ap[, 4] <- pheno[, 5]
pheno_ap[, 5] <- pheno[, 6] / pheno[, 7]
pheno_ap[, 6] <- pheno[, 8]
pheno_ap[, 7] <- pheno[, 9]
pheno_ap[, 8] <- pheno[, 10] / pheno[, 11]
pheno_ap[, 9] <- pheno[, 12]
pheno_ap[, 10] <- pheno[, 13]
pheno_ap[, 11] <- pheno[, 14]
pheno_ap[, 12] <- pheno[, 15]
pheno_ap[, 13] <- ifelse(pheno[, 16]<0, NA, pheno[, 16])
pheno_ap[, 14] <- pheno[, 17]
pheno_ap[, 15] <- pheno[, 18]
pheno_ap[, 16] <- ifelse(pheno[, 19]<0, NA, pheno[, 19])
age <- 2018 - pheno[cnd, 20]
pheno_ap <- pheno_ap[cnd, ]

###############
###sqc data
############
sqc.i <- sqc.i[cnd, ]
sqc.i[, 12] <- ifelse(sqc.i[, 12] == "M", 0, 1)

# ##################
# ###sample
# #################
sample_idx <- cbind(c(1: length(cnd)), c(1: length(cnd)))[cnd, ]
write.table(sample_idx, file = "/net/mulan/disk2/yasheng/phenotype_file/idx_tot.txt",
            row.names = F, col.names = F, quote = F)

###########
###residual and label
##########
seed_str <- c(522, 529, 1111, 723, 93)
fold <- 5
PC <- sqc.i[, c(27:46)]
sex <- sqc.i[, 12]
covVar_v2 <- as.matrix(cbind(PC[, c(1:10)], sex))
covVar_v3 <- as.matrix(cbind(PC, sex, age, age^2, sex*age, sex*age^2))

tot_pheno <- matrix(NA, nrow = 337198, ncol = 16)
for (i in 1: 16){  
  na_idx <- ifelse(is.na(pheno_ap[, i]), T, F)
  pheno_na <- ifelse(na_idx, NA, pheno_ap[, i])
  pheno_scale <- scale(pheno_na)
  resid_pheno <- vector("numeric", length(pheno_scale))
  resid <- lm(pheno_scale[!na_idx] ~ covVar_v2[!na_idx, ])$residual
  resid_pheno[!na_idx] <- qqnorm(resid, plot.it = F)$x
  resid_pheno[na_idx] <- -9
  tot_pheno[, i] <- resid_pheno
  cat(i, "\n")
}
write.table(tot_pheno, file = "/net/mulan/disk2/yasheng/phenotype_file/pheno_tot_v2.txt",
            quote = F, row.names = F, col.names = F)
write.table(pheno_ap, file = "/net/mulan/disk2/yasheng/phenotype_file/pheno_tot.txt",
            quote = F, row.names = F, col.names = F)
write.table(age, file = "/net/mulan/disk2/yasheng/phenotype_file/age.txt",
            quote = F, row.names = F, col.names = F)


# library(rbgen)
# neale <- fread("/net/mulan/home/yasheng/Biobank/data/QC/neale_pheno1_chr21.txt", 
#               header = T) %>>% {data.frame(.)}
# snp_position <- neale[c(113496: 113695), 1] %>>% 
#              {str_split(., ":")} %>>%
#              {laply(., function(a) a[2])} 
# snp_code <- neale[c(113496: 113695), 2] 
# 
# summ <- function(snp_pos, snp_rs){
#   
#   snp_info <- data.frame(chromosome = "21", start = as.numeric(snp_pos), end = as.numeric(snp_pos))
#   snp = bgen.load("/net/mulan/Biobank/rawdata/EGAD00010001225/001/ukb_imp_chr21_v2.bgen", snp_info)
#   snp_exp <- snp$data[snp_rs,,][, 2]+snp$data[snp_rs,,][, 3]+snp$data[snp_rs,,][, 3]
#   na_idx <- which(is.na(resid_pheno))
#   fit1 <- lm(resid_pheno[-na_idx] ~ snp_exp[-na_idx])
#   res1 <- coef(summary(fit1))[2, ]
#   fit2 <- lm(pheno_scale ~ snp_exp + as.matrix(sqc.i[, c(12, 27:36)]))
#   res2 <- coef(summary(fit2))[2, ]
#   dosage <- sum(snp_exp[-na_idx]*resid_pheno[-na_idx])
#   res <- c(dosage, res1, res2)
#   return(res)
# }

