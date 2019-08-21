library(data.table)
## input phenotype
# pheno <- data.frame(fread("/net/mulan/disk2/yasheng/phenotype_file/pheno_tot_v3.txt"))
# sampleIdx <- read.table("/net/mulan/disk2/yasheng/phenotype_file/idx_tot.txt")
# freq <- apply(pheno, 1, function(a) sum(a != -9))
# set.seed(2019)
# idx <- sort(sample(which(freq == 16), 1000))
# 
# ## output index and phenotype
# write.table(pheno[idx, ], file = "/net/mulan/disk2/yasheng/plink_file/subsample/pheno_sub_v2.txt", 
#             col.names = F, row.names = F, quote = F)
# write.table(sampleIdx[idx, ], file = "/net/mulan/disk2/yasheng/plink_file/subsample/idx_sub.txt", 
#             col.names = F, row.names = F, quote = F)

pheno <- data.frame(fread("/net/mulan/disk2/yasheng/phenotype_file/pheno_tot_v2.txt"))
sampleIdx <- read.table("/net/mulan/disk2/yasheng/phenotype_file/idx_tot.txt")
sampleSub <- read.table("/net/mulan/disk2/yasheng/plink_file/subsample/idx_sub.txt")
idx <- which(sampleIdx[, 1] %in% sampleSub[, 1])

write.table(pheno[idx, ], file = "/net/mulan/disk2/yasheng/plink_file/subsample/pheno_sub_v2.txtx",
            col.names = F, row.names = F, quote = F)

## idx of na and subset
idx_na <- apply(pheno, 2, function(a) ifelse(a == -9, 1, 0))
idx_sub <- ifelse(sampleIdx[, 1] %in% sampleIdx[idx, 1], 1, 0)
pheno[idx, ] <- -9

fold <- 5
seed_str <- c(522, 529, 1111, 723, 93)
for (i in 13){
  
  pheno_train <- matrix(-9, nrow = nrow(pheno), ncol = fold)
  pheno_test <- matrix(-9, nrow = nrow(pheno), ncol = fold)
  cnd <- idx_na[, i] == 0 & idx_sub == 0
  n_train <- vector("numeric", fold)
  for(j in 1: fold){
    label_cross <- vector("numeric", nrow(pheno))
    set.seed(seed_str[j])
    label_cross[cnd == 1] <- rbinom(sampleIdx[cnd, 1], 1, c(1 - 1/fold)) 
    label_cross[label_cross == 0 & pheno[, i] != -9] <- -1
    pheno_train[label_cross == 1, j] <- pheno[label_cross == 1, i]
    pheno_test[label_cross == -1, j] <- pheno[label_cross == -1, i]
    n_train[j] <- sum(label_cross == 1)
    n <- sum(label_cross == 1) + sum(label_cross == -1)
    cat ("pheno: ", i, " sample size: ", n, "\n")
    testIdx <- sampleIdx[label_cross == -1, ]
    write.table(testIdx, file = paste0("/net/mulan/disk2/yasheng/phenotype_file/testIdx/phenoIdx_", i, "_cross", j, ".txtx"),
                quote = F, row.names = F, col.names = F)
  }
  write.table(n_train, file = paste0("/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_", i, "_n_train.txtx"),
              quote = F, row.names = F, col.names = F)
  write.table(pheno_train, file = paste0("/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_", i, "_train.txtx"),
              quote = F, row.names = F, col.names = F)
  write.table(pheno_test, file = paste0("/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_", i, "_test.txtx"),
              quote = F, row.names = F, col.names = F)
} 

