library(data.table)
tot_idx <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/genotype/chr1.fam"))[, c(1, 2)]
set.seed(20170529)
idx <- tot_idx[sample(c(1: nrow(tot_idx))), ]

##
n_train <- 10000
n_test <- 1000
n_val <- 1000
n_sample <- dim(tot_idx)[1]

seq_train <- c(1: n_train)
seq_test <- c(c(n_train + 1): c(n_train + n_test))
seq_val <- c(c(n_sample - n_val + 1): n_sample)

setwd("/net/mulan/disk2/yasheng/plink_file/simulation")
write.table(idx[seq_train, ], file = "train_idx2.txt", col.names = F, 
            row.names = F, quote = F)
write.table(idx[seq_test, ], file = "test_idx2.txt", col.names = F, 
            row.names = F, quote = F)
write.table(idx[seq_val, ], file = "ref_idx2.txt", col.names = F, 
            row.names = F, quote = F)
write.table(idx[c(seq_train, seq_test, seq_val), ], 
            file = "idx2.txt", col.names = F, 
            row.names = F, quote = F)

