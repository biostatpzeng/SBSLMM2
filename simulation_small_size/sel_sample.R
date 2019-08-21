setwd("/net/mulan/disk2/yasheng/plink_file/simulation")
idx <- read.table("idx.txt")

##
n_train <- 2000
n_test <- 200
n_val <- 200
n_sample <- dim(idx)[1]

seq_train <- c(1: n_train)
seq_test <- c(c(n_train + 1): c(n_train + n_test))
seq_val <- c(c(n_sample - n_val + 1): n_sample)

write.table(idx[seq_train, ], file = "train_bslmm_idx.txt", col.names = F, 
            row.names = F, quote = F)
write.table(idx[seq_test, ], file = "test_bslmm_idx.txt", col.names = F, 
            row.names = F, quote = F)
write.table(idx[seq_val, ], file = "ref_bslmm_idx.txt", col.names = F, 
            row.names = F, quote = F)
write.table(idx[c(seq_train, seq_test, seq_val), ], 
            file = "bslmm_idx.txt", col.names = F, 
            row.names = F, quote = F)

