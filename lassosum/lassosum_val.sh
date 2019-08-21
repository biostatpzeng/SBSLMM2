#!/bin/bash

#SBATCH --partition=nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=lassosum_p
#SBATCH --mem=3G

#SBATCH --array=1-11
#SBATCH --output=/net/mulan/disk2/yasheng/out/lassosum_val1_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/lassosum_val1_%a.err

bash 
let k=0
plink2=/net/mulan/home/yasheng/plink2_linux_x86_64/plink2
for PHENO in 1; do
for chr in 1 ;do
# let k=${k}+1
# if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
for cross in 2 3 4 5; do
idx=/net/mulan/disk2/yasheng/phenotype_file/testIdx/phenoIdx_${PHENO}_cross${cross}.txt
est=/net/mulan/disk2/yasheng/pheno${PHENO}/lassosum_ukb/esteff_cross${cross}
bfile=/net/mulan/disk2/yasheng/plink_file/genotype/chr${chr}
phenoPred=/net/mulan/disk2/yasheng/pheno${PHENO}/lassosum_ukb/pheno_cross${cross}_chr${chr}
plink-1.9 --bfile ${bfile} --score ${est}_chr${chr}.txt 1 2 3 sum --keep ${idx} --out ${phenoPred}
rm ${phenoPred}.log
rm ${phenoPred}.nopred
done
# fi
done
done
