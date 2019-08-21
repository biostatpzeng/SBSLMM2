#!/bin/bash

#SBATCH --partition=nomosix,mulan
#SBATCH --time=1-00:00:00
#SBATCH --job-name=lassosum
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=2

#SBATCH --array=1-3
#SBATCH --output=/net/mulan/disk2/yasheng/out/lassosum%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/lassosum%a.err

lassosum=/net/mulan/home/yasheng/Biobank/code/lassosum/lassosum_real.R
ref=/net/mulan/disk2/yasheng/sampleEUR/frq/merge
val=/net/mulan/disk2/yasheng/plink_file/subsample/frq/merge
thread=4
for PHENO in 1 2 6;do
for cross in 1 ; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
summ=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}
cat ${summ}_chr*.assoc.txt > ${summ}.assoc.txt
beta=/net/mulan/disk2/yasheng/pheno${PHENO}/lassosum_ukb/esteff_cross${cross}_chr
Rscript ${lassosum} --summ ${summ}.assoc.txt --ref ${ref} --valid ${val} --thread ${thread} --pheno ${PHENO} --beta ${beta}
rm ${summ}.assoc.txt
fi
done
done
