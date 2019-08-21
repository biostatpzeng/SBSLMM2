#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=PT_t
#SBATCH --mem=2G

#SBATCH --array=1-12
#SBATCH --output=/net/mulan/disk2/yasheng/out/PT_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/PT_%a.err
bash
let k=0

toplinkf=/net/mulan/home/yasheng/Biobank/code/PT/toplinkf.R
clumpf=/net/mulan/home/yasheng/Biobank/code/PT/toclumpf.R
sumpred=/net/mulan/home/yasheng/Biobank/code/PT/sumPred.R
max=/net/mulan/home/yasheng/Biobank/code/PT/max.R
res=/net/mulan/home/yasheng/Biobank/code/PT/res.R

for PHENO in 8 ;do
for cross in 5 ;do

###################################
#### find best p value 
###################################
summf=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}
idx=/net/mulan/disk2/yasheng/phenotype_file/testIdx/phenoIdx_${PHENO}_cross${cross}.txt

# best p value
pbestf=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/pbest_cross${cross}
pbest=`cat ${pbestf}.txt`

### different cutoff values
clump=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/summary_cross${cross}
pred=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/pheno_cross${cross}
r2=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/r2_cross${cross}
for ((chr=11;chr<23;chr++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
## clumping
ref=/net/mulan/disk2/yasheng/plink_file/subsample/frq/chr${chr}
Rscript ${clumpf} --gemma ${summf}_chr${chr}.assoc.txt --plink ${summf}_chr${chr}.plink.txt --ref ${ref} --pth ${pbest} \
--clump ${clump}_p${pbest}_chr${chr} --pred ${pred}_p${pbest}_chr${chr}
rm ${clump}_p${pbest}_chr${chr}.clumped
rm ${clump}_p${pbest}_chr${chr}.log
rm ${pred}_p${pbest}_chr${chr}.log
rm ${pred}_p${pbest}_chr${chr}.nopred
rm ${pred}_p${pbest}_chr${chr}.profile
bfile=/net/mulan/disk2/yasheng/plink_file/genotype/chr${chr}
mv ${clump}_p${pbest}_chr${chr}.txt ${clump}_best_chr${chr}.txt
plink-1.9 --bfile ${bfile} --score ${clump}_best_chr${chr}.txt 1 2 3 sum --keep ${idx} --out ${pred}_chr${chr}
rm ${pred}_chr${chr}.log
fi
done
done
done
