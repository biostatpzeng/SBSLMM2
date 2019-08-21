#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=PT_t
#SBATCH --mem=2G

#SBATCH --array=1-5
#SBATCH --output=/net/mulan/disk2/yasheng/out/PT_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/PT_%a.err
bash
let k=0

toplinkf=/net/mulan/home/yasheng/Biobank/code/PT/toplinkf.R
clumpf=/net/mulan/home/yasheng/Biobank/code/PT/toclumpf.R
sumpred=/net/mulan/home/yasheng/Biobank/code/PT/sumPred.R
max=/net/mulan/home/yasheng/Biobank/code/PT/max.R
res=/net/mulan/home/yasheng/Biobank/code/PT/res.R

for PHENO in 16 ;do
for cross in 1 2 3 4 5 ;do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

###################################
#### find best p value 
###################################
#### use subset
## to plink format
summf=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}
for ((chr=1;chr<23;chr++));do
Rscript ${toplinkf} --gemma ${summf}_chr${chr}.assoc.txt --plink ${summf}_chr${chr}.plink.txt
done

### different cutoff values
clump=/net/mulan/disk2/yasheng/pheno${PHENO}/PT/summary_cross${cross}
pred=/net/mulan/disk2/yasheng/pheno${PHENO}/PT/pheno_cross${cross}
r2=/net/mulan/disk2/yasheng/pheno${PHENO}/PT/r2_cross${cross}
for p in 5e-8 1e-6 1e-4 1e-3 1e-2 5e-2 1e-1 2e-1 5e-1 1.0 ; do 
for ((chr=1;chr<23;chr++)); do
## clumping
ref=/net/mulan/disk2/yasheng/plink_file/subsample/frq/chr${chr}
Rscript ${clumpf} --gemma ${summf}_chr${chr}.assoc.txt --plink ${summf}_chr${chr}.plink.txt --ref ${ref} --pth ${p} \
--clump ${clump}_p${p}_chr${chr} --pred ${pred}_p${p}_chr${chr}
rm ${clump}_p${p}_chr${chr}.clumped
rm ${clump}_p${p}_chr${chr}.log
rm ${pred}_chr${chr}.log
rm ${pred}_chr${chr}.nopred
done

## predict for each threshold
Rscript ${sumpred} --pred ${pred}_p${p}_chr --pheno ${PHENO} --r2 ${r2}_p${p}.txt
rm ${pred}_p${p}_chr*.profile
rm ${pred}_p${p}_chr*.log
done

# best p value
pbestf=/net/mulan/disk2/yasheng/pheno${PHENO}/PT/pbest_cross${cross}
Rscript ${max} --r2 ${r2}_p --pbest ${pbestf}.txt
pbest=`cat ${pbestf}.txt`

## output SNP effect
idx=/net/mulan/disk2/yasheng/phenotype_file/testIdx/phenoIdx_${PHENO}_cross${cross}.txt
for ((chr=1;chr<23;chr++)); do
bfile=/net/mulan/disk2/yasheng/plink_file/genotype/chr${chr}
mv ${clump}_p${pbest}_chr${chr}.txt ${clump}_best_chr${chr}.txt
plink-1.9 --bfile ${bfile} --score ${clump}_best_chr${chr}.txt 1 2 3 sum --keep ${idx} --out ${pred}_chr${chr}
rm ${pred}_chr${chr}.log
done
rm ${clump}_p*
fi
done
done
