#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=PT
#SBATCH --mem=2G

#SBATCH --array=1-5
#SBATCH --output=/net/mulan/disk2/yasheng/out/PT_pheno_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/PT_pheno_%a.err
bash
let k=0

toplinkf=/net/mulan/home/yasheng/Biobank/code/PT/toplinkf.R
clumpf=/net/mulan/home/yasheng/Biobank/code/PT/toclumpf.R
sumpred=/net/mulan/home/yasheng/Biobank/code/PT/sumPred.R
max=/net/mulan/home/yasheng/Biobank/code/PT/max.R
res=/net/mulan/home/yasheng/Biobank/code/PT/res.R

let PHENO=1
let s=200
for cross in 1 2 3 4 5;do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

###################################
#### find best p value 
###################################
#### use subset
## to plink format
summf=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}
# for ((chr=1;chr<23;chr++));do
# Rscript ${toplinkf} --gemma ${summf}_chr${chr}.assoc.txt --plink ${summf}_chr${chr}.plink.txt
# done

### different cutoff values
for p in 5e-8 1e-6 1e-4 1e-3 1e-2 5e-2 1e-1 2e-1 5e-1 1.0; do 
for ((chr=1;chr<23;chr++)); do
## clumping
clump=/net/mulan/disk2/yasheng/sample${s}/PT/summary_cross${cross}_chr${chr}
ref=/net/mulan/disk2/yasheng/plink_file/subsample/frq/chr${chr}
pred=/net/mulan/disk2/yasheng/sample${s}/PT/pheno_cross${cross}
Rscript ${clumpf} --gemma ${summf}_chr${chr}.assoc.txt --plink ${summf}_chr${chr}.plink.txt --ref ${ref} --pth ${p} --clump ${clump} --pred ${pred}_chr${chr}
rm ${clump}*
rm ${pred}_chr${chr}.log
rm ${pred}_chr${chr}.nopred
done

## predict for each threshold
r2=/net/mulan/disk2/yasheng/sample${s}/PT/r2_cross${cross}
Rscript ${sumpred} --pred ${pred}_chr --pheno ${PHENO} --r2 ${r2}_${p}.txt
rm ${pred}_chr*.profile
done

## best p value
pbestf=/net/mulan/disk2/yasheng/sample${s}/PT/pbest_cross${cross}
Rscript ${max} --r2 ${r2} --pbest ${pbestf}.txt
pbest=`cat ${pbestf}.txt`
rm ${r2}_*.txt

## output SNP effect
for ((chr=1;chr<23;chr++)); do
clump=/net/mulan/disk2/yasheng/sample${s}/PT/summary_cross${cross}_chr${chr}
ref=/net/mulan/disk2/yasheng/sample${s}/frq/chr${chr}
eff=/net/mulan/disk2/yasheng/sample${s}/PT/esteff_cross${cross}_chr${chr}.txt
bfile=/net/mulan/disk2/yasheng/plink_file/genotype/chr${chr}
pred=/net/mulan/disk2/yasheng/sample${s}/PT/pheno_cross${cross}_chr${chr}
Rscript ${res} --gemma ${summf}_chr${chr}.assoc.txt --plink ${summf}_chr${chr}.plink.txt --ref ${ref} --pth ${pbest} --clump ${clump} --eff ${eff}
plink-1.9 --bfile ${bfile} --score ${eff} 2 6 9 sum --out ${pred} 
rm ${clump}*
rm ${pred}.log
done
fi
done

