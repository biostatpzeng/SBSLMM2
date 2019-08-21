#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=2-00:00:00
#SBATCH --job-name=sblup
#SBATCH --mem=24G
#SBATCH --cpus-per-task=6

#SBATCH --array=1-2
#SBATCH --output=/net/mulan/disk2/yasheng/out/sblup_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/sblup_%a.err

bash
let k=0
let PHENO=1
let s=500

gcta=/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64
plink2=/net/mulan/home/yasheng/plink2_linux_x86_64/plink2

for cross in 1 2 3 4 5; do

herit=/net/mulan/disk2/yasheng/sample${s}/heritability/h2_cross${cross}.log

## heritability
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`

## snp number  
merge=/net/mulan/disk2/yasheng/sample${s}/frq/merge
m=`cat ${merge}.bim | wc -l`
cojo=$(echo "${m}*(1/${h2}-1)" | bc -l)

idx=/net/mulan/disk2/yasheng/phenotype_file/testIdx/phenoIdx_${PHENO}_cross${cross}.txt

for ((chr=1;chr<23;chr++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
## est
summ=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}_chr${chr}
ref=/net/mulan/disk2/yasheng/sample${s}/frq/chr${chr}
est=/net/mulan/disk2/yasheng/sample${s}/sblup/aesteff_cross${cross}_chr${chr} 
${gcta} --bfile ${ref} --chr ${chr} --cojo-file ${summ}.ma --cojo-sblup ${cojo} --cojo-wind 200 --thread-num 12 --out ${est} 

## valid
bfile=/net/mulan/disk2/yasheng/plink_file/genotype/chr${chr}
pheno=/net/mulan/disk2/yasheng/sample${s}/sblup/apheno_cross${cross}_chr${chr}
plink-1.9 --bfile ${bfile} --score ${est}.sblup.cojo 1 2 4 sum --keep ${idx} --out ${pheno}
## remove file
rm ${est}.log
rm ${est}*badsnps
# # # rm ${est}.sblup.cojo
rm ${pheno}.log
fi
done
done
