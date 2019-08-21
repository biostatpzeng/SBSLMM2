#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=sblup
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

#SBATCH --array=1-1540%30
#SBATCH --output=/net/mulan/disk2/yasheng/out/sblup_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/sblup_%a.err

bash
let k=0
let ldr=200
let s=200

gcta=/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64

for ((PHENO=3;PHENO<17;PHENO++));do
for cross in 1 2 3 4 5; do
for ((chr=1;chr<23;chr++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
herit=/net/mulan/disk2/yasheng/pheno${PHENO}/heritability/h2_cross${cross}_ukb.log

## heritability
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`

## snp number  
merge=/net/mulan/disk2/yasheng/sample${s}/frq/merge
m=`cat ${merge}.bim | wc -l`
cojo=$(echo "${m}*(1/${h2}-1)" | bc -l)

## est
summ=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}_chr${chr}
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summ}.assoc.txt > ${summ}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summ}.ma
ref=/net/mulan/disk2/yasheng/sample${s}/frq/chr${chr}
est=/net/mulan/disk2/yasheng/pheno${PHENO}/sblup_ukb/esteff_cross${cross}_ldr${ldr}_chr${chr}
${gcta} --bfile ${ref} --chr ${chr} --cojo-file ${summ}.ma --cojo-sblup ${cojo} --cojo-wind ${ldr} --thread-num 8 --out ${est} 

## remove file
rm ${est}.log
rm ${est}.*badsnps
rm ${pheno}.log
rm ${summ}.ma
fi
done
done
done
