#!/bin/bash

#SBATCH --partition=nomosix,mulan
#SBATCH --time=1-00:00:00
#SBATCH --job-name=summ
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

#SBATCH --array=1-44
#SBATCH --output=/net/mulan/disk2/yasheng/out/summ%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/summ%a.err

bash
let k=0
gemma=/net/mulan/home/yasheng/Biobank/program/GEMMA/gemma-0.98-linux-static

for p in 13 16;do
for((chr=1;chr<23;chr++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
let antichr=23-${chr}
bfile=/net/mulan/disk2/yasheng/plink_file/genotype/chr${antichr}
summ=summary_pheno${p}_chr${antichr}
cd /net/mulan/disk2/yasheng/
${gemma} -bfile ${bfile} -notsnp -lm 1 -n ${p} -o ${summ}
sed -i'1d'/net/mulan/disk2/yasheng/output/${summ}.assoc.txt
rm /net/mulan/disk2/yasheng/output/${summ}.log.txt
fi
done
done
