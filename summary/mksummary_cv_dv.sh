#!/bin/bash

#SBATCH --partition=nomosix,mulan
#SBATCH --time=1-00:00:00
#SBATCH --job-name=summ
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

#SBATCH --array=1
#SBATCH --output=/net/mulan/disk2/yasheng/out/summ%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/summ%a.err

bash
let k=0
gemma=/net/mulan/home/yasheng/Biobank/program/GEMMA/gemma-0.98-linux-static

for p in 15; do
for chr in 12; do
	for cross in 4; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            let antichr=23-${chr}
			let col=(${p}-1)*5+${cross}
            bfile=/net/mulan/disk2/yasheng/tot_plink_file/chr${antichr}
            summ=summary_crosss${cross}_chr${antichr}
            cd /net/mulan/disk2/yasheng/pheno${p}
            ${gemma} -bfile ${bfile} -notsnp -lm 1 -n ${col} -o ${summ}
			sed -i '1d' /net/mulan/disk2/yasheng/pheno${p}/output/${summ}.assoc.txt
			rm /net/mulan/disk2/yasheng/pheno${p}/output/${summ}.log.txt
        fi
	done
done
done
