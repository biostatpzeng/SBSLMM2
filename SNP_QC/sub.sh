#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=sub
#SBATCH --mem=2G

#SBATCH --array=1-22
#SBATCH --output=/net/mulan/disk2/yasheng/out/sub_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/sub_%a.err

bash
let k=0

subIdx=/net/mulan/disk2/yasheng/plink_file/subsample/idx_sub.txt

for ((chr=1;chr<23;chr++));do
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
    bfile1=/net/mulan/disk2/yasheng/plink_file/genotype/chr${chr}
	bfile2=/net/mulan/disk2/yasheng/plink_file/subsample/frq/chr${chr}
	## change the allele order
	plink-1.9 --bfile ${bfile1} --keep ${subIdx} --make-bed --out ${bfile2}

	fi
done
