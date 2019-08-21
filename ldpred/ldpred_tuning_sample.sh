#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=ldpred
#SBATCH --mem=30G

#SBATCH --array=1-4
#SBATCH --output=/net/mulan/disk2/yasheng/out/ldpred_sample_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/ldpred_sample_%a.err

toldpred=/net/mulan/home/yasheng/Biobank/code/ldpred/toldpred.R
py=/net/mulan/home/yasheng/py3/bin/python
ldpred=/net/mulan/home/yasheng/Biobank/program/ldpred/LDpred.py
split=/net/mulan/home/yasheng/Biobank/code/ldpred/split_res2.R
calcr2=/net/mulan/home/yasheng/Biobank/code/ldpred/r2.R
max=/net/mulan/home/yasheng/Biobank/code/ldpred/max.R

# bash
# let k=0

let PHENO=1
for s in 1000;do
for cross in 4; do
# let k=${k}+1
# if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
summ=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}
herit=/net/mulan/disk2/yasheng/sample${s}/heritability/h2_cross${cross}.log
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`
n_file=/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_${PHENO}_n_train.txt
n=`sed -n "${cross}, 1p" ${n_file}`

ref=/net/mulan/disk2/yasheng/EUR/genotype/merge
coord=/net/mulan/disk2/yasheng/sample${s}/ldpred/summary_ref_cross${cross}.HDF5
ldest=/net/mulan/disk2/yasheng/sample${s}/ldpred/ld_cross${cross}
infest=/net/mulan/disk2/yasheng/sample${s}/ldpred/esteff_cross${cross}
${py} ${ldpred} coord --gf ${ref} --ssf ${summ}_LDpred.sumstat --out ${coord} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
${py} ${ldpred} gibbs --cf ${coord} --ldr 200 --ldf ${ldest} --out ${infest} --N ${n} --h2 ${h2} --f 1.0

Rscript ${split} --tot ${infest}_LDpred-inf.txt --sp ${infest}_inf_chr
Rscript ${split} --tot ${infest}_LDpred_p1.0000e+00.txt --sp ${infest}_pbest_chr

rm ${coord}
rm ${ldest}*.pkl.gz
# rm ${infest}_LDpred*.txt
# fi
done
done

