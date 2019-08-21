#!/bin/bash

#SBATCH --partition=mulan
#SBATCH --time=1-00:00:00
#SBATCH --job-name=ldpred
#SBATCH --mem=4G

#SBATCH --array=1-15
#SBATCH --output=/net/mulan/disk2/yasheng/out/ldpred_est_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/ldpred_est_%a.err

toldpred=/net/mulan/home/yasheng/Biobank/code/ldpred/toldpred.R
py=/net/mulan/home/yasheng/py3/bin/python
ldpred=/net/mulan/home/yasheng/Biobank/program/ldpred/LDpred.py
split=/net/mulan/home/yasheng/Biobank/code/ldpred/split_res2.R
calcr2=/net/mulan/home/yasheng/Biobank/code/ldpred/r2.R
max=/net/mulan/home/yasheng/Biobank/code/ldpred/max.R
PLINK=/net/mulan/home/yasheng/plink2_linux_x86_64/plink2

bash
let k=0

let ldr=200
let s=500
for ((PHENO=2;PHENO<17;PHENO++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
for cross in 1 2 3 4 5 ; do

###################################
#### find best p value 
###################################
#### use subset
## to ldpred format
summ=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}
cat ${summ}_chr*.assoc.txt > ${summ}.assoc.txt
Rscript ${toldpred} --gemma ${summ}.assoc.txt  --ldpred ${summ}_LDpred.sumstat

## coord data
n_file=/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_${PHENO}_n_train.txt
n=`sed -n "${cross}, 1p" ${n_file}`
ref1=/net/mulan/disk2/yasheng/plink_file/subsample/frq/merge
coord1=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/summary_cv_cross${cross}.HDF5
${py} ${ldpred} coord --gf ${ref1} --ssf ${summ}_LDpred.sumstat --out ${coord1} --N ${n} --ssf-format STANDARD

## ldpred
ldest=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/ld_cross${cross}_ldr${ldr}
infest=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/esteff_cross${cross}_ldr${ldr}
${py} ${ldpred} gibbs --cf ${coord1} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --f 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 

## get best p value
ref1=/net/mulan/disk2/yasheng/plink_file/subsample/frq/merge
pred=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/pheno_cross${cross}_ldr${ldr}
r2=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/r2_cross${cross}_ldr${ldr}
for p in 1.0000e+00 3.0000e-01 1.0000e-01 3.0000e-02 1.0000e-02 3.0000e-03 1.0000e-03 3.0000e-04 1.0000e-04;do
plink-1.9 --bfile ${ref1} --score ${infest}_LDpred_p${p}.txt 3 4 7 header sum --out ${pred}_p${p}
Rscript ${calcr2} --pred ${pred}_p${p}.profile --pheno ${PHENO} --r2 ${r2}_p${p}.txt
rm ${pred}_p${p}.log
rm ${pred}_p${p}.nopred
rm ${pred}_p${p}.profile
done
pbestf=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/r2_cross${cross}_ldr${ldr}_best.txt
Rscript ${max} --r2 ${r2}_p --pbest ${pbestf} 
pbest=`cat ${pbestf}`

## remove
rm ${ldest}*.pkl.gz
rm ${infest}_LDpred*.txt

###################################
#### use best p value to precition
###################################
#### use extra sample
## ldpred
ref2=/net/mulan/disk2/yasheng/sample${s}/frq/merge
coord2=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/summary_ref_cross${cross}_ldr${ldr}.HDF5
${py} ${ldpred} coord --gf ${ref2} --ssf ${summ}_LDpred.sumstat --out ${coord2} --N ${n} --ssf-format STANDARD
${py} ${ldpred} gibbs --cf ${coord2} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --f ${pbest}

## split to chr
Rscript ${split} --tot ${infest}_LDpred-inf.txt --sp ${infest}_inf_chr
Rscript ${split} --tot ${infest}_LDpred_p${pbest}.txt --sp ${infest}_pbest_chr

## remove
rm ${coord2}
rm ${ldest}*.pkl.gz
rm ${infest}_LDpred*.txt

done
fi
done
