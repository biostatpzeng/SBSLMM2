#!/bin/bash

#SBATCH --partition=mulan
#SBATCH --time=10-00:00:00
#SBATCH --job-name=ldpred
#SBATCH --mem=40G

#SBATCH --array=1-2
#SBATCH --output=/net/mulan/disk2/yasheng/out/ldpred500_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/ldpred500_%a.err

toldpred=/net/mulan/home/yasheng/Biobank/code/ldpred/toldpred.R
py=/net/mulan/home/yasheng/py3/bin/python
ldpred=/net/mulan/home/yasheng/Biobank/program/ldpred/LDpred.py
split=/net/mulan/home/yasheng/Biobank/code/ldpred/split_res2.R
calcr2=/net/mulan/home/yasheng/Biobank/code/ldpred/r2.R
max=/net/mulan/home/yasheng/Biobank/code/ldpred/max.R

bash
let k=0

let s=500
let ldr=200

for PHENO in 14;do
for cross in 1 2; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
## summary 
summ=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}
cat ${summ}_chr*.assoc.txt > ${summ}.assoc.txt
Rscript ${toldpred} --gemma ${summ}.assoc.txt  --ldpred ${summ}_LDpred.sumstat

## h2
herit=/net/mulan/disk2/yasheng/pheno${PHENO}/heritability/h2_cross${cross}_ukb.log
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`
echo herit: ${h2}

## ldpred (cross validation)
n_file=/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_${PHENO}_n_train.txt
n=`sed -n "${cross}, 1p" ${n_file}`
sub=/net/mulan/disk2/yasheng/plink_file/subsample/frq/merge
coord1=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/summary_cv_cross${cross}.HDF5
${py} ${ldpred} coord --gf ${sub} --ssf ${summ}_LDpred.sumstat --out ${coord1} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
ldest=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/ld_cross${cross}_ldr${ldr}
infest=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/esteff_cross${cross}_ldr${ldr}
${py} ${ldpred} gibbs --cf ${coord1} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${h2} --f 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 

## cross validation
pred=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/pheno_cross${cross}
r2=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/r2_cross${cross}
for p in 1.0000e+00 3.0000e-01 1.0000e-01 3.0000e-02 1.0000e-02 3.0000e-03 1.0000e-03 3.0000e-04 1.0000e-04;do
plink-1.9 --bfile ${sub} --score ${infest}_LDpred_p${p}.txt 3 4 7 header sum --out ${pred}_p${p}
Rscript ${calcr2} --pred ${pred}_p${p}.profile --pheno ${PHENO} --r2 ${r2}_p${p}.txt
rm ${pred}_p${p}.log
rm ${pred}_p${p}.nopred
rm ${pred}_p${p}.profile
done
pbestf=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/r2_cross${cross}_pbest.txt
Rscript ${max} --r2 ${r2}_p --pbest ${pbestf} 
pbest=`cat ${pbestf}`

## delete files
rm ${ldest}*.pkl.gz
rm ${infest}_LDpred*.txt

## ldpred (reference data)
ref2=/net/mulan/disk2/yasheng/sample${s}/frq/merge
coord2=/net/mulan/disk2/yasheng/pheno${PHENO}/ldpred_ukb/summary_ref_cross${cross}_ldr${ldr}.HDF5
${py} ${ldpred} coord --gf ${ref2} --ssf ${summ}_LDpred.sumstat --out ${coord2} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
${py} ${ldpred} gibbs --cf ${coord2} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${h2} --f ${pbest}

## split to chr
Rscript ${split} --tot ${infest}_LDpred-inf.txt --sp ${infest}_inf_chr
Rscript ${split} --tot ${infest}_LDpred_p${pbest}.txt --sp ${infest}_pbest_chr

## remove
rm ${coord2}
rm ${ldest}*.pkl.gz
rm ${infest}_LDpred*.txt
rm ${summ}.assoc.txt
rm ${summ}_LDpred.sumstat
fi
done
done
