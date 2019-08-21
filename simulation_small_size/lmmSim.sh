#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=bslmmLMM
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

#SBATCH --array=1-30%15
#SBATCH --output=/net/mulan/disk2/yasheng/out/simLMM_bslmm_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/simLMM_bslmm_%a.err

bash
let k=0

let ldr=16
let block=5
let thread=2

py=/net/mulan/home/yasheng/py3/bin/python
gemma=/net/mulan/home/yasheng/Biobank/program/GEMMA/gemma-0.98-linux-static
gcta=/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64
sbslmm=/net/mulan/home/yasheng/Biobank/code/slmm2.1/sbslmm
PLINK=/net/mulan/home/yasheng/plink2_linux_x86_64/plink2
ldsc=/net/mulan/home/yasheng/Biobank/program/ldsc/ldsc.py
ldpred=/net/mulan/home/yasheng/Biobank/program/ldpred/LDpred.py
lasso=/net/mulan/home/yasheng/Biobank/code/simulation_bslmm/biglasso.R
lassosum=/net/mulan/home/yasheng/Biobank/code/lassosum/lassosum.R

lmmPheno=/net/mulan/home/yasheng/Biobank/code/simulation_bslmm/lmmPheno.R
toplinkf=/net/mulan/home/yasheng/Biobank/code/PT/toplinkf.R
clumpf1=/net/mulan/home/yasheng/Biobank/code/PT/toclumpf.R
simpred=/net/mulan/home/yasheng/Biobank/code/PT/simPred.R
max1=/net/mulan/home/yasheng/Biobank/code/PT/max.R
mkldsc=/net/mulan/home/yasheng/Biobank/code/heritability/mkldsc.R
toldpred=/net/mulan/home/yasheng/Biobank/code/ldpred/toldpred.R
max2=/net/mulan/home/yasheng/Biobank/code/ldpred/max.R
bslmmc=/net/mulan/home/yasheng/Biobank/code/simulation_bslmm/bslmm.R

blockf=/net/mulan/home/yasheng/Biobank/data/LDblock/chr1.txt
vaIdx=/net/mulan/disk2/yasheng/plink_file/simulation/ref_bslmm_idx.txt
trIdx=/net/mulan/disk2/yasheng/plink_file/simulation/train_bslmm_idx.txt
teIdx=/net/mulan/disk2/yasheng/plink_file/simulation/test_bslmm_idx.txt

for dist in 1 ; do
for herit in 0.1 0.2 0.5; do
for ((cross=1;cross<11;cross++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
## simulate genotype data and phenotype data
path=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/data/
Rscript ${lmmPheno} --path ${path} --block ${block} --herit ${herit} --cross ${cross} --type 0 --dist ${dist}

## training and test set
bfile=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/data/geno_block${block}_her${herit}_cross${cross}_dist${dist}
bfiletr=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/data/train_block${block}_her${herit}_cross${cross}_dist${dist}
bfilete=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/data/test_block${block}_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfile} --keep ${trIdx} --make-bed --out ${bfiletr}
plink-1.9 --bfile ${bfile} --keep ${teIdx} --make-bed --out ${bfilete}
rm ${bfiletr}.log
rm ${bfilete}.log

## reference data
bfileva=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/data/val_block${block}_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfile} --keep ${vaIdx} --make-bed --out ${bfileva}
rm ${bfileva}.log
rm ${bfile}*

## summary statistics (GEMMA)
cd /net/mulan/disk2/yasheng/simulation_bslmm/LMM
summ=summary_block${block}_her${herit}_cross${cross}_dist${dist}
${gemma} -bfile ${bfiletr} -notsnp -lm 1 -n 1 -o ${summ}
sed -i '1d' /net/mulan/disk2/yasheng/simulation_bslmm/LMM/output/${summ}.assoc.txt
rm /net/mulan/disk2/yasheng/simulation_bslmm/LMM/output/${summ}.log.txt

## PT 
summf=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/output/summary_block${block}_her${herit}_cross${cross}_dist${dist}
Rscript ${toplinkf} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt
clump1=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/PT/summary_block${block}_her${herit}_cross${cross}_dist${dist}
pred1=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/PT/pheno_block${block}_her${herit}_cross${cross}_dist${dist}
for p in 5e-8 1e-6 1e-6 1e-4 1e-3 1e-2 5e-2 1e-1 2e-1 5e-1 1.0; do 
Rscript ${clumpf1} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt --ref ${bfileva} --pth ${p} \
--clump ${clump1}_p${p} --pred ${pred1}
rm ${clump1}_p${p}.clumped
rm ${clump1}_p${p}.log
rm ${pred1}.log
rm ${pred1}.nopred
r21=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/PT/r2_block${block}_her${herit}_cross${cross}_dist${dist}
Rscript ${simpred} --pheno ${pred1}.profile --r2 ${r21}_p${p}.txt
rm ${pred1}.profile
done
pbestf1=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/PT/pbest1_block${block}_her${herit}_cross${cross}_dist${dist}
pbestf2=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/PT/pbest2_block${block}_her${herit}_cross${cross}_dist${dist}
Rscript ${max1} --r2 ${r21}_p --pbest1 ${pbestf1}.txt --pbest2 ${pbestf2}.txt
pbest1=`cat ${pbestf1}.txt`
pbest2=`cat ${pbestf2}.txt`
plink-1.9 --bfile ${bfilete} --score ${clump1}_p${pbest1}.txt 1 2 3 sum --out ${pred1}
mv ${clump1}_p${pbest1}.txt ${clump1}_best.txt
rm ${r21}_p*
rm ${pred1}.log
rm ${pbestf1}.txt
rm ${pbestf2}.txt
rm ${clump1}*

m=`cat ${summf}.assoc.txt | wc -l`
n=`cat ${bfiletr}.fam | wc -l`
refld=/net/mulan/disk2/yasheng/sample500/frq/chr1

## SBSLMM
fix=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/slmm/summary_fix_block${block}_her${herit}_cross${cross}_dist${dist}
ran=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/slmm/summary_ran_block${block}_her${herit}_cross${cross}_dist${dist}
Rscript ${clumpf2} --summ ${summf} --ref ${refld} --p 0.001 --max 1000 --effl ${fix} --effs ${ran}
snp=/net/mulan/disk2/yasheng/simulation2/LMM/slmm/eff_block${block}_her${herit}_cross${cross}_dist${dist}
${sbslmm} -s ${ran}.txt -l ${fix}.txt -r ${refld} -nsnp ${m} -n ${n} -b ${blockf} -h ${herit} \
-maxMaf 0.2 -t ${thread} -eff ${snp}
pred2=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/slmm/pheno_block${block}_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfilete} --score ${snp}.txt 1 2 4 sum --out ${pred2}
rm ${pred2}.log
rm ${fix}.txt
rm ${ran}.txt

## BSLMM
cd /net/mulan/disk2/yasheng/simulation_bslmm/LMM/bslmm
snpeff1=snpeff_block${block}_her${herit}_cross${cross}_dist${dist}
${gemma} -bfile ${bfiletr} -bslmm 1 -o ${snpeff1} -w 6000 -s 2000 -rpace 1000
snpeff2=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/bslmm/output/snpeff_block${block}_her${herit}_cross${cross}_dist${dist}
snpeff3=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/bslmm/snpeff_block${block}_her${herit}_cross${cross}_dist${dist}
Rscript ${bslmmc} --bim ${bfiletr}.bim --eff ${snpeff2}.param.txt --effc ${snpeff3}.txt
predt=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/bslmm/pheno_block${block}_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfilete} --score ${snpeff3}.txt 1 2 3 sum --out ${predt}
rm ${snpeff2}*
rm ${snpeff3}*
rm ${predt}.log

## lasso
esttr3=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/lasso/train_block${block}_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfiletr} --lasso ${herit} --out ${esttr3}
predt=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/lasso/pheno_block${block}_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfilete} --score ${esttr3}.lasso 2 3 4 header sum --out ${predt}
rm ${predt}.log
rm ${esttr3}*

## SBLUP
cojo=$(echo "${m}*(1/${herit}-1)" | bc -l)
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summf}.assoc.txt > ${summf}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summf}.ma
est200=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/sblup/esteff_block${block}_her${herit}_cross${cross}_dist${dist}_200
est1000=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/sblup/esteff_block${block}_her${herit}_cross${cross}_dist${dist}_1000
${gcta} --bfile ${refld} --chr 1 --cojo-file ${summf}.ma --cojo-sblup ${cojo} --cojo-wind 200 --thread-num ${thread} --out ${est200}
${gcta} --bfile ${refld} --chr 1 --cojo-file ${summf}.ma --cojo-sblup ${cojo} --cojo-wind 1000 --thread-num ${thread} --out ${est1000}
pred3200=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/sblup/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_200
pred31000=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/sblup/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_1000
plink-1.9 --bfile ${bfilete} --score ${est200}.sblup.cojo 1 2 4 sum --out ${pred3200}
plink-1.9 --bfile ${bfilete} --score ${est1000}.sblup.cojo 1 2 4 sum --out ${pred31000}
rm ${est200}.log
rm ${est1000}.log
rm ${est200}*badsnps
rm ${est1000}*badsnps
rm ${est200}.sblup.cojo
rm ${est1000}.sblup.cojo
rm ${pred3200}.log
rm ${pred31000}.log

## ldpred
Rscript ${toldpred} --gemma ${summf}.assoc.txt  --ldpred ${summf}_LDpred.sumstat
coord1=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/summary_cv_block${block}_her${herit}_cross${cross}_dist${dist}.HDF5
${py} ${ldpred} coord --gf ${bfileva} --ssf ${summf}_LDpred.sumstat --out ${coord1} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
ldest=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/ld_block${block}_her${herit}_cross${cross}_dist${dist}
infest=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/esteff_block${block}_her${herit}_cross${cross}_dist${dist}
${py} ${ldpred} gibbs --cf ${coord1} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} --f 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 

pred41=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/pheno_block${block}_her${herit}_cross${cross}_dist${dist}
r24=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/r2_block${block}_her${herit}_cross${cross}_dist${dist}
for p in 1.0000e+00 3.0000e-01 1.0000e-01 3.0000e-02 1.0000e-02 3.0000e-03 1.0000e-03 3.0000e-04 1.0000e-04;do
plink-1.9 --bfile ${bfileva} --score ${infest}_LDpred_p${p}.txt 3 4 7 header sum --out ${pred41}_p${p}
Rscript ${simpred} --pheno ${pred41}_p${p}.profile --r2 ${r24}_p${p}.txt
rm ${pred41}_p${p}.log
rm ${pred41}_p${p}.nopred
rm ${pred41}_p${p}.profile
done
pbestf4=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/r2_block${block}_her${herit}_cross${cross}_dist${dist}_best.txt
Rscript ${max2} --r2 ${r24}_p --pbest ${pbestf4} 
pbest4=`cat ${pbestf4}`
rm ${r24}_p*

coord2=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/summary_ref_block${block}_her${herit}_cross${cross}_dist${dist}.HDF5
${py} ${ldpred} coord --gf ${refld} --ssf ${summf}_LDpred.sumstat --out ${coord2} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
let ldr1=ldr+1
${py} ${ldpred} gibbs --cf ${coord2} --ldr ${ldr1} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} --f ${pbest4}
predInf4=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_inf
predp4=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/ldpred/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_pbest
plink-1.9 --bfile ${bfilete} --score ${infest}_LDpred-inf.txt 3 4 7 sum --out ${predInf4}
plink-1.9 --bfile ${bfilete} --score ${infest}_LDpred_p${pbest4}.txt 3 4 7 sum --out ${predp4}
rm ${coord1}
rm ${coord2}
rm ${ldest}*.pkl.gz
rm ${infest}_LDpred*.txt
rm ${predInf4}.nopred
rm ${predInf4}.log
rm ${predp4}.nopred
rm ${predp4}.log


## lassosum
est5=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/lassosum/esteff_block${block}_her${herit}_cross${cross}_dist${dist}
Rscript ${lassosum} --summ ${summf}.assoc.txt --ref ${refld} --test ${bfileva} --beta ${est5}.txt
pred5=/net/mulan/disk2/yasheng/simulation_bslmm/LMM/lassosum/pheno_block${block}_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfilete} --score ${est5}.txt 1 2 3 sum --out ${pred5}

fi
done
done
done
