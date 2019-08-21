#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=2-00:00:00
#SBATCH --job-name=simSLM
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

#SBATCH --array=1-60%15
#SBATCH --output=/net/mulan/disk2/yasheng/out/simSLM10_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/simSLM10_%a.err

bash
let k=0

let block=10
let ldr=30
let thread=2

py=/net/mulan/home/yasheng/py3/bin/python
gemma=/net/mulan/home/yasheng/Biobank/program/GEMMA/gemma-0.98-linux-static
gcta=/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64
sbslmm=/net/mulan/home/yasheng/Biobank/code/slmm2.1/sbslmm
PLINK=/net/mulan/home/yasheng/plink2_linux_x86_64/plink2
ldsc=/net/mulan/home/yasheng/Biobank/program/ldsc/ldsc.py
ldpred=/net/mulan/home/yasheng/Biobank/program/ldpred/LDpred.py
lassosum=/net/mulan/home/yasheng/Biobank/code/lassosum/lassosum.R

slmPheno=/net/mulan/home/yasheng/Biobank/code/simulation/slmPheno.R
toplinkf=/net/mulan/home/yasheng/Biobank/code/PT/toplinkf.R
clumpf1=/net/mulan/home/yasheng/Biobank/code/PT/toclumpf.R
simpred=/net/mulan/home/yasheng/Biobank/code/PT/simPred.R
max1=/net/mulan/home/yasheng/Biobank/code/PT/max.R
mkldsc=/net/mulan/home/yasheng/Biobank/code/heritability/mkldsc.R
clumpf2=/net/mulan/home/yasheng/Biobank/code/slmm2/toclumpf.R
toldpred=/net/mulan/home/yasheng/Biobank/code/ldpred/toldpred.R
max2=/net/mulan/home/yasheng/Biobank/code/ldpred/max.R

blockf=/net/mulan/home/yasheng/Biobank/data/LDblock/chr1.txt
valIdx=/net/mulan/disk2/yasheng/plink_file/simulation/ref_idx2.txt
trIdx=/net/mulan/disk2/yasheng/plink_file/simulation/train_idx2.txt
teIdx=/net/mulan/disk2/yasheng/plink_file/simulation/test_idx2.txt
ps=0.001

for dist in 2 3 ; do
for herit in 0.1 0.2 0.5 ; do
for ((cross=1;cross<11;cross++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
## simulate genotype data and phenotype data
path=/net/mulan/disk2/yasheng/simulation2/SLM/data/
Rscript ${slmPheno} --path ${path} --block ${block} --herit ${herit} --cross ${cross} --type 0 --dist ${dist} --p ${ps}

## training and test set
bfile=/net/mulan/disk2/yasheng/simulation2/SLM/data/geno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
bfiletr=/net/mulan/disk2/yasheng/simulation2/SLM/data/train_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
bfilete=/net/mulan/disk2/yasheng/simulation2/SLM/data/test_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
plink-1.9 --bfile ${bfile} --keep ${trIdx} --make-bed --out ${bfiletr}
plink-1.9 --bfile ${bfile} --keep ${teIdx} --make-bed --out ${bfilete}
rm ${bfiletr}.log
rm ${bfilete}.log

## reference data
bfileval=/net/mulan/disk2/yasheng/simulation2/SLM/data/val_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
plink-1.9 --bfile ${bfile} --keep ${valIdx} --make-bed --out ${bfileval}
rm ${bfileval}.log
rm ${bfile}*

## summary statistics (GEMMA)
cd /net/mulan/disk2/yasheng/simulation2/SLM
summ=summary_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
${gemma} -bfile ${bfiletr} -notsnp -lm 1 -n 1 -o ${summ}
sed -i '1d' /net/mulan/disk2/yasheng/simulation2/SLM/output/${summ}.assoc.txt
rm /net/mulan/disk2/yasheng/simulation2/SLM/output/${summ}.log.txt

## PT 
summf=/net/mulan/disk2/yasheng/simulation2/SLM/output/summary_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
Rscript ${toplinkf} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt
clump1=/net/mulan/disk2/yasheng/simulation2/SLM/PT/summary_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
pred1=/net/mulan/disk2/yasheng/simulation2/SLM/PT/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
for p in 5e-8 1e-6 1e-4 1e-3 1e-2 5e-2 1e-1 2e-1 5e-1 1.0; do 
Rscript ${clumpf1} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt --ref ${bfileval} --pth ${p} \
--clump ${clump1}_p${p} --pred ${pred1}
rm ${clump1}_p${p}.clumped
rm ${clump1}_p${p}.log
rm ${pred1}.log
rm ${pred1}.nopred
r21=/net/mulan/disk2/yasheng/simulation2/SLM/PT/r2_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
Rscript ${simpred} --pheno ${pred1}.profile --r2 ${r21}_p${p}.txt
rm ${pred1}.profile
done
pbestf1=/net/mulan/disk2/yasheng/simulation2/SLM/PT/pbest1_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
pbestf2=/net/mulan/disk2/yasheng/simulation2/SLM/PT/pbest2_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
Rscript ${max1} --r2 ${r21}_p --pbest1 ${pbestf1}.txt --pbest2 ${pbestf2}.txt
pbest1=`cat ${pbestf1}.txt`
pbest2=`cat ${pbestf2}.txt`
plink-1.9 --bfile ${bfilete} --score ${clump1}_p${pbest1}.txt 1 2 3 sum --out ${pred1}
mv ${clump1}_p${pbest1}.txt ${clump1}_best.txt
rm ${r21}_p*
rm ${pred1}.log
rm ${pred1}.txt
rm ${pred2}.txt
rm ${clump1}*

m=`cat ${summf}.assoc.txt | wc -l`
n=`cat ${bfiletr}.fam | wc -l`
refld=/net/mulan/disk2/yasheng/sample500/frq/chr1

## SLMM
fix=/net/mulan/disk2/yasheng/simulation2/SLM/slmm/summary_fix_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
ran=/net/mulan/disk2/yasheng/simulation2/SLM/slmm/summary_ran_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
Rscript ${clumpf2} --summ ${summf} --ref ${refld} --p 0.001 --max 1000 --effl ${fix} --effs ${ran}
snp=/net/mulan/disk2/yasheng/simulation2/SLM/slmm/eff_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
${sbslmm} -s ${ran}.txt -l ${fix}.txt -r ${refld} -nsnp ${m} -n ${n} -b ${blockf} -h ${herit} -maxMaf 0.2 -t ${thread} -eff ${snp}
pred2=/net/mulan/disk2/yasheng/simulation2/SLM/slmm/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
plink-1.9 --bfile ${bfilete} --score ${snp}.txt 1 2 4 sum --out ${pred2}
rm ${pred2}.log
rm ${fix}.txt
rm ${ran}.txt

## SBLUP
cojo=$(echo "${m}*(1/${herit}-1)" | bc -l)
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summf}.assoc.txt > ${summf}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summf}.ma
est200=/net/mulan/disk2/yasheng/simulation2/SLM/sblup/esteff_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_200
est1000=/net/mulan/disk2/yasheng/simulation2/SLM/sblup/esteff_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_1000
${gcta} --bfile ${refld} --chr 1 --cojo-file ${summf}.ma --cojo-sblup ${cojo} --cojo-wind 200 --thread-num ${thread} --out ${est200}
${gcta} --bfile ${refld} --chr 1 --cojo-file ${summf}.ma --cojo-sblup ${cojo} --cojo-wind 1000 --thread-num ${thread} --out ${est1000}
pred3200=/net/mulan/disk2/yasheng/simulation2/SLM/sblup/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_200
pred31000=/net/mulan/disk2/yasheng/simulation2/SLM/sblup/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_1000
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
coord1=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/summary_cv_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}.HDF5
${py} ${ldpred} coord --gf ${bfileval} --ssf ${summf}_LDpred.sumstat --out ${coord1} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
ldest=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/ld_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
infest=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/esteff_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
${py} ${ldpred} gibbs --cf ${coord1} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} \
--f 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 

pred41=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
r24=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/r2_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
for p in 1.0000e+00 3.0000e-01 1.0000e-01 3.0000e-02 1.0000e-02 3.0000e-03 1.0000e-03 3.0000e-04 1.0000e-04;do
plink-1.9 --bfile ${bfileval} --score ${infest}_LDpred_p${p}.txt 3 4 7 header sum --out ${pred41}_p${p}
Rscript ${simpred} --pheno ${pred41}_p${p}.profile --r2 ${r24}_p${p}.txt
rm ${pred41}_p${p}.log
rm ${pred41}_p${p}.nopred
rm ${pred41}_p${p}.profile
done
pbestf4=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/r2_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_best.txt
Rscript ${max2} --r2 ${r24}_p --pbest ${pbestf4} 
pbest4=`cat ${pbestf4}`

coord2=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/summary_ref_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}.HDF5
${py} ${ldpred} coord --gf ${refld} --ssf ${summf}_LDpred.sumstat --out ${coord2} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
let ldr1=ldr+1
${py} ${ldpred} gibbs --cf ${coord2} --ldr ${ldr1} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} --f ${pbest4}
predInf4=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_inf
predp4=/net/mulan/disk2/yasheng/simulation2/SLM/ldpred/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_pbest
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
rm ${r24}*

## lassosum
est5=/net/mulan/disk2/yasheng/simulation2/SLM/lassosum/esteff_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
Rscript ${lassosum} --summ ${summf}.assoc.txt --ref ${refld} --test ${bfileval} --beta ${est5}.txt
pred5=/net/mulan/disk2/yasheng/simulation2/SLM/lassosum/pheno_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}
plink-1.9 --bfile ${bfilete} --score ${est5}.txt 1 2 3 sum --out ${pred5}

fi
done
done
done
