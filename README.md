## SBSLMM
"SBSLMM", refering to Summary statistics for Bayesian Sparse Linear Mixed Model (SBSLMM), which wrapped the Plink and Sbslmm code. SBSLMM both defines the indepnedent large effect SNPs and re-estimate the merginal effect with the reference LD panel and block information. The aim of SBSLMM can handle large-scale of SNPs with low memory and high running speed, expecailly UK Biobank scale (~9 million SNPs). 

### Installation
- PLINK <br>
SBSLMM uses PLINK to select independent large effect SNPs. The link of Plink-1.9: https://www.cog-genomics.org/plink/1.9/.
- sbslmm <br>
SBSLMM uses sbslmm to re-restiamte the marginal effect (z-score). The executable file of sbslmm is included in the file folder SBSLMM.
- block file <br>
SBSLMM uses the block information. Reference: [Berisa and Pickrell (2015)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium)
- Linkage Disequilibrium SCore regression (LDSC) <br>
SBSLMM uses the estimated heritability. You can use any sofeware to estimate the heritability, such as [LDSC](https://github.com/bulik/ldsc), [GEMMA](https://github.com/genetics-statistics/GEMMA) and [GCTB](http://cnsgenomics.com/software/gctb/#SummaryBayesianAlphabet). <br>
In the paper, we use LDSC to estimate heritability. Here is the code: 
````shell
py=/net/mulan/home/yasheng/py3/bin/python
ldsc=/net/mulan/home/yasheng/Biobank/program/ldsc/ldsc.py
summ=/net/mulan/disk2/yasheng/GWAS/summary_pheno${PHENO}
ref=/net/mulan/disk2/yasheng/sample500/ldsc/
herit=/net/mulan/disk2/yasheng/GWAS/heritability/h2_pheno${PHENO}
## LDSC
${py} ${ldsc} --h2 ${summ}.ldsc.gz --ref-ld-chr ${ref} --w-ld-chr ${ref} --out ${herit}
rm ${summ}.ldsc.gz
## get herit
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`
````
- SBSLMM requires the following R packages: data.table, optparse. Install them by: <br>
````R
install.packages(c("data.table", "optparse"), dependencies=TRUE)
````
### Input file format
- summary statistics (GEMMA format) <br>
The separate is tab (\t).
````
chr rs ps n_mis n_obs allele1 allele0 af beta se p_wald
1	rs554763599	41975332	0	10000	A	G	0.476	1.768782e-03	1.731524e-02	9.186383e-01
1	rs678110	41976017	0	10000	C	T	0.409	2.573095e-02	1.756221e-02	1.429165e-01
1	rs72669005	41976217	110	9890	A	G	0.040	-7.861154e-02	4.424035e-02	7.561186e-02
1	rs140839222	41976345	88	9912	C	T	0.027	-1.061644e-01	5.334030e-02	4.658253e-02
1	rs377343544	41976500	0	10000	A	G	0.058	2.293323e-02	3.689809e-02	5.342658e-01
1	rs11809423	41976529	0	10000	T	C	0.037	-6.884283e-02	4.559762e-02	1.311286e-01
````
- block information (.bed)
The block information is download from the website: https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/

### Parameter Setting and example code
You can directly use the Rscript `/SBSLMM/SBSLMM.R`. You can get the explaination of each parameter: 
````R
Rscript SBSLMM.R -h
Rscript SBSLMM.R --help
````
The details is: 
````
--summary=CHARACTER
		INPUT: the summary statistics (gemma output format)
--tmpPath=CHARACTER
		INPUT: the temptorary path
--outPath=CHARACTER
		INPUT: the output path
--plink=CHARACTER
		INPUT: the perfix of Plink software
--ref=CHARACTER
		INPUT: the perfix of reference panel
--r2=CHARACTER
		INPUT: the cutoff of SNPs clumping (default:0.1)
--pv=CHARACTER
		INPUT: the cutoff of SNPs pruning (default:1e-6)
--sbslmm=CHARACTER
		INPUT: the perfix of sbslmm software
--mafMax=CHARACTER
		INPUT: the maximium of the difference between reference panel and summary data
--nsnp=CHARACTER
		INPUT: the number of SNPs in whole genome
--block=CHARACTER
		INPUT: the block information (Berisa and Pickrell 2015)
--h2=CHARACTER
		INPUT: the heritability of trait
--thread=CHARACTER
		INPUT: the number of threads
````
### Output file format
The example of output file is: 
````
rs559141705 G -0.0686326 -0.220596 1
rs116836856 T 0.0825543 0.160458 1
rs139791561 T -0.0339878 -0.0541334 1
rs80205519 A 0.0170505 0.057553 1
rs114747222 A -0.106633 -0.399451 1
rs1357783 A 0.114697 0.163234 1
rs7543904 C 0.0425442 0.0673737 1
rs1202804 A 0.0691282 0.111027 1
rs1749778 C 0.10026 0.144773 1
rs12024338 C 0.105702 0.150094 1
rs991302 C 0.0612944 0.132898 1
````
The first column is SNP ID. The second column is allele code. The third code is scaled effect size. The forth is non-scaled effect size. Here, we use the MAF from the summary data. You can also use the testing data to transfer it. The fifth column is the index of whether it is large effect or not (1: large effect, 0: small effect). You can directly use the output file to plink-1.9. The code is: 
````shell
bfilete=/path/test
est=/path/est_SBSLMM.txt
pred=/path/pheno_pred.txt
## plink 1.9
plink-1.9 --bfile ${bfilete} --score ${est}.txt 1 2 4 sum --out ${pred}
## plink 2
plink-1.9 --bfile ${bfilete} --score ${est}.txt 1 2 4 cols=+scoresums --out ${pred}
````
### Simulation
- P+T <br>
To select the best cutoff, we ues a validate data to do P+T. We use the plink-1.9. The `toplinkf`, `clumpf`, `simPred` and `max` functions are in `PT` folder. Here, the code is as following: 
````shell
toplinkf=~/PT/toplinkf.R
clumpf=~/PT/toclumpf.R
simpred=/PT/simPred.R
max=~/PT/max.R

bfileval=/path/validate
bfilete=/path/test
## transfer gemma format to plink format
summf=/path/summary_statistics_from_gemma
Rscript ${toplinkf} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt
## 
clump1=/path/summary
pred1=/path/pheno
for p in 5e-8 1e-6 1e-4 1e-3 1e-2 5e-2 1e-1 2e-1 5e-1 1.0; do 
Rscript ${clumpf1} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt --ref ${bfileval} --pth ${p} \
--clump ${clump1}_p${p} --pred ${pred1}
rm ${clump1}_p${p}.clumped
rm ${clump1}_p${p}.log
rm ${pred1}.log
rm ${pred1}.nopred
r21=/net/mulan/disk2/yasheng/simulation2/SLMM/PT/r2_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_prop${prop}
Rscript ${simpred} --pheno ${pred1}.profile --r2 ${r21}_p${p}.txt
rm ${pred1}.profile
done
pbestf1=/path/pbest1
pbestf2=/path/pbest2
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
````
- SBLUP <br>
SBLUP is performed by GCTA. The code is similar to the example of GCTA. 
````shell
bfiletr=/path/train
summf=/path/summary_statistics_from_gemma
refld=/path/ref
herit=0.1
m=`cat ${summf}.assoc.txt | wc -l`
n=`cat ${bfiletr}.fam | wc -l`
cojo=$(echo "${m}*(1/${herit}-1)" | bc -l)
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summf}.assoc.txt > ${summf}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summf}.ma
est3=/path/esteff
${gcta} --bfile ${refld} --chr 1 --cojo-file ${summf}.ma --cojo-sblup ${cojo} --cojo-wind 200 --thread-num ${thread} --out ${est3}
pred3=/path/pheno
plink-1.9 --bfile ${bfilete} --score ${est}.sblup.cojo 1 2 4 sum --out ${pred3}
rm ${est3}.log
rm ${est3}*badsnps
rm ${est3}.sblup.cojo
rm ${pred3}.log
````
- LDpred
LDpred is performed by manual of LDpred, including the cutoff and radius. 
````shell
bfiletr=/path/train
bfileval=/path/validate
bfilete=/path/test
n=`cat ${bfiletr}.fam | wc -l`
herit=0.1
## transfer gemma format to LDpred format
Rscript ${toldpred} --gemma ${summf}.assoc.txt  --ldpred ${summf}_LDpred.sumstat
## validate
coord1=/path/summary_cv.HDF5
${py} ${ldpred} coord --gf ${bfileval} --ssf ${summf}_LDpred.sumstat --out ${coord1} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
ldest=/path/ld
infest=/path/esteff
${py} ${ldpred} gibbs --cf ${coord1} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} \
--f 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 
pred41=/path/pheno
r24=/path/r2
for p in 1.0000e+00 3.0000e-01 1.0000e-01 3.0000e-02 1.0000e-02 3.0000e-03 1.0000e-03 3.0000e-04 1.0000e-04;do
plink-1.9 --bfile ${bfileval} --score ${infest}_LDpred_p${p}.txt 3 4 7 header sum --out ${pred41}_p${p}
Rscript ${simpred} --pheno ${pred41}_p${p}.profile --r2 ${r24}_p${p}.txt
rm ${pred41}_p${p}.log
rm ${pred41}_p${p}.nopred
rm ${pred41}_p${p}.profile
done
pbestf4=/path/r2_best.txt
Rscript ${max2} --r2 ${r24}_p --pbest ${pbestf4} 
pbest4=`cat ${pbestf4}`

## test
coord2=/path/summary_ref.HDF5
${py} ${ldpred} coord --gf ${refld} --ssf ${summf}_LDpred.sumstat --out ${coord2} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
${py} ${ldpred} gibbs --cf ${coord2} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} --f ${pbest4}
predInf4=/path/pheno_inf
predp4=/path/pheno_pbest
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
````
- lassosum <br>
We use the lassosum by the R package `lassosum`. The 'lassosum' function is in the `lassosum` file folder. 
````shell
res5=/net/mulan/disk2/yasheng/simulation_bslmm/SLMM/lassosum/res_block${block}_her${herit}_cross${cross}_dist${dist}_ps${ps}_prop${prop}
Rscript ${lassosum} --summ ${summf}.assoc.txt --ref ${refld} --valid ${bfileva} --test ${bfilete} --res ${res5}.txt
````
- lasso (sample size = 2,000)<br>
We use the lasso algorithm in Plink-1.9. The code is as following: 
````shell
## lasso
esttr3=/path/train
plink-1.9 --bfile ${bfiletr} --lasso ${herit} --out ${esttr3}
predt=/path/pheno
plink-1.9 --bfile ${bfilete} --score ${esttr3}.lasso 2 3 4 header sum --out ${predt}
rm ${predt}.log
rm ${esttr3}*
````
- BSLMM (sample size = 2,000)<br>
We use the BSLMM in the `GEMMA` software. 'bslmm.R' is in the 'simulation_small_size' file folder. 
````shell
## BSLMM
gemma=/path/gemma-0.98-linux-static
bslmmc=/bslmm/bslmm.R
cd /path
snpeff1=snpeff
${gemma} -bfile ${bfiletr} -bslmm 1 -o ${snpeff1} -w 6000 -s 2000 -rpace 1000
snpeff2=/path/output/snpeff
snpeff3=/path/snpeff
Rscript ${bslmmc} --bim ${bfiletr}.bim --eff ${snpeff2}.param.txt --effc ${snpeff3}.txt
predt=/path/pheno
plink-1.9 --bfile ${bfilete} --score ${snpeff3}.txt 1 2 3 sum --out ${predt}
rm ${snpeff2}*
rm ${snpeff3}*
rm ${predt}.log
````
All the simulation result are included in the manuscript.

### UK Biobank data process
For SBSLMM, we also process the UKB data. Following the Neale lab and official file of UKB, we perform qualtiy control of samples and SNPs. 
#### SNP QC
We use `snp_qc.R` and `toplink.sh` to perform the SNP QC. The details are shown in the manuscript.
- MAF <br>
We select SNPs with MAF>0.2. 
- INFO score<br>
We select the SNPs with INFO score>0.8. 
- Duplicate<br>
We delete the duplicated SNPs. 
- HWE<br>
We delete the SNPs with HWE<1e-7.
- Missing rate<br>
We delete the SNPS with Pm>0.05.
#### Sample QC
We use `mkpheno.R` to select the samples. The details are shown in the manuscript.<br>
