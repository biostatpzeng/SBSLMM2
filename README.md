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

### Simulation
- P+T
- SBLUP
- LDpred
- lasso (sample size = 2,000)
- BSLMM (sample size = 2,000)
figure1 and figure 2

### UK Biobank data process
#### SNP QC
#### Sample QC
figure of process

