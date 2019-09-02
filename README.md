## SBSLMM
"SBSLMM", refering to Summary statistics for Bayesian Sparse Linear Mixed Model (SBSLMM), which wrapped the Plink and Sbslmm code. SBSLMM both defines the indepnedent large effect SNPs and re-estimate the merginal effect with the reference LD panel and block information. The aim of SBSLMM can handle large-scale of SNPs with low memory and high running speed, expecailly UK Biobank scale (~9 million SNPs). 

### Installaiton
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

### Run

