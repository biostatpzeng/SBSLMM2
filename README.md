## SBSLMM
"SBSLMM", refering to Summary statistics for Bayesian Sparse Linear Mixed Model (SBSLMM), which wrapped the Plink and Sbslmm code. SBSLMM both defines the indepnedent large effect SNPs and re-estimate the merginal effect with the reference LD panel and block information. The aim of SBSLMM can handle large-scale of SNPs with low memory and high running speed, expecailly UK Biobank scale (~9 million SNPs). 

### Installaiton
- PLINK <br>
SBSLMM uses PLINK to select independent large effect SNPs. The link of Plink-1.9: https://www.cog-genomics.org/plink/1.9/.
- sbslmm <br>
SBSLMM uses sbslmm to re-restiamte the marginal effect (z-score). The executable file of sbslmm is included in the file folder SBSLMM.
- block file <br>

- LDSC <br>

