#!/bin/bash
#SBATCH --partition=nomosix,mulan
#SBATCH --time=1-00:00:00
#SBATCH --job-name=PT
#SBATCH --mem-per-cpu=1G

#SBATCH --array=1
#SBATCH --output=/net/mulan/disk2/yasheng/out/PT_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/PT_%a.err


mkPT=/net/mulan/home/yasheng/Biobank/code/PT/PT.R
for PHENO in 1 ;do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	for cross in 1 2 3 4 5;do
		for ((chr=1;chr<23;chr++)); do
			summ_gemma=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}_chr${chr}.assoc.txt
			plinkf=/net/mulan/disk2/yasheng/pheno${PHENO}/summ/summary_cross${cross}_chr${chr}.assoc.plink
			# refeur=/net/mulan/disk2/yasheng/EUR/genotype/chr${chr}
			clump=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/summary_cross${cross}_chr${chr}
			# effeur=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_eur/esteff_cross${cross}_chr${chr}.txt
			# fixeur=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_eur/fix_cross${cross}_chr${chr}.assoc.txt
			# timeeur=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_eur/time_cross${cross}_chr${chr}.assoc.txt
			ref2000=/net/mulan/disk2/yasheng/sample2000/frq/chr${chr}
			eff2000=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/esteff_cross${cross}_chr${chr}.txt
			fix2000=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/fix_cross${cross}_chr${chr}.assoc.txt
			time2000=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/time_cross${cross}_chr${chr}.assoc.txt
			
			#### PT #####
			# Rscript ${mkPT} --summ ${summ_gemma} --assoc ${plinkf} --ref ${refeur} --clump ${clump} --eff ${effeur} --fix ${fixeur} --time ${timeeur}
			Rscript ${mkPT} --summ ${summ_gemma} --assoc ${plinkf} --ref ${ref2000} --clump ${clump} --eff ${eff2000} --fix ${fix2000} --time ${time2000}
			rm ${plinkf} 
			rm ${clump}*
			rm ${eff2000}

			#### PT validation #####
			bfile=/net/mulan/disk2/yasheng/tot_plink_file/chr${chr}
			# phenoeur=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_eur/pheno_cross${cross}_chr${chr}.txt
			pheno2000=/net/mulan/disk2/yasheng/pheno${PHENO}/PT_2000/pheno_cross${cross}_chr${chr}
			# ${valid} -bfile ${bfile} -pfile ${pfile} -ncol ${cross} -thread 5 -scaled T -eff ${effeur} -predPheno ${phenoeur} 
			plink-1.9 --bfile ${bfile} --score ${fix2000} 2 6 9 sum --out ${pheno2000}
			# rm ${effeur} 
			rm ${pheno2000}.log
		done
	done
	fi
done

