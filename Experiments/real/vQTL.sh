#!/bin/bash

cd ./QUAIL; #https://github.com/qlu-lab/QUAIL

/s/bin/Rscript  Step1_QUAIL_rank_score.R \
--pheno ./pheno/Total.lab.y.txt \
--covar ./covar/Total.lab.covar.txt \
--output ./results/rs/Total.lab.y_rankscore.txt \
--num_levels 500 \
--num_cores 10

/s/bin/Rscript  Step1_QUAIL_rank_score.R \
--pheno ./pheno/Total.unlab.y_hat.txt \
--covar ./covar/Total.unlab.covar.txt \
--output ./results/rs/Total.unlab.y_hat_rankscore.txt \
--num_levels 500 \
--num_cores 10

/s/bin/Rscript  Step1_QUAIL_rank_score.R \
--pheno ./pheno/Total.lab.y_hat.txt \
--covar ./covar/Total.lab.covar.txt \
--output ./results/rs/Total.lab.y_hat_rankscore.txt \
--num_levels 500 \
--num_cores 10

for chr in {1..22}; do
pheno_path="./results/rs/Total.lab.y_rankscore.txt" 
covar_path="./data/covar/Total.lab.covar.txt" 
dir_out="./results/gwas/Total.lab.y_chr.${chr}"  

geno_path="./genotype/chr${chr}.EUR.idupdated"
./plink2 \ # https://www.cog-genomics.org/plink/2.0/
--bfile $geno_path \
--pheno $pheno_path \
--covar $covar_path \
--glm 'qt-residualize' hide-covar \
--out $dir_out

pheno_path="./results/rs/Total.lab.y_hat_rankscore.txt" 
covar_path="./data/covar/Total.lab.covar.txt" 
dir_out="./results/gwas/Total.lab.y_hat_chr.${chr}"   

geno_path="./genotype/chr${chr}.EUR.idupdated"
./plink2 \
--bfile $geno_path \
--pheno $pheno_path \
--covar $covar_path \
--glm 'qt-residualize' hide-covar \
--out $dir_out

pheno_path="./results/rs/Total.unlab.y_hat_rankscore.txt" 
covar_path="./data/covar/Total.unlab.covar.txt" 
dir_out="./results/gwas/Total.unlab.y_hat_chr.${chr}"   

geno_path="./genotype/chr${chr}.EUR.idupdated"

./plink2 \
--bfile $geno_path \
--pheno $pheno_path \
--covar $covar_path \
--glm 'qt-residualize' hide-covar \
--out $dir_out
done
