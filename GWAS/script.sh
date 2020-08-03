#!/bin/bash

#filter
plink --vcf wheat_genome.vcf.gz --maf 0.05 --geno 0.2 --recode vcf-iid --out filter_wheat_geno --allow-extra-chr
#pca=5
plink --vcf filter_wheat_geno.vcf --pca 5 --out pca_wheat_geno --allow-extra-chr
#kinship
/data/user/wangxb/software/tassel/run_pipeline.pl -Xmx200g -importGuess filter_wheat_geno.vcf -KinshipPlugin -method Centered_IBS -en
dPlugin -export kin_wheat_genome.txt -exportType SqrMatrix
#cmlm
/data/user/wangxb/software/tassel/run_pipeline.pl -Xmx300G -fork1 -vcf filter_wheat_genome.vcf.gz -fork2 -r pheno.txt -fork3 -r pca_wheat_geno.txt -excludeLastTrait -fork4 -k kin_wheat_genome.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export mlm_geno -runfork1 -runfork2 -runfork3 -runfork4
