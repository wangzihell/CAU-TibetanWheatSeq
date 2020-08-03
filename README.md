# The Tibetan wheat sequencing project 
## Performed by Wheat Genetics and Genomics Center (WGGC) at China Argiculture University (CAU)

This repo contains customized codes used in Tibetan wheat sequencing project. Codes are grouped into multiple folders according to their usage. Appropriate comments and user manuals have been added within these scripts. A README file is provided including the descriptions of each folder.

## System requirements
Samtools (v1.3.1)
Bamtools (v2.4.1)
Bedtools (v2.27.1)
GTAK (v3.8)
BWA (v0.7.15)
Trimmomatic (v0.36)
Bcftools (v1.8-41)
VCFtools (v0.1.13)
R (v3.6)
Admixture (v1.3.0)
Tassel (5.2.50)
Plink (v1.90)
smartpca (v6.1.4)

## Installation guide
All softwares listed in system requirements could be downloaded and installed following instructions on their offical websites.

## Demos and corresponding instructions
Demos and corresponding instructions are in each foler.

## Brief description of code in each folder
Zang1817_assembly:  
Construction of pseudomolecules of Zang1817 genome.  
variant_calling:  
SNPs/INDELs calling, filtering and annotating of samples used this study from raw fastq files to VCF files.  
PAV:  
PAV calling between Chinese Spring IWGSC Refseq assembly and Zang1817 genome assembled in this study.  
NJ_tree:  
Neighbor joining tree construction and ploting from variants stored in VCF format.  
PCA:  
PCA computation and ploting from VCF files.  
Admixtue  
doing Admixtue analyze from VCF files.  
GWAS:  
Genome-wide scanning of candidate regions for rachis brittleness during wheat de-domestication process using GWAS analysis.  
Fst:  
Genome-wide Fst scaning between HA(high altitude), LA(low altitude) and between DE(de-domesticated), DO(domesticated) wheat accessions.  
pi-ratio:  
Genome-wide Pi-ration(PiDO/PiDE) scaning between DE(de-domesticated) and DO(domesticated) wheat accessions.  
CNVindex:  
Genome-wide CNVindex analyze between DE(de-domesticated) and DO(domesticated) wheat accessions.  
