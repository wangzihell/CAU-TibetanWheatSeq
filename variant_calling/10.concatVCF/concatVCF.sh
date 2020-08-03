#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

source /WORK/app/osenv/ln1/set2.sh

WP="/WORK/pp192/"
CHR=$1

tmp1=`mktemp -d /WORK/pp192/tmp/XXXXXXX`
tmp2=`mktemp -d /WORK/pp192/tmp/XXXXXXX`

bcftools concat `ls ../09.filterVCF/${CHR}.1.*.filter.final.vcf.gz` -O b --threads 4 |\
    bcftools sort -T $tmp1 -O b -o ${CHR}.1.sort.bcf.gz &

bcftools concat `ls ../09.filterVCF/${CHR}.2.*.filter.final.vcf.gz` -O b --threads 4 |\
    bcftools sort -T $tmp2 -O b -o ${CHR}.2.sort.bcf.gz &
#
wait
#
(bcftools view -h ${CHR}.1.sort.bcf.gz | grep -v "##contig" | sed "/^##reference/ r header_length.txt";
bcftools view --threads 4 -H ${CHR}.1.sort.bcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);print}'; 
bcftools view --threads 4 -H ${CHR}.2.sort.bcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);$2=$2+400000000;print}' )|\
bcftools view --threads 4 -Ob -o ${CHR}.bcf.gz

bcftools index ${CHR}.bcf.gz