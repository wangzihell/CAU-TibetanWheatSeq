#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start BWA"

WP="/WORK/pp192/"

DIR=`dirname $PWD`
SM=`basename $DIR`
#
ID=$1
R1="../02.*/${ID}_R1.clean.fq.gz"
R2="../02.*/${ID}_R2.clean.fq.gz"
#
bwa mem -t 20 \
  -R '@RG\tID:'${SM}'\tLB:'${SM}'\tPL:ILLUMINA\tSM:'${SM} \
  ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
  ${R1} ${R2} | \
  samtools view -Sbh - > ${ID}.bam
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] To sort bam file"
#
samtools sort -@ 20 -m 48G -o ${ID}.sort.bam ${ID}.bam
#
samtools index ${ID}.sort.bam
#
samtools flagstat ${ID}.sort.bam > ${ID}.flagstat
#
rm ${ID}.bam*
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish BWA"
#
