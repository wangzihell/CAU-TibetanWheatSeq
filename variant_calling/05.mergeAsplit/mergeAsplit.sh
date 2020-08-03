#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start mergeAsplit 1"

WP="/WORK/pp192/"

for CHR in `samtools view -H ../04.*/aa.filter.bam | grep '^@SQ' | gawk '{print substr($2, 4)}' | grep '\.1'` Mt; do
  (echo $CHR ;
  #
   ( samtools view -H ../04.*/aa.filter.bam;
  for BAM in ../04.*/??.filter.bam ; do samtools view ${BAM} $CHR; done ) | samtools view -Sb - > ${CHR}.bam
  #
  samtools sort  -o ${CHR}.sort.bam ${CHR}.bam; rm ${CHR}.bam
  samtools rmdup -S ${CHR}.sort.bam ${CHR}.dedup.bam
  samtools index ${CHR}.dedup.bam
  samtools flagstat ${CHR}.dedup.bam > ${CHR}.flagstat
  rm ${CHR}.sort.bam) &
done
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish mergeAsplit 1"
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start mergeAsplit 2"

for CHR in `samtools view -H ../04.*/aa.filter.bam | grep '^@SQ' | gawk '{print substr($2, 4)}' | grep '\.2'` Pt; do
  (echo $CHR ;
  #
   ( samtools view -H ../04.*/aa.filter.bam;
  for BAM in ../04.*/??.filter.bam ; do samtools view ${BAM} $CHR; done ) | samtools view -Sb - > ${CHR}.bam
  #
  samtools sort  -o ${CHR}.sort.bam ${CHR}.bam; rm ${CHR}.bam
  samtools rmdup -S ${CHR}.sort.bam ${CHR}.dedup.bam
  samtools index ${CHR}.dedup.bam
  samtools flagstat ${CHR}.dedup.bam > ${CHR}.flagstat
  rm ${CHR}.sort.bam) &
done
#
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish mergeAsplit 2"
#
rm ../04.*/*.filter.bam*
#
cd ../06.*/
sh 06.*.sh


