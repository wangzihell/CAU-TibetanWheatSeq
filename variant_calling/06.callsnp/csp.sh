#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Callsnp"

WP="/WORK/pp192/"
#
CHR=$1
#
samtools view -h ../05.*/${CHR}.dedup.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${CHR}.uniq.bam
samtools index ${CHR}.uniq.bam
samtools flagstat ${CHR}.uniq.bam > ${CHR}.flagstat
#
echo $CHR
#  -jdk_deflater -jdk_inflater are new features in 3.8 but not support in GZ
java -Xmx60g \
  -Djava.io.tmpdir=${WP}/tmp \
  -jar ${WP}/Install/GenomeAnalysisTK.jar \
  -jdk_deflater \
  -jdk_inflater \
  -nct 24 \
  -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
  -T HaplotypeCaller -ERC GVCF -L ${CHR} \
  -I ../05.*/${CHR}.bam \
  -o ${CHR}.g.vcf.gz
#
wait
#
rm ${CHR}.uniq.bam*
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Callsnp"

