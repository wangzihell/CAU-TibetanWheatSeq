#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19
#
source /WORK/app/osenv/ln1/set2.sh
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start filter"
#
WP="/WORK/pp192/"
#
for BAM in ../03.*/??.sort.bam; do
  echo $BAM
  ID=`basename ${BAM} | sed s/.sort.bam//g`
  (bamtools filter \
    -in ${BAM}  \
    -out ${ID}.filter.bam \
    -forceCompression -script ../../filter.list;
   samtools index ${ID}.filter.bam;
   samtools flagstat ${ID}.filter.bam > ${ID}.flagstat
   rm ${BAM})&
done
#
wait
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Filter"
#
cd ../05*/
sh ./05.mergeAsplit.sh
#
