#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19
#
source /WORK/app/osenv/ln1/set2.sh
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start FlagStat"
#
WP="/WORK/pp192/"
#
#
for BAM in *.bam; do
   ID=`echo $BAM | cut -d"." -f1` 
   samtools flagstat $BAM > ${ID}.flagstat &
done
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish FlagStat"
#

