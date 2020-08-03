#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Split"

WP="/WORK/pp192/"

#for ID in `for F in ../01*/R1_??.gz; do basename $F; done | cut -c4-5 | tr '\n' ' '`; do
for ID in `while read line; do echo $line|tr '\n' " "; done < $1`; do 
  echo $ID
  R1="R1_"${ID}.gz
  R2="R2_"${ID}.gz
  (java \
     -jar ${WP}/Install/Trimmomatic-0.36/trimmomatic-0.36.jar \
     PE -phred33 -threads 1 \
     ../01.split/${R1} ../01.split/${R2} \
     ${ID}_R1.clean.fq.gz ${ID}_R1.unpaired.fq.gz ${ID}_R2.clean.fq.gz ${ID}_R2.unpaired.fq.gz \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;
   zcat ../01.split/${R1} | gawk '{T++;} END{print T/4;}' > ../01.split/${ID}.count ;
   rm ../01.split/${R1} ../01.split/${R2} ;
   zcat ${ID}_R1.clean.fq.gz | gawk '{T++;} END{print T/4;}' > ${ID}.count )&
done

wait

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Split"

cd ../03.*/
sh 03.*.sh

