#!/usr/bin/env bash

WP="/WORK/pp192/"

#for ID in `for F in ../01*/R1_??.gz; do basename $F; done | cut -c4-5 | tr '\n' ' '`; do
for ID in aa; do 
  echo $ID
  R1="R1_"${ID}.gz
  R2="R2_"${ID}.gz
  (java \
     -jar /home/wangzh/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
     PE -phred33 -threads 1 \
     ../01.split/${R1} ../01.split/${R2} \
     ${ID}_R1.clean.fq.gz ${ID}_R1.unpaired.fq.gz ${ID}_R2.clean.fq.gz ${ID}_R2.unpaired.fq.gz \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;
   zcat ../01.split/${R1} | gawk '{T++;} END{print T/4;}' > ../01.split/${ID}.count ;
   zcat ${ID}_R1.clean.fq.gz | gawk '{T++;} END{print T/4;}' > ${ID}.count )&
done