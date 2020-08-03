#!/usr/bin/env bash
# Guo, Weilong; guoweilong@126.com; 2017-10-24

for BAM in *.uniq.bam; do
  CHR=`basename $BAM | sed s/.uniq.bam//g`
  echo $CHR
  yhbatch -n 1 -e ${CHR}.e -o ${CHR}.o ./justcall.sh $CHR
done