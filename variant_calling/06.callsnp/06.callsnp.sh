#!/usr/bin/env bash
# Guo, Weilong; guoweilong@126.com; 2017-10-24

for BAM in ../05.mergeAsplit/*.dedup.bam; do
  CHR=`basename $BAM | sed s/.dedup.bam//g`
  echo $CHR
  wait_yhq
  yhbatch -n 1 -e ${CHR}.e -o ${CHR}.o ./csp.sh $CHR
done