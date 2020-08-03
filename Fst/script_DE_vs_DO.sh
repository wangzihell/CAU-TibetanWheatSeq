#!/bin/bash

PAIR=$(basename $(pwd))
S1=DE
S2=DO
size=200

# do computation

cut -f1 /data2/rawdata2/GRP/final/${S1}.grp > ${S1}.grp
cut -f1 /data2/rawdata2/GRP/final/${S2}.grp > ${S2}.grp

for BCF in chr[1-7][ABD].ann.bcf.gz; do
  CHR=`basename $BCF | sed s/.ann.bcf.gz//g`
  echo $CHR
  vcftools --bcf $BCF --weir-fst-pop ${S1}.grp --weir-fst-pop ${S2}.grp  --fst-window-size 200000 --fst-window-step 200000 --out ${CHR}.${PAIR}.${size}k &
  wait_all
done

wait
