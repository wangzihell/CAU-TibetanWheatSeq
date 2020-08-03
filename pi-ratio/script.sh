#!/bin/bash
set -euxo pipefail

PAIR=DE_vs_DO
S1=`echo $PAIR | cut -d"_" -f1`
S2=`echo $PAIR | cut -d"_" -f3`
size=100

for F in chr[1-7][ABD].ann.bcf.gz; do
  CHR=`basename $F | sed s/.ann.bcf.gz//g`
  (for Code in ${S1} ${S2}; do
    vcftools --bcf $F --keep ${Code}.grp --window-pi ${size}000 --out ${CHR}.${Code}.${size}k &
    vcftools --bcf $F --keep ${Code}.grp --window-pi ${size}000 --out ${CHR}.${Code}.${size}k &
  done
  wait
  gawk -vOFS="\t" 'ARGIND==1{A[$1":"$2]=$5;} ARGIND==2&&(FNR>1){if($1":"$2 in A){print $1, $2, $3, A[$1":"$2], $5, A[$1":"$2]/$5, $5/A[$1":"$2];}}' ${CHR}.${S1}.${size}k.windowed.pi ${CHR}.${S2}.${size}k.windowed.pi > ${CHR}.${PAIR}.${size}k.windowed.pi2 ) &
  wait_all
done

wait

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -F"\t" -vOFS="\t" '{print $1,$2,$3,$7}' ${CHR}.${PAIR}.${size}k.windowed.pi2 > ${CHR}.${PAIR}.${size}k.piratio
  sed -i '1i CHROM\tSTART\tEND\tPiDO/PiDE' ${CHR}.${PAIR}.${size}k.piratio # add header 200521
done

cat chr?A.${PAIR}.${size}k.piratio | gawk 'NR==1 || /^chr/' > chrA.${PAIR}.${size}k.piratio
cat chr?B.${PAIR}.${size}k.piratio | gawk 'NR==1 || /^chr/' > chrB.${PAIR}.${size}k.piratio
cat chr?D.${PAIR}.${size}k.piratio | gawk 'NR==1 || /^chr/' > chrD.${PAIR}.${size}k.piratio
gawk -vOFS="\t" '{$4=log($4+1)/log(10);print}' chrA.${PAIR}.${size}k.piratio > chrA.${PAIR}.${size}k.logpiratio
gawk -vOFS="\t" '{$4=log($4+1)/log(10);print}' chrB.${PAIR}.${size}k.piratio > chrB.${PAIR}.${size}k.logpiratio
gawk -vOFS="\t" '{$4=log($4+1)/log(10);print}' chrD.${PAIR}.${size}k.piratio > chrD.${PAIR}.${size}k.logpiratio
