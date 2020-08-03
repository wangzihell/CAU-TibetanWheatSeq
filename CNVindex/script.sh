#!/usr/bin/bash
set -euxo pipefail
# --------

PAIR=DE_vs_DO
S1=`echo $PAIR | cut -d"_" -f1`
S2=`echo $PAIR | cut -d"_" -f3`
size=100

if true;then
for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  Rscript CNVindexScan.R -A ${S2}.grp -B ${S1}.grp -c ${CHR} -d
  #sleep 1
done
wait

cat chr?A.CNV-index.xls | gawk 'NR==1 || /^chr/' | cut -f1,2,3,4 > chrA.${PAIR}.CNV-index.txt
cat chr?B.CNV-index.xls | gawk 'NR==1 || /^chr/' | cut -f1,2,3,4 > chrB.${PAIR}.CNV-index.txt
cat chr?D.CNV-index.xls | gawk 'NR==1 || /^chr/' | cut -f1,2,3,4 > chrD.${PAIR}.CNV-index.txt
fi
