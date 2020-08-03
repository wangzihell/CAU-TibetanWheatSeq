#!/bin/bash
set -x

WD=.
for CHR in DD AABB AABBDD;do
smartpca.perl \
    -i ${WD}/${CHR}.bed \
    -a ${WD}/${CHR}.bim \
    -b ${WD}/${CHR}.fam \
    -o ${WD}/${CHR}.pca \
    -p ${WD}/${CHR}.plot \
    -e ${WD}/${CHR}.eval \
    -l ${WD}/${CHR}.log &
done

wait

for CHR in AABB AABBDD DD;do
  cat ${CHR}.pca.evec|tr -s " " "\t"|sed 's/#//'|cut -f1,13 --complement|cut -d":" -f2|sed '1s/^/exp/;1i Sample\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10' > ../../dev/${CHR}.pca
done
