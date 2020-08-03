#!/usr/bin/env bash
set -euo pipefail

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Split"

batch=/home/wangzh/Project/reseq/haplo_map
for i in ${batch}/*_1*;do
    SM=`basename $i|cut -d '_' -f1`
    mkdir $SM;
    cd $SM;
    zcat $i | split -l 100000000 - R1_ &
    zcat ${batch}/${SM}_2* | split -l 100000000 - R2_ &
    wait
    for F in R?_??; do
        gzip $F &
        wait_all
    done
    wait
    cd ..
    GZrsync ${SM} &
done
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Split"