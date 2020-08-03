#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

CHR=$1
SMlst=$2

WP="/WORK/pp192/ReseqPJ/mergeTest"
source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start SplitGVCF"

for SM in `echo ${SMlst}|tr "," " "`;do
    (cd ${WP}/$SM/07.splitGVCF;
    ln -s ../../07.splitGVCF/* .
    python SplitGvcfByBin.py -i ../06.callsnp/${CHR}*.vcf.gz -B 20000000;
    for ID in ${CHR}.??.g.vcf;do
        bgzip ${ID};
        tabix -p vcf ${ID}.gz
    done
    cd ../../ )&
done
#
wait

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish SplitGVCF"