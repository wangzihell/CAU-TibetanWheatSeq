#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

source /WORK/app/osenv/ln1/set2.sh

WP="/WORK/pp192/"

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start filter and merge GVCF"

CHRlst=$1

for CHR in `echo ${CHRlst} | tr "," " "`; do
    echo ${CHR}
    echo '#!/usr/bin/env sh' > filterGVCF_${CHR}.sh
    echo 'set -euxo pipefail' >> filterGVCF_${CHR}.sh
    echo '#SBATCH -N 1' >> filterGVCF_${CHR}.sh
    echo '#SBATCH -n 1' >> filterGVCF_${CHR}.sh
    echo '#SBATCH -c 24' >> filterGVCF_${CHR}.sh
    echo '#SBATCH -t 2400' >> filterGVCF_${CHR}.sh
    echo 'WP="/WORK/pp192/"' >> filterGVCF_${CHR}.sh
    for GVCF in `ls ../08.mergeGVCF/${CHR}.??.raw.vcf.gz`;do
        ID=`basename $GVCF|cut -d. -f3`
        echo ${ID}
        #
        sed "s/\$CHR/$CHR/g;s/\$ID/$ID/g" VCFfilter.sh > VCFfilter_${CHR}_${ID}.sh
        cat VCFfilter_${CHR}_${ID}.sh >> filterGVCF_${CHR}.sh
        rm VCFfilter_${CHR}_${ID}.sh
    done
    echo "wait"  >> filterGVCF_${CHR}.sh
    yhbatch -n 1 -e ${CHR}.e -o ${CHR}.o ./filterGVCF_${CHR}.sh
done

echo "["`date +%Y-%m-%d,%H:%M:%S`"] End filter and merge GVCF"