#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# there should be changed to samples you want to merge
CHRlst=$1
SMlst=$(cat ../sample.txt|xargs|tr " " ",")
WP="/WORK/pp192/ReseqPJ/mergeTest"

source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge GVCF"

FSM=`echo ${SMlst}|cut -d, -f1`

for CHR in `echo ${CHRlst} | tr "," " "`; do
    # make a new file 
    echo '#!/usr/bin/env sh' > mergeGVCF_${CHR}.sh
    for GVCF in `ls ${WP}/${FSM}/07.splitGVCF/${CHR}.??.g.vcf.gz`;do
        ID=`basename $GVCF|cut -d. -f3`
        echo ${ID}
        #
        (echo "java -Xmx3g \\";
        echo "  -Djava.io.tmpdir=/WORK/pp192/tmp \\";
        echo "  -jar /WORK/pp192/Install/GenomeAnalysisTK.jar \\";
        echo "  -R /WORK/pp192/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \\";
        echo "  -T GenotypeGVCFs \\";) >> mergeGVCF_${CHR}.sh
        #
        (for SM in `echo ${SMlst} | tr "," " "`; do
            echo "  --variant ${WP}/${SM}/07.splitGVCF/${CHR}.${ID}.g.vcf.gz \\";
        done) >> mergeGVCF_${CHR}.sh
        #
        echo "  -o ${CHR}.${ID}.raw.vcf.gz &"  >> mergeGVCF_${CHR}.sh
    done
    echo "wait"  >> mergeGVCF_${CHR}.sh
    yhbatch -n 1 -e ${CHR}.e -o ${CHR}.o ./mergeGVCF_${CHR}.sh
done

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge GVCF"