#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

SMlst=S0,S1

source /WORK/app/osenv/ln1/set2.sh

WP="/WORK/pp192/"

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge GVCF"

FSM=`echo ${SMlst}|cut -d, -f1`
FCHR='chr1A.1'

# for CHR in chr1A.1;do
# for CHR in `ls ${WP}/ReseqPJ/${FSM}/07.splitGVCF/*.vcf.gz | cut -d/ -f7 | cut -d. -f1,2 | uniq;do`
for CHR in `zcat ${WP}/ReseqPJ/${FSM}/07.splitGVCF/${FCHR}.aa.g.vcf.gz|grep contig|gawk -F'=|,' '{print $3}'`;do
    echo ${CHR}
    declare -i ID_num=0
    for GVCF in `ls ${WP}/ReseqPJ/${FSM}/07.splitGVCF/${CHR}.??.g.vcf.gz`;do
        ID_num+=1
        if [[ "ID_num" -le 10 ]]; then parts=1; else parts=2; fi
        ID=`basename $GVCF|cut -d. -f3`
        echo ${ID}
        #
        (echo '#!/usr/bin/env sh';
        echo "java  \\";
        echo "  -Djava.io.tmpdir=${WP}/tmp \\";
        echo "  -jar ${WP}/Install/GenomeAnalysisTK.jar \\";
        echo "  -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \\";
        echo "  -T GenotypeGVCFs \\";) >> mergeGVCF_${CHR}_${parts}.sh
        #
        (for SM in `echo ${SMlst} | tr "," " "`; do
            echo "  --variant ${WP}/ReseqPJ/${SM}/07.splitGVCF/${CHR}.${ID}.g.vcf.gz \\";
        done) >> mergeGVCF_${CHR}_${parts}.sh
        #
        echo "  -o ${CHR}.${ID}.raw.vcf.gz &"  >> mergeGVCF_${CHR}_${parts}.sh
    done
    for parts in 1 2;do
        echo "wait"  >> mergeGVCF_${CHR}_${parts}.sh
        yhbatch -n 1 -e ${CHR}_${parts}.e -o ${CHR}_${parts}.o ./mergeGVCF_${CHR}_${parts}.sh
    done
done

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge GVCF"