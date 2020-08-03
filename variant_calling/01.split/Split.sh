#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Split"

WP="/WORK/pp192/"

#/WORK/pp192/Project/ReSeqPip/S7/01.split
DIR=`dirname $PWD`
SM=`basename $DIR`

#mkdir $SM; cd $SM

#mkdir 01.split; cd 01.split;

zcat ${WP}/Rawdata/${SM}/*_1*.fq.gz | split -l 100000000 - R1_ &
zcat ${WP}/Rawdata/${SM}/*_2*.fq.gz | split -l 100000000 - R2_ &

wait

rm ${WP}/Rawdata/${SM}/*_[12]*.fq.gz

for F in R?_??; do
  gzip $F &
done

wait

lst=`ls R1_?? | sed s/R1_//g | tr '\n' ',' | sed s/,$//g`

#cd ..


echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Split"

cd ../02.*/
sh ./02.trimmomatic.sh

