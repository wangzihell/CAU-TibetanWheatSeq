#!/usr/bin/env bash
# set -x
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

source /WORK/app/osenv/ln1/set2.sh

WP="/WORK/pp192/"

# retry failed samples
CHRlst=$1
sfile=${2:-sample.txt}
nsm=$(wc -l ../${sfile} |cut -f 1 -d" ");

# From help let:

# Exit Status:
# If the last ARG evaluates to 0, let returns 1; let returns 0 otherwise..
let nbatch=$nsm/20

for CHR in `echo ${CHRlst}|tr "," " "`;do
    for (( i = 0; i <= $nbatch; i++ ));do
        batch=$(tail -n +$[i*20+1] ../${sfile}|head -n 20|xargs|tr " " ",")
        echo ${batch}
        echo ${CHR}
        wait_yhq
        yhbatch -n 1 -e splitBatch${i}.e -o splitBatch${i}.o ./splitScript.sh ${CHR} ${batch}
    done
done