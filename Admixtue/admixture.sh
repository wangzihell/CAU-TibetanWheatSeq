#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

CHR=$1
K=$2
SEED=$3
mkdir -p S${SEED}; cd S${SEED};
admixture -j16 --seed=${SEED} ../${CHR}.bed ${K}
