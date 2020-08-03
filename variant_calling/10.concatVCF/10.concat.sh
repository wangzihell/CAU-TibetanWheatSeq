#!/usr/bin/env bash

# without parts, ie. chr1A,chr1B,chr1D

CHRl=$1
for CHR in `echo $CHRl|tr "," " "`;do
    yhbatch -n 1 -e concat_${CHR}.e -o concat_${CHR}.o ./concatVCF.sh ${CHR}
done