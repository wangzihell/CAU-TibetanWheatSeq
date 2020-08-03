#!/usr/bin/env bash

for BAM in ../05.mergeAsplit/*.dedup.bam; do
    CHR=`basename $BAM | sed s/.dedup.bam//g`
    echo $CHR
    (samtools view -h ../05.*/${CHR}.dedup.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${CHR}.uniq.bam
    samtools index ${CHR}.uniq.bam
    samtools flagstat ${CHR}.uniq.bam > ${CHR}.flagstat)&
done