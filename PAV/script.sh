#!/bin/bash

# map CS and Zang1817 reads to CS and Zang1817 genome seprately.
for SM in CS Zang1817;do
bwa mem -t 20 \
  -R '@RG\tID:'${SM}'\tLB:'${SM}'\tPL:ILLUMINA\tSM:'${SM} \
  IWGSCv1_pesudo_bwa \
  ${SM}_1.fq ${SM}_2.fq | \
  samtools view -Sbh - > ${SM}_to_CS.bam
done

for SM in CS Zang1817;do
bwa mem -t 20 \
  -R '@RG\tID:'${SM}'\tLB:'${SM}'\tPL:ILLUMINA\tSM:'${SM} \
  Zang1817_bwa \
  ${SM}_1.fq ${SM}_2.fq | \
  samtools view -Sbh - > ${SM}_to_Zang1817.bam
done

# count depth
for SM in CS Zang1817;do
for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  bedtools coverage -a CS_${CHR}.bed -b ${CHR}.${SM}_to_CS.bam -counts > ${CHR}.${SM}_to_CS.DP
  bedtools coverage -a Zang1817_${CHR}.bed -b ${CHR}.${SM}_to_CS.bam -counts > ${CHR}.${SM}_to_CS.DP
done
done

# normlize
for SM in CS_to_CS CS_to_Zang1817 Zang1817_to_CS Zang1817_to_Zang1817;do

sum=$(gawk -vc=$DPcol '{sum+=$c} END{print sum}' ${SM}/*.DP)
line=$(wc -l ${SM}/*.DP|tail -n 1|gawk '{print $1}')
ave=$(perl -e "print $sum/$line")

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -vOFS="\t" -vave=$ave -vc=$DPcol '{print $c/ave}' ${SM}/${CHR}.1.DP > ${SM}/${CHR}.norm
  gawk -vOFS="\t" -vave=$ave -vc=$DPcol '{print $c/ave}' ${SM}/${CHR}.2.DP >> ${SM}/${CHR}.norm

done
done

# call PAV
for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  # to formated bed
  paste ../CS_to_CS/${CHR}.${bin}.norm ../Zang1817_to_CS/${CHR}.${bin}.norm | gawk -vOFS="\t" -vc=$CHR -vb=$bin1 '$1>=0.8 && ($1==0 || ($1-$2)/$1>=0.8){print c,(NR-1)*b-2501,NR*b+2501}' > ${CHR}.${bin}.Zang1817_absbin_ext2k5.bed
  sed -i 's/-.*\t/0\t/' ${CHR}.${bin}.Zang1817_absbin_ext2k5.bed
  bedtools merge -i ${CHR}.${bin}.Zang1817_absbin_ext2k5.bed |gawk -vOFS="\t" '{$4=$3-$2;if($4>10002){$2=$2+2501;$3=$3-2501;$4=$4-5002;print}}' > ${CHR}.${bin}.Zang1817_absbin_atleasttwobin.bed
done

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  # to formated bed
  paste ../Zang1817_to_Zang1817/${CHR}.${bin}.norm ../CS_to_Zang1817/${CHR}.${bin}.norm | gawk -vOFS="\t" -vc=$CHR -vb=$bin1 '$1>=0.8 && ($1==0 || ($1-$2)/$1>=0.8){print c,(NR-1)*b-2501,NR*b+2501}' > ${CHR}.${bin}.CS_absbin_ext2k5.bed
  sed -i 's/-.*\t/0\t/' ${CHR}.${bin}.CS_absbin_ext2k5.bed
  bedtools merge -i ${CHR}.${bin}.CS_absbin_ext2k5.bed |gawk -vOFS="\t" '{$4=$3-$2;if($4>10002){$2=$2+2501;$3=$3-2501;$4=$4-5002;print}}' > ${CHR}.${bin}.CS_absbin_atleasttwobin.bed
done
