#!/usr/bin/env bash
# Guo, Weilong; guoweilong@126.com; 2018-04-02

# combine sequence of each scaffold according to the order and orientation defined in Previous step.
for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D; do
  echo $CHR
  (cut -f1,4 ../../Dotplot_Zang1817_IWGSCv1/scf_chr_median_dir.${CHR} > scf_dir.${CHR}
  FastaOrderByList -i ../Zang1817.fa -c scf_dir.${CHR} -o ${CHR}.order.fa -l 100
  ThreadFasta -i ${CHR}.order.fa -c `echo $CHR | sed s/chr//g`  -g 300 -M 1000000000 -d ${CHR}.contig.dict -o ${CHR}.fa
  sed 's/_1//g' ${CHR}.fa | gzip > ${CHR}.fa.gz
  rm ${CHR}.order.fa ${CHR}.fa
  ) &
  sleep 20
done

for F in *.dict; do
  echo $F
  sed -i s/_1//g $F
done

# extract unplaced scaffols
cat ../../Dotplot_Zang1817_IWGSCv1/scf_chr_median_dir.chr[1-7][ABD] | cut -f1 > chrABD_scf
grep '^>' ../Zang1817.fa | tr -d '>'  > Zang1817_scf
gawk 'ARGIND==1{A[$1];} ARGIND==2{if($1 in A==0){print;}}' chrABD_scf Zang1817_scf | gawk '{print $1"\t+;"}' > chrUn_scf

# combine sequence of unplaced scaffols to form chrUn.
FastaOrderByList -i ../Zang1817.fa -c scf_dir.Un -o chrUn.order.fa -l 100
ThreadFasta -i chrUn.order.fa -c 1A -g 300 -M 1000000000 -d chrUn.contig.dict -o chrUn.fa
sed 's/_1//g' chrUn.fa | gzip > chrUn.fa.gz
cat chrUn.fa | sed s/chr1A_1/chrUn/g | gzip > chrUn.fa.gz
rm chrUn.order.fa chrUn.fa

# combine sequence of each pseudo-chromosome.
zcat chr1A.fa.gz chr1B.fa.gz chr1D.fa.gz chr2A.fa.gz chr2B.fa.gz chr2D.fa.gz chr3A.fa.gz chr3B.fa.gz chr3D.fa.gz chr4A.fa.gz chr4B.fa.gz chr4D.fa.gz chr5A.fa.gz chr5B.fa.gz chr5D.fa.gz chr6A.fa.gz chr6B.fa.gz chr6D.fa.gz chr7A.fa.gz chr7B.fa.gz chr7D.fa.gz chrUn.fa.gz | gzip > Zang1817_PseudoGenome.fa.gz
