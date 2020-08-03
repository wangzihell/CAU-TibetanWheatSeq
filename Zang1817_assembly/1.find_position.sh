#!/usr/bin/evn sh
# Guo, Weilong; guoweilong@126.com; 2018-03-15

# Generate Reads
GenerateTilingReads_PE -i IWGSCv1.fa -l 500 -o IWGSCv1Simu -I 500 -s 1000 &
GenerateTilingReads_PE -i Zang1817.fa -l 500 -o Zang1817Simu -I 500 -s 1000 &
#

wait

# build bwa index for two genomes
bwa index IWGSCv1.fa -p IWGSCv1_bwa
bwa index Zang1817.fa -p Zang1817_bwa

# map simulated reads to each other
SM="Zang1817Simu"
bwa mem -t 20 \
  -R '@RG\tID:'${SM}'\tLB:'${SM}'\tPL:ILLUMINA\tSM:'${SM} \
  IWGSCv1_bwa \  
  ${SM}_1.fa | \
  samtools view -Sbh - > ${SM}.bam &


SM="IWGSCv1Simu"
bwa mem -t 20 \
  -R '@RG\tID:'${SM}'\tLB:'${SM}'\tPL:ILLUMINA\tSM:'${SM} \
  Zang1817_bwa \
  ${SM}_1.fa | \
  samtools view -Sbh - > ${SM}.bam &


wait

# The 1st genome refer to contigs
samtools view Zang1817Simu.bam | cut -f1,3,4 > MapRes_Zang1817.xls   &
# The 2nd genome shall be chr-level
samtools view IWGSCv1Simu.bam  | cut -f1,3,4 > MapRes_IWGSCv1.xls    &

wait

# -- Mapping Statistics --

# -- Find the chr for each scaffold --
gawk -F"_" '{print $2"\t"$4;}' MapRes_Zang1817.xls | gawk '($3~/^[Cc]hr/)' > Zang1817_pos_IWGSCv1_pos &
gawk -F"_" '($2~/^[Cc]hr/){print $2"\t"$4;}' MapRes_IWGSCv1.xls | gawk '$3!="*"' > IWGSCv1_pos_Zang1817_pos &
wait

#
gawk -vOFS="\t" 'ARGIND==1{A[$1]} ARGIND==2&&($3 in A){print $3, $4, $1, $2;}' Zang1817_pos_IWGSCv1_pos IWGSCv1_pos_Zang1817_pos > scf_pos_chr_pos.xls

#  ===

cut -f1,3 scf_pos_chr_pos.xls | sort | uniq -c | gawk -vOFS="\t" '{print $2, $3, $1;}' | sort -k3nr > scaffold_chr_count.xls
gawk '!/\*/{if($1 in A){if($3==B[$1]){print $1"\t"$2"\t"$3;}}else{print $1"\t"$2"\t"$3; A[$1]=$2; B[$1]=$3;}}' scaffold_chr_count.xls > scaffold_chr_maxcount.xls

# when using threshold count>=22, the scaffold and chr pair is unique
gawk '$3>=22' scaffold_chr_maxcount.xls | cut -f1 | sort | uniq -c | gawk '{print $1}' | uniq

gawk 'ARGIND==1&&($3>=22){A[$1]=$2; if($1 in A){B[$1]}} ARGIND==2{if(($1 in A) &&($1 in B==0) &&A[$1]==$3){T[$1]+=$4; C[$1]++;}} END{for(scf in A){print scf"\t"A[scf];}}' scaffold_chr_maxcount.xls scf_pos_chr_pos.xls > scaffold_chr_matched ; wc -l scaffold_chr_matched

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7; do
  for S in A B D; do
    gawk -vCHR=${CHR} 'ARGIND==1&&($2==CHR){A[$1];} ARGIND==2&&($1 in A)&&($3==CHR){print;}' scaffold_chr_matched scf_pos_chr_pos.xls | Sort_chr_pos -c 1 -p 2 > scaffold_chr_pos.${CHR}
    for SCF in `cut -f1 scaffold_chr_pos.${CHR} | sort -u`; do
      Median=`gawk -vSCF=$SCF '$1==SCF{print $4;}' scaffold_chr_pos.${CHR} | GetMedian.R`
      Direction=`gawk -vSCF=$SCF '$1==SCF{print $2"\t"$4;}' scaffold_chr_pos.${CHR} | ./CheckDirection.R | gawk '{if($1==1){print "+"}else{print "-";}}'`
      echo ${SCF}";$CHR;"${Median}";"${Direction}
    done | tr ';' '\t' | sort -k3g > scf_chr_median_dir.${CHR}
    
    gawk -vOFS="\t" 'BEGIN{Start=0;} ARGIND==1{L[$1]=$2;} ARGIND==2&&($1 in L){print $1, Start, Start+L[$1], L[$1], $2, $3, $4; Start=Start+L[$1];}' ../NRGeneStat/scanffold_len_pct.xls scf_chr_median_dir.${CHR} > scf_start_end_len_chr_median_dir.${CHR}
    
    
    gawk -vOFS="\t" 'ARGIND==1{CS[$1]=$2; CE[$1]=$3; D[$1]=$7;} ARGIND==2&&($1 in CS){if(D[$1]=="+"){print CS[$1]+$2, $4;}else{print CE[$1]-$2, $4;}}' scf_start_end_len_chr_median_dir.${CHR} scaffold_chr_pos.${CHR} > scf_chr_dotplot.${CHR}
    
    gawk -vCHR=${CHR} '$1==CHR&&(!/gap/){print $1"\t"$2;}' /data/genome/wheat/CS_IWGSC/v1/161010_Chinese_Spring_v1.0_pseudomolecules_to_scaffolds.bed > chr_scf_pos.${CHR}
  done
done
