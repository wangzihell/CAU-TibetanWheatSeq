#!/bin/bash

# prepare VCF files
for CHR in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  python 4DTv.py TLTS.${CHR}.vcf TLTS.${CHR}.4dtv.thin.vcf 161010_Chinese_Spring_v1.0_pseudomolecules.fasta &
done
wait

bcftools concat TLTS.chr?D.4dtv.thin.vcf -o DD_snp.4dtv.thin.vcf

# generate SFS files
python easySFS.py -i DD_snp.4dtv.thin.vcf -p samplelist.txt -a --proj 30,40 -f -o DD_snp_4dtv_thin

# calculate using dadi
for i in {1..1000};do
for j in 2epoch_m.py 2epoch.py 3epoch_m1.py 3epoch_m2.py 3epoch_m.py 3epoch.py 3epoch_m_single.py 3epoch_m_single2.py;do
  python ${j} -i DD_snp_4dtv_thin/dadi/TL-TS.sfs -n ${i} -f "yes" &
done
done

# summary results
for j in 2epoch_m 2epoch 3epoch_m1 3epoch_m2 3epoch_m 3epoch 3epoch_m_single 3epoch_m_single2;do
  cat ${j} | gawk -F"," '{N=NF-1;print $N"\t"$0}' | sort -k1,1nr |cut -f2-> ${j}.summary.csv
done

# estimate theta
for CHR in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  bcftools view ${CHR}.bcf.gz|grep -v "#"|cut -f 8|cut -f 1 -d";"|cut -f 2 -d"=" > ${CHR}.AC.txt
  sort ${CHR}.AC.txt|uniq -c |sed 's/^\ \+//'|sort -k2,2n > tmp
  csvtk join -H -f "1;2" -k -d " " --fill 0 index tmp | sponge tmp
  (echo "617 unfold \"TSTL\""; cat tmp|cut -f 1 -d","|tr "\n" " ") > ${CHR}.unfilter.sfs
  n=$(grep ${CHR} chr_length.txt|cut -f 2)
  theta=$(python esti_theta_from_sfs.py ${CHR}.unfilter.sfs)
  gawk -vn=$n -vtheta=$theta -vOFS="\t" 'BEGIN{print theta,n,theta/n/4/6.5*10^9}' >> summary.txt
done

gawk '{sum+=$1;n+=$2} END{print sum/n/4/6.5*10^9}' summary.txt
