#!/usr/bin/env bash
set -euxo pipefail

# do admixture analyze multi rounds
for s in {1..10};do
for i in {2..6};do
  for CHR in DD AABBDD AABB;do
    if [ ! -f S${s}/${CHR}.nocultivar.${i}.P ];then
      bash ./admixture.sh ${CHR}.nocultivar ${i} ${s} &
      wait_all
    fi
  done
done
done

# select round with highest likelihood
for CHR in AABBDD AABB DD;do
for k in {2..6};do
  s=$(tail -n 1 rawdata_log/${CHR}_${k}_likelihood_order.txt | cut -d"." -f 1| cut -d"_" -f3)
  cp S${s}/${CHR}.*${k}.Q select_${CHR}/
done
done

# plot 
for CHR in AABBDD AABB DD;do
  Rscript admixture.R -d ../01.bed/select_${CHR}
done
