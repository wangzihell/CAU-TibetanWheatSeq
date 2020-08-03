#!/usr/bin/env sh
# Guo, Weilong; guoweilong@126.com; 2017-10-24

n=$(wc -l ../01.split/ID.list | cut -f 1 -d" ")
if (( $n < 25 ));then
    yhbatch -n 1 -e trim.e -o trim.o ./trimmomatic.sh ../01.split/ID.list
elif [[ $n < 50 ]];then
    yhbatch -n 1 -e trim1.e -o trim1.o ./trimmomatic.sh ../01.split/ID.list.part1
    yhbatch -n 1 -e trim2.e -o trim2.o ./trimmomatic.sh ../01.split/ID.list.part2
else 
    yhbatch -n 1 -e trim1.e -o trim1.o ./trimmomatic.sh ../01.split/ID.list.part1
    yhbatch -n 1 -e trim2.e -o trim2.o ./trimmomatic.sh ../01.split/ID.list.part2
    yhbatch -n 1 -e trim3.e -o trim3.o ./trimmomatic.sh ../01.split/ID.list.part3
fi