#! /bin/bash

# script to download bam and bai files from SNP-Seek servers

cult=$1
mkdir -p ${cult}
cd ${cult}

wget https://s3.amazonaws.com/3kricegenome/Nipponbare/${cult}.realigned.bam
wget https://s3.amazonaws.com/3kricegenome/Nipponbare/${cult}.realigned.bam.bai
wget https://s3.amazonaws.com/3kricegenome/Nipponbare/${cult}.realigned.bam.stats

# this part requires samtools to be available 
# module load samtools/1.3
[ -f ${cult}.realigned.bam ] && ! [ -f ${cult}.realigned.bam.bai ] && samtools index ${cult}.realigned.bam

cd ..
