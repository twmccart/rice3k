#! /bin/bash

# script to download bam and bai files from SNP-Seek servers
# these save us some time because we do not have to map the reads
# alternatively, we could get the raw reads using the xgetreads script from Tregga
# https://github.com/BrendelGroup/TRegGA/tree/master/reads

cult=$1
#mkdir -p ${cult}
#pushd ${cult}

wget https://s3.amazonaws.com/3kricegenome/Nipponbare/${cult}.realigned.bam

# this part requires samtools to be available 
# module load samtools/1.3
samtools index ${cult}.realigned.bam

#popd

