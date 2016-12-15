#!/bin/bash


GATK="java -jar ../../GenomeAnalysisTK.jar"
REFERENCE="../reference/all.chrs.fix.fasta" 
#CULT=$1
LOCI="chr12:17290000-17320000"
TARGET="NbSWEET13"

module load java
module load tabix
module load bcftools/1.3

ROOT="/N/u/muroztur/Karst/bulkdata/SNP-Seek/"
cd ${ROOT};
while read CULT; do
cd ${ROOT};
cd ${CULT} &&
[ ! -f ${TARGET}.vcf.gz ] &&
${GATK} -L ${LOCI} -R ${REFERENCE} -T HaplotypeCaller  -I ${CULT}.realigned.bam -ERC GVCF -o ${TARGET}.g.vcf &&
${GATK} -L ${LOCI} -R ${REFERENCE} -T GenotypeGVCFs  --variant ${TARGET}.g.vcf -allSites -o ${TARGET}.vcf &&

bgzip -f ${TARGET}.vcf &&
tabix -f ${TARGET}.vcf.gz ;
done  < cultivars

cd ${ROOT};
bcftools merge */${TARGET}.vcf.gz -o ${TARGET}.vcf
