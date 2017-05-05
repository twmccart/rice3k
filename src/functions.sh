#!/bin/bash


# function that will generate the gvcf file for a given cultivar

function call_variants() {
	which java
	echo $GATK

	CULT=$1

	$GATK -T HaplotypeCaller \
      		-R $REFERENCE \
      		-I $MAPS/${CULT}.realigned.bam \
      		-ERC GVCF \
      		-o $CALLS/${CULT}.g.vcf \
    		-nct 8
}

function genotype() {
	which java
	echo $GATK

	CULT=$1

	$GATK -T GenotypeGVCFs \
      		-R $REFERENCE \
      		-V $CALLS/${CULT}.g.vcf \
      		-o $CALLS/${CULT}.vcf
}

function full_genotype() {
	which java
	echo $GATK

	CULT=$1

	$GATK -T GenotypeGVCFs \
      		-R $REFERENCE \
      		-V $CALLS/${CULT}.g.vcf \
            -allSites \
      		-o $CALLS/${CULT}.full.vcf \
	    	-nt 24
}


# same as above, but only for a given loci rather than a full chromosome
# loci can be in the gatk format chrname:begin-end
# or it can be a file with a list of such formatted loci

function call_variants_range() {
	which java
	echo $GATK

	CULT=$1
	LOCI=$2

	$GATK -T HaplotypeCaller \
		-L $LOCI \
      		-R $REFERENCE \
      		-I $MAPS/${CULT}.realigned.bam \
      		-ERC GVCF \
      		-o $CALLS/${CULT}-${LOCI}.g.vcf
}


function clean_vcf() {
	cultivar=$1
	bcftools view --exclude-uncalled --exclude-types 'indels' -O z -o ${cultivar}.cleaned.vcf.gz ${cultivar}.full.vcf
	tabix ${cultivar}.cleaned.vcf.gz
}
