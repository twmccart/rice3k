#!/bin/bash


# function that will generate the gvcf file for a given cultivar

function call_variants() {
	which java
	echo $gatk

	CULT=$1

	$gatk -T HaplotypeCaller \
      		-R $reference \
      		-I $maps/${CULT}.realigned.bam \
      		-ERC GVCF \
      		-o $calls/${CULT}.g.vcf \
    		-nct 8
}

function genotype() {
	which java
	echo $gatk

	CULT=$1

	$gatk -T GenotypeGVCFs \
      		-R $reference \
      		-V $calls/${CULT}.g.vcf \
      		-o $calls/${CULT}.vcf
}

function full_genotype() {
	which java
	echo $gatk

	CULT=$1

	$gatk -T GenotypeGVCFs \
      		-R $reference \
      		-V $calls/${CULT}.g.vcf \
            -allSites \
      		-o $calls/${CULT}.full.vcf \
	    	-nt 24
}


# same as above, but only for a given loci rather than a full chromosome
# loci can be in the gatk format chrname:begin-end
# or it can be a file with a list of such formatted loci

function call_variants_range() {
	which java
	echo $gatk

	CULT=$1
	LOCI=$2

	$gatk -T HaplotypeCaller \
		-L $LOCI \
      		-R $reference \
      		-I $maps/${CULT}.realigned.bam \
      		-ERC GVCF \
      		-o $calls/${CULT}-${LOCI}.g.vcf
}


function clean_vcf() {
	cultivar=$1
	bcftools view --exclude-uncalled --exclude-types 'indels' -O z -o ${cultivar}.cleaned.vcf.gz ${cultivar}.full.vcf
	tabix ${cultivar}.cleaned.vcf.gz
}
