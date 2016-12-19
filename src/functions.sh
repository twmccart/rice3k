#!/bin/bash




function call_variants_range {
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


function call_variants {
	which java
	echo $GATK

	CULT=$1

	$GATK -T HaplotypeCaller \
      		-R $REFERENCE \
      		-I $MAPS/${CULT}.realigned.bam \
      		-ERC GVCF \
      		-o $CALLS/${CULT}-${LOCI}.g.vcf
}

