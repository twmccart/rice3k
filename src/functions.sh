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
		-o $calls/${CULT}.vcf \
		-nt 24
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
	vcf=$1
	bgzip ${vcf}
	tabix ${vcf}.gz
	bcftools view --exclude-uncalled --exclude-types 'indels' --genotype ^het -O v ${vcf}.gz | awk ' /^#/ {print} length($4) == 1 {print} ' > ${vcf%.vcf}.noindels_hets_mnps.vcf
	#tabix ${cultivar}.cleaned.vcf.gz
}

function clean_and_split_vcf() {
	cultivar=$1
	bgzip ${cultivar}.full.vcf
	tabix -f ${cultivar}.full.vcf.gz
	for chromosome in chr{01,02,03,04,05,06,07,08,09,10,11,12}; do (
		bcftools view --exclude-uncalled --exclude-types 'indels' --genotype ^het -r ${chromosome} -O v ${cultivar}.full.vcf.gz | awk ' /^#/ {print} length($4) == 1 {print} ' | bgzip -c > ../split/${cultivar}.${chromosome}.noindels_hets_mnps.vcf.gz; tabix ../split/${cultivar}.${chromosome}.noindels_hets_mnps.vcf.gz) &
	done
	wait
	echo "SPLIT IS FINISHED"
}

function merge() {
	chromosome=$1
	vcf-merge *${chromosome}*.vcf.gz > ../merges/${chromosome}.merge.vcf
	# The awk command filters any Multiple Nucleotide Polymorphisms, which are apparently a thing
	< ../merges/${chromosome}.merge.vcf bcftools view --exclude-uncalled --exclude-types 'indels' --min-ac 1 --max-af 0.99 --genotype ^miss -O v | awk ' /^#/ {print} length($4) == 1 {print} ' > ../merges/${chromosome}.merge.cleaned.vcf
}

function refilter_merged() {
	file=$1
	< $file bcftools view --min-ac 2 -O v > ${file%.vcf}.min2.vcf
	< $file bcftools view --min-ac 3 -O v > ${file%.vcf}.min3.vcf
	< $file bcftools view --min-ac 4 -O v > ${file%.vcf}.min4.vcf
	< $file bcftools view --min-ac 5 -O v > ${file%.vcf}.min5.vcf
	< $file bcftools view --min-ac 8 -O v > ${file%.vcf}.min8.vcf
}
