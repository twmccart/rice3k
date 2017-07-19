## Genotype Variation in Rice
This is ....

### Requirements
The following software needs to be loaded on your machine:
samtools/1.3
java/1.8.0_40
tabix
bcftools/1.3
vcftools/0.1.13
Python/2.7.x
RAxML

### Installation Instructions
In order to use this workflow, first you need to clone the following repository from GitHub into your working
directory:
`git clone https://github.com/twmccart/rice3k`
`cd ~/rice3k`

Then you need to run this one-time script using the command:
` something `

This script will take care of few things for you. It will install some files necessary to this workflow; like the GATK
tool, the reference genome used in the 3k rice genome project, and it will create some necessary directories.

### How To Use
The workflow has two main steps. In the first step, both the HaplotypeCaller and GenotypeGVCFs tools will be called to
generate a `.vcf` file from the `.bam` files. This simply done by running the `generateVCF` script provided with a
cultivar name. An example of job submission is:

`qsub -l nodes=1:ppn=8,walltime=48:00:00,vmem=30gb -N $CULT.callvariant src/generateVCF`

where `$CULT` is a cultivar name, e.g. IRIS_313-10295.
