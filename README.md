## Genotype Variation in Rice
This workflow is designed to study the relationship among cultivars obtained from the 3k Rice Genome project. It utilizes the GATK 
tools developed by the Broad Institute to generate VCF files containing all of the called sites in a given cultivar. It also 
generate VCF files for each chromosome in the cultivar with indels, Multiple Nucleotide Polymorphisms, and any uncalled or 
heterozygous sites removed. The output of this workflow is a tree showing the relationship among all processed cultivars.  

### Requirements
The following modules need to be loaded on your machine:
samtools/1.3, java/1.8.0_40, tabix/0.2.6, bcftools/1.3, vcftools/0.1.13, and python/2.7.3

### Installation Instructions
In order to use this workflow, first you need to clone the following repository from GitHub into your working 
directory: 
`git clone https://github.com/twmccart/rice3k`
`cd ~/rice3k`

Then you need to run this one-time script:
`bash INSTALL`

This script will take care of few things for you. It will install some files necessary to this workflow; like the GATK 
tool, the reference genome used in the 3k rice genome project, and it will create some necessary directories. 

### How To Use
The workflow has two main steps. In the first step, both the HaplotypeCaller and GenotypeGVCFs tools will be called to 
generate a `.vcf` file from the `.bam` files. This simply done by running the `generateVCF` script provided with a 
cultivar name. An example of job submission is:

`qsub -l nodes=1:ppn=8,walltime=48:00:00,vmem=30gb -N $CULT.genotype src/generatevcf`

where `$CULT` is a cultivar name, e.g. IRIS_313-10295.

In the second step, the script will take all the files generated using the `generatevcf` script as input files, and merge the 
chromosomes in parallel. It will also strip any site that doesn't have SNPs, and concatenate them into a new VCF file containing 
only the sites with SNPs. A subset of these SNPs are chosen at random and used as an alignment from which a maximum-likelihood tree 
is generated.  

So for this step you need to make a list containing thee names of all cultivars that need to be processed, then call the script 
`generatetree` on that list. For example:

`qsub -l nodes=1:ppn=1,walltime=10:00:00,vmem=10gb -N $CULT.tree src/generatetree`

**Note**: The walltime assigned was based on processing 20 cultivars. If more cultivars being processed, then more time might be 
needed. 

At this step you should end up with a single file ready to be viewed by a tree viewer program. For instance, we used the **FigTree 
v1.4.3** to visualize out tree.
