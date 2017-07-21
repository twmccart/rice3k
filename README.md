## Genotype Variation in Rice
This workflow is designed to study the relationship among cultivars obtained from the 3k Rice Genome project. It utilizes the GATK 
tools developed by the Broad Institute to generate VCF files containing all of the called sites in a given cultivar. It also 
generate VCF files for each chromosome in the cultivar with indels, Multiple Nucleotide Polymorphisms, and any uncalled or 
heterozygous sites removed. The output of this workflow is a tree file showing the relationship among all processed cultivars.  

### Requirements
The following software needs to be installed on your machine. The versions listed have been tested successfully:
Samtools version 1.3.x
Oracle Java version 1.8.0_131 OR OpenJDK Java version 1.8.0_40 
Tabix version 0.2.6
bcftools version 1.3.* 
vcftools version 0.1.13 (version 0.1.14 does not work)
Python version 2.7.x
Biopython version 1.70
PyVCF version 0.6.8
RAxML version 8.2.10

### Installation Instructions
In order to use this workflow, first you need to clone the following repository from GitHub into your working
directory:
`git clone https://github.com/twmccart/rice3k`
`cd ~/rice3k`

Then you need to run one-time this script:
`bash INSTALL`

This script will take care of few things for you. It will install some files necessary to this workflow; like the GATK 
tools themselves, the reference genome used in the 3k rice genome project, and it will create some necessary directories. 

### How To Use
The workflow has two main steps. In the first step, both the HaplotypeCaller and GenotypeGVCFs tools will be called to
generate a `.vcf` file from the `.bam` files. This simply done by running the `generatevcf` script provided with a
cultivar name. An example of job submission is:

`bash generatevcf cultivar`

where `cultivar` is any cultivar name, e.g. IRIS_313-15900. This step can very time-consuming; it might take up to 48 hours when run on an 8-cpu computer with a 
30GB of ram.

In the second step, the script will take all the files generated using the `generatevcf` script as input files, and merge the 
chromosomes in parallel. It will also strip any site that doesn't have SNPs, and concatenate them into a new VCF file containing 
only the sites with SNPs. A subset of these SNPs are chosen at random and used as an alignment from which a maximum-likelihood tree 
is generated.  

For this step you need to make a list containing the names of all cultivars that need to be processed, then call the script 
`generatetree` on that list, as follows:

`bash generatetree cultivarlist`

**Note**: Here it took about 10 hours to process 20 cultivars. If more cultivars being processed, then more time might be needed.

At this step you should end up with a single file ready to be viewed by a tree viewer program. For instance, we used the **FigTree 
v1.4.3** to visualize our tree.
