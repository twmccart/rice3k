## Description

This workflow pulls mapped read data from [SNP-seek](http://snp-seek.irri.org/) database and calls variants using the more recent and recommended HaplotypeCaller tool and best practices from [GATK](https://software.broadinstitute.org/gatk/) rather than the obsolete UnifiedGenotyper tool that was used to build the SNP-Seek database. 

## Requirements

  * Mapped reads for 2999 cultivars require around 15 TB of space
  * Variant calling requires Picard and GATK
