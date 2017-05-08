#!/bin/bash
set -euo pipefail

#find . -name "*.full.vcf" print0 | xargs -0 -I{} bgzip -c {} > {}.gz
#find . -name "*.full.vcf.gz" print0 | xargs -0 tabix 

vcf-merge *.vcf.gz > 6test.merge.vcf
vcftools --remove-indels --non-ref-af-any 0.5 --max-missing 1 --recode -out 6test.merge.intersect

#/home/twmccart/bcftools-1.4/bcftools merge -m both *full.vcf.gz -O v > 6test.bcftoolsmerge.vcf
# -m2 is minimum of 2 alleles listed. 
#/home/twmccart/bcftools-1.4/bcftools view --trim-alt-alleles --exclude-uncalled --exclude-types 'indels' -m 2 -O v 6test.bcftoolsmerge.vcf > 6test.bcftoolsview.intersect.vcf
