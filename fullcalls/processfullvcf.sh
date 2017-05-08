#!/bin/bash
set -euo pipefail

#find . -name "*.full.vcf" print0 | xargs -0 -I{} bgzip -c {} > {}.gz
#find . -name "*.full.vcf.gz" print0 | xargs -0 tabix 

#/home/twmccart/bcftools-1.4/bcftools merge -m both *.vcf.gz -O v -o 6test.merge.vcf
vcf-merge *.vcf.gz > 6test.merge.vcf
bcftools view --exclude-uncalled --exclude-types 'indels' --min-ac 1 -O v -o 6test.merge.cleaned.vcf 6test.merge.vcf
#vcftools --remove-indels --non-ref-af-any 0.5 --max-missing 1 --recode -out 6test.merge.intersect
#vcftools --vcf ${accession}.full.vcf --remove-indels --recode --stdout | bgzip -c > ${accession}.full.noindels.vcf.gz && \
#tabix -p vcf ${accession}.full.noindels.vcf.gz
