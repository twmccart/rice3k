#!/bin/bash
set -euo pipefail

#find . -name "*.full.vcf" print0 \
#	| parallel -0 { \
#	bcftools view --exclude-uncalled --exclude-types 'indels' -O z -o {.}.cleaned.vcf.gz {} ; \
#	tabix {.}.cleaned.vcf.gz \
#	}

bcftools merge -m both *.noindels.nouncalled.vcf.gz -O v > 6test.bcftoolsmerge.vcf
# -m2 is minimum of 2 alleles listed. 
bcftools view --trim-alt-alleles -m 2 -O v 6test.bcftoolsmerge.vcf > 6test.bcftoolsview.intersect.vcf
