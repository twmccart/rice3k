#!/bin/bash
set -euo pipefail

#find . -name "*.full.vcf" print0 | xargs -0 -I{} bgzip -c {} > {}.gz
#find . -name "*.full.vcf.gz" print0 | xargs -0 tabix 

#vcf-merge *.vcf.gz | tee 6test.merge.vcf | bcftools view --exclude-uncalled --exclude-types 'indels' --min-ac 1 --genotype ^miss -O u | bcftools view --genotype ^het -O v | awk ' /^#/ {print} length($4) == 1 {print} ' > 6test.merge.cleaned.short.awk.vcf

bcftools view --exclude-uncalled --exclude-types 'indels' --min-ac 1 --genotype ^miss -O u $1 | bcftools view --genotype ^het -O v | awk ' /^#/ {print} length($4) == 1 {print} ' > ${1}.cleaned
