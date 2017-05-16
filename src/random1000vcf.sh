#!/bin/bash
set -euo pipefail

file=$1

header_line_number=$(grep --line-number "^#CHROM" ${file} | cut -d":" -f1)
let vcf_line_number=$header_line_number+1

# The header region of the VCF file needs to be preserved
head -n "${header_line_number}" ${file} > ${file}.random1000
# The 1000 random lines should only come from legitimate sites
shuf -n 1000 <(tail --lines=+${vcf_line_number} ${file}) >> ${file}.random1000
