#!/bin/bash
set -euo pipefail

# Takes a vcf file and an integer, outputs a VCF file with $2 random sites from the original

file=$1
lines=$2

header_line_number=$(grep --line-number "^#CHROM" ${file} | cut -d":" -f1)
let vcf_line_number=$header_line_number+1

# The header region of the VCF file needs to be preserved
head -n "${header_line_number}" ${file} #> ${file}.random1000
# The 1000 random lines should only come from legitimate sites
shuf -n $lines <(tail --lines=+${vcf_line_number} ${file}) #>> ${file}.random1000
