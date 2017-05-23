#!/bin/bash
set -euo pipefail

file=$1

for gene in gene{1..32914}; do
	echo "$gene	$(grep -c "$gene[,	]" $file)" >> ${file%.vcf}.geneSNPcount.tsv
	done

