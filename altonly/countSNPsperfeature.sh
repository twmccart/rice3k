#!/bin/bash
set -euo pipefail

file=$1
rm -f ${file%.vcf}.promoterSNPcount.tsv
rm -f ${file%.vcf}.geneSNPcount.tsv

for promoter in promoter{0..32914}; do
    echo "$promoter $(grep -c "$promoter[,;	]" $file)" >> ${file%.vcf}.promoterSNPcount.tsv
    done

for gene in gene{0..32914}; do
	echo "$gene	$(grep -c "$gene[,;	]" $file)" >> ${file%.vcf}.geneSNPcount.tsv
	done
