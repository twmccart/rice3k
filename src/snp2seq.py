#!/bin/python
# This script originates with Murat Ozturk.


import vcf

a=vcf.Reader(open('out.recode.vcf'))
seq = dict()
for sample in a.samples :
    seq[sample]=''

for record in a :
    for sample in record.samples :
        seq[sample.sample] += sample.gt_bases[0]



from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

bunch = [ SeqRecord(Seq(seq[cult]), cult) for cult in seq]

SeqIO.write(bunch , 'rice.fasta', 'fasta')
