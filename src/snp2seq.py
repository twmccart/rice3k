#!/usr/bin/env python
# This script originates with Murat Ozturk.

from sys import argv
import vcf
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main(arguments):
    """Converts a multi-sample VCF file into a FASTA alignment

    Args:
        arguments: a list that should contain 2 strings, the name of the invoked
            program, and the name of input VCF file
    Returns:
        Writes out a FASTA file containing sequences for each sample in the input
        VCF file. These function as an alignment for the purposes of treeing.
    """

    input_file_name=arguments[1]
    with open(input_file_name) as input_file:
        vcf=vcf.Reader(input_file)
        sequences={}
        for sample in vcf.samples:
            seq[sample]=''

        for record in vcf:
            for sample in record.samples:
                seq[sample.sample] += sample.gt_bases[0]



        bunch = [ SeqRecord(Seq(seq[cult]), cult) for cult in seq]

        SeqIO.write(bunch , 'rice.fasta', 'fasta')

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
