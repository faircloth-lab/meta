#!/usr/bin/env python
# encoding: utf-8
"""
File: fasta_qual_to_fastq.py
Author: Brant Faircloth

Created by Brant Faircloth on 18 October 2013 11:10 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import argparse
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Convert fasta+qual to fastq""")
    parser.add_argument(
            "--fasta",
            required=True,
            help="""The fasta file to convery"""
        )
    parser.add_argument(
            "--qual",
            required=True,
            help="""The quality file to convert"""
        )
    parser.add_argument(
            "--fastq",
            required=True,
            help="""The output fastq to generate"""
        )
    return parser.parse_args()

args = get_args()
with open(args.fastq, "w") as outf:
    records = PairedFastaQualIterator(open(args.fasta), open(args.qual))
    count = SeqIO.write(records, outf, "fastq")
print "Converted {} records".format(count)
