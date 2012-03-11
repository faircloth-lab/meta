"""
File: get_reads_from_voucher_file.py
Author: Brant Faircloth

Created by Brant Faircloth on 17 December 2011 12:12 PST (-0800)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""


import os
import re
import string
import argparse
import ConfigParser

from seqtools.sequence import fasta
from seqtools.fs.args import FullPaths, is_dir


import pdb

def get_args():
    parser = argparse.ArgumentParser(description="""Parse a voucher species file
        into smaller file(s)""")
    parser.add_argument("input", help = """The input database to match
        against""", action = FullPaths)
    parser.add_argument("conf", help = """The path to a config file containing
        the species in a given group to subset""", action = FullPaths)
    parser.add_argument("output", help = """The directory in which to store the
        output""", action = FullPaths, type = is_dir)
    #parser.add_argument('--raw', help = 'Output the raw psl',
    #        action = 'store_true')
    return parser.parse_args()

def get_search_string(species):
    #search = '.*Alseis[.,!?;\s-]blackiana.*'
    search = ''
    for k,sp in species:
        search += ".*" + sp.replace(" ", "[.,!?;\s-]") + ".*" + "|"
    #strip off extra |
    search = search.rstrip('|')
    return re.compile(search, re.IGNORECASE)

def main():
    args = get_args()
    config = ConfigParser.ConfigParser()
    config.read(args.conf)
    for section in config.sections():
        outf = fasta.FastaWriter(
            os.path.join(
                args.output, section.replace(' ','-')
                ) + ".fasta"
            )
        species = config.items(section)
        regex = get_search_string(species)
        for record in fasta.FastaReader(args.input):
            if regex.search(record.identifier):
                outf.write(record)
        outf.close()

if __name__ == '__main__':
    main()
