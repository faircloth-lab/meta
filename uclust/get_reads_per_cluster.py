"""
File: get_reads_per_cluster.py
Author: Brant Faircloth

Created by Brant Faircloth on 15 December 2011 14:12 PST (-0800)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import argparse
from collections import defaultdict
from seqtools.fs.args import FullPaths, is_dir

import pdb

def get_args():
    parser = argparse.ArgumentParser(description="Get counts of reads in UCLUST clusters")
    parser.add_argument('input', help = """The input directory containing fastas
            to cluster""", action = FullPaths, type = is_dir)
    parser.add_argument('--sorted', help = """Sort by counts rather than clusters""", 
        action = "store_true")
    return parser.parse_args()


def main():
    args = get_args()
    all_counts = {}
    mx = 0
    for file in glob.glob(os.path.join(args.input, "*.results.uc")):
        name = os.path.basename(file)
        counts = defaultdict(set)
        for line in open(file, "rU"):
            if not line.startswith("#"):
                ls = line.strip().split('\t')
                counts[int(ls[1])].add(ls[8].split(' ')[0])
        if max(counts.keys()) > mx:
            mx = max(counts.keys())
        all_counts[name] = counts
    for file in sorted(all_counts.keys()):
        #use None here instead of '', so will sort correctly
        spacer = [None] * (mx + 1)
        for slot in sorted(all_counts[file].keys()):
            spacer[slot] = len(all_counts[file][slot])
        if args.sorted:
            spacer = [str(s) for s in sorted(spacer, reverse = True)]
        else:
            spacer = [str(s) for s in spacer]
        # replace None with "" to => pretty
        print "{0},{1}".format(file,','.join(spacer).replace("None", ""))
            

if __name__ == '__main__':
    main()
    