#!/usr/bin/env python
# encoding: utf-8
"""
File: get_fasta_from_clusters.py
Author: Brant Faircloth

Created by Brant Faircloth on 10 March 2012 17:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Build fasta files (for assembly) from reads that
cluster together.

"""

import os
import sys
import glob
import uclust
import argparse
from collections import defaultdict
from seqtools.sequence import fasta

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Build individual fasta files for each UCLUST cluster""")
    parser.add_argument(
            "uclust",
            help="""The input directory containing UCLUST results"""
        )
    parser.add_argument(
            "fasta",
            help="""The directory containing the fasta reads clustered"""
        )
    parser.add_argument(
            "outdir",
            help="""The outut directory to hold our fastas"""
        )
    return parser.parse_args()


def map_fastas_to_cluster(infile, fastas, base):
    """ from matching results, determine which reads
    from a given fasta go into which cluster """
    # build a dictionary of which reads go with which
    # cluster
    clusters = {}
    reads = defaultdict(list)
    for result in uclust.Reader(infile):
        clusters[result.query_label] = result.cluster_number
    fasta_file = os.path.join(
            fastas,
            "{0}.fasta".format(base)
        )
    quality_file = os.path.join(
            fastas,
            "{0}.qual".format(base)
        )
    for read in fasta.FastaQualityReader(fasta_file, quality_file):
        # get cluster id
        cluster_id = clusters[read.identifier.strip('>')]
        reads[cluster_id].append(read)
    return reads


def write_results(reads, outdir):
    """ given our read map and an output directory
    write out results out"""
    for k, v in reads.iteritems():
        outf = "Cluster-{}.fasta".format(k)
        outq = "Cluster-{}.qual".format(k)
        outfp = os.path.join(outdir, outf)
        outfq = os.path.join(outdir, outq)
        fw = fasta.FastaWriter(outfp, outfq)
        sys.stdout.write("Writing Cluster {}".format(k))
        sys.stdout.flush()
        for read in v:
            sys.stdout.write(".")
            sys.stdout.flush()
            fw.write(read, qual=True)
        fw.close
        print ""


def main():
    args = get_args()
    for infile in glob.glob(os.path.join(args.uclust, "*.uc")):
        base = os.path.basename(infile).split(os.extsep)[0]
        print base
        print '=' * (len(base) + 2)
        reads = map_fastas_to_cluster(infile, args.fasta, base)
        outdir = os.path.join(args.outdir, base)
        # TODO: add a little user interaction here
        os.makedirs(outdir)
        write_results(reads, outdir)
        print "\n"


if __name__ == '__main__':
    main()
