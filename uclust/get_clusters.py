"""
File: get_cluster_counts.py
Author: Brant Faircloth

Created by Brant Faircloth on 14 December 2011 21:12 PST (-0800)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: Given a fasta file of sequences, cluster those reads by
user-specific parameter and write results to outdir.

"""

import os
import re
import glob
import argparse
import subprocess
from seqtools.fs.args import FullPaths, is_dir


import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Match UCE probes to assembled
            contigs and store the data""")
    parser.add_argument(
            'indir',
            action=FullPaths,
            type=is_dir,
            help="""The input directory containing fastas
            to cluster"""
        )
    parser.add_argument(
            'outdir',
            action=FullPaths,
            type=is_dir,
            help="""The output directory to hold results""",
        )
    parser.add_argument(
            '--id',
            default=0.9,
            help="""The % identity to use"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    regex = re.compile('([0-9]+)\sclusters')
    for fasta in glob.glob(os.path.join(args.indir, "*.fasta")):
        name = os.path.splitext(fasta)[0]
        outname = os.path.join(args.outdir, os.path.split(name)[1])
        sort_cmd = [
                "uclust",
                "--sort",
                fasta,
                "--output",
                "{0}.sorted.fasta".format(outname)
            ]
        sort = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sort_stdout, sort_stderr = sort.communicate()
        cluster_cmd = [
                "uclust",
                "--input",
                "{0}.sorted.fasta".format(outname),
                "--uc",
                "{0}.sorted.results.uc".format(outname),
                "--id",
                args.id
            ]
        cluster = subprocess.Popen(cluster_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cluster_stdout, cluster_stderr = cluster.communicate()
        clusters = regex.search(cluster_stderr)
        cluster_count = int(clusters.groups()[0])
        print os.path.split(fasta)[1], cluster_count


if __name__ == '__main__':
    main()
