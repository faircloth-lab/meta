"""
File: get_cluster_counts.py
Author: Brant Faircloth

Created by Brant Faircloth on 14 December 2011 21:12 PST (-0800)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import sys
import glob
import shutil
import argparse
from seqtools.fs.args import FullPaths, is_dir
from subprocess import Popen, PIPE

import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Match UCE probes to assembled
            contigs and store the data""")
    parser.add_argument(
            'input',
            action=FullPaths,
            type=is_dir,
            help="""The input directory containing fastas
            to cluster"""
        )
    parser.add_argument(
            '--output',
            default=os.getcwd(),
            action=FullPaths,
            type=is_dir,
            help="""The output directory to hold results""",
        )
    parser.add_argument(
            '--id',
            help="""The % identity to use""",
            default=0.9
        )
    return parser.parse_args()


def main():
    args = get_args()
    #pdb.set_trace()
    # remove prior results
    regex = re.compile('([0-9]+)\sclusters')
    for file in glob.glob(os.path.join(args.input, "*.fasta")):
        name = os.path.splitext(file)[0]
        outname = os.path.split(name)[1]
        sort = "uclust --sort {0}.fasta --output {1}.sorted.fasta 2> /dev/null"
        cluster = "uclust --input {0}.sorted.fasta --uc {0}.sorted.results.uc --id {1}"
        uc_sort = Popen(sort.format(name, os.path.join(args.output, outname)), shell=True, stdout=PIPE)
        sort_stdout, sort_stderr = uc_sort.communicate()
        uc_cluster = Popen(cluster.format(os.path.join(args.output, outname), args.id), shell=True,
                stdout=PIPE, stderr=PIPE)
        cluster_stdout, cluster_stderr = uc_cluster.communicate()
        clusters = regex.search(cluster_stderr)
        cluster_count = int(clusters.groups()[0])
        print os.path.split(file)[1], cluster_count
        #pdb.set_trace()


if __name__ == '__main__':
    main()
