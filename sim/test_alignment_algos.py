#!/usr/bin/env python
# encoding: utf-8
"""
File: unnamed_file.py
Author: Brant Faircloth

Created by Brant Faircloth on 24 October 2012 10:10 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import re
import sys
import copy
import tempfile
import argparse
import subprocess
import multiprocessing
from Bio import SeqIO
from collections import defaultdict

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Test difference alignment programs""")
    parser.add_argument(
            "fasta",
            help="""The fasta file to use"""
        )
    parser.add_argument(
            "db",
            help="""The database to use"""
        )
    parser.add_argument(
            "algo",
            help="""The algorithm to pass"""
        )
    parser.add_argument(
            "output",
            help="""The output file name"""
        )
    parser.add_argument(
            "--taxa",
            type=int,
            help="""The number of taxa to screen""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            help="""The number of compute cores to use""",
        )
    parser.add_argument(
            "--five-prime",
            action="store_true",
            default=False,
            help="""Trim sequence from 5' => 3' instead of reverse""",
        )
    return parser.parse_args()


def get_swipe_cmd(file, db, output):
    # ~/src/swipe-2.0.4/swipe -d bci.fasta -i test.fasta -p 0 -M BLOSUM80 -G 100 > out.results
    cmd = [
            '/home/bcf/src/swipe-2.0.4/swipe',
            '-d', db,
            '-i', file,
            '-p', '0',
            '-M', 'BLOSUM80',
            '-G', '100',
            '-o', output
        ]
    return cmd


def get_blast_cmd(file, db, output):
    # blastn -db bci.fasta -query abarema-macradenia-0.fasta -task megablast -evalue 1e-50 > abarema-out.txt
    cmd = [
            'blastn',
            '-db', db,
            '-query', file,
            '-task', 'megablast',
            '-evalue', '1e-50',
            '-gapopen', '3',
            '-gapextend', '1',
            '-out', output
        ]
    return cmd


def get_blat_cmd(file, db, output):
    # blastn -db bci.fasta -query abarema-macradenia-0.fasta -task megablast -evalue 1e-50 > abarema-out.txt
    cmd = [
            'blat',
            db,
            file,
            '-t=dna',
            '-q=dna',
            '-fastMap',
            '-out=blast',
            output
        ]
    return cmd


def get_cmd_line(algo, file, db, output):
    # run alignment of sequence against database
    if algo == 'swipe':
        cmd = get_swipe_cmd(file, db, output)
    elif algo == 'blast':
        cmd = get_blast_cmd(file, db, output)
    elif algo == 'blat':
        cmd = get_blat_cmd(file, db, output)
    return cmd


def parse_alignment_results(regex, output):
    container = defaultdict(list)
    for line in open(output, 'rU'):
        if line.startswith('>'):
            break
        else:
            result = regex.search(line)
            if result:
                grp = result.groups()
                container[int(grp[1])].append(grp)
    return container


def run_algo(work):
    rv = None
    sequence, regex, algo, five_prime, db = work
    if five_prime:
        slcs = range(0, len(sequence.seq), 10)
    else:
        slcs = range(10, len(sequence.seq), 10)
        slcs.append(len(sequence.seq))
        slcs = reversed(slcs)
    for slc in slcs:
        # write sequence out to tempfile
        handle, file = tempfile.mkstemp(suffix='.fasta')
        if five_prime:
            os.write(handle, sequence[slc:].format('fasta'))
        else:
            os.write(handle, sequence[:slc].format('fasta'))
        os.close(handle)
        # create output tempfile
        handle, output = tempfile.mkstemp(suffix='.swipe')
        os.close(handle)
        # get the appropriate cmd
        cmd = get_cmd_line(algo, file, db, output)
        process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
        stdout, stderr = process.communicate()
        # parse the alignment results
        container = parse_alignment_results(regex, output)
        # remove the temp files
        for f in [file, output]:
            os.remove(f)
        if container:
            # get max bit score results - these are "best" matches with highest bit score
            best_hits = container[max(container.keys())]
            best_taxa = set(['_'.join(taxon[0].split('_')[:2]) for taxon in best_hits])
            best_taxa_with_ids = set([taxon[0] for taxon in best_hits])
            if five_prime and slc == 0:
                first_taxa = copy.deepcopy(best_taxa)
            elif not five_prime and slc == len(sequence.seq):
                first_taxa = copy.deepcopy(best_taxa)
        # if we reach the maximum, or the results change, or we don't get results
        if five_prime and (slc == max(slcs) or (best_taxa != first_taxa) or not container):
            rv = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
                    sequence.id,
                    len(sequence),
                    len(sequence) - previous_slc,
                    previous_max, previous_runner_up,
                    ', '.join(previous_taxa),
                    ', '.join(best_taxa_with_ids)
                )
            break
        # if we reach the minimum, or the results change, or we don't get results
        elif not five_prime and (slc == 10 or (best_taxa != first_taxa) or not container):
            rv = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
                    sequence.id,
                    len(sequence),
                    previous_slc,
                    previous_max,
                    previous_runner_up,
                    ', '.join(previous_taxa),
                    ', '.join(best_taxa_with_ids)
                )
            break
        previous_taxa = copy.deepcopy(best_taxa_with_ids)
        previous_slc = slc
        if len(container) > 1:
            previous_max, previous_runner_up = sorted(container.keys())[::-1][:2]
        else:
            previous_max, previous_runner_up = max(container.keys()), None
    # give some indication of progress
    sys.stdout.write('.')
    sys.stdout.flush()
    if rv:
        return rv
    else:
        # if all else fails
        return "{0}\t{1}\tNone\tNone\tNone\n".format(sequence.id, len(sequence))


def get_algo_regex(algo):
    if algo == 'swipe':
        # compile regular expression to pick up summary info line once
        regex = re.compile("^gnl\|(?:.*)\|[0-9]+\s([A-Za-z0-9_]+)\s+(?:.*)\s+([0-9]+)\s+(.*)$")
    elif algo == 'blast':
        regex = re.compile("^\s+([A-Za-z0-9_]+)\s+([0-9]+)\s+(.*)$")
    elif algo == 'blat':
        regex = re.compile("^([A-Za-z0-9_]+)\s+([0-9]+)\s+(.*)$")
    return regex


def main():
    """main loop"""
    args = get_args()
    regex = get_algo_regex(args.algo)
    work = [(seq, regex, args.algo, args.five_prime, args.db) for seq in SeqIO.parse(open(args.fasta, 'rU'), 'fasta')]
    if args.taxa:
        work = work[:args.taxa]
    sys.stdout.write("Running")
    sys.stdout.flush()
    if args.cores > 1:
        pool = multiprocessing.Pool(args.cores)
        results = pool.map(run_algo, work)
    else:
        results = map(run_algo, work)
    outp = open(args.output, 'w')
    for result in results:
        outp.write(result)
    outp.close()


if __name__ == '__main__':
    main()
