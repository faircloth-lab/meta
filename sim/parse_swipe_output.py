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
import subprocess
import multiprocessing
from Bio import SeqIO
from collections import defaultdict

import pdb

def run_swipe(work):
    rv = None
    sequence, regex = work
    slcs = range(10, len(sequence.seq), 10)
    slcs.append(len(sequence.seq))
    for slc in reversed(slcs):
        # write sequence out to tempfile
        handle, file = tempfile.mkstemp(suffix='.fasta')
        os.write(handle, sequence[:slc].format('fasta'))
        os.close(handle)
        # create output tempfile
        handle, output = tempfile.mkstemp(suffix='.swipe')
        os.close(handle)
        # run alignment of sequence against database
        # ~/src/swipe-2.0.4/swipe -d bci.fasta -i test.fasta -p 0 -M BLOSUM80 -G 100 > out.results
        cmd = [
            '/home/bcf/src/swipe-2.0.4/swipe',
            '-d', 'bci.fasta',
            '-i', file,
            '-p', '0',
            '-M', 'BLOSUM80',
            '-G', '100',
            '-o', output]
        process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
        stdout, stderr = process.communicate()
        # parse the alignment results
        container = defaultdict(list)
        for line in open(output, 'rU'):
            if line.startswith('>'):
                break
            else:
                result = regex.search(line)
                if result:
                    grp = result.groups()
                    container[int(grp[2])].append(grp)
        # remove the temp files
        for f in [file, output]:
            os.remove(f)
        if container:
            # get max bit score results - these are "best" matches with highest bit score
            best_hits = container[max(container.keys())]
            best_taxa = set(['_'.join(taxon[0].split('_')[:2]) for taxon in best_hits])
            best_taxa_with_ids = set([taxon[0] for taxon in best_hits])
            if slc == len(sequence.seq):
                first_taxa = copy.deepcopy(best_taxa)
        # if we reach the minimum, or the results change, or we don't get results
        if slc == 10 or (best_taxa != first_taxa) or not container:
            rv = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(sequence.id, len(sequence), previous_slc, previous_max, previous_runner_up, ', '.join(previous_taxa))
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

def main():
    """main loop"""
    # compile regular expression to pick up summary info line once
    regex = re.compile("^gnl\|(?:.*)\|[0-9]+\s([A-Za-z0-9_]+)\s+(.*)\s+([0-9]+)\s+(.*)$")
    work = [(seq, regex) for seq in SeqIO.parse(open('bci.fasta', 'rU'), 'fasta')]
    pool = multiprocessing.Pool(12)
    sys.stdout.write("Running")
    sys.stdout.flush()
    results = pool.map(run_swipe, work)
    outp = open('matches.tdt', 'w')
    for result in results:
        outp.write(result)
    outp.close()
    

if __name__ == '__main__':
    main()