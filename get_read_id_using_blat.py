import os
import sys
import glob
import shutil
import argparse
from operator import attrgetter, itemgetter
from tools.fs.args import FullPaths
from tools.align.psl import score, percent_id
from subprocess import Popen, PIPE

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('db', help = 'The input database to match against',
        action = FullPaths)
    parser.add_argument('query', help = 'The input file containing the' + \
        ' fasta read', action = FullPaths)
    parser.add_argument('--top-five', dest = 'top_five', help = 'Output the top five matches',
            action = 'store_true')
    parser.add_argument('--raw', help = 'Output the raw psl',
            action = 'store_true')
    return parser.parse_args()

def main():
    args = get_args()
    blat_args = "blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0" + \
        " -noHead {0} {1} stdout"
    blat_args = blat_args.format(args.db, args.query)
    #pdb.set_trace()
    blat = Popen(blat_args, shell=True, stdout=PIPE)
    stdout, stderr = blat.communicate()
    stdout = stdout.strip()
    stdout = [[int(item) if item.isdigit() else item for item in line.split('\t')] for line in stdout.split('\n')]
    # get best result in list
    #pdb.set_trace()
    s = sorted(stdout, key = itemgetter(0), reverse = True)
    best = sorted(s, key = itemgetter(1))
    fieldnames = ['match',
                'mismatch',
                'repmatch',
                'N',
                'q_gap_count',
                'q_gap_bases',
                't_gap_count',
                't_gap_bases',
                'strand',
                'q_name',
                'q_size',
                'q_start',
                'q_end',
                't_name',
                't_size',
                't_start',
                't_end',
                'block_count',
                'block_sizes',
                'qstarts',
                't_starts']
    if not args.top_five:
        psl_dict = dict(zip(fieldnames, best[0]))
        print psl_dict['t_name'].split('|')[0], score(psl_dict), percent_id(psl_dict)
    else:
        if not args.raw:
            for match in best[:5]:
                psl_dict = dict(zip(fieldnames, match))
                name =  psl_dict['t_name'].split('|')[0]
                print "{} {} {}\t matches = {}\t mismatches = {} gaps = {}".format(name, score(psl_dict), 
                    percent_id(psl_dict), psl_dict['match'],
                    psl_dict['mismatch'], psl_dict['t_gap_count'])
        else:
            print "Remember that --raw drop block columns"
            print '\t'.join(fieldnames[:-4])
            for match in best[:5]:
                print '\t'.join([str(m) for m in match[:-4]])
    #pdb.set_trace()
    

if __name__ == '__main__':
    main()


