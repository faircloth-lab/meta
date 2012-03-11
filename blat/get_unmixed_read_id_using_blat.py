import os
import sys
import glob
import shutil
import argparse
import tempfile

from collections import defaultdict
from seqtools.align import psl
from operator import attrgetter, itemgetter
from seqtools.fs.args import FullPaths
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

def run_blat(db, query, top_five = False, raw = False):
    fd, temp_psl = tempfile.mkstemp(suffix='.blat')
    os.close(fd) #kludge
    blat_args = "blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0" + \
        " -noHead {0} {1} {2}"
    blat_args = blat_args.format(db, query, temp_psl)
    blat = Popen(blat_args, shell=True, stdout=PIPE)
    stdout, stderr = blat.communicate()
    p = psl.PslReader(temp_psl)
    # delete temp file
    os.remove(temp_psl)
    # get best result in list
    results = defaultdict(list)
    for match in p:
        # severely punish mismatches
        score = psl.score(match) - 15 * match.mismatch
        results[score].append(match)
    best = sorted(results.keys(), reverse = True)
    try:
        if not top_five:
            top = results[best[0]]
            if len(top) == 1:
                winner = top[0]
                print os.path.basename(query), winner.t_name.split('|')[0], winner[0], psl.percent_id(winner)
            else:
                for winner in top:
                    print "\t", os.path.basename(query), winner.t_name.split('|')[0], winner[0], psl.percent_id(winner)


        else:
            if not raw:
                for score in best[:5]:
                    winner = results[score][0]
                    name = winner.t_name.split('|')[0]
                    print "{} {} {} {}\t matches = {}\t mismatches = {} gaps = {}".format(
                            os.path.basename(query),
                            name,
                            score,
                            psl.percent_id(winner),
                            winner.match,
                            winner.mismatch,
                            winner.t_gap_count)
            else:
                print "Remember that --raw drop block columns"
                #pdb.set_trace()
                for k,score in enumerate(best[:5]):
                    if k == 0:
                        print '\t'.join(results[score][0]._asdict().keys())
                    for item in results[score]:
                        print '\t'.join([str(i) for i in item._asdict().values()])
    except IndexError:
        print "{} No contig".format(os.path.basename(query))
    


def main():
    args = get_args()
    if not os.path.isdir(args.query):
        run_blat(args.db, args.query, args.top_five, args.raw)
    else:
        for file in glob.glob(os.path.join(args.query, "*.fasta")):
            run_blat(args.db, file, args.top_five, args.raw)
        #pdb.set_trace()
    

if __name__ == '__main__':
    main()


