import os
import sys
import glob
import math
import shutil
import tempfile
import argparse
import ConfigParser
from collections import defaultdict, Counter

from blat import Blat
from psl import PslReader

import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class CreateDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        d = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(d):
            answer = raw_input("[WARNING] Output directory exists, REMOVE [Y/n]? ")
            if answer == "Y":
                shutil.rmtree(d)
            else:
                print "[QUIT]"
                sys.exit()
        # create the new directory
        os.makedirs(d)
        # return the full path
        setattr(namespace, self.dest, d)


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

class Core:
    def __init__(self, name):
        self.name = name
        self.read = defaultdict(lambda: defaultdict(list))

    def add_read(self, psl):
        identifier = psl.q_name
        self.read[identifier][psl.match].append(psl)


def get_args():
    parser = argparse.ArgumentParser(description="""Match UCE probes to assembled contigs and store the data""")
    parser.add_argument("db", help = """The input database to match against""",
        action = FullPaths)
    parser.add_argument("query", help = """The input file containing the
        fasta read""", action = FullPaths, type = is_dir)
    parser.add_argument("--filter-length", dest = "length", help = """The shortest reads to
            accept""", type = int, default = 50)
    parser.add_argument("--raw", help = """Output the raw psl""",
            action = "store_true")
    parser.add_argument("--scale", help = """Output the raw psl""",
            action = "store_true")
    parser.add_argument("--conf", help = """Input desired order and output as
        csv""")
    parser.add_argument("--section", help = """Section of --conf to use as
        dict""")

    return parser.parse_args()


def main():
    args = get_args()
    if args.conf:
        config = ConfigParser.ConfigParser()
        config.read(args.conf)
        expected_order = config.items(args.section)
        expected_dict = {v:None for k,v in expected_order}
        print ',' + ','.join(["{}".format(e[1]) for e in expected_order])
    for file in glob.glob(os.path.join(args.query, '*.fasta')):
        # get core name
        core = Core(os.path.basename(file))
        # just store raw psl in tempfile
        fd, temp_psl = tempfile.mkstemp(suffix='.blat')
        os.close(fd) #kludge
        blat = Blat(args.db, file, temp_psl)
        blat.run()
        psl = PslReader(temp_psl)
        if not args.conf:
            # construct psl result library for given core and read
            sys.stdout.write("\nScanning PSL file {}".format(core.name))
        for count, row in enumerate(psl):
            core.add_read(row)
            if not args.conf:
                if count % 1000 == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
        if not args.conf:
            print "\n"
        os.remove(temp_psl)
        candidates = []
        for id in core.read.keys():
            # get max score
            mx = max(core.read[id].keys())
            if args.scale:
                mx = mx - math.floor(0.1 * mx)
            #pdb.set_trace()
            [candidates.extend(core.read[id][score])
                    for score in core.read[id].keys() \
                    if score >= mx
                ]
            names = [cand.t_name.split('|')[0] for cand in candidates \
                        if cand.q_size >= args.length
                        ]
        species_counts = Counter(names)
        if not args.conf:
            print "{0},{1}".format(core.name, species_counts)
        else:
            # merge dicts
            td = dict(expected_dict.items() + dict(species_counts).items())
            outlist = ["{0}".format(td[k])
                    if td[k] is not None else ''
                    for v,k in expected_order
                    ]
            print "{0},{1}".format(core.name, ','.join(outlist))
            #pdb.set_trace()

if __name__ == '__main__':
    main()
