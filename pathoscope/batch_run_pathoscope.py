#!/usr/bin/env python
# encoding: utf-8
"""
File: batch_run_pathoscope.py
Author: Brant Faircloth

Created by Brant Faircloth on 18 October 2013 15:10 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import sys
import glob
import shutil
import logging
import argparse
import subprocess
import ConfigParser
import PathoID
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "--fastas",
            required=True,
            type=is_dir,
            action=FullPaths,
            help="""The directory of input fastas (containing sample-specific dirs/*(.fasta|qual) files."""
        )
    parser.add_argument(
            "--output",
            required=True,
            action=CreateDir,
            help="""The path to the output directory"""
        )
    parser.add_argument(
            "--reference",
            required=True,
            type=is_file,
            action=FullPaths,
            help="""The path to the reference database fasta file."""
        )
    parser.add_argument(
            "--primer-conf",
            type=is_file,
            action=FullPaths,
            help="""The path to the primers to trim."""
        )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    return parser.parse_args()


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


def setup_logging(args):
    import __main__ as main
    my_name = os.path.basename(os.path.splitext(main.__file__)[0])
    log = logging.getLogger(my_name)
    console = logging.StreamHandler(sys.stdout)
    if args.log_path is not None:
        logfile = logging.FileHandler(os.path.join(args.log_path, "{}.log".format(my_name)))
    else:
        logfile = logging.FileHandler("{}.log".format(my_name))
    if args.verbosity == "INFO":
        log.setLevel(logging.INFO)
        console.setLevel(logging.INFO)
        logfile.setLevel(logging.INFO)
    if args.verbosity == "WARN":
        log.setLevel(logging.WARN)
        console.setLevel(logging.WARN)
        logfile.setLevel(logging.WARN)
    if args.verbosity == "CRITICAL":
        log.setLevel(logging.CRITICAL)
        console.setLevel(logging.CRITICAL)
        logfile.setLevel(logging.CRITICAL)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logfile.setFormatter(formatter)
    log.addHandler(console)
    log.addHandler(logfile)
    text = " Starting {} ".format(my_name)
    log.info(text.center(65, "="))
    for arg, value in sorted(vars(args).items()):
        log.info("Argument --%s: %r", arg, value)
    return log, my_name


def convert_fasta_qual_to_fastq(log, directory, sample, output):
    # current fasta and qual files are
    fasta = os.path.join(directory, "{}.fasta".format(sample))
    quality = os.path.join(directory, "{}.qual".format(sample))
    # make a similar directory in args.output
    outdir = os.path.join(output, sample)
    os.mkdir(outdir)
    # convert current sequences to fastq, writing them out
    fastq = os.path.join(outdir, "{}.fastq".format(sample))
    with open(fastq, "w") as outf:
        records = PairedFastaQualIterator(open(fasta), open(quality))
        count = SeqIO.write(records, outf, "fastq")
    log.info("Converted {} fasta+qual records to fastq".format(count))
    return outdir, fastq


def create_trim_file(log, output, primer_conf):
    log.info("Creating primers.fasta file for trimming")
    config = ConfigParser.ConfigParser()
    # case sensitive
    config.optionxform=str
    config.read(primer_conf)
    primers = dict(config.items('Primers'))
    outname = os.path.join(output, "primers.fasta")
    with open(outname, "w") as outf:
        for name, primer in primers.iteritems():
            outf.write(">{}\n{}\n".format(
                name,
                primer
            ))
    return outname


def trim_primers(log, primers, sample, fastq):
    log.info("Trimming primers")
    pth, name = os.path.split(fastq)
    name = os.path.splitext(name)[0]
    new_fastq = "{}.trimmed.fastq".format(name)
    new_fastq_pth = os.path.join(pth, new_fastq)
    cmd = [
        "fastq-mcf",
        "-q",
        "5",
        primers,
        fastq,
        "-o",
        new_fastq_pth
    ]
    with open(os.path.join(pth, "fastq-mcf.stdout"), 'w') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
        stdout, stderr = proc.communicate()
    if stderr:
        raise IOError("[fastq-mcf] {}".format(stderr))
    return new_fastq_pth


def create_bowtie_dict(log, output, reference):
    good = False
    dict_name = os.path.splitext(reference)[0]
    pth = os.path.split(dict_name)[0]
    for ext in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        if os.path.exists("{}{}".format(dict_name, ext)):
            good = True
        else:
            good = False
            break
    if not good:
        log.info("Creating bowtie dictionary")
        # make a new bowtie ref db
        cmd = [
            "bowtie2-build",
            reference,
            dict_name
        ]
        with open(os.path.join(pth, "bowtie2-build.stdout"), 'w') as out:
            proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
            stdout, stderr = proc.communicate()
        if stderr:
            raise IOError("[bowtie2-build] {}".format(stderr))
    else:
        log.info("Bowtie2 dictionary exists.  Skipping dictionary creation.")
    return dict_name


def run_bowtie(log, ref_dict_name, fastq):
    log.info("Running bowtie2")
    pth, name = os.path.split(fastq)
    name = name.split(".")[0]
    out_sam = os.path.join(pth, "{}.sam".format(name))
    cmd = [
        "bowtie2",
        "-k",
        "100",
        "-x",
        ref_dict_name,
        "-U",
        fastq,
        "-S",
        out_sam
    ]
    with open(os.path.join(pth, "bowtie2.stdout"), 'w') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
        stdout, stderr = proc.communicate()
    if stderr:
        raise IOError("[bowtie2] {}".format(stderr))
    return out_sam


def run_pathoscope(log, sam, sample, outdir):
    log.info("Running pathoscope")
    out_matrix = False
    verbose = False
    score_cutoff = 0.01
    exp_tag = sample
    ali_format = "sam"
    emEpsilon = 1e-7
    maxIter = 50
    PathoID.pathoscope_reassign(
        out_matrix,
        verbose,
        score_cutoff,
        exp_tag,
        ali_format,
        sam,
        outdir,
        emEpsilon,
        maxIter,
        True
    )


def main():
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    sorted_dirs = sorted(glob.glob(os.path.join(args.fastas, "*")))
    # check for dict
    ref_dict_name = create_bowtie_dict(log, args.output, args.reference)
    # create_trim_file if we're trimming
    if args.primer_conf:
        primers = create_trim_file(log, args.output, args.primer_conf)
    else:
        primers = False
    for cnt, directory in enumerate(sorted_dirs):
        sample = os.path.basename(directory)
        text = " Running sample {} ".format(sample)
        log.info(text.center(65, "-"))
        # make fastq file
        sample_outdir, fastq = convert_fasta_qual_to_fastq(log, directory, sample, args.output)
        if args.primer_conf:
            # trim primers from fastq
            fastq = trim_primers(log, primers, sample, fastq)
        # run bowtie
        sam = run_bowtie(log, ref_dict_name, fastq)
        # run pathoscope
        run_pathoscope(log, sam, sample, sample_outdir)
        if cnt == 10:
            break
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))



if __name__ == '__main__':
    main()
