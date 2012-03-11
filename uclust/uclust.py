#!/usr/bin/env python
# encoding: utf-8
"""
File: uclust.py
Author: Brant Faircloth

Created by Brant Faircloth on 10 March 2012 17:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Parse a file of UCLUST results.

"""

import pdb


class UclustResult():
    def __init__(self):
        self.type = None
        self.cluster_number = None
        self.seq_length = None
        self.pct_id = None
        self.strand = None
        self.query_start = None
        self.seed_start = None
        self.alignment = None
        self.query_label = None
        self.target_label = None

    def from_dict(self, d):
        for k, v in d.iteritems():
            self.__dict__[k] = v


class Reader():
    def __init__(self, infile):
        self.uclust_file = open(infile)

    def close(self):
        """close files"""
        self.uclust_file.close()

    def _convert(self, value):
        try:
            value = int(value)
        except ValueError:
            value = float(value)
        finally:
            return value
        return value

    def next(self):
        """ parse a UCLUST file """
        line = self.uclust_file.readline()
        while line and line.startswith('#'):
            line = self.uclust_file.readline()
        if not line or line == '':
            raise StopIteration
        ls = line.strip().split('\t')
        assert len(ls) == 10, "Fields are too short"
        ls = [self._convert(i) for i in ls]
        fieldnames = (
                "type",
                "cluster_number",
                "seq_length",
                "pct_id",
                "strand",
                "query_start",
                "seed_start",
                "alignment",
                "query_label",
                "target_label"
            )
        fields = dict(zip(fieldnames, ls))
        uc_result = UclustResult()
        uc_result.from_dict(fields)
        return uc_result

    def __iter__(self):
        """iterator"""
        while True:
            yield self.next()
