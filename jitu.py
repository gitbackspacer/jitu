#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module is helper utilities for an alternate use for Fasta parsing on Python2.7 and test on Python >=3.6
Using duration for timing section of code 

"""

import sys
import itertools
import re
from collections import defaultdict
import pprint
from itertools import groupby
from contextlib import contextmanager
from time import time
from sys import stdout


# from dottree import dotree as dot

@contextmanager
def duration_in_seconds(outfile=stdout):
    """ use with duration to time the code sections  """
    start = time()
    yield
    end = time()
    outfile.write(str(end - start) + "\n")


def getEasy(f):
    """ Pure python fasta reader  """
    tag = None
    for header, group in itertools.groupby(f, key=lambda x: x.startswith(">")):
        if header:
            line = next(group)
            tag = line[1:].strip().split()
        else:
            sequence = ''.join(line.strip() for line in group)
            if not tag:
                sys.exit("Error parsing the Fasta file make sure the Fasta starts with '>' and first line is not blank!")
            yield tag, sequence.strip()


def getSeqD(fastaFiler="ReferenceGenome.fasta"):
    """
    Load multiple tuple sequence we are losing the yield advantage or use getEasy
    """
    seqD = {}
    tube = []
    with open(fastaFiler, "r") as inpFast:
        for h, seq in getEasy(inpFast):
            seqD[tuple(h)] = seq
            tube.append(tuple(h))
    return tube, seqD


def getTubeD(fastaFiler="ReferenceGenome.fasta"):
    """
    Load all fasta record at once and get concat headers  into a single string
    """
    seqD = {}
    tube = []
    with open(fastaFiler, "r") as inpFast:
        for h, seq in getEasy(inpFast):
            tag = "".join(h)
            seqD[tag] = seq
            tube.append(tag)
    return tube, seqD
