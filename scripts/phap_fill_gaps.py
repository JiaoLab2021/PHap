#!/usr/bin/env python
'''
time: 2024-09-27
author: pxxiao
version: 1.0
description: 基于 ONT reads 比对，填充 Gap
----
'''

import argparse
import pickle
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import os

import glob

from util import *




def parse_arguments():
    parser = argparse.ArgumentParser(description="Haplotype clustering of autopolyploid genome")

    parser.add_argument('--p_utg', required=True, type=str, help='Path to p_utg file')
    parser.add_argument('--mT2T', required=True, type=str, help='Path to mT2T file or reference genome')
    parser.add_argument('--contig_type', required=True, type=str, help='Path to contig type from dosage analysis')
    parser.add_argument('--threads', type=int, default=10, help='The number of threads [10]')

    find_longest = parser.add_argument_group('>>> Parameters for find longest subsequences')
    find_longest.add_argument('--min_align_length', type=int, default=200, help="Minimum alignment length [200]")
    find_longest.add_argument("--min_unitig_length", type=int, default=1000, help="Minimum unitig length [1000]")
    find_longest.add_argument('--min_alignment_distance', type=int, default=500000, help='The minimum distance between two alignments [500000]')
    find_longest.add_argument('--min_match_ratio', type=float, default=0.05, help='Minimum match ratio [0.05]')
    find_longest.add_argument('--min_lis_size', type=int, default=5, help='Minimum number of alignments in a LIS to be considered for mergeing [5]')
    find_longest.add_argument('--min_lis_length', type=int, default=1000000, help='Minimum length of alignments in a LIS to be considered for mergeing [1000000]')
    find_longest.add_argument('--max_lis_distance', type=int, default=3000000, help='Maximum distance between LIS to be considered for merge [3000000]')

    allelic_table = parser.add_argument_group('>>> Allelic table generated')
    allelic_table.add_argument("--bin_size", type=int, default=100000, help="Bin size [100000]")
    allelic_table.add_argument('--chr_num', type=int, default=12, help='The number of chromosomes [12]')
    allelic_table.add_argument('--top_n', type=int, default=4, help='The number of haplotypes, for example, tetraploid is 4 [4]')
    allelic_table.add_argument('--search_range', type=int, default=20, help='Search range [20]')

    cluster = parser.add_argument_group('>>> Cluster based on Hi-C links')
    cluster.add_argument('--full_links', type=str, required=True, help='Path to full links pickle')
    cluster.add_argument('--flank', type=int, help='Flank')

    recluster = parser.add_argument_group('>>> Recluster based on Hi-C links')
    recluster.add_argument('--clm', type=str, required=True, help='Path to clm, from parsed hi-c links bam')

    args = parser.parse_args()
    return args



def main():
    args = parse_arguments()

    script_realpath = os.path.dirname(os.path.realpath(__file__))
    utils_realpath = os.path.join(script_realpath, '..', 'utils')

    cwd = os.getcwd()       # current working dir