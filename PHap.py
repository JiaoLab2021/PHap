#!/usr/bin/env python
'''
time: 2024-08-07
author: pxxiao
version: 1.0
description: 该流程用于组装多倍体基因组
'''


import argparse
import os
import subprocess
import sys

# Constants
# __version__ = "1.0.0"
# __update_time__ = "2024-08-08"

__version__ = "1.1.0"
__update_time__ = "2025-04-22"


HELP_DOC = """
  Usage: phap [command] <parameters>
  
  Command       Description
  --------      ------------------------------------------------------------------------
  mt2t          Generate a mosaic telomere-to-telomere genome using the p_ctg genome as
                the reference. This is used for creating an allelic contig table to solve
                the allelic conflict problem.
             
  cluster       Cluster contigs. First, extract contigs corresponding to each chromosome 
                based on mT2T alignment. Then, cluster the contigs for each haplotype of
                each chromosome using Hi-C signals. Third, cluster the contigs not on mT2T
                chromosomes to clustered group.
  
  phase_reads   Based on the clustering results, the corresponding reads of each haplotype
                are extracted, then assembled and anchored separately.

  Use phap [command] --help/-h to see detailed help for each individual command.
  """


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Telomere-to-telomere Autopolyploid assembly pipeline')
    subparsers = parser.add_subparsers(dest='function', metavar='', help='Available functions')

    ### mT2T  subcommand
    mT2T = subparsers.add_parser('mt2t', help='Use p_ctg genome to assemble a mosaic telomere-to-telomere '
                                              'genome as a reference genome for subsequent production of '
                                              'allelic contig tables')
    mT2T.add_argument('--p_ctg', required=True, type=str, help='Path to p_ctg genome')
    mT2T.add_argument('--min_distance', required=False, type=float, default=0.15, help='Minimum distance between two contigs [0.15]')
    mT2T.add_argument('--match_ratio', required=False, type=float, default=0.2, help='Minimum match ratio between contigs [0.2]')
    mT2T.add_argument('--min_contig_length', required=False, type=int, default=100000, help='Minimum contig length [100000]')
    mT2T.add_argument('--min_alignment_length', required=False, type=int, default=200, help='Minimum alignment length [200]')

    ### cluster
    cluster = subparsers.add_parser('cluster', help='Use the comparison information of p_utg and mT2T to create '
                                                    'an allelic contigs table, and combine the Hi-C signal to cluster '
                                                    'the contigs of the same haplotype.')
    cluster.add_argument('--p_utg', required=True, type=str, help='Path to p_utg genome')

    ### version subcommand
    version_parser = subparsers.add_parser('version', help='Print version information')
    version_parser.set_defaults(function='version')  # 设置默认的 function 值

    args = parser.parse_args()
    return args


def execute_script(script_name):
    script_realpath = os.path.dirname(os.path.realpath(__file__))
    sub_program = os.path.join(script_realpath, 'scripts', script_name)
    commands = [sys.executable, sub_program] + sys.argv[2:]

    try:
        subprocess.run(commands, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing {script_name}: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    # args = parse_arguments()

    # Map each function to its corresponding script
    params = {
        'mt2t': 'phap_mT2T.py',
        'cluster': 'phap_cluster.py',
        'phase_reads': 'phap_phase_reads.py',
        'fill_gap': 'phap_fill_gaps.py'
    }

    if len(sys.argv) > 1 and sys.argv[1] in params:
        execute_script(params[sys.argv[1]])

    elif len(sys.argv) > 1 and sys.argv[1] in ['version', '-v', '--version', '-version']:
        print(f"Version: {__version__}\nUpdate Time: {__update_time__}")
    elif len(sys.argv) < 2:
        print(HELP_DOC)

    else:
        print(HELP_DOC)


if __name__ == '__main__':
    main()