#!/usr/bin/env python
'''
time: 2024-07-29
author: pxxiao
version: 1.0
description:
根据 HiFi reads 比对的 bam 文件；contig type；unitigs groups：
    根据 contig 的类型，将 collapsed unitigs HiFi reads 随机分为 2 份、3 份、四份；
    最后得到每个 Group 包含的 HiFi reads；
    然后，根据 HiFi reads，对每个 Group 分别进行 hifiasm 组装。
----
time: 2024-07-30
author: pxxiao
version: 2.0
description:
增加 ONT reads，Hi-C reads 的 group-specific reads 划分；
然后，使用 HiFi+ONT 策略进行组装；
haphic 利用 Hi-C reads 进行挂载（还得看一下效果如何）
'''


import argparse
import os
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description="Haplotype assembly and scaffolding of autopolyploid genome")
    parser.add_argument('--bam_hifi', type=str, required=True, help='Path to HiFi bam file')
    parser.add_argument('--bam_hic', type=str, required=True, help='Path to Hi-C bam file')
    parser.add_argument('--bam_ont', type=str, required=True, help='Path to ONT bam file')
    parser.add_argument('--contig_type', required=True, type=str, help='Path to contig type from dosage analysis')
    parser.add_argument('--group', type=str, help='Path to merge group file from cluster, default is /02.cluster/05.rescue/merge.group.reassignment.cluster.txt')
    parser.add_argument('--hifi', type=str, required=True, help='Path to HiFi reads file')
    parser.add_argument('--ont', type=str, required=True, help='Path to ONT reads file')
    parser.add_argument('--hic1', type=str, required=True, help='Path to Hi-C forward reads file')
    parser.add_argument('--hic2', type=str, required=True, help='Path to Hi-C reverse reads file')
    parser.add_argument('--ont_length', type=int, default=1, help='Length of ONT reads file for hifiasm asm [1]')
    parser.add_argument('--ont_quality', type=int, default=0, help='Quality of ONT [0]')
    parser.add_argument('--threads', type=int, default=10, help='The number of threads [10]')
    parser.add_argument('--process', type=int, default=4, help='The number of processes [4]')
    parser.add_argument('--seed', type=int, default=100, help='Random seed for reproducibility [100]')

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    script_realpath = os.path.dirname(os.path.realpath(__file__))
    utils_realpath = os.path.join(script_realpath, '..', 'utils')

    cwd = os.getcwd()       # current working dir

    ### Step 1: p_utg vs mT2T
    step1_dir = os.path.join(cwd, '03.phase_reads')
    os.makedirs(step1_dir, exist_ok=True)
    os.chdir(step1_dir)
    cmd = (f'python {utils_realpath}/phase_reads_assemble_anchor.py '
           f'--bam_hifi {args.bam_hifi} '
           f'--bam_hic {args.bam_hic} '
           f'--bam_ont {args.bam_ont} '
           f'--contig_type {args.contig_type} '
           f'--group {cwd}/02.cluster/05.rescue/merge.group.reassignment.cluster.txt '
           f'--hifi {args.hifi} '
           f'--ont {args.ont} '
           f'--ont_length {args.ont_length} '
           f'--ont_quality {args.ont_quality} '
           f'--hic1 {args.hic1} '
           f'--hic2 {args.hic2} '
           f'--threads {args.threads} '
           f'--process {args.process} '
           f'--seed {args.seed} '
           f'> log_phase_reads_assemble_anchor_out 2> log_phase_reads_assemble_anchor_err')
    subprocess.run(cmd, shell=True, check=True)


if __name__ == '__main__':
    main()