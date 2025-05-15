#!/usr/bin/env python
'''
time: 2024-06-24
author: pxxiao
version: 1.0
----
description: 解析 Hi-C bam 文件，获得 unitig 之间的 Hi-C links

'''
import argparse
import gc
import pickle

import pysam
from math import ceil
from collections import defaultdict
from array import array


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, type=str, help='Hi-C bam file')
    parser.add_argument('--fasta', required=True, type=str, help='The fasta file')
    parser.add_argument('--flank', required=False, type=int, default=500000, help='The flank of the unitig [500000]')
    parser.add_argument('--threads', required=False, type=int, default=10, help='The number of threads [10]')
    args = parser.parse_args()
    return args


## 解析 FASTA: ctg: [seq, seq length, RE site's count]
def parse_fasta(fasta, RE='GATC'):
    """Parse input FASTA file and save sequences, lengths, and RE site counts of contigs into a dict."""
    ## 对序列计算 RE 个数
    def count_RE_sites(sequence, RE):
        ''' Count the number of restriction enzyme (RE) sites in a sequence. '''
        return sequence.count(RE)

    fa_dict = {}
    with open(fasta) as f:
        for line in f:
            if not line.strip(): continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                fa_dict[ctg] = []
            else:
                fa_dict[ctg].append(line.strip())
    for ctg, seq_list in fa_dict.items():
        # Joining list is faster than concatenating strings
        seq = ''.join(seq_list)
        # Add pseudo-count of 1 to prevent division by zero (as what ALLHiC does)
        RE_sites = count_RE_sites(seq, RE) + 1
        fa_dict[ctg] = [seq, len(seq), RE_sites]
    return fa_dict


def bam_generator(bam, threads, format_options):
    # just a wrapper of pysam.AlignmentFile used to keep the style of Hi-C link parsing consistent
    with pysam.AlignmentFile(bam, mode='rb', threads=threads, format_options=format_options) as f:
        for aln in f:
            yield aln.reference_name, aln.next_reference_name, aln.reference_start, aln.next_reference_start


def is_flank(coord, length, flank):

    """determine whether the coord is inside the flanking regions of a given fragment"""
    if flank and (coord <= flank or coord > length - flank):
        return True
    elif not flank:
        return True
    else:
        return False


def update_clm_dict(clm_dict, ctg_name_pair, len_i, len_j, coord_i_0, coord_j_0):

    clm_dict[ctg_name_pair].extend((
        len_i - coord_i_0 + coord_j_0,
        len_i - coord_i_0 + len_j - coord_j_0,
        coord_i_0 + coord_j_0,
        coord_i_0 + len_j - coord_j_0))


def output_pickle(dict_, to):

    with open(to, 'wb') as fpkl:
        pickle.dump(dict_, fpkl)


def load_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data


def output_clm(clm_dict):

    ori_tuple = (('+', '+'), ('+', '-'), ('-', '+'), ('-', '-'))

    with open('paired_links.clm', 'w') as fout:
        for ctg_name_pair, list_ in clm_dict.items():
            # minLinks == 3, links = len(list_)/2
            if len(list_) < 8:
                continue
            for n in range(4):
                new_list = ['{0} {0}'.format(v) for v in sorted(list_[n::4])]
                fout.write('{}{} {}{}\t{}\t{}\n'.format(
                    ctg_name_pair[0], ori_tuple[n][0],
                    ctg_name_pair[1], ori_tuple[n][1],
                    len(new_list)*2, ' '.join(new_list)))


def divide_into_bins(fa_dict, bin_size):
    """Divide contigs into bins of specified size."""
    bins = {}
    for ctg, data in fa_dict.items():
        length = data[1]
        num_bins = ceil(length / bin_size)
        for bin_index in range(num_bins):
            bin_key = f"{ctg}_bin{bin_index+1}"
            bins[bin_key] = {
                'start': bin_index * bin_size,
                'end': min((bin_index + 1) * bin_size, length),
                'length': min(bin_size, length - bin_index * bin_size)
            }
    return bins


def calculate_bin_hic_density(full_link_dict, bins):
    """Calculate Hi-C link density for each bin."""
    bin_link_counts = defaultdict(int)

    # Accumulate links per bin
    for (ctg_i, ctg_j), count in full_link_dict.items():
        for bin_key, bin_info in bins.items():
            ctg, bin_num = bin_key.rsplit('_bin', 1)
            if ctg == ctg_i or ctg == ctg_j:
                bin_link_counts[bin_key] += count

    # Calculate density based on bin size and count
    bin_density = {}
    for bin_key, count in bin_link_counts.items():
        bin_length_kb = bins[bin_key]['length'] / 1000  # Convert length to kilobases
        bin_density[bin_key] = count / bin_length_kb if bin_length_kb else 0

    return bin_density


def output_hic_density(hic_density, filename):
    """Output Hi-C density to a text file."""
    with open(filename, 'w') as f:
        for bin_key, density in sorted(hic_density.items()):
            f.write(f"{bin_key}\t{density:.4f}\n")


def update_HT_link_dict(HT_link_dict, ctg_i, ctg_j, len_i, len_j, coord_i, coord_j):

    def add_suffix(ctg, ctg_len, coord):

        if coord * 2 > ctg_len:
            return ctg + '_T'
        else:
            return ctg + '_H'

    ctg_name_HT_i = add_suffix(ctg_i, len_i, coord_i)
    ctg_name_HT_j = add_suffix(ctg_j, len_j, coord_j)

    HT_link_dict[(ctg_name_HT_i, ctg_name_HT_j)] += 1


def parse_hic_data(fa_dict, args):

    # 初始化字典
    full_link_dict = defaultdict(int)
    flank_link_dict = defaultdict(int)
    HT_link_dict = defaultdict(int)
    clm_dict = defaultdict(lambda: array('i'))

    # 解析 bam 文件
    format_options = [b'filter=flag.read1']
    alignments = bam_generator(args.bam, args.threads, format_options)
    print(f'Debugs {alignments}')

    for ref, mref, pos, mpos in alignments:
        if ref not in fa_dict or mref not in fa_dict:
            print(f'Debugs: {alignments} {ref} {mref} {pos} {mpos}')
            continue

        if ref == mref:     # 过滤 contig 内部的 links
            continue

        (ctg_i, coord_i), (ctg_j, coord_j) = sorted(((ref, pos + 1), (mref, mpos + 1)))
        ctg_name_pair = (ctg_i, ctg_j)

        if is_flank(coord_i, fa_dict[ctg_i][1], args.flank) and is_flank(coord_j, fa_dict[ctg_j][1], args.flank):
            flank_link_dict[ctg_name_pair] += 1

        full_link_dict[ctg_name_pair] += 1

        # 假设已定义辅助函数 update_clm_dict
        update_clm_dict(clm_dict, ctg_name_pair, fa_dict[ctg_i][1], fa_dict[ctg_j][1], coord_i - 1, coord_j - 1)

        # HT link pickle file
        update_HT_link_dict(HT_link_dict, ctg_i, ctg_j, fa_dict[ctg_i][1], fa_dict[ctg_j][1], coord_i, coord_j)
    return full_link_dict, flank_link_dict, clm_dict, HT_link_dict


def calculate_density(full_link_dict, fa_dict):
    ''' 计算每个 ctg 的 Hi-C links density '''
    link_count_dict = defaultdict(int)
    density_dict = {}

    # 累计每个 contig 的 links 数量
    for (ctg_i, ctg_j), link_count in full_link_dict.items():
        # if ctg_i != ctg_j:
        link_count_dict[ctg_i] += link_count
        if ctg_i != ctg_j:      # 避免重复计算自链接
            link_count_dict[ctg_j] += link_count

    # 计算每个 contig 的链接密度
    for ctg, total_links in link_count_dict.items():
        RE_count = fa_dict[ctg][2]      # 获得酶切位点个数
        density = total_links / RE_count if RE_count else 0 # 计算密度
        density_dict[ctg] = density
    return density_dict


def output_density(density_dict):
    with open('hic_density.txt', 'w') as out:
        for ctg, density in density_dict.items():
            out.write(f'{ctg}\t{density}\n')


def main():
    args = parse_args()

    ### 解析 FASTA
    fa_dict = parse_fasta(args.fasta)

    full_link_dict, flank_link_dict, clm_dict, HT_link_dict = parse_hic_data(fa_dict, args)
    #
    output_pickle(full_link_dict, 'full_links.pkl')
    output_pickle(flank_link_dict, 'flank_links.pkl')
    output_pickle(HT_link_dict, 'HT_links.pkl')

    # clm
    output_clm(clm_dict)

    # full links
    out = open('full.links.txt', 'w')
    for key, value in full_link_dict.items():
        out.write(str(key) + '\t' + str(value) + '\n')
    out.close()

    # flank links
    out = open('flank.links.txt', 'w')
    for key, value in flank_link_dict.items():
        out.write(str(key) + '\t' + str(value) + '\n')
    out.close()
    # del flank_link_dict
    # gc.collect()

    # HT links
    out  = open('HT_links.txt', 'w')
    for key, value in HT_link_dict.items():
        out.write(str(key) + '\t' + str(value) + '\n')
    out.close()

    # full links density
    density_dict = calculate_density(full_link_dict, fa_dict)
    output_density(density_dict)

    # Calculate Hi-C density for bins
    # bins = divide_into_bins(fa_dict, 10000)  # Define bin size, e.g., 10kb
    # bin_hic_density = calculate_bin_hic_density(full_link_dict, bins)
    # output_hic_density(bin_hic_density, 'bin_hic_density.txt')


if __name__ == '__main__':
    main()