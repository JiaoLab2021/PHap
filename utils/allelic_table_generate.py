#!/usr/bin/env python
'''
time: 2024-05-25
author: pxxiao
version: 1.0
----
description: 通过划分染色体窗口的形式，找到 allelic contig
'''

import re
from collections import defaultdict
import os
import numpy as np
import argparse

'''
alignment filter:
    1. Skip supplementary alignments
    2. 每一个 contig，只保留最长匹配的 reference
    # 3. 过滤离群值
'''


## 解析 contig type
def parse_contig_type_based_on_dosage(file_contig_type):
    dic_dosage_contig_type = {}
    for lines in open(file_contig_type, 'r'):
        if lines.startswith('contig_ID'):
            continue
        else:
            line = lines.strip().split()
            dic_dosage_contig_type[line[0]] = line[2]
    return dic_dosage_contig_type


def parse_paf_line(line, min_align_length, min_unitig_length):
    parts = line.strip().split('\t')
    query_id = parts[0]
    query_len = int(parts[1])
    query_start = int(parts[2])
    query_end = int(parts[3])
    target_id = parts[5]
    target_len = int(parts[6])
    target_start = int(parts[7])
    target_end = int(parts[8])
    match_len = int(parts[9])

    ## filter: primary alignment or supplementary alignment
    align_type = parts[16][-1]
    if align_type == 'S':
        return None

    ## filter: min alignment length
    if match_len < min_align_length:
        return None

    ## filter: min unitig length
    if query_len < min_unitig_length:
        return None

    return {
        'query_id': query_id,
        'query_len': query_len,
        'query_start': query_start,
        'query_end': query_end,
        'target_id': target_id,
        'target_len': target_len,
        'target_start': target_start,
        'target_end': target_end,
        'match_len': match_len
    }


def read_paf(file_paf, min_align_length, min_unitig_length, chr_num):
    alignments = []
    dic_scaffold_length = {}        # 存储 Scaffolds 长度
    with open(file_paf, 'r') as f:
        for line in f:
            alignment = parse_paf_line(line, min_align_length, min_unitig_length)
            if alignment:                       # Skip supplementary alignment
                alignments.append(alignment)
                dic_scaffold_length[alignment['target_id']] = alignment['target_len']

    ## 获得染色体匹配数量的 scaffold id
    scffolds_chrs = sorted(dic_scaffold_length.items(), key=lambda x: x[1], reverse=True)[:chr_num]
    scffolds_chr = [scaffold_id for scaffold_id, length in scffolds_chrs]
    print(scffolds_chr)
    return alignments, scffolds_chr


def filter_longest_reference_per_contig(alignments):
    ''' 每一个 contig，只保留最长匹配的 reference '''
    contig_to_refs = defaultdict(lambda: defaultdict(int))

    # Calculate total match length per reference for each contig
    for aln in alignments:
        contig_to_refs[aln['query_id']][aln['target_id']] += aln['match_len']

    filtered_alignments = []
    for query_id, refs in contig_to_refs.items():
        # Find the reference with the longest total match length
        best_ref = max(refs, key=refs.get)
        filtered_alignments.extend(
            aln for aln in alignments if aln['query_id'] == query_id and aln['target_id'] == best_ref
        )
    # for i in filtered_alignments:
    #     print(i)
    return filtered_alignments


def mad(arr):
    median = np.median(arr)
    return np.median(np.abs(arr - median))


def filter_outliers(alignments):
    if not alignments:
        return alignments

    target_starts = [aln['target_start'] for aln in alignments]
    target_ends = [aln['target_end'] for aln in alignments]

    start_median = np.median(target_starts)
    start_mad = mad(target_starts)

    end_median = np.median(target_ends)
    end_mad = mad(target_ends)

    # for aln in alignments:
    #     if aln['query_id'] == 'utg000018l' and aln['target_id'] == 'scaffold_10':
    #         print(start_median, start_mad, end_median, end_mad)

    filtered_alignments = [aln for aln in alignments if abs(aln['target_start'] - start_median) <= 3 * start_mad and abs(aln['target_end'] - end_median) <= 3 * end_mad]

    return filtered_alignments


def split_alignment(aln, bin_size):
    target_start = aln['target_start']
    target_end = aln['target_end']
    split_alignments = []

    while target_start < target_end:
        bin_start = (target_start // bin_size) * bin_size
        bin_end = bin_start + bin_size
        aln_copy = aln.copy()
        aln_copy['target_start'] = max(target_start, bin_start)
        aln_copy['target_end'] = min(target_end, bin_end)
        aln_copy['match_len'] = aln_copy['target_end'] - aln_copy['target_start']
        split_alignments.append(aln_copy)
        target_start = bin_end

    return split_alignments


def assign_to_bins_and_filter(alignments, bin_size):
    bins = defaultdict(list)
    for aln in alignments:
        split_aligns = split_alignment(aln, bin_size)
        for split_aln in split_aligns:
            bin_start = (split_aln['target_start'] // bin_size) * bin_size
            bin_end = bin_start + bin_size
            bin_key = (split_aln['target_id'], bin_start, bin_end)
            bins[bin_key].append(split_aln)

    return bins


def merge_contig_alignments1(bins, dic_contig_type):
    merged_bins = defaultdict(lambda: defaultdict(int))
    for bin_key, alignments in bins.items():
        for aln in alignments:
            query_id = aln['query_id']
            match_len = aln['match_len']
            unitig_type = dic_contig_type.get(query_id, 'haplotig')

            # 根据 unitig 类型调整 match_len
            if unitig_type == 'diplotig':
                match_len *= 2
            elif unitig_type == 'triplotig':
                match_len *= 3
            elif unitig_type == 'tetraplotig':
                match_len *= 4

            merged_bins[bin_key][aln['query_id']] += aln['match_len']
    return merged_bins


def get_top_n_contigs_per_bin1(merged_bins, dic_contig_type, top_n):
    top_contigs_per_bin = {}
    for bin_key, contigs in merged_bins.items():
        sorted_contigs = sorted(contigs.items(), key=lambda x: x[1], reverse=True)

        adjusted_top_n = top_n
        selected_contigs = []

        for query_id, match_len in sorted_contigs:
            unitig_type = dic_contig_type.get(query_id, 'haplotig')

            # 调整 top_n 根据 unitig 类型
            if unitig_type == 'diplotig':
                if adjusted_top_n > 0:
                    selected_contigs.append((query_id, match_len))
                    adjusted_top_n -= 2
            elif unitig_type == 'triplotig':
                if adjusted_top_n > 1:
                    selected_contigs.append((query_id, match_len))
                    adjusted_top_n -= 3
            elif unitig_type == 'tetraplotig':
                if adjusted_top_n > 2:
                    selected_contigs.append((query_id, match_len))
                    adjusted_top_n -= 4
            else:  # 其他类型视为 haplotig
                if adjusted_top_n > 0:
                    selected_contigs.append((query_id, match_len))
                    adjusted_top_n -= 1

            if adjusted_top_n <= 0:
                break

        top_contigs_per_bin[bin_key] = selected_contigs
    return top_contigs_per_bin


def merge_contig_alignments(bins, dic_contig_type):
    merged_bins = defaultdict(lambda: defaultdict(int))
    for bin_key, alignments in bins.items():
        for aln in alignments:
            merged_bins[bin_key][aln['query_id']] += aln['match_len']
    return merged_bins


def get_top_n_contigs_per_bin(merged_bins, dic_contig_type, top_n):
    top_contigs_per_bin = {}
    for bin_key, contigs in merged_bins.items():
        sorted_contigs = sorted(contigs.items(), key=lambda x: x[1], reverse=True)
        top_contigs_per_bin[bin_key] = sorted_contigs[:top_n]
    return top_contigs_per_bin


def write_top_contigs_per_bin(file_paf, min_align_length, min_unitig_length, bin_size, top_n, chr_num, output_file, out_allelic_table, file_contig_type):

    dic_contig_type = parse_contig_type_based_on_dosage(file_contig_type)

    alignments, scffolds_chr = read_paf(file_paf, min_align_length, min_unitig_length, chr_num)

    filtered_alignments = filter_longest_reference_per_contig(alignments)

    bins = assign_to_bins_and_filter(filtered_alignments, bin_size)

    merged_bins = merge_contig_alignments(bins, dic_contig_type)
    top_contigs_per_bin = get_top_n_contigs_per_bin(merged_bins, dic_contig_type, top_n)

    out_allelic = open(out_allelic_table, 'w')
    with open(output_file, 'w') as out:
        out.write("#target_id\tbin_start\tbin_end\tquery_id\tmatch_len\n")
        for bin_key, contigs in top_contigs_per_bin.items():
            target_id, bin_start, bin_end = bin_key
            if target_id not in scffolds_chr: continue
            out_allelic.write(f'{target_id}\t{bin_start}\t{bin_end}\t')
            for query_id, match_len in contigs:
                out.write(f"{target_id}\t{bin_start}\t{bin_end}\t{query_id}\t{match_len}\n")
                out_allelic.write(query_id + '\t')
            out_allelic.write('\n')
    out_allelic.close()


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf_file", type=str, required=True, help="PAF file")
    parser.add_argument('--min_align_length', type=int, default=1000, help="Minimum alignment length [1000]")
    parser.add_argument("--min_unitig_length", type=int, default=100000, help="Minimum unitig length [100000]")
    parser.add_argument("--bin_size", type=int, required=False, default=100000,
                        help="Bin size [100000]")
    parser.add_argument("--top_n", type=int, required=False, default=4,
                        help="The maximum number of contigs corresponding to each bin."
                             " This corresponds to the number of haplotypes, for example, tetraploid is 4 [4]")
    parser.add_argument('--out_top_contigs_per_bin', type=str, required=False, default='top_contigs_per_bin.txt',
                        help="Output file [top_contigs_per_bin.txt]")
    parser.add_argument('--out_allelic_table', type=str, required=False, default='allelic.ctg.table',
                        help="Output file [allelic.ctg.table]")
    parser.add_argument('--chr_num', type=int, default=12, help='number of chromosomes [12]')
    parser.add_argument('--contig_type', type=str, required=True, help='contig type based on dosage analysis')
    args = parser.parse_args()
    return args


def main():

    args = parse_arguments()
    # 使用示例
    file_paf = args.paf_file
    min_align_length = args.min_align_length
    min_unitig_length = args.min_unitig_length
    bin_size = args.bin_size
    top_n = args.top_n
    out_top_contigs_per_bin = args.out_top_contigs_per_bin
    out_allelic_table = args.out_allelic_table
    chr_num = args.chr_num
    file_contig_type = args.contig_type

    write_top_contigs_per_bin(file_paf, min_align_length, min_unitig_length, bin_size, top_n, chr_num, out_top_contigs_per_bin, out_allelic_table, file_contig_type)


if __name__ == '__main__':
    main()