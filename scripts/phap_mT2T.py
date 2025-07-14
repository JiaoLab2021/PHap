#!/usr/bin/env python
'''
time: 2024-05-16
author: pxxiao
version: 1.0
----
description: 生成 consensus 序列
----
time: 2024-05-22
author: pxxiao
version: 2.0
description:
    1. all_vs_all，判断 contig 之间的距离；
    2. 距离较近的 contigs 进行两两比对；
    3. 合并 PAF，对结果进行过滤
        ① 过滤 query 长度小于 ？ 的 alignment
        ② 过滤 length 长度小于 ？ 的 alignment
        ③ 计算 query 的比对率：match length / query length * 100，大于 ？ 则认为是冗余的序列，可以考虑去除
    4. 数据结构：(contig1, contig2): [contig1 length, contig2 length, contig1 match length, contig2 match length, match ratio]
----
time: 2025-04-23
author: pxxiao
version: 3.0
description:
    fix bugs: remove_redundancy()
'''


import argparse
from util import *


def parse_args():
    parser = argparse.ArgumentParser('Get mosaic T2T (mT2T) reference from primary contig assembly (p_ctg).')
    # Parameters for mT2T
    mT2T = parser.add_argument_group('>>> mT2T reference')
    mT2T.add_argument('--p_ctg', required=True,
                                        help='path to p_ctg file')
    mT2T.add_argument('--min_distance', required=False, type=float, default=0.15,
                                        help='minimum distance between two contigs [0.15]')
    mT2T.add_argument('--match_ratio_target', required=False, type=float, default=0.05,
                                        help='minimum match ratio between contigs [0.05]')
    mT2T.add_argument('--match_ratio_query', required=False, type=float, default=0.1,
                                        help='minimum match ratio between contigs [0.1]')
    mT2T.add_argument('--r2_threshold', required=False, type=float, default=0.5,
                      help='minimum R² value to consider alignment as collinear [0.5]')
    mT2T.add_argument('--min_contig_length', required=False, type=int, default=100000,
                                        help='minimum contig length [100000]')
    mT2T.add_argument('--min_alignment_length', required=False, type=int, default=200,
                                        help='minimum alignment length [200]')
    mT2T.add_argument('--min_chr_length', required=False, type=int, default=10000000,
                                        help='minimum chromosomes length [10000000]')
    parser.add_argument('--internal_margin_ratio', type=float, default=0.05,
                        help='Proportional margin to define internal alignments (default: 0.05)')
    parser.add_argument('--internal_ratio_threshold', type=float, default=0.75,
                        help='Ratio of internal alignments to consider a contig contained (default: 0.75)')

    # global
    parser.add_argument('--threads', required=False, type=int, default=1, help='number of threads [1]')
    parser.add_argument('--process', required=False, type=int, default=50, help='number of processes [50]')

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    cwd = os.getcwd()

    ### step0:
    dic_contig = fasta_read_filter(args.p_ctg, args.min_contig_length)

    ### step1: 拆分 contig.fa -- 一个序列为一个文件
    os.makedirs(f'{cwd}/01.mT2T/01.split.fa', exist_ok=True)
    file_split_list = split_fa(dic_contig, args.min_contig_length, f'{cwd}/01.mT2T/01.split.fa')     # > min_contig_length
    # print(file_split_list)

    ### step2: 计算 contig 之间的距离, 筛选出相似性高的 contig pairs
    os.makedirs(f'{cwd}/01.mT2T/02.distance', exist_ok=True)
    cal_distance(file_split_list, f'{cwd}/01.mT2T/02.distance', args)

    ### step3: contig pairs 进行比对, 得到 03.alignment/merge.fa
    os.makedirs(f'{cwd}/01.mT2T/03.alignment', exist_ok=True)
    contig_pair_aln(f'{cwd}/01.mT2T/02.distance/similar_contig_pairs.txt', f'{cwd}/01.mT2T/03.alignment', args)

    ### step5: 去除冗余的 contig
    os.makedirs(f'{cwd}/01.mT2T/04.remove.redundancy', exist_ok=True)
    remove_redundancy(f'{cwd}/01.mT2T/03.alignment/merge.paf', f'{cwd}/01.mT2T/04.remove.redundancy', dic_contig, args)

    # ### step6: 获得 consensus 序列
    # os.makedirs(f'{cwd}/01.mT2T/05.consensus', exist_ok=True)


'''
### to do list:
1. 只处理 Nx 的 contig
    · 保留下来的 contig，还是有很多比较短的序列，例如：100 kb 以下的 contig 数量有~420，长度为~26Mb；大于 100 kb 的 contig 数量有~170，长度为~802Mb。
2. 
'''


if __name__ == '__main__':
    main()