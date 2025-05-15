#!/usr/bin/env python
'''
time: 2024-06-06
author: pxxiao
version: 1.0
description: 我们已经获得了 allelic table，dosage contig type，和 Hi-C links；利用这些信息，对 allelic unitig 分型，聚类为 4 个 group。
----
time: 2024-06-13
author: pxxiao
version: 2.0
description: 优化 Hi-C links 的分配
----
time: 2024-06-15
author: pxxiao
version: 3.0
description: 解决 allelic table 不准确的问题
            genotype 的时候，会发现 allelic table 不准确的现象，即 allelic unitig 之间存在重叠，但是却属于一个 haplotype
            针对这种情况，我们需要判断 allelic unitig 之间的信号与grouped unitig之间的信号强弱，
            如果 allelic unitig 之间的信号强于与 grouped unitig之间的信号强弱，则将allelic unitig 归为一类。
----
time: 2024-06-16
author: pxxiao
version: 4.0
description: 解决 collapsed 的问题
            目前，并没有考虑 collapsed unitig，现在需要把这个信息考虑进去。
----
time: 2024-06-18
author: pxxiao
version: 5.0
description: 解决 collapsed 的问题
            对于 allelic table 不规则的 group，即：如果有两个 group 为空，new_unitig为 triplotig，这时候不能正确分配 new_unitig。修改代码，完成这个功能。
---
version7.0:
    chr1 genotype 亲本验证结果：
        g1/log_g1_trioevale_out:N	1230	357877	0.003425
        g2/log_g2_trioevale_out:N	1418120	26237	0.018165
        g3/log_g3_trioevale_out:N	24258	219496	0.099518
        g4/log_g4_trioevale_out:N	829353	58010	0.065373
----
time: 2024-06-18
author: pxxiao
version7.0.v2
description: 进行优化
    1. 增加 flank or full 模式选择参数，如果是 flank，则 RE 酶切位点个数为两端 flank 区域的 RE 个数，而不是全长的个数。
---

To do:
    1. 有的多倍体基因组不需要 collapsed 恢复，但是其他挂载软件，例如：HapHiC 又不能解决这些基因组的挂载问题。
        这时候需要添加一个参数，跳过 collapsed 的步骤。

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


## 解析 FASTA: ctg: [seq, seq length, RE site's count]
def parse_fasta(fasta, flank, RE='GATC'):
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

        if flank is not None:
            # Ensure flank does not exceed sequence length
            flank_size = min(flank, len(seq) // 2)
            seq_segment = seq[:flank_size] + seq[-flank_size:]
        else:
            seq_segment = seq

        # Count RE sites in the determined sequence segment
        RE_sites = count_RE_sites(seq_segment, RE) + 1  # Add pseudo-count of 1 to prevent division by zero
        fa_dict[ctg] = [seq, len(seq), RE_sites]

        # Add pseudo-count of 1 to prevent division by zero (as what ALLHiC does)
        # RE_sites = count_RE_sites(seq, RE) + 1
        # fa_dict[ctg] = [seq, len(seq), RE_sites]
    return fa_dict


def adjust_hic_signals(hic_data, dic_contig_type):
    """
    Adjust Hi-C interaction signals based on the types of unitigs involved.

    Args:
    hic_data (dict): Dictionary with keys as (ctg1, ctg2) and values as interaction counts.
    unitig_types (dict): Dictionary with unitig names as keys and their types as values.

    Returns:
    dict: Adjusted Hi-C interaction data.
    """
    adjusted_hic_data = {}
    weights = {
        'haplotig': 1.0,
        'diplotig': 0.5,
        'triplotig': 0.333,
        'tetraplotig': 0.25
    }

    for (ctg1, ctg2), count in hic_data.items():
        weight1 = weights.get(dic_contig_type.get(ctg1, 'haplotig'), 1)
        weight2 = weights.get(dic_contig_type.get(ctg2, 'haplotig'), 1)
        adjusted_count = count * weight1 * weight2
        adjusted_hic_data[(ctg1, ctg2)] = adjusted_count

        print(f"Adjusted Hi-C signal between {ctg1} and {ctg2}: Original = {count}, Adjusted = {adjusted_count}")

    return adjusted_hic_data


def parse_pickle(fa_dict, pickle_file, dic_contig_type):
    ''' 从 pickle_file 获得 full_links '''
    with open(pickle_file, 'rb') as f:
        full_link_dict = pickle.load(f)

    adjust_hic_signals(full_link_dict, dic_contig_type)

    # sort by contig length
    sorted_ctg_list = sorted([(ctg, fa_dict[ctg][1]) for ctg in fa_dict], key=lambda x: x[1], reverse=True)
    # RE site
    RE_site_dict = {ctg: ctg_info[2] for ctg, ctg_info in fa_dict.items()}

    return full_link_dict, sorted_ctg_list, RE_site_dict


def load_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data


def get_link_list(dic_pair_hic):
    ''' 将 contig pair 的 Hi-C 信号，转换为 contig 的 Hi-C 信号 '''
    dic_contig_hic = defaultdict(list)
    ## 处理 contig pair links dic
    for (contig1, contig2), links in dic_pair_hic.items():
        dic_contig_hic[contig1].append((contig2, links))
        if contig1 != contig2:
            dic_contig_hic[contig2].append((contig1, links))

    # 对每个 contig 的 links 按照从大到小排序
    for contig, links in dic_contig_hic.items():
        links.sort(key=lambda x: x[1], reverse=True)

    return dic_contig_hic


def parse_contig_type(file_contig_type):
    df = pd.read_csv(file_contig_type, delim_whitespace=True)
    dic_contig_type = dict(zip(df['contig_ID'], df['contig_type']))
    return dic_contig_type


def unitig_overlap_ratio(file_allelic_table):
    ''' 根据 allelic table 获得 unitig 之间的 overlap ratio '''
    def read_allelic_table(file_path):
        '''
        读取 allelic table 文件并解析数据，返回每个 unitig 的 bin 覆盖情况。
        :param file_path:
        :return: 包含每个 unitig 覆盖 bin 的字典
        '''
        unitig_bins = {}

        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split()
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                unitigs = columns[3:]

                for unitig in unitigs:
                    if unitig not in unitig_bins:
                        unitig_bins[unitig] = set()
                    unitig_bins[unitig].add((start, end))

        return unitig_bins

    def calculate_overlap_ratio(unitig_bins, unitig1, unitig2):
        '''
        计算两个 unitig 之间的 overlap ratio
        :param unitig_bins: 包含每个 unitig 覆盖 bin 的字典
        :param unitig1: 第一个 unitig
        :param unitig2: 第二个 unitig
        :return: 两个 overlap ratio，(ratio1, ratio2)
        '''
        bins1 = unitig_bins.get(unitig1, set())
        bins2 = unitig_bins.get(unitig2, set())

        if not bins1 or not bins2:
            return 0, 0

        overlap_bins = bins1 & bins2
        overlap_count = len(overlap_bins)

        ratio1 = overlap_count / len(bins1) if bins1 else 0
        ratio2 = overlap_count / len(bins2) if bins2 else 0

        return ratio1, ratio2

    unitig_bins = read_allelic_table(file_allelic_table)

    overlap_ratios = {}
    unitigs = list(unitig_bins.keys())

    for i in range(len(unitigs)):
        for j in range(i + 1, len(unitigs)):
            unitig1 = unitigs[i]
            unitig2 = unitigs[j]
            ratio1, ratio2 = calculate_overlap_ratio(unitig_bins, unitig1, unitig2)
            overlap_ratios[(unitig1, unitig2)] = ratio1
            overlap_ratios[(unitig2, unitig1)] = ratio2

    return overlap_ratios


def genotype_allelic_table(file_allelic_table, dic_contig_type, dic_contig_hic, dic_overlap_ratios, RE_site_dict):
    def get_contig_types(contigs, dic_contig_type):
        return [dic_contig_type.get(contig, 'unknown') for contig in contigs]

    def find_group_for_new_unitig1(new_unitig, groups, dic_contig_hic, dic_contig_type, dic_RE_site):
        ''' 根据 Hi-C links 的密度找到适合的 group 分配 new_unitig '''

        ## 获取 new_unitig 的类型
        contig_type = dic_contig_type.get(new_unitig, 'haplotig')
        # 所需要的 group 数量
        required_group_count = {'haplotig': 1, 'diplotig': 2, 'triplotig': 3, 'tetraplotig': 4}.get(contig_type, 1)

        ## 计算每个 group 的 Hi-C links 总数和密度
        group_density = []
        for idx, group in enumerate(groups):
            total_hic_links = 0
            total_re_sites = 0

            # 计算该 group 的 Hi-C links 总数
            for unitig in group:
                hic_links = dic_contig_hic.get(unitig, [])
                hic_count = sum(count for u, count in hic_links if u == new_unitig)
                total_hic_links += hic_count
                total_re_sites += dic_RE_site.get(unitig, 1)  # 默认每个 unitig 至少有一个酶切位点，避免除零错误

            # 计算 Hi-C links 密度
            if total_hic_links > 0:
                density = total_hic_links / total_re_sites
            else:
                density = 0

            group_density.append((idx, density))

        # 筛选出密度最高的 groups
        group_density.sort(key=lambda x: x[1], reverse=True)

        # debugs -- 输出 new_unitig 以及 候选的 group links density
        print(f'Debugs: new_unitig group_density {new_unitig} {group_density}')

        selected_groups = [idx for idx, density in group_density[:required_group_count]]

        # 如果找到合适的 groups，则返回它们的索引
        if selected_groups:
            return selected_groups
        else:
            print(f"No suitable group found for {new_unitig} based on Hi-C link density.")
            return None

    def find_group_for_new_unitig(new_unitig, groups, dic_contig_hic, dic_contig_type, threshold=50):

        ## new_unitig 的类型
        contig_type = dic_contig_type.get(new_unitig, 'haplotig')

        ## 找到每个 group 中最后一个 unitig
        # last_unitigs = [group[-1] if group else None for group in groups]

        valid_last_unitigs = []
        for group in groups:
            last_unitig = group[-1] if group else None
            # 检查并处理 collapsed unitig
            while last_unitig and dic_contig_hic.get(last_unitig) in  ['diplotig', 'triplotig', 'tetraplotig']:
                # 如果当前 last_unitig 是 collapsed unitig，则回溯找前一个 unitig
                last_index = group.index(last_unitig)
                last_unitig = group[last_index - 1] if last_index > 0 else None
            valid_last_unitigs.append(last_unitig)

        ## 计算 new_unitig 与 grouped_last_unitig 的 Hi-C links
        last_link_counts = []
        for last_unitig in valid_last_unitigs:
            if last_unitig:
                links = dic_contig_hic.get(last_unitig, [])     # last unitig 的 Hi-C links 列表
                link_count = next((count for u, count in links if u == new_unitig), 0)      # last unitig 与 new_unitig 的 Hi-C link
                last_link_counts.append(link_count)
            else:
                last_link_counts.append(0)

        ## 计算 new_unitig 与所有 group 的 link 数量
        total_links_per_group = {}
        # 遍历所有的 groups
        for idx, group in enumerate(groups):
            total_links = 0
            # 遍历 group 中的每个 unitig
            for unitig in group:
                # 获取 unitig 的 Hi-C 链接数据
                links = dic_contig_hic.get(unitig, [])
                # 累加与 new_unitig 的 links
                total_links += sum(count for u, count in links if u == new_unitig)
            # 将计算结果存储到字典中
            total_links_per_group[idx] = total_links

        ## 计算 new_unitig 与所有 group 的 links density


        # debugs:
        if new_unitig == 'utg000961l':
            print(f'Debugs contig_type {contig_type}')
            print(f'Debugs xpx last_link_counts {last_link_counts}')    # [14, 0, 0]

        ## 如果 group 为空，并且满足 collapsed unitig 对应的数量，则直接返回 empty group 的索引；否则，返回 empty group 的索引和 Hi-C 匹配的 group 索引
        empty_groups = [i for i, group in enumerate(groups) if not group]
        print(f'Debugs xpx empty_groups {empty_groups}')
        required_empty_count = {
            'haplotig': 1,
            'diplotig': 2,
            'triplotig': 3,
            'tetraplotig': 4
        }.get(contig_type, 1)

        result_indices = []

        if len(empty_groups) >= required_empty_count:
            return empty_groups[:required_empty_count]      # tetraplotig 的时候，是不是直接返回四个group？
        else:
            result_indices.extend(empty_groups)

        print(f'Debugs result_indices-1 {new_unitig} {result_indices}')

        ## 根据 unitig 类型确定需要的 group 数量
        required_group_count = {
            'haplotig': 1,
            'diplotig': 2,
            'triplotig': 3,
            'tetraplotig': 4
        }.get(contig_type, 1)

        ## 检查 new_unitig 与 group 的所有 Hi-C 数量差异
        ## 找到链接数最多的group
        sorted_groups = sorted(range(len(total_links_per_group)), key=lambda i: total_links_per_group[i], reverse=True)
        best_group_indices = sorted_groups[:required_group_count]

        ## 如果 group 为空，并且满足 collapsed unitig 对应的数量，则直接返回 empty group 的索引；否则，返回 empty group 的索引和 Hi-C 匹配的 group 索引
        empty_groups = [i for i, group in enumerate(groups) if not group]
        if len(empty_groups) >= required_group_count:
            return empty_groups[:required_group_count]

        ## 如果最高链接数大于阈值，返回链接数最高的group，否则返回None
        if all(total_links_per_group[idx] > threshold for idx in best_group_indices):
            return best_group_indices
        else:
            print(
                f"Debug: No Hi-C signal found for {new_unitig} meeting the threshold with any group. Not adding to any group.")
            return None

    def get_remaining_to_original_index(groups, processed_groups):
        """
        为 remaining_groups 中的每个组生成唯一标识符，并映射到原始组的索引。
        """
        remaining_groups = [groups[i] for i in range(len(groups)) if i not in processed_groups]
        remaining_to_original_index = {}
        used_indices = set()

        for i, group in enumerate(remaining_groups):
            for j, g in enumerate(groups):
                if g == group and j not in used_indices:
                    remaining_to_original_index[i] = j
                    used_indices.add(j)
                    break  # 只匹配第一个未使用的组

        return remaining_groups, remaining_to_original_index

    # 创建四个 list，用来存储 haplotype‘s unitig
    g1, g2, g3, g4 = [], [], [], []
    groups = [g1, g2, g3, g4]
    first_line = True  # 标志位，用于标记第一行
    for lines in open(file_allelic_table, 'r'):
        line = lines.strip().split()
        ## 染色体的开头
        if first_line:
            unitigs = line[3:]
            unitig_types = get_contig_types(unitigs, dic_contig_type)
            print(unitig_types)
            ## 一个 diplotig，两个 haplotig
            if unitig_types.count('diplotig') == 1 and len(unitigs) == 3:
                diplotig_index = unitig_types.index('diplotig')
                haplotig_indices = [i for i, t in enumerate(unitig_types) if t == 'haplotig']
                g1.append(unitigs[diplotig_index])
                g2.append(unitigs[diplotig_index])
                g3.append(unitigs[haplotig_indices[0]])
                g4.append(unitigs[haplotig_indices[1]])

            ## 一个  triplotig，一个 haplotig
            elif unitig_types.count('triplotig') == 1 and len(unitigs) == 2:
                triplotig_index = unitig_types.index('triplotig')
                haplotig_index = unitig_types.index('haplotig')
                g1.append(unitigs[triplotig_index])
                g2.append(unitigs[triplotig_index])
                g3.append(unitigs[triplotig_index])
                g4.append(unitigs[haplotig_index])          # 不是 g3, 是 g4! 这里竟然写错了！

            ## 一个 tetraplotig
            elif unitig_types.count('tetraplotig') == 1 and len(unitigs) == 1:
                tetraplotig_index = unitig_types.index('tetraplotig')
                g1.append(unitigs[tetraplotig_index])
                g2.append(unitigs[tetraplotig_index])
                g3.append(unitigs[tetraplotig_index])
                g4.append(unitigs[tetraplotig_index])

            ## 四个 unitig 全部为 haplotig
            elif len(unitigs) == 4:
                g1.append(unitigs[0])
                g2.append(unitigs[1])
                g3.append(unitigs[2])
                g4.append(unitigs[3])
            elif unitig_types.count('diplotig') ==0 and unitig_types.count('triplotig') == 0 and unitig_types.count('tetraplotig') == 0:
                # 只有 haplotig，但数量不够四个
                for i, unitig in enumerate(unitigs):
                    groups[i % 4].append(unitig)
            else:
                print('Please check your input')
            first_line = False  # 标志第一行已经处理完毕
        else:
            ##
            # 处理染色体后端的情况，不做处理或进行其他操作
            unitigs = line[3:]
            unitig_types = get_contig_types(unitigs, dic_contig_type)       # unitig type [list]: haplotig or diplotig or triplotig or tetraplotig
            print(f'Debugs: unitig_types {unitig_types}')

            processed_groups = set()  # 跟踪在这一行已经处理过的 groups
            processed_unitigs = set()   # 已经处理过的 unitigs

            ## 首先，检查所有 unitig 是否已经存在某个 group 中
            for unitig in unitigs:      # 遍历 unitigs
                for group in groups:
                    if unitig in group and groups.index(group) not in processed_groups:
                        group.append(unitig)
                        processed_unitigs.add(unitig)
                        processed_groups.add(groups.index(group))
                        # break     # 这里运算速度变慢了

            ## 将每个 unitig 分配给已有的 group
            for unitig in unitigs:
                if unitig in processed_unitigs: continue

                for group in groups:
                    # if unitig in group and groups.index(group) not in processed_groups:
                    if unitig in group:
                        print(f'Debugs: {unitig}')
                        group.append(unitig)
                        processed_unitigs.add(unitig)
                        processed_groups.add(groups.index(group))
                        # break     # # 这里运算速度变慢了
            print(f'Debugs processed_groups {processed_groups}')

            ##
            for i, unitig in enumerate(unitigs):
                if unitig in processed_unitigs:
                    print('unitig', i, unitig, unitigs)
                    continue

                # ## 针对 haplotype 开头的 unitig 进行处理
                # remaining_groups = [groups[i] for i in range(len(groups)) if i not in processed_groups]     # 待处理的 groups: [[u1, u2...], [u11, u22, ...], ]
                # 假设每个 remaining_groups 中的组在原始 groups 中的索引是已知的，可以通过一个字典进行映射
                # remaining_to_original_index = {i: groups.index(group) for i, group in enumerate(remaining_groups)}

                remaining_groups, remaining_to_original_index = get_remaining_to_original_index(groups,
                                                                                                processed_groups)

                print(f'Debugs processed_groups remaining_groups {unitig} {processed_groups} {remaining_groups}')
                print(f'Debugs: remaining_to_original_index {remaining_to_original_index}')
                print(len(g1), g1)
                print(len(g2), g2)
                print(len(g3), g3)
                print(len(g4), g4)

                ## allelic unitig 之间的信号值 -- （这里没有考虑 contig type）
                allelic_unitig_signals = [
                    next((count for u, count in dic_contig_hic.get(unitig, []) if u == other_unitig), 0)
                    for other_unitig in unitigs if other_unitig != unitig
                ]
                max_allelic_signal = max(allelic_unitig_signals, default=0)

                ## allelic other unitig
                allelic_other_unitigs = [u for u in unitigs if u != unitig]

                ## unitig 与 group unitig 之间的信号值
                group_signals = []
                for group in groups:
                    group_signal = max(
                        (next((count for u, count in dic_contig_hic.get(u, []) if u == unitig), 0)
                         for u in group), default=0
                    )
                    group_signals.append(group_signal)
                max_group_signal = max(group_signals, default=0)

                print(f'Debugs: unitig, unitigs, max_allelic_signal, max_group_signal, remaining_groups {unitig} {unitigs} {max_allelic_signal} {max_group_signal} {remaining_groups}')

                ## 如果 allelic 之间的信号值更大，并且两者之间的 overlap 比例较低，则将 allelic unitig 归为一个 haplotype
                if max_allelic_signal >= max_group_signal > 100 and dic_overlap_ratios[(unitig, allelic_other_unitigs[allelic_unitig_signals.index(max_allelic_signal)])] < 0.2 and dic_overlap_ratios[(allelic_other_unitigs[allelic_unitig_signals.index(max_allelic_signal)], unitig)] < 0.2:
                    # allelic unitig 之间的信号更强，将 allelic unitig 归为一类
                    max_signal_index = allelic_unitig_signals.index(max_allelic_signal)
                    allelic_unitig = unitigs[max_signal_index]

                    # 找到 allelic unitig 所在的 group 并添加当前 unitig
                    for group in groups:
                        if allelic_unitig in group:
                            group.append(unitig)
                            processed_unitigs.add(unitig)
                            # break
                else:
                    ## 按照原始逻辑分配 unitig
                    # target_group_indices = find_group_for_new_unitig(unitig, remaining_groups, dic_contig_hic, dic_contig_type)
                    target_group_indices = find_group_for_new_unitig1(unitig, remaining_groups, dic_contig_hic, dic_contig_type, RE_site_dict)
                    print(f'Debugs target_group_indices {unitig} {target_group_indices}')

                    if target_group_indices is not None:
                        for target_group_index in target_group_indices:
                            target_group = remaining_groups[target_group_index]
                            # original_group_index = remaining_to_original_index[target_group_index]
                            original_group_index = remaining_to_original_index[target_group_index]
                            print(f'Debugs: target_group {target_group} {target_group_index}')      # Debugs: target_group [] 2
                            print(f'Debugs: original_group_index {original_group_index}')

                            # 如果 group 中只有一个元素且与 new_unitig 都没有 Hi-C 信号，则替换
                            if len(target_group) == 1 or len(set(target_group)) == 1:
                                last_unitig = target_group[-1]
                                links = dic_contig_hic.get(last_unitig, [])
                                link_count = next((count for u, count in links if u == unitig), 0)
                                # 开头的 unitig 没有 Hi-C 信号，不管是一个 unitig，还是多个相同的 unitig
                                if link_count == 0 and last_unitig not in dic_contig_hic:
                                    print(f'Debugs: Replacing {last_unitig} with {unitig} in group {target_group_index + 1}.')
                                    # 替换
                                    for i, u in enumerate(target_group):
                                        target_group[i] = unitig
                                else:
                                    target_group.append(unitig)
                            else:
                                groups[original_group_index].append(unitig)

                            processed_groups.add(original_group_index)
                    else:
                        #
                        if all(len(group) > 0 for group in groups):
                            print(f'Debugs: No valid Hi-C signal for {unitig}, discarding.')
                        else:
                            # 当 group 为空时，添加 unitig
                            empty_group = next(group for group in groups if len(group) == 0)
                            empty_group.append(unitig)

    return g1, g2, g3, g4


def list_to_file(list_1, file_name):
    out = open(file_name, 'w')
    list_1 = list(set(list_1))
    for i in list_1:
        out.write(i + '\n')
    out.close()


def list_to_cluster_file(list_1, list_name, file='group.cluster.txt'):
    out = open(file, 'a+')
    list_1 = list(set(list_1))
    out.write(list_name + '\t' + str(len(list_1)) + '\t' + ' '.join(list_1) + '\n')
    out.close()


def write_combined_genotypes_to_file(g1, g2, g3, g4, filename):
    """ 将 g1, g2, g3, g4 写入同一个文件，每一列是一个 genotype """
    # 找到最长的 list 的长度
    max_length = max(len(g1), len(g2), len(g3), len(g4))

    # 填充较短的 list 使其长度相同
    g1.extend([''] * (max_length - len(g1)))
    g2.extend([''] * (max_length - len(g2)))
    g3.extend([''] * (max_length - len(g3)))
    g4.extend([''] * (max_length - len(g4)))

    # 将组合后的内容写入文件
    with open(filename, 'w') as f:
        f.write('g1' + '\t' + 'g2' + '\t' + 'g3' + '\t' + 'g4' + '\n')
        for row in zip(g1, g2, g3, g4):
            f.write('\t'.join(row) + '\n')


def output_fa_from_list(list_1, dic_fasta, file_name):
    out = open(file_name, 'w')
    list_1 = list(set(list_1))
    for i in list_1:
        if i in dic_fasta:
            out.write('>' + i + '\n' + dic_fasta[i] + '\n')
        else:
            print(f"Warning: {i} not found in the fasta dictionary.")
    # for i in list_1:
    #     out.write('>' + i + '\n' + dic_fasta[i] + '\n')
    out.close()


def rename_id(g_fasta, merge):
    '''
    Mark 单倍型 ID，添加到 merge.fa 后面
    例如：g1.fa，ID 更新为>g1_unitigID
    '''
    bn = os.path.basename(g_fasta).replace('.fa', '')
    out = open(merge, 'a+')
    for lines in open(g_fasta, 'r'):
        if lines.startswith('>'):
            out.write('>' + bn + '_' + lines.strip().replace('>', '') + '\n')
        else:
            out.write(lines)


def false_genotype_unitig(g_list, dic_fasta, dic_contig_hic):
    ''' 被错误分型的 unitig '''
    corrected_list = []
    for u in g_list:
        if u not in dic_contig_hic:
            print(f'Debugs: No hi-c link, false genotype {u}')
        else:
            corrected_list.append(u)
    return corrected_list


def cluster(fasta, full_links, flank, contig_type, allelic_table_file, wd):
    dic_fasta = fasta_read(fasta)
    fa_dict = parse_fasta(fasta, flank)
    dic_contig_type = parse_contig_type(contig_type)
    full_link_dict, sorted_ctg_list, RE_site_dict = parse_pickle(fa_dict, full_links, dic_contig_type)
    out = open('full.links.txt', 'w')
    for key, value in full_link_dict.items():
        out.write(str(key) + '\t' + str(value) + '\n')
    out.close()

    ### step1: 解析获得 unitig flank Hi-C links
    dic_pair_hic = load_pickle_file(full_links)
    out = open('flank.links.txt', 'w')
    for key, value in dic_pair_hic.items():
        out.write(str(key) + '\t' + str(value) + '\n')
        # print(key, value)
    out.close()

    dic_contig_hic = get_link_list(dic_pair_hic)
    dic_overlap_ratios = unitig_overlap_ratio(allelic_table_file)
    g1, g2, g3, g4 = genotype_allelic_table(allelic_table_file,
                                            dic_contig_type,
                                            dic_contig_hic,
                                            dic_overlap_ratios,
                                            RE_site_dict)
    ### 错误分型的 unitig
    g1 = false_genotype_unitig(g1, dic_fasta, dic_contig_hic)
    g2 = false_genotype_unitig(g2, dic_fasta, dic_contig_hic)
    g3 = false_genotype_unitig(g3, dic_fasta, dic_contig_hic)
    g4 = false_genotype_unitig(g4, dic_fasta, dic_contig_hic)
    print(f'Debugs: g1 {g1}')
    print(f'Debugs: g2 {g2}')
    print(f'Debugs: g3 {g3}')
    print(f'Debugs: g4 {g4}')
    list_to_file(g1, f'{wd}/g1.txt')
    list_to_file(g2, f'{wd}/g2.txt')
    list_to_file(g3, f'{wd}/g3.txt')
    list_to_file(g4, f'{wd}/g4.txt')
    if os.path.exists('group.cluster.txt'):
        os.remove('group.cluster.txt')
    list_to_cluster_file(g1, 'group1')
    list_to_cluster_file(g2, 'group2')
    list_to_cluster_file(g3, 'group3')
    list_to_cluster_file(g4, 'group4')
    output_fa_from_list(g1, dic_fasta, f'{wd}/g1.fa')
    output_fa_from_list(g2, dic_fasta, f'{wd}/g2.fa')
    output_fa_from_list(g3, dic_fasta, f'{wd}/g3.fa')
    output_fa_from_list(g4, dic_fasta, f'{wd}/g4.fa')
    # note: 这里的g1, g2, g3, g4已经更改了，末尾添加了' '
    write_combined_genotypes_to_file(g1, g2, g3, g4, f'{wd}/g1g2g3g4.txt')
    ## merge g1+g2+g3+g4
    if os.path.exists('g1g2g3g4.fa'):
        os.remove('g1g2g3g4.fa')
    rename_id(f'{wd}/g1.fa', f'{wd}/g1g2g3g4.fa')
    rename_id(f'{wd}/g2.fa', f'{wd}/g1g2g3g4.fa')
    rename_id(f'{wd}/g3.fa', f'{wd}/g1g2g3g4.fa')
    rename_id(f'{wd}/g4.fa', f'{wd}/g1g2g3g4.fa')


def parse_arguments():
    parser = argparse.ArgumentParser("Haplotype clustering of autopolyploid genome.")

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

    ### Step 1: p_utg vs mT2T
    step1_dir = os.path.join(cwd, '02.cluster', '01.putg_vs_mT2T')
    os.makedirs(step1_dir, exist_ok=True)

    paf_file = os.path.join(step1_dir, 'p_utg_vs_mT2T.paf')
    sorted_paf_file = os.path.join(step1_dir, 'p_utg_vs_mT2T.sort.paf')
    best_paf_file = os.path.join(step1_dir, 'putg_vs_mT2T.best.paf')
    allelic_table_file = os.path.join(step1_dir, 'allelic.ctg.table')
    allelic_table_sorted = os.path.join(step1_dir, 'allelic.ctg.table.sort')
    top_contigs_file = os.path.join(step1_dir, 'top_contigs_per_bin.txt')

    try:
        # Minimap2 alignment
        if not os.path.exists(paf_file):
            subprocess.run( ' '.join([
                'minimap2', '-cx', 'asm5', '-t', str(args.threads), args.mT2T, args.p_utg,
                '>', paf_file
            ]), shell=True, check=True)

        # Sort PAF file
        subprocess.run(' '.join([
            'sort', '-k1,1', '-k6,6', '-k8,8n', paf_file,
            '>', sorted_paf_file
        ]), shell=True, check=True)

        # Find the longest subsequence
        subprocess.run(' '.join([
            'nohup', 'time', '-v', os.path.join(utils_realpath, 'find_longest_subsequence.py'),
            '--paf', sorted_paf_file,
            '--min_alignment_distance', str(args.min_alignment_distance),
            '--min_align_length', str(args.min_align_length),
            '--min_unitig_length', str(args.min_unitig_length),
            '--best_lis_output', best_paf_file,
            '--min_match_ratio', str(args.min_match_ratio),
            '--min_lis_size', str(args.min_lis_size),
            '--max_lis_distance', str(args.max_lis_distance)
        ]), shell=True, check=True)

        # Generate allelic table
        subprocess.run(' '.join([
            'nohup', 'time', '-v', os.path.join(utils_realpath, 'allelic_table_generate.py'),
            '--paf_file', best_paf_file,
            '--bin_size', str(args.bin_size),
            '--chr_num', str(args.chr_num),
            '--contig_type', args.contig_type,
            '--out_top_contigs_per_bin', top_contigs_file,
            '--out_allelic_table', allelic_table_file
        ]), shell=True, check=True)

        # Sort allelic table
        subprocess.run(' '.join([
            'sort', '-k1,1', '-k2,2n', allelic_table_file,
            '>', allelic_table_sorted
        ]), shell=True, check=True)

        # Refresh allelic table
        subprocess.run(' '.join([
            'nohup', 'time', '-v', os.path.join(utils_realpath, 'allelic_table_refresh.py'),
            '--wd', step1_dir,
            '--allelic_table', allelic_table_sorted,
            '--fasta', args.p_utg,
            '--contig_type', args.contig_type,
            '--top_n', str(args.top_n),
            '--search_range', str(args.search_range)
        ]), shell=True, check=True)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    ### Step 2: Extract chromosome sequences from p_utg
    step2_dir = os.path.join(cwd, '02.cluster', '02.chr_seq')
    os.makedirs(step2_dir, exist_ok=True)
    try:
        cmd = ' '.join([
            os.path.join(utils_realpath, 'extract_chr_from_putg.py'),
            '--p_utg', args.p_utg,
            '--paf', best_paf_file,
            '--wd', step2_dir,
            '--chr_num', str(args.chr_num)
        ])
        subprocess.run(cmd, shell=True, check=True)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    ### step3: Cluster based on hi-c links
    step3_dir = os.path.join(cwd, '02.cluster', '03.cluster')
    os.makedirs(step3_dir, exist_ok=True)
    # 待聚类的染色体文件列表
    file_list = glob.glob(os.path.join(step2_dir, '*putg.fa'))
    print(f'Debugs: file_list {file_list}')

    for file in file_list:
        file_bn = os.path.basename(file)
        chr = file_bn.replace('.putg.fa', '')
        # 聚类
        cluster_chr_dir = os.path.join(step3_dir, chr)
        os.makedirs(cluster_chr_dir, exist_ok=True)
        os.chdir(cluster_chr_dir)
        cmd = f'grep {chr} {cwd}/02.cluster/01.putg_vs_mT2T/corrected_allelic_table.txt > {chr}.corrected_allelic_table.txt'
        subprocess.run(cmd, shell=True, check=True)
        allelic_table_chr = f'{chr}.corrected_allelic_table.txt'
        cluster(file, args.full_links, args.flank, args.contig_type, allelic_table_chr, cluster_chr_dir)
        # cluster(args.p_utg, args.full_links, args.flank, args.contig_type, allelic_table_chr, cluster_chr_dir)
    os.chdir(cwd)

    ### step4: Re-cluster unclustered unitig
    commands = []
    step4_dir = os.path.join(cwd, '02.cluster', '04.recluster')
    os.makedirs(step4_dir, exist_ok=True)
    for file in file_list:
        file_bn = os.path.basename(file)
        chr = file_bn.replace('.putg.fa', '')
        recluster_chr_dir = os.path.join(step4_dir, chr)
        os.makedirs(recluster_chr_dir, exist_ok=True)
        # os.chdir(recluster_chr_dir)
        cmd = (f'cd {recluster_chr_dir}; {utils_realpath}/chr_uncluster_recluster.py --fasta {file} --contig_type {args.contig_type} '
               f'--full_links {args.full_links} --clusters_file {step3_dir}/{chr}/group.cluster.txt --clm_file {args.clm} > log_re_out 2> log_re_err')
        commands.append(cmd)
    run_in_parallel(commands, 12)
    os.chdir(cwd)

    ### step5: rescue: 与 mT2T 未比对的 unitigs，被拯救
    step5_dir = os.path.join(cwd, '02.cluster', '05.rescue')
    os.makedirs(step5_dir, exist_ok=True)
    os.chdir(step5_dir)
    # 生成总的 merge.group.reassignment.cluster.txt
    recluster_cluster_file_list = glob.glob(os.path.join(step4_dir, 'chr*', 'group.reassignment.cluster.txt'))
    merged_cluster_file = os.path.join(step5_dir, 'merge.group.reassignment.cluster.txt')
    # 使用 cat 命令合并文件
    cmd_cat = 'cat ' + ' '.join(recluster_cluster_file_list) + ' > ' + merged_cluster_file
    subprocess.run(cmd_cat, shell=True, check=True)
    # 生成总的 un_chr.fa
    file_un_chr_fasta = os.path.join(step2_dir, 'un_chr.fa')
    # rescue
    cmd_rescue = (f'nohup time -v {utils_realpath}/unchr_recluster.py --draft_fasta {args.p_utg} --fasta {file_un_chr_fasta} '
                 f'--contig_type {args.contig_type} --full_links {args.full_links} --clusters_file {merged_cluster_file} --clm_file {args.clm} '
                 f'> log_rescue_out 2> log_rescue_err')
    subprocess.run(cmd_rescue, shell=True, check=True)


if __name__ == '__main__':
    main()