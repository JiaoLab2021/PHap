#!/usr/bin/env python
'''
time: 2024-06-22
author: pxxiao
version: 1.0
----
description: 对 un_grouped unitig reassignment

To do list:
    1. hi-c links density
'''


import argparse
import pickle
from util import *
from collections import defaultdict


def load_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data


def parse_contig_type(file_contig_type):
    df = pd.read_csv(file_contig_type, delim_whitespace=True)
    dic_contig_type = dict(zip(df['contig_ID'], df['contig_type']))
    return dic_contig_type


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


def get_RE_site_dict_full_length(fasta_file, RE):
    ''' RE_site_dict '''
    ## dic_fa
    fa_dict = parse_fasta(fasta_file, RE)
    RE_site_dict = {}
    for ctg, (seq, ctg_len, RE_sites) in fa_dict.items():
        RE_site_dict[ctg] = RE_sites
    return RE_site_dict


def parse_pickle(fa_dict, pickle_file):
    ''' 从 pickle_file 获得 full_links '''
    with open(pickle_file, 'rb') as f:
        full_link_dict = pickle.load(f)

    # sort by contig length
    sorted_ctg_list = sorted([(ctg, fa_dict[ctg][1]) for ctg in fa_dict], key=lambda x: x[1], reverse=True)
    # RE site
    RE_site_dict = {ctg: ctg_info[2] for ctg, ctg_info in fa_dict.items()}

    return full_link_dict, sorted_ctg_list, RE_site_dict


def parse_clusters(clusters_file, RE_site_dict):

    ctg_group_dict = dict()
    group_RE_dict = dict()

    with open(clusters_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            group = cols[0]
            # add pseudo-count of 1
            if group not in group_RE_dict:
                group_RE_dict[group] = 1  # Initialize with pseudo-count
            for ctg in cols[2:]:
                if ctg not in ctg_group_dict:
                    ctg_group_dict[ctg] = set()
                ctg_group_dict[ctg].add(group)
                # minus the pseudo-count of 1 for each contig
                group_RE_dict[group] += RE_site_dict.get(ctg, 1) - 1  # Use get() to handle missing keys gracefully
            # add pseudo-count of 1
            # group_RE_dict[group] = 1
            # for ctg in cols[2:]:
            #     ctg_group_dict[ctg] = group
            #     minus the pseudo-count of 1 for each contig
            #     group_RE_dict[group] += RE_site_dict[ctg] - 1

    return ctg_group_dict, group_RE_dict


def parse_link_dict(link_dict, ctg_group_dict, normalize_by_nlinks=False):

    def add_ctg_group(ctg, groups, links):
        for group in groups:
            if group != 'ungrouped':
                if group in ctg_group_link_dict[ctg]:
                    ctg_group_link_dict[ctg][group] += links
                else:
                    ctg_group_link_dict[ctg][group] = links

    ctg_group_link_dict = defaultdict(lambda: defaultdict(int))         # 这个 contig 对应所有 group 的 links 数量：{ctg: {group1: 10, group2: 100, group3: 10000, group4: 200}, ...}
    linked_ctg_dict = defaultdict(set)
    ctg_link_dict = defaultdict(int)

    # if normalize_by_nlinks:         # 归一化：计算每个 contig 的总链接数 以及 链接的总和
    #     total_links = 0
    #     normalized_total_links = 0
    #     for (ctg_i, ctg_j), links in link_dict.items():
    #         ctg_link_dict[ctg_i] += links
    #         ctg_link_dict[ctg_j] += links
    #         total_links += links
    if normalize_by_nlinks:
        total_links = sum(link_dict.values())
        for (ctg_i, ctg_j), links in link_dict.items():
            ctg_link_dict[ctg_i] += links
            ctg_link_dict[ctg_j] += links
    for (ctg_i, ctg_j), links in link_dict.items():
        if normalize_by_nlinks:
            links /= (ctg_link_dict[ctg_i] * ctg_link_dict[ctg_j]) ** 0.5
            link_dict[(ctg_i, ctg_j)] = links

        groups_i = ctg_group_dict.get(ctg_i, set())
        groups_j = ctg_group_dict.get(ctg_j, set())
        add_ctg_group(ctg_i, groups_j, links)
        add_ctg_group(ctg_j, groups_i, links)
        linked_ctg_dict[ctg_i].add(ctg_j)
        linked_ctg_dict[ctg_j].add(ctg_i)

    if normalize_by_nlinks:
        scale_factor = total_links / sum(link_dict.values())
        for (ctg_i, ctg_j), links in link_dict.items():
            link_dict[(ctg_i, ctg_j)] *= scale_factor
            add_ctg_group(ctg_i, ctg_group_dict.get(ctg_j, set()), links)
            add_ctg_group(ctg_j, ctg_group_dict.get(ctg_i, set()), links)

    # for (ctg_i, ctg_j), links in link_dict.items():
    #     if ctg_i in ctg_group_dict and ctg_j in ctg_group_dict:
    #         if normalize_by_nlinks:
                # normalize by geometric mean
                # links /= (ctg_link_dict[ctg_i] * ctg_link_dict[ctg_j]) ** 0.5
                # link_dict[(ctg_i, ctg_j)] = links
                # normalized_total_links += links
            # else:
            #     group_i, group_j = ctg_group_dict[ctg_i], ctg_group_dict[ctg_j]
            #     add_ctg_group(ctg_i, group_j, links)
            #     add_ctg_group(ctg_j, group_i, links)
            #     linked_ctg_dict[ctg_i].add(ctg_j)
            #     linked_ctg_dict[ctg_j].add(ctg_i)

    # make sure the total links are the same
    # if normalize_by_nlinks:
    #     scale_factor = total_links / normalized_total_links
    #     for (ctg_i, ctg_j), links in link_dict.items():
    #         links *= scale_factor
    #         link_dict[(ctg_i, ctg_j)] = links
    #         group_i, group_j = ctg_group_dict[ctg_i], ctg_group_dict[ctg_j]
    #         add_ctg_group(ctg_i, group_j, links)
    #         add_ctg_group(ctg_j, group_i, links)
    #         linked_ctg_dict[ctg_i].add(ctg_j)
    #         linked_ctg_dict[ctg_j].add(ctg_i)

    return ctg_group_link_dict, linked_ctg_dict


def run_reassignment(sorted_ctg_list, ctg_group_link_dict, ctg_group_dict, full_link_dict, linked_ctg_dict, fa_dict, RE_site_dict, group_RE_dict,
                     min_RE_sites=1, min_links=0, nround=1, ambiguous_cutoff=0.6, min_link_density=0.0001):

    ''' 对 unitig 进行重新分配 '''

    def update(ctg, max_group):

        # reassign new group (max_group) to the ctg
        former_group = ctg_group_dict[ctg]
        ctg_group_dict[ctg] = max_group
        # update ctg_group_link_dict
        for each_ctg in linked_ctg_dict[ctg]:
            group_links = ctg_group_link_dict[each_ctg]
            ctg_name_pair = tuple(sorted([ctg, each_ctg]))
            links = full_link_dict[ctg_name_pair]
            if former_group != 'ungrouped':
                group_links[former_group] -= links
            if max_group in group_links:
                group_links[max_group] += links
            elif max_group != 'ungrouped':
                group_links[max_group] = links

    def cal_link_density(ctg, group, former_group, links):

        group_RE_sites = group_RE_dict[group]
        if group == former_group:
            return links / group_RE_sites
        else:
            return links / (group_RE_sites + RE_site_dict[ctg] - 1)

    result_dict = defaultdict(int)
    ## 遍历 contig list
    for ctg, ctg_len in sorted_ctg_list:
        former_group = ctg_group_dict.get(ctg, 'ungrouped')
        group_links = ctg_group_link_dict[ctg]                  # 用来验证已经聚类的 contig 是不是真正属于这个 group
        print(f'Debugs group_links {group_links}')

        # filter RE_site
        if (RE_site_dict[ctg] - 1 < min_RE_sites ) or not group_links:
            max_group, max_links = 'ungrouped', 0
            result_dict['not_rescued'] += 1
            filtered = True
        else:
            sorted_group_links = sorted(group_links.items(), key=lambda x: x[1], reverse=True)
            max_group, max_links = sorted_group_links[0]
            if len(sorted_group_links) > 1:
                second_links = sorted_group_links[1][1]
            else:
                second_links = 0
            # filter max_links
            if max_links < min_links :
                max_group, max_links = 'ungrouped', 0
                result_dict['not_rescued'] += 1
                filtered = True
            # we define contigs with second_links / max_links >= ambiguous_cutoff as ambiguous,
            # these contigs will NOT be reassigned, but only rescued in the additional round of rescue
            elif nround and second_links / max_links >= ambiguous_cutoff:
                max_group, max_links = 'ungrouped', 0
                result_dict['not_rescued'] += 1
                filtered = True
            else:
                # calculate link densities
                max_link_density = cal_link_density(ctg, max_group, former_group, max_links)
                # filter max_link_density
                if max_link_density < min_link_density:
                    max_group, max_links = 'ungrouped', 0
                    result_dict['not_rescued'] += 1
                    filtered = True
                else:
                    # calculate the average link density to other nonbest groups
                    # other_group_density_sum = sum([cal_link_density(ctg, group, former_group, links) for group, links in sorted_group_links[1:]])

                    # a patch for gfa file

                    other_group_density_sum = sum(
                        [cal_link_density(ctg, group, former_group, links) for group, links in
                         sorted_group_links[1:]])

                    if other_group_density_sum:
                        avg_other_group_density = other_group_density_sum / (len(group_RE_dict) - 1)
                    else:
                        # assign a big number to prevent division by zero problem
                        avg_other_group_density = 1000000000


def assign_unitigs_by_density(sorted_ctg_list, ctg_group_dict, dic_contig_type, full_link_dict, RE_site_dict, group_RE_dict):
    ''' 分配 unitigs 到 groups，基于 Hi-C links density '''

    list_unreassigned_unitig = []

    ## 遍历 unitigs
    for ctg, ctg_len in sorted_ctg_list:
        # 存储每个 group 对每个 contig 的链接密度
        group_density_dict = defaultdict(float)

        # 跳过 grouped unitig
        if ctg in ctg_group_dict.keys():
            print(f'Debugs: {ctg} has grouped, pass!')
            continue
        print(f'Debugs: ctg ctg_len {ctg} {ctg_len}')

        # 获得当前 unitig type
        contig_type = dic_contig_type.get(ctg, 'haplotig')
        num_groups_to_assign = {
            'haplotig': 1,
            'diplotig': 2,
            'triplotig': 3,
            'tetraplotig': 4
        }.get(contig_type, 1)

        links_list = []

        # 计算与所有已分组 unitig 的 Hi-C links density
        for (ctg1, ctg2), links in full_link_dict.items():
            if ctg1 == ctg or ctg2 == ctg:
                other_ctg = ctg2 if ctg1 == ctg else ctg1
                if other_ctg in ctg_group_dict:
                    for assigned_group in ctg_group_dict[other_ctg]:
                        links_list.append(links)
                        group_RE_sites = group_RE_dict[assigned_group]
                        density = links / group_RE_sites if group_RE_sites else 0
                        print(f'Debugs: new_unitig assigned_group density {ctg}-{assigned_group} {links} {group_RE_sites} {density}')
                        group_density_dict[assigned_group] += density

        if len(links_list) > 0 and max(links_list) < 5:     # 对 Hi-C links 数量进行过滤，如果 Hi-C links 的数量小于 5，那么可以认为是没有 Hi-C 信号的，考虑跳过。
            list_unreassigned_unitig.append(ctg)
            continue

        # 选择链接密度最高的几个 group
        sorted_groups_by_density = sorted(group_density_dict.items(), key=lambda x: x[1], reverse=True)
        top_groups = [group for group, _ in sorted_groups_by_density[:num_groups_to_assign]]

        if not top_groups:
            list_unreassigned_unitig.append(ctg)
            print(f'Debugs: No Hi-C signal, {ctg} has not been assigned to any group.')
        else:
            print(f'Debugs: top_groups sorted_groups_by_density {top_groups} {sorted_groups_by_density}')
            # 分配到最适合的 group
            ctg_group_dict[ctg] = set()
            for group in top_groups:
                ctg_group_dict[ctg].add(group)
                print(f'Debugs: {ctg} reassigned to {group}')
                # 更新 group 的酶切位点计数
                group_RE_dict[group] += RE_site_dict[ctg]
            # for group in top_groups:
            #     if ctg not in ctg_group_dict:
            #         ctg_group_dict[ctg] = set()
            #     ctg_group_dict[ctg].add(group)  # 假设 ctg_group_dict[ctg] 是一个 set
            #     print(f'Debugs: {ctg} reassigned to {group}')
            #     更新 group 的酶切位点计数
                # group_RE_dict[group] += RE_site_dict[ctg]  # 假设我们有每个 ctg 的 RE_sites_dict
    return ctg_group_dict, list_unreassigned_unitig


def write_group_ids_and_sequences(ctg_group_dict, list_unreassigned_unitig, fa_dict):
    ''' 将 grouped unitig ID and sequences 写入文件 '''

    # 初始化存储每个 group 的 unitig ID 和序列
    group_ids = {'group1': [], 'group2': [], 'group3': [], 'group4': []}
    group_seqs = {'group1': [], 'group2': [], 'group3': [], 'group4': []}
    combined_seqs = []  # 用于存储所有group的序列，包含group ID前缀

    # 按 group 分类 unitig IDs 和 sequences
    for ctg, groups in ctg_group_dict.items():
        print(ctg, groups)
        for group in groups:
            group_suffix = group[-1]        # 最后一个是编号，1，2，3，4
            group_ids[f'group{group_suffix}'].append(ctg)
            group_seqs[f'group{group_suffix}'].append(f'>{ctg}\n{fa_dict[ctg][0]}')
            combined_seqs.append(f'>g{group_suffix}_{ctg}\n{fa_dict[ctg][0]}')
            print(f'Debugs: {ctg} is added to group g{group_suffix}')

    # 写入每个 group 的 ID 和 sequences 到文件
    for group in group_ids:
        with open(f'{group}.reassignment.txt', 'w') as f:
            f.write('\n'.join(group_ids[group]) + '\n')

        with open(f'{group}.reassignment.fa', 'w') as f:
            f.write('\n'.join(group_seqs[group]) + '\n')

        with open(f'{group}.txt', 'w') as f:
            f.write('#Contig' + '\t' + 'RECounts' + '\t' + 'Length' + '\n')
            for ctg in group_ids[group]:
                f.write(ctg + '\t' + str(fa_dict[ctg][2]) + '\t' + str(fa_dict[ctg][1]) + '\n')

    # reassignment group cluster: group.reassignment.cluster.txt
    with open('group.reassignment.cluster.txt', 'w') as f:
        for group in group_ids:
            f.write(group + '\t' + ' '.join(group_ids[group]) + '\n')

    # 将所有序列合并到一个文件
    with open('g1g2g3g4.reassignment.fa', 'w') as f:
        f.write('\n'.join(combined_seqs) + '\n')

    # 处理未分配的 unitig
    if list_unreassigned_unitig:
        with open('unassigned_unitigs.txt', 'w') as f:
            f.write('\n'.join(list_unreassigned_unitig) + '\n')
            print(f'Debugs: Unassigned unitigs written to unassigned_unitigs.txt')


def split_clm_file(clm_file, group_ctg_dict, ctg_group_dict):

    # logger.info('Splitting clm file into subfiles by group...')

    # make directory for final groups
    # final_dir = 'final_groups'
    # os.mkdir(final_dir)

    # create symbolic links for final groups
    # if subdir == 'reassigned_groups':
    #     prefix = 'reassigned'
    # else:
    #     assert subdir == 'hc_groups'
    #     prefix = 'hc'

    # for group in group_ctg_dict:
    #     # group files
    #     os.symlink('../{0}/{1}_{2}.txt'.format(subdir, prefix, group),
    #                '{0}/{1}.txt'.format(final_dir, group))
    #
    # # clusters file
    # os.symlink('../{0}/{1}_clusters.txt'.format(subdir, prefix),
    #            '{0}/final_clusters.txt'.format(final_dir))

    # make directory for clm splitting
    subdir = 'split_clms'
    os.makedirs(subdir, exist_ok=True)

    fp_dict = dict()

    for group in group_ctg_dict:
        fp_dict[group] = open('{}/{}.clm'.format(subdir, group), 'w')

    with open(clm_file) as f:
        for line in f:
            cols = line.split()
            ctg_1, ctg_2 = cols[0][:-1], cols[1][:-1]
            if ctg_1 in ctg_group_dict and ctg_2 in ctg_group_dict:
                common_groups = ctg_group_dict[ctg_1].intersection(ctg_group_dict[ctg_2])
                if common_groups:  # 检查是否存在共同的组
                    for group in common_groups:  # 为每个共同的组写入
                        fp_dict[group].write(line)

    for group, fp in fp_dict.items():
        fp.close()

    # fp_dict = dict()
    #
    # for group in group_ctg_dict:
    #     fp_dict[group] = open('{}/{}.clm'.format(subdir, group), 'w')
    #
    # with open(clm_file) as f:
    #     for line in f:
    #         cols = line.split()
    #         ctg_1, ctg_2 = cols[0][:-1], cols[1][:-1]
    #         if ctg_1 == ctg_2: continue
    #         print(f'Debugs: {ctg_1} {ctg_2}')
    #         if ctg_1 in ctg_group_dict and ctg_2 in ctg_group_dict and ctg_group_dict[ctg_1] == ctg_group_dict[ctg_2]:
    #             print(f'{ctg_group_dict[ctg_1]}')
    #             for i in fp_dict[ctg_group_dict[ctg_1]]:
    #                 fp_dict[i].write(line)
    #
    # for group, fp in fp_dict.items():
    #     fp.close()


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True, help='fasta file')
    parser.add_argument('--contig_type', required=True, help='contig type based on dosage analysis')
    parser.add_argument('--RE', required=False, default='GATC', help='RE')
    parser.add_argument('--full_links', required=True, help='full links from haphic cluster')
    parser.add_argument('--clusters_file', required=True, help='the cluster file')
    parser.add_argument('--clm_file', type=str, required=True, help='clm file')
    parser.add_argument('--normalize_by_nlinks', required=False, action='store_true',
                        help='normalize inter-contig and inter-group Hi-C links by the number of links to other contigs or groups, default: %(default)s')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()


    def reassign(fasta, contig_type, full_links, clusters_file, normalize_by_nlinks, clm_file):
        ### step1: hi-c links density
        ## fa dict
        fa_dict = parse_fasta(fasta)

        ## contig type
        dic_contig_type = parse_contig_type(contig_type)
        print(f'Debugs: dic_contig_type {dic_contig_type}')

        ## RE site
        full_link_dict, sorted_ctg_list, RE_site_dict = parse_pickle(fa_dict, full_links)
        print(f'Debugs: full_link_dict {full_link_dict}')
        print(f'Debugs: sorted_ctg_list {len(sorted_ctg_list)} {sorted_ctg_list}')
        print(f'Debugs: RE_site_dict {RE_site_dict}')
        out = open('full.links.txt', 'w')
        for key, value in full_link_dict.items():
            out.write(str(key) + '\t' + str(value) + '\n')
        out.close()

        ## parse_clusters
        ctg_group_dict, group_RE_dict = parse_clusters(clusters_file, RE_site_dict)
        print(f'Debugs: ctg_group_dict {ctg_group_dict} ')
        print(f'Debugs: group_RE_dict {group_RE_dict}')

        ## cal link density
        ctg_group_link_dict, linked_ctg_dict = parse_link_dict(full_link_dict, ctg_group_dict,
                                                               normalize_by_nlinks=normalize_by_nlinks)
        print(f'Debugs ctg_group_link_dict {ctg_group_link_dict}')
        print(f'Debugs link_dict {linked_ctg_dict}')

        ctg_group_dict_new, list_unreassigned_unitig = assign_unitigs_by_density(sorted_ctg_list, ctg_group_dict,
                                                                                 dic_contig_type, full_link_dict,
                                                                                 RE_site_dict, group_RE_dict)
        write_group_ids_and_sequences(ctg_group_dict_new, list_unreassigned_unitig, fa_dict)

        ## split clm
        split_clm_file(clm_file, group_RE_dict, ctg_group_dict)


    ### step1: hi-c links density
    ## fa dict
    fa_dict = parse_fasta(args.fasta)

    ## contig type
    dic_contig_type = parse_contig_type(args.contig_type)
    print(f'Debugs: dic_contig_type {dic_contig_type}')

    ## RE site
    full_link_dict, sorted_ctg_list, RE_site_dict = parse_pickle(fa_dict, args.full_links)
    print(f'Debugs: full_link_dict {full_link_dict}')
    print(f'Debugs: sorted_ctg_list {len(sorted_ctg_list)} {sorted_ctg_list}')
    print(f'Debugs: RE_site_dict {RE_site_dict}')
    out = open('full.links.txt', 'w')
    for key, value in full_link_dict.items():
        out.write(str(key) + '\t' + str(value) + '\n')
    out.close()

    ## parse_clusters
    ctg_group_dict, group_RE_dict = parse_clusters(args.clusters_file, RE_site_dict)
    print(f'Debugs: ctg_group_dict {ctg_group_dict} ')
    print(f'Debugs: group_RE_dict {group_RE_dict}')

    ## cal link density
    ctg_group_link_dict, linked_ctg_dict = parse_link_dict(full_link_dict, ctg_group_dict, normalize_by_nlinks=args.normalize_by_nlinks)
    print(f'Debugs ctg_group_link_dict {ctg_group_link_dict}')
    print(f'Debugs link_dict {linked_ctg_dict}')

    ctg_group_dict_new, list_unreassigned_unitig = assign_unitigs_by_density(sorted_ctg_list, ctg_group_dict, dic_contig_type, full_link_dict, RE_site_dict, group_RE_dict)
    write_group_ids_and_sequences(ctg_group_dict_new, list_unreassigned_unitig, fa_dict)

    ## split clm
    split_clm_file(args.clm_file, group_RE_dict, ctg_group_dict)

    # run_reassignment(sorted_ctg_list, ctg_group_link_dict, ctg_group_dict, full_link_dict, linked_ctg_dict, fa_dict, RE_site_dict, group_RE_dict)


if __name__ == '__main__':
    main()