#!/usr/bin/env python
'''
time: 2024-05-23
author: pxxiao
version: 1.0
----
description: 找到 p_utg 与 mT2T 之间的最佳匹配
'''


import argparse
import re
from collections import defaultdict
import os
import subprocess
import pandas as pd
import pysam
import pickle


def fasta_read(file_input):
    '''
    解析 FASTA 文件，存入字典 d
    :param file_input: .fasta
    :return: d: id as key, seq as value
    '''
    d = {}
    for lines in open(file_input, "r"):
        if lines.startswith(">"):
            id = lines.strip().replace(">", "")
            d[id] = []
        else:
            d[id].append(lines.strip().upper())
    for key, values in d.items():
        d[key] = "".join(values)
    return d


def genome_aln(file_ref, file_qry, wd, args):
    ''' p_utg vs mT2T '''
    if not os.path.exists(f'{wd}/putg_vs_mT2T.paf'):
        cmd_minimap2 = (f'minimap2 -cx asm5 -t {args.threads} {file_ref} {file_qry} > {wd}/putg_vs_mT2T.paf')
        subprocess.run(cmd_minimap2, shell=True, close_fds=True)


def contig_pair_hic_links(bam_file_name, prefix, dic_fasta):

    def get_contig_pair(contig1_id, contig2_id, len_contig1, len_contig2):
        '''Returns contig pair based on contig length'''
        return (contig1_id, contig2_id) if len_contig1 > len_contig2 else (contig2_id, contig1_id)

    def output_pickle(dict, file):
        with open(file, 'wb') as fpkl:
            pickle.dump(dict, fpkl)

    if os.path.exists(f'{prefix}.contig.pair.links.pkl'):
        return None

    dic_full_links = defaultdict(int)

    ### output file
    out_inter_contig = open(f'{prefix}.inter-contig.txt', 'w')

    ### 是否处理过这个 read-pair
    processed_read_ids = set()

    ### 处理 bam 文件
    bam = pysam.AlignmentFile(bam_file_name, 'rb')

    ### 定义染色体名称和索引的映射字典
    chromosome_names = dict(enumerate(bam.references))

    for read in bam.fetch():
        # 检查当前 read 是否处理过
        if read.query_name not in processed_read_ids:
            # reference id
            ref1 = chromosome_names.get(read.reference_id, 'Unknown')
            ref2 = chromosome_names.get(read.next_reference_id, 'Unknown')

            # 检查当前 read 和 它的 mate 是否比对到同一条染色体上
            if ref2 == ref1:     # 同一条染色体：intral-contig
                continue

            elif ref2 != ref1:       # 不同一条染色体：inter-contig
                # inter-contig
                out_inter_contig.write(f'{read.query_name}\t{ref1}\t{ref2}\n')
                contig_pair = get_contig_pair(ref1, ref2, len(dic_fasta[ref1]), len(dic_fasta[ref2]))
                dic_full_links[contig_pair] += 1

            else:
                continue
            processed_read_ids.add(read.query_name)
    bam.close()
    out_inter_contig.close()
    output_pickle(dic_full_links, f'{prefix}.contig.pair.links.pkl')
    return None


def calculate_contig_match_ratio(file_paf, wd, args):

    def convert_to_matrix(contig_match_ratios, wd):

        def export_matrix_to_csv(matrix, output_file):
            matrix.to_csv(output_file, sep='\t', index=True)

        def export_matrix_to_excel(matrix, output_file):
            matrix.to_excel(output_file, index=True)

        # Prepare data for DataFrame
        data = []
        putg_ids = set()
        scaffold_ids = set()

        for putg_id, scaffolds in contig_match_ratios.items():
            putg_ids.add(putg_id)
            for scaffold_id, values in scaffolds.items():
                scaffold_ids.add(scaffold_id)
                data.append([putg_id, scaffold_id, values[2]])  # values[2] is the match ratio
        # Create DataFrame
        df = pd.DataFrame(data, columns=['putg_id', 'scaffold_id', 'match_ratio'])
        # Pivot DataFrame to create matrix
        matrix = df.pivot(index='putg_id', columns='scaffold_id', values='match_ratio').fillna(0)
        # output
        export_matrix_to_csv(matrix, f'{wd}/contig_match_ratios.csv')
        export_matrix_to_excel(matrix, f'{wd}/contig_match_ratios.xlsx')
        return matrix


    def find_max_match_ratios(contig_match_ratios, scffolds_chr, wd):
        ''' 找到每个 contig 最佳的匹配 '''
        max_match_ratios = {}
        out = open(os.path.join(wd+'/putg_vs_mT2T.best.match.txt'), 'w')
        out.write('#contig_id\tscaffold_id\n')
        for contig_id, scaffolds in contig_match_ratios.items():
            max_scaffold = None
            max_ratio = 0
            for scaffold_id, values in scaffolds.items():
                if values[2] > max_ratio:
                    max_ratio = values[2]
                    max_scaffold = scaffold_id
            max_match_ratios[contig_id] = (max_scaffold, max_ratio)

        ## 去除染色体之外的匹配并更新max_match_ratios
        dic_putg_match_mT2T = {}
        dic_chr_posses = defaultdict(list)
        for contig_id, scaffold_id_max_ratio in max_match_ratios.items():
            max_scaffold, max_ratio = scaffold_id_max_ratio[0], scaffold_id_max_ratio[1]
            if scaffold_id_max_ratio[0] in scffolds_chr:
                out.write(contig_id + '\t' + max_scaffold + '\n')
                dic_putg_match_mT2T[contig_id] = max_scaffold
                dic_chr_posses[max_scaffold].append(contig_id)
        return dic_putg_match_mT2T, dic_chr_posses


    def find_more_than_cutoff(contig_match_ratios, scffolds_chr, wd, ratio_cutoff):
        ''' 找到每个 contig 所有匹配率大于给定阈值的 Scaffold '''
        dic_matches_above_cutoff = defaultdict(list)
        out = open(os.path.join(wd, 'putg_vs_mT2T.match_above_cutoff.txt'), 'w')
        out.write('#contig_id\tscaffold_id\tmatch_ratio\n')

        for contig_id, scaffolds in contig_match_ratios.items():
            # matching_scaffolds = []
            for scaffold_id, values in scaffolds.items():
                match_ratio = values[2]
                if match_ratio > ratio_cutoff and scaffold_id in scffolds_chr:
                    dic_matches_above_cutoff[contig_id].append((scaffold_id, match_ratio))
                    out.write(f'{contig_id}\t{scaffold_id}\t{match_ratio}\n')
        out.close()
        return dic_matches_above_cutoff

    chr_num = args.chr_num
    contig_match_ratios = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    dic_scaffold_length = {}        # 存储 Scaffolds 长度
    cigar_pattern = re.compile(r'(\d+)([A-Z])')

    for lines in open(file_paf, 'r'):
        line = lines.strip().split()

        # 跳过次要比对
        if line[16][-1] == 'S':
            continue
        # match count
        cigar = line[-1][5:]
        matches = cigar_pattern.findall(cigar)
        match_count = sum(int(length) for length, op in matches if op == 'M')

        putg_id = line[0]
        putg_length = int(line[1])
        scaffold_length = int(line[6])
        scaffold_id = line[5]
        # putg_length, match count, match ratio
        contig_match_ratios[putg_id][scaffold_id][0] = putg_length
        contig_match_ratios[putg_id][scaffold_id][1] += match_count
        contig_match_ratios[putg_id][scaffold_id][2] = contig_match_ratios[putg_id][scaffold_id][1] / putg_length
        # scaffold length
        dic_scaffold_length[scaffold_id] = scaffold_length

    ## contig_match_ratios 写入文件
    matrix = convert_to_matrix(contig_match_ratios, wd)

    ## 获得染色体匹配数量的 scaffold id
    scffolds_chrs = sorted(dic_scaffold_length.items(), key=lambda x: x[1], reverse=True)[:chr_num]
    scffolds_chr = [scaffold_id for scaffold_id, length in scffolds_chrs]
    print(scffolds_chr)

    ## putg vs mT2T 最佳匹配
    dic_putg_match_mT2T, dic_chr_posses = find_max_match_ratios(contig_match_ratios, scffolds_chr, wd)
    print(f'Debugs: dic_putg_match_mT2T {dic_putg_match_mT2T}')
    print(f'Debugs: dic_chr_posses {dic_chr_posses}')

    dic_matches_above_cutoff = find_more_than_cutoff(contig_match_ratios, scffolds_chr, wd, 0.2)
    print(f'Debugs: dic_matches_above_cutoff {dic_matches_above_cutoff}')

    return dic_putg_match_mT2T, dic_chr_posses, dic_matches_above_cutoff


def reassignment_contig_to_chromosome(dic_matches_above_cutoff, dic_putg_match_mT2T, prefix, min_links_threshold, top_n=3):

    ## 读取 pickle 文件
    def load_pickle_file(file_path):
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
        return data

    ## 获得与给定 contig links 数量前三个的 contig
    def get_top_hic_links(dic_contig_pair_hic_links, query_contig):
        dic_link_counts = defaultdict(dict)
        for (contig1, contig2), links_n in dic_contig_pair_hic_links.items():
            if contig1 == query_contig:
                dic_link_counts[query_contig][contig2] = links_n
            elif contig2 == query_contig:
                dic_link_counts[query_contig][contig1] = links_n
            else:
                continue
        # 获取 query_contig 对应的链接字典
        link_counts = dic_link_counts[query_contig]
        # 对链接数量进行排序并取前 top_n 个
        sorted_links = sorted(link_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
        # print(f'DEBUGS link_counts', sorted_links)  # [('utg000107l', 1), ('utg000653l', 1)]
        return sorted_links

    ## Hi-C links 信息
    dic_contig_pair_links = load_pickle_file(f'{prefix}.contig.pair.links.pkl')
    print(f'Debugs: dic_contig_pair_links {dic_contig_pair_links}')

    # for key, values in dic_contig_pair_links.items():
    #     print(key, values)

    contig_assignments = {}
    for contig_id, matches in dic_matches_above_cutoff.items():
        ## 利用 mT2T 比对，不能确定 contig 的归属
        if len(matches) > 1:
            print(contig_id, matches)
            # 与不确定 contig，Hi-C links 数量最多的前 top_n 个 contig
            sorted_links = get_top_hic_links(dic_contig_pair_links, contig_id)
            print(f'Debugs -- mT2T:', contig_id, 'Hi-C links', sorted_links)
            if len(sorted_links) == 0:      # 利用 Hi-C 数据，不能确定与 不确定 contig 相邻的 contig
                continue
            dic_ = defaultdict(list)
            for i in sorted_links:
                mate_contig_by_hic = i[0]       # hic 支持的 mate contig
                links = i[1]                    # Hi-C reads links
                if links < min_links_threshold:
                    continue
                mate_contig_belong_by_mT2T = dic_putg_match_mT2T[mate_contig_by_hic]        # hic 支持的 mate contig 属于哪个染色体（mT2T 证据）
                # mT2T比对鉴定到比对到多个染色体的contig，利用Hi-C数据这个contig的邻近contig，邻近contig的mT2T归属染色体
                print(f'Debugs -- mT2T:', contig_id, mate_contig_by_hic, mate_contig_belong_by_mT2T)
                dic_[contig_id].append((mate_contig_by_hic, links, mate_contig_belong_by_mT2T))

            # 候选的染色体如果都是同一个，那么就认为这个 contig 属于这个染色体；否则，定一个 links 的截取值，按照 link 的多少，选择最多的那个。
            print(dic_)
            candidate_chromosomes = [value[2] for key, values in dic_.items() for value in values if len(values) > 0]     # 候选染色体
            print('candidate_chromosomes', candidate_chromosomes)       # 候选的染色体
            if len(set(candidate_chromosomes)) == 1:         # 候选的染色体是同一个染色体
                print('Debugs: xpx')
                contig_assignments[contig_id] = candidate_chromosomes[0]
            # elif len(set(candidate_chromosomes)) == 0:      # 没有候选的染色体，也就是说利用 Hi-C 数据，不能确定与 不确定 contig 相邻的 contig

            else:
                top_chromosome, top_links = None, 0
                for contig_id, values in dic_.items():
                    for value in values:
                        if value[1] > top_links and value[1] > min_links_threshold:
                            top_chromosome = value[2]
                            top_links = value[1]
                        else:
                            continue
                if top_chromosome:
                    contig_assignments[contig_id] = top_chromosome
                else:
                    contig_assignments[contig_id] = 'Uncertain'
        else:
            continue
            # contig_assignments[contig_id] = matches[0][0]
    return contig_assignments


def extract_chr_seq(wd, dic_putg, dic_chr_posses):
    for chrom, unitigs in dic_chr_posses.items():
        output_file = os.path.join(wd, f'{chrom}.putg.fa')
        with open(output_file, 'w') as f:
            for unitig in unitigs:
                f.write('>' + unitig + '\n' + dic_putg[unitig] + '\n')
    return None


def extract_unchr_seq(wd, dic_putg, dic_chr_posses):
    list_unitig_on_chr = []
    output_file = os.path.join(wd, f'un_chr.fa')
    for chrom, unitigs in dic_chr_posses.items():
        for unitig in unitigs:
            list_unitig_on_chr.append(unitig)
    with open(output_file, 'w') as f:
        for unitig, seq in dic_putg.items():
            if unitig not in list_unitig_on_chr:
                f.write('>' + unitig + '\n' + seq + '\n')
    return None


def parse_arguments():
    parser = argparse.ArgumentParser('Find best match between mT2T and p_utg')
    parser.add_argument('--wd', required=True, help='Workding directory for this scripts')
    parser.add_argument('--paf', required=True, type=str, help='path to putg vs mT2T paf file')
    parser.add_argument('--chr_num', type=int, default=12, help='number of chromosomes [12]')
    parser.add_argument('--p_utg', required=True, type=str, help='path to putg')
    args = parser.parse_args()
    return args


def main():

    args = parse_arguments()

    # prefix = 'test.links'

    dic_putg = fasta_read(args.p_utg)
    paf_file = args.paf
    wd = args.wd

    ### step1: p_utg vs mT2T
    # os.makedirs('01.putg_vs_mT2T', exist_ok=True)
    # genome_aln(args.mT2T, args.putg, '01.putg_vs_mT2T', args)

    ### step2: get best match between p_utg and mT2T
    dic_putg_match_mT2T, dic_chr_posses, dic_matches_above_cutoff = calculate_contig_match_ratio(paf_file, wd, args)

    ### step3: Output the sequence of each chromosome
    extract_chr_seq(wd, dic_putg, dic_chr_posses)

    ### step4: Output sequences that are not on chromosome
    extract_unchr_seq(wd, dic_putg, dic_chr_posses)

    # ### step3: 利用 Hi-C reads 的信息，确定多个比对的 contig 属于哪个染色体
    # contig_pair_hic_links(args.bam, prefix, dic_putg)
    # contig_assignments = reassignment_contig_to_chromosome(dic_matches_above_cutoff, dic_putg_match_mT2T, prefix, args.min_links_threshold)
    # print(contig_assignments)
    #
    # ## output: mT2T 比对之后，不能确定的染色体归属的 contig，再次利用 Hi-C 数据，能够确定的 contig 输出到一个新的文件当中。
    # output_file = open('01.putg_vs_mT2T/putg_vs_mT2T_Hi-C.reassignments.contigs.txt', 'w')
    # for contig, value in contig_assignments.items():
    #     output_file.write(contig + '\t' + value + '\n')
    #     print(contig, value)
    # output_file.close()


if __name__ == '__main__':
    main()