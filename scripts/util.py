#!/usr/bin/env python

import os
import subprocess
import multiprocessing
import time
import pandas as pd
from collections import defaultdict
import re
from scipy.stats import linregress
from bisect import bisect_left

### 制作文件夹，先判断文件夹是否存在？
def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return None


def execute_command(command):
    '''
    执行命令
    :param command:
    :return:
    '''
    process = subprocess.Popen(command, shell=True)
    process.wait()
    return None


def run_in_parallel(commands, num_processes):
    '''
    多线程运行命令
    :param commands:
    :param num_processes:
    :return:
    '''
    # processes = []
    pool = multiprocessing.Pool(processes=num_processes)
    for command in commands:
        pool.apply_async(execute_command, args=(command, ))
    pool.close()
    pool.join()
    return None


def check_file_in_path(file, cmd):
    if os.path.exists(file):
        print(f'[info] {file} exists, CMD: {cmd}; PASS!')
    else:
        subprocess.run(cmd, shell=True, close_fds=True)
    return None


### 判断文件是否存在，然后执行某个函数
def check_file_in_path_fun(file, fun):
    '''
    文件不存在，执行某个函数；存在则跳过
    :param file:
    :param fun:
    :return:
    '''
    if os.path.exists(file):
        print('File exits')
    else:
        print('File not exits')
        # fun()
    return None


### fasta to dic
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


def fasta_read_filter(file_input, min_len):
    '''
    解析 FASTA 文件，存入字典 d
    过滤长度小于 min_len 的 contig
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
        if len("".join(values)) >= min_len:
            d[key] = "".join(values)
    return d


def delta_filter_same_seq(file_delta, file_filtered_delta):
    '''
    对于short to short比对，需要去除比对到自身的部分
    :param file_delta: 待过滤的 delta 文件
    :param file_filtered_delta: 过滤相同 seq 匹配之后的 delta 文件
    :return:
    '''
    out = open(file_filtered_delta, 'w')
    for lines in open(file_delta, 'r'):
        line = lines.strip().split()
        if lines.startswith('/'):
            out.write(lines)
            continue
        if lines.startswith('NUCMER'):
            out.write(lines)
            continue
        ## >ptg000011l ptg000028l 41964882 2810910
        if lines.startswith('>'):
            if line[0][1:] != line[1]:  # 不是自身比对
                block = 1
                out.write(lines)
            else:
                block = 0               # 自身比对
        ## 具体的比对信息
        if len(line) == 1 or len(line) == 7:
            if block == 0:
                continue
            else:
                out.write(lines)
    return None


def cal_aln_rate(paf_file):
    def interval_merge(list):
        # list: [a,b], [c,d], [e,f]...
        merged = []
        # 对list进行排序
        sorted_list = sorted(list, key=lambda x: x[0])
        for interval in sorted_list:
            if not merged or interval[0] > merged[-1][1]:  # 判断是不是第一个区间，或者当前区间的开始是否大于已有区间的结尾
                merged.append(interval)
            else:  # 当前区间的开始小于已有区间的结尾，进行合并
                merged[-1][1] = max(merged[-1][1], interval[1])
        return merged

    def sum_interval(list):
        su = 0
        for i in list:
            su += (i[1] - i[0])
        return su

    output_file = paf_file + 'aln.ratio.txt'
    output_file = open(output_file, 'w')
    dic_paf = parse_paf(paf_file)
    print('qryID\trefID\tqryLen\trefLen\taln_direction\tratio_qry\tratio_ref\tratio_qry2ref')
    output_file.write('#qryID\trefID\tqryLen\trefLen\taln_direction\tratio_qry\tratio_ref\tratio_qry2ref\n')
    for key, values in dic_paf.items():
        qry_len = values[0][0]
        ref_len = values[1][0]
        # print(qry_len, ref_len)
        ref_list = interval_merge(values[3])
        qry_list = interval_merge(values[2])
        positive_qry_list = interval_merge(values[4])
        negative_qry_list = interval_merge(values[5])
        sum_ref = sum_interval(ref_list)  # 比对上的长度：ref
        sum_qry = sum_interval(qry_list)  # 比对上的长度：qry
        sum_positive_qry = sum_interval(positive_qry_list)  # 比对上的，正向的长度
        sum_negative_qry = sum_interval(negative_qry_list)  # 比对上的，负向的长度

        # 各自比对上的比例
        ratio_ref = sum_ref / ref_len
        ratio_qry = sum_qry / qry_len
        # qry / ref 的比例
        ratio_qry2ref = qry_len / ref_len
        # 比对的方向，positive or negative
        if sum_positive_qry > sum_negative_qry:
            aln_direction = '+'
        else:
            aln_direction = '-'

        ref_ID = key.split('_to_')[1]
        qry_ID = key.split('_to_')[0]
        print(f'{qry_ID}\t{ref_ID}\t{qry_len}\t{ref_len}\t{aln_direction}\t{ratio_qry}\t{ratio_ref}\t{ratio_qry2ref}')
        output_file.write(f'{qry_ID}\t{ref_ID}\t{qry_len}\t{ref_len}\t{aln_direction}\t{ratio_qry}\t{ratio_ref}\t{ratio_qry2ref}\n')
    return dic_paf


def parse_paf(paf_file):
    ''' 解析 paf 文件 '''
    # 字典：存储长度信息
    dic_paf = {}
    for lines in open(paf_file, "r"):
        line = lines.strip().split()
        qry_ID, qry_len, qry_start, qry_end, qry_strand = line[0], int(line[1]), int(line[2]), int(line[3]), line[4]
        tag_ID, tag_len, tag_start, tag_end = line[5], int(line[6]), int(line[7]), int(line[8])
        match_count, block_length, mapping_quality = line[9], line[10], line[11]

        qry_to_tag = qry_ID + "_to_" + tag_ID
        tag_to_qry = tag_ID + "_to_" + qry_ID
        if tag_to_qry not in dic_paf:       # 20240416 fix bugs: ctg1_VS_ctg2 与 ctg2_VS_ctg1只保留一个比对
            if qry_to_tag not in dic_paf.keys():
                # qry len, tag len, qry coordinate interval, tag coordinate interval, sum_qry_forward, sum_qry_reverse
                dic_paf[qry_to_tag] = [[], [], [], [], [], []]
                dic_paf[qry_to_tag][0].append(qry_len)
                dic_paf[qry_to_tag][1].append(tag_len)
                dic_paf[qry_to_tag][2].append([qry_start, qry_end])
                dic_paf[qry_to_tag][3].append([tag_start, tag_end])
            else:
                dic_paf[qry_to_tag][2].append([qry_start, qry_end])
                dic_paf[qry_to_tag][3].append([tag_start, tag_end])
            if qry_strand == '+':
                # 链一致
                dic_paf[qry_to_tag][4].append([qry_start, qry_end])
            elif qry_strand == '-':
                # 链不一致
                dic_paf[qry_to_tag][5].append([qry_start, qry_end])
    print('20240416 fix bugs: 过滤之前，互相比对的染色体对数', len(dic_paf))
    return dic_paf


def parse_delta(file_delta):

    ### delta's basename
    bn_delta = os.path.basename(file_delta).replace('.delta', '')

    ### step1: 过滤 DELTA 文件
    file_filtered_delta = bn_delta + '.filtered.same.seq.delta'
    print('pctg_self_align.filtered.same.seq.delta')
    # delta_filter_same_seq(file_delta, file_filtered_delta)
    check_file_in_path_fun(file_filtered_delta, lambda: delta_filter_same_seq(file_delta, file_filtered_delta))

    ### step2: 将过滤之后的 DELTA 文件转为 PAF 文件
    file_filtered_delta_paf = bn_delta + '.filtered.same.seq.delta.paf'
    cmd_delta2paf = f'delta2paf {file_filtered_delta} > {file_filtered_delta_paf}'
    check_file_in_path(file_filtered_delta_paf, cmd_delta2paf)

    ### step3: 计算比对率
    cal_aln_rate(f'{file_filtered_delta_paf}')

    ### step4: 去除冗余
    list_remove = []
    # for lines in open('pctg_self_align.filtered.same.seq.delta.pafaln.ratio.txt', 'r'):
    for lines in open('pctg_self_align_maxmatch.filtered.same.seq.delta.pafaln.ratio.txt', 'r'):
        if lines.startswith('#'):
            continue
        line = lines.strip().split()
        len_qry = int(line[2])
        len_ref = int(line[3])
        ratio_qry = float(line[5])
        ratio_ref = float(line[6])

        ## 判断 reference 与 query 的长短
        if len_qry < len_ref:
            if ratio_qry > 0.1 and line[0] not in list_remove:
                print(f'Ref length long: {line[0]}')
                print(line)
                list_remove.append(line[0])
        elif len_ref < len_qry:
            if ratio_ref > 0.1 and line[1] not in list_remove:
                print(f'Qry length long: {line[1]}')
                print(line)
                list_remove.append(line[1])

    for i in list_remove:
        print(i)
    return None


def process_delta(file_delta, args):

    def delta2paf(file_delta, file_paf):
        ''' convert delta to PAF '''
        cmd_delta2paf = f'delta2paf {file_delta} > {file_paf}'
        check_file_in_path(file_paf, cmd_delta2paf)

    def get_contig_pair(contig1, contig2, len_contig1, len_contig2):
        ''' contig pair based on contig length '''
        if len_contig1 > len_contig2:
            return (contig1, contig2)
        else:
            return (contig2, contig1)

    def filter_paf(file_paf, file_paf_filter, length_alignment, length_query):
        ''' filter PAF '''
        dic_paf_info = defaultdict(lambda: [0, 0, 0, 0, 0.00, 0.00])
        processed_contig_pair = set()
        for lines in open(file_paf, 'r'):
            line = lines.strip().split()
            ## same seq
            if line[0] == line[5]: continue
            ## alignment length
            if int(line[9]) < length_alignment: continue
            ## query length
            if int(line[1]) < length_query: continue
            ## 判断是否处理过？ctg1, ctg2
            if (line[5], line[0]) not in processed_contig_pair:
                processed_contig_pair.add((line[0], line[5]))
            else:
                continue
            ## contig pair
            contig_pair = get_contig_pair(line[0], line[5], int(line[1]), int(line[6]))
            dic_paf_info[contig_pair][0] = int(line[1])     # # ctg1_len, ctg2_len, +=alignment_qry_len, +=alignment_ref_len
            dic_paf_info[contig_pair][1] = int(line[6])
            dic_paf_info[contig_pair][2] += int(line[9])
            dic_paf_info[contig_pair][3] += int(line[10])
            dic_paf_info[contig_pair][4] = dic_paf_info[contig_pair][2] / dic_paf_info[contig_pair][0]
            dic_paf_info[contig_pair][5] = dic_paf_info[contig_pair][3] / dic_paf_info[contig_pair][1]
        return dic_paf_info

    ## delta2paf
    basename = os.path.basename(file_delta).replace('.delta', '')
    file_paf = basename + '.paf'
    delta2paf(file_delta, file_paf)

    ## filter paf
    file_paf_filter = basename + '.filter.paf'
    dic_paf_info = filter_paf(file_paf, file_paf_filter, args.min_alginment_length, args.min_contig_length)

    ## Remove redundancy
    remove_contig = set()
    out = open('removed_contigs_based_on_mummer.id.txt', 'w')
    for contg_pair, values in dic_paf_info.items():
        if values[4] > args.match_ratio or values[5] > args.match_ratio:
            contig_long = contg_pair[0]     # 长的 contig
            contig_short = contg_pair[1]    # 短的 contig
            if contig_short not in remove_contig:
                remove_contig.add(contig_short)
            # print(contg_pair, values)
            # print(contig_short)
    for i in remove_contig:
        out.write(i + '\n')
    return dic_paf_info


## split fa
def split_fa(dic_fasta, min_contig_length, wd):
    ''' split fasta file one sequence one file '''
    file_split_list = set()                     # 存储拆分文件的路径
    # 确保工作目录存在
    os.makedirs(wd, exist_ok=True)
    for ctgid, ctgseq in dic_fasta.items():
        if len(ctgseq) < min_contig_length:
            continue
        elif len(ctgseq) > min_contig_length:
            # 构建输出文件的路径
            output_path = os.path.join(wd, ctgid + '.fa')
            file_split_list.add(output_path)
            # 如果文件已经存在，则跳过
            if os.path.exists(output_path):
                continue
            with open(output_path, 'w') as out:
                out.write('>' + ctgid + '\n' + ctgseq + '\n')
                out.close()
    return file_split_list


def cal_distance(file_split_set, wd, args):
    ''' contig all_vs_all distance calculate '''
    def run_mash_sketch(file_split_set, wd):
        commands = []
        for file in file_split_set:
            basename = os.path.basename(file)
            cmd_mash = f'mash sketch -o {wd}/{basename}.mash {file}'
            commands.append(cmd_mash)
        return commands
    if not os.path.exists(f'{wd}/distances.txt'):
        ## step1: mash sketch
        os.makedirs(wd, exist_ok=True)
        commands = run_mash_sketch(file_split_set, wd)
        run_in_parallel(commands, args.process)
        ## step2: mach dist
        commands = []
        file_split_list = list(file_split_set)
        for i in range(len(file_split_list)):
            for j in range(i + 1, len(file_split_list)):
                file1 = file_split_list[i]
                file2 = file_split_list[j]
                cmd_mash = f'mash dist {file1} {file2} >> {wd}/distances.txt'
                commands.append(cmd_mash)
        run_in_parallel(commands, args.process)

        ## step3: filter dist
    if not os.path.exists(f'{wd}/similar_contig_pairs.txt'):
        df = pd.read_csv(f'{wd}/distances.txt', sep='\t', header=None)  # 读取 Mash 距离文件
        df.columns  = ['seq1', 'seq2', 'dist', 'pvalue', 'shared_hashes']
        similar_contig_pairs = df[df['dist'] < args.min_distance]                            # 筛选距离小于阈值的 contig pairs
        similar_contig_pairs[['seq1', 'seq2']].to_csv(f'{wd}/similar_contig_pairs.txt', sep='\t', index=False)


def contig_pair_aln(file_similar_contig_pairs, wd, args):
    ''' contig pair alignment '''
    if not os.path.exists(f'{wd}/merge.paf') or os.path.getsize(f'{wd}/merge.paf') == 0:
        commands = []
        for lines in open(file_similar_contig_pairs, 'r'):
            if not lines.startswith('seq1'):
                line = lines.strip().split()
                file1 = line[0]
                file2 = line[1]
                file1_bn = os.path.basename(file1)
                file2_bn = os.path.basename(file2)
                cmd_minimap2 = f'minimap2 -cx asm5 -t {args.threads} {file1} {file2} > {wd}/{file1_bn}_{file2_bn}.paf'
                commands.append(cmd_minimap2)
        run_in_parallel(commands, args.process)
        ## merge paf
        # subprocess.run(f'cat {wd}/*.paf > {wd}/merge.paf', shell=True, close_fds=True)    # Too many files will cause an error
        subprocess.run(f'find {wd}/*.paf -name "*.paf" -print0 | xargs -0 cat > {wd}/merge.paf')


def find_identity_contig_pair(file_paf):
    cigar_pattern = re.compile(r'(\d+)([A-Z])')
    dic_paf_info = defaultdict(lambda: [0, 0, 0, 0, 0])

    for lines in open(file_paf, 'r'):
        line = lines.strip().split()

        # 跳过次要比对
        if line[16][-1] == 'S':
            continue

        cigar = line[-1][5:]
        matches = cigar_pattern.findall(cigar)
        match_count = sum(int(length) for length, op in matches if op == 'M')

        contig1, contig2 = line[0], line[5]
        contig_pair = (contig1, contig2)

        contig1_length = int(line[1])
        contig2_length = int(line[6])

        dic_paf_info[contig_pair][0] = contig1_length     # ctg1_len, ctg2_len, match count, ctg1_match_ratio, ctg2_match_ratio
        dic_paf_info[contig_pair][1] = contig2_length
        dic_paf_info[contig_pair][2] += match_count
        dic_paf_info[contig_pair][3] = dic_paf_info[contig_pair][2] / contig1_length
        dic_paf_info[contig_pair][4] = dic_paf_info[contig_pair][2] / contig2_length
    return dic_paf_info


def remove_redundancy(merge_paf, wd, dic_contig, args):
    def get_contig_pair(contig1, contig2, len1, len2):
        ''' contig pair based on contig length, Return (longer_contig, shorter_contig) '''
        if len1 > len2:
            return (contig1, contig2)
        else:
            return (contig2, contig1)

    def get_contig_coord(len1, len2, start1, start2):
        ''' contig coord based on contig length, Return (longer_contig start, shorter_contig start) '''
        if len1 > len2:
            return (start1, start2)
        else:
            return (start2, start1)

    def filter_paf(file_paf, length_alignment, length_query):
        ''' filter PAF '''
        cigar_pattern = re.compile(r'(\d+)([A-Z])')

        # Dictionary to store contig pair statistics
        dic_paf_info = defaultdict(lambda: [0, 0, 0, 0, 0.00, 0.00])
        # Dictionary to store alignment coordinates for R² calculation
        coord_dict = defaultdict(list)
        # Dictionary to store alignment coordinates in target sequence
        alignment_position = defaultdict(list)

        # processed_contig_pair = set()
        for lines in open(file_paf, 'r'):
            line = lines.strip().split()

            ## jumping supplementary alignment
            if line[16][-1] == 'S':
                continue

            ## Filter: alignment length
            if int(line[9]) < length_alignment: continue

            ## Filter: query length
            if int(line[1]) < length_query: continue

            cigar = line[-1][5:]
            matches = cigar_pattern.findall(cigar)
            match_count = sum(int(length) for length, op in matches if op == 'M')

            # time: 2025-04-23 22:20:36 -- match count 可能太严谨了，造成Haplotype之间的比对率特别低，大概只有13-15%，这时候在去除冗余的时候，阈值太高，则会过多保留；阈值太低，则会过多去除。
            #   可以尝试使用paf 文件第十列和十一列的信息来计算。
            #   尝试之后，没有太大区别
            match_count = int(line[9])

            ## contig pair
            qname, tname = line[0], line[5]
            qlen, tlen = int(line[1]), int(line[6])
            qstart, tstart = int(line[2]), int(line[7])
            qend, tend = int(line[3]), int(line[8])

            contig_pair = get_contig_pair(qname, tname, qlen, tlen)     # contig_pair：（long contig, short contig）

            if match_count > 5000:
                # 关键修复：alignment 坐标必须为 longer contig 上的
                if qlen > tlen:
                    long_start, long_end, long_len = qstart, qend, qlen
                else:
                    long_start, long_end, long_len = tstart, tend, tlen
                alignment_position[contig_pair].append((min(long_start, long_end), max(long_start, long_end), long_len))

                coord_pair = get_contig_coord(qlen, tlen, qstart, tstart)
                # Record coordinates for R² evaluation
                coord_dict[contig_pair].append(coord_pair)

            dic_paf_info[contig_pair][0] = max(qlen, tlen)     # dic_paf_info: key: (long ctg, short ctg); value: [long len, short len, long match count, short match count, long match ratio, short match ratio]
            dic_paf_info[contig_pair][1] = min(qlen, tlen)
            dic_paf_info[contig_pair][2] += match_count
            dic_paf_info[contig_pair][3] += match_count

        for pair, values in dic_paf_info.items():
            values[4] = values[2] / values[0] if values[0] else 0.0
            values[5] = values[3] / values[1] if values[1] else 0.0
            # dic_paf_info[contig_pair][4] = dic_paf_info[contig_pair][2] / dic_paf_info[contig_pair][0]      # ctg1 match ratio / long conitg
            # dic_paf_info[contig_pair][5] = dic_paf_info[contig_pair][3] / dic_paf_info[contig_pair][1]      # ctg2 match ratio / short contig
        return dic_paf_info, coord_dict, alignment_position

    def compute_r2(coord_dict, threshold):
        """Compute R² (collinearity) for each contig pair"""
        r2_flag = {}
        for pair, coords in coord_dict.items():
            if len(coords) < 2:
                r2_flag[pair] = False
                continue
            qstarts, tstarts = zip(*coords)
            # print(pair, coords)
            r_squared = linregress(qstarts, tstarts).rvalue ** 2
            print(pair, coords, r_squared)
            r2_flag[pair] = r_squared >= threshold
        return r2_flag

    def is_contained(contig_pair, alns, tlen, margin_ratio=0.05, ratio_threshold=0.8):
        """
        判断一个 contig 是否被包含在另一个 contig 中。
        - alns: [(qstart, qend), ...]
        - tlen: target contig 长度
        - margin_ratio: 边缘比例，如 0.05 表示前5%和后5%
        - ratio_threshold: 比对落在内部的比例阈值
        """
        if not alns or tlen <= 0:
            return False
        margin = int(tlen * margin_ratio)
        total = len(alns)
        inside = sum(1 for start, end in alns if start > margin and end < tlen - margin)
        # print(alns, inside / total)
        a = (inside / total) >= ratio_threshold
        print('is contained?', contig_pair, {inside/total}, a)
        return (inside / total) >= ratio_threshold

    ### step1: parse PAF and calculate match ratios and coordinate data
    file_paf_filter = os.path.join(wd, 'merge.filter.paf')
    dic_paf_info, coord_dict, aln_pos = filter_paf(merge_paf, args.min_alignment_length, args.min_contig_length)

    ### step2: compute collinearity (R²) for each contig pair
    r2_result = compute_r2(coord_dict, threshold=args.r2_threshold)

    ### step3: write debug info
    with open(os.path.join(wd, 'removed_contigs_match.ratio.txt'), 'w') as f:
        for key, values in dic_paf_info.items():
            r2_status = r2_result.get(key, False)
            str_values = "\t".join(map(str, values))
            f.write(f'{key[0]}\t{key[1]}\t{str_values}\tR2_pass={r2_status}\n')

    ### step4: determine which contigs to remove
    remove_contig = set()
    # processed_contig = set()      # 存储处理过的 contig -- 目的：对于 contig1-contig2，如果 contig1 已经被移除，那么对于 contig1-contig3，contig3 较短，这时候按照规则需要去除 contig3，但是 contig1 已经去除，contig3 需要保留。
    for contig_pair, values in dic_paf_info.items():
        # match_ratio_pass = values[4] > args.match_ratio or values[5] > args.match_ratio
        match_ratio_pass = values[4] > args.match_ratio_target and values[5] > args.match_ratio_query        # long contig and short contig 都超过阈值; 冗余比较多
        collinear_pass = r2_result.get(contig_pair, False)

        alns_full = aln_pos.get(contig_pair, [])
        if alns_full:
            tlen = alns_full[0][2]
            alns = [(s, e) for s, e, _ in alns_full]
            contained_pass = is_contained(contig_pair, alns, tlen, margin_ratio=args.internal_margin_ratio,
                                          ratio_threshold=args.internal_ratio_threshold)
        else:
            contained_pass = False

        # if contig_pair[1] not in remove_contig:
        #     if (match_ratio_pass and collinear_pass) or contained_pass:
        #         remove_contig.add(contig_pair[1])

        if contig_pair[1] not in remove_contig:
            # if (match_ratio_pass and collinear_pass) or (contained_pass and collinear_pass):
            if (match_ratio_pass and contained_pass and collinear_pass):
                remove_contig.add(contig_pair[1])

    ### Step 5: write output FASTA files
    with open(f'{wd}/removed.fa', 'w') as out_remove, open(f'{wd}/retained.fa', 'w') as out_retain:
        for id, seq in dic_contig.items():
            if len(seq) >= args.min_chr_length:
                output = out_remove if id in remove_contig else out_retain
                output.write(f'>{id}\n{seq}\n')

        # if values[4] > args.match_ratio or values[5] > args.match_ratio:
        #     contig_long = contg_pair[0]     # 长的 contig
        #     contig_short = contg_pair[1]    # 短的 contig
        #     if contig_short not in remove_contig:
        #         # if contig_long not in remove_contig:
        #             remove_contig.add(contig_short)     # 需要提升：一个长的对应多个短的
    #
    # ## output: output removed.fa and retained.fa
    # with open(f'{wd}/removed.fa', 'w') as out_remove, open(f'{wd}/retained.fa', 'w') as out_retain:
    #     for id, seq in dic_contig.items():
    #         if len(seq) > 1000000:
    #             output = out_remove if id in remove_contig else out_retain
    #             output.write(f'>{id}\n{seq}\n')

    return None



import re

def remove_redundancy_v2(merge_paf, wd, dic_contig, args):
    def get_contig_pair(contig1, contig2, len1, len2):
        if len1 > len2:
            return (contig1, contig2)
        else:
            return (contig2, contig1)

    def get_contig_coord(len1, len2, start1, start2):
        if len1 > len2:
            return (start1, start2)
        else:
            return (start2, start1)

    def filter_paf(file_paf, length_alignment, length_query):
        cigar_pattern = re.compile(r'(\d+)([A-Z])')
        dic_paf_info = defaultdict(lambda: [0, 0, 0, 0, 0.00, 0.00])
        coord_dict = defaultdict(list)
        alignment_position = defaultdict(list)

        with open(file_paf, 'r') as f:
            for line in f:
                line = line.strip().split()
                if line[16][-1] == 'S':
                    continue
                if int(line[9]) < length_alignment: continue
                if int(line[1]) < length_query: continue

                match_count = int(line[9])
                qname, tname = line[0], line[5]
                qlen, tlen = int(line[1]), int(line[6])
                qstart, tstart = int(line[2]), int(line[7])
                qend, tend = int(line[3]), int(line[8])

                contig_pair = get_contig_pair(qname, tname, qlen, tlen)

                if match_count > 5000:
                    if qlen > tlen:
                        long_start, long_end, long_len = qstart, qend, qlen
                    else:
                        long_start, long_end, long_len = tstart, tend, tlen
                    alignment_position[contig_pair].append(
                        (min(long_start, long_end), max(long_start, long_end), long_len))

                    coord_pair = get_contig_coord(qlen, tlen, qstart, tstart)
                    coord_dict[contig_pair].append(coord_pair)

                dic_paf_info[contig_pair][0] = max(qlen, tlen)
                dic_paf_info[contig_pair][1] = min(qlen, tlen)
                dic_paf_info[contig_pair][2] += match_count
                dic_paf_info[contig_pair][3] += match_count

        for pair, values in dic_paf_info.items():
            values[4] = values[2] / values[0] if values[0] else 0.0
            values[5] = values[3] / values[1] if values[1] else 0.0
        return dic_paf_info, coord_dict, alignment_position

    def compute_r2(coord_dict, threshold):
        r2_flag = {}
        for pair, coords in coord_dict.items():
            if len(coords) < 2:
                r2_flag[pair] = False
                continue
            qstarts, tstarts = zip(*coords)
            r_squared = linregress(qstarts, tstarts).rvalue ** 2
            r2_flag[pair] = r_squared >= threshold
        return r2_flag

    def is_contained(contig_pair, alns, tlen, margin_ratio=0.05, ratio_threshold=0.8):
        if not alns or tlen <= 0:
            return False
        margin = int(tlen * margin_ratio)
        total = len(alns)
        inside = sum(1 for start, end in alns if start > margin and end < tlen - margin)
        return (inside / total) >= ratio_threshold

    def try_merge_lis_lds(paf_lines, qname, tname, qseq, tseq, r2_thresh, cov_thresh, dic_contig):

        # paf_subset = [line.strip().split('\t') for line in paf_lines
        #               if (line.startswith(qname + '\t') and f'\t{tname}\t' in line)
        #               or (line.startswith(tname + '\t') and f'\t{qname}\t' in line)]
        paf_subset = []
        for line in paf_lines:
            cols = line.strip().split('\t')
            qn, tn = cols[0], cols[5]
            qlen, tlen = int(cols[1]), int(cols[6])

            # 判断是否是我们关心的 qname/tname 配对
            if set([qn, tn]) != set([qname, tname]):
                continue

            # 如果当前 query（cols[0]）比 target 短，则交换
            if qn == qname and qlen < tlen:
                # 保持顺序，qname 是 query，tname 是 target
                paf_subset.append(cols)
            elif tn == qname and tlen < qlen:
                # 当前是反方向，但 target 更长，应交换为 query 和 target
                cols_swapped = cols[:]
                cols_swapped[0], cols_swapped[5] = cols[5], cols[0]
                cols_swapped[1], cols_swapped[6] = cols[6], cols[1]
                cols_swapped[2], cols_swapped[7] = cols[7], cols[2]
                cols_swapped[3], cols_swapped[8] = cols[8], cols[3]
                # Strand 不变，这里只管位置
                paf_subset.append(cols_swapped)
            elif qlen >= tlen:
                # qname 已经是 query，并且比 target 更长，交换角色
                cols_swapped = cols[:]
                cols_swapped[0], cols_swapped[5] = cols[5], cols[0]
                cols_swapped[1], cols_swapped[6] = cols[6], cols[1]
                cols_swapped[2], cols_swapped[7] = cols[7], cols[2]
                cols_swapped[3], cols_swapped[8] = cols[8], cols[3]
                paf_subset.append(cols_swapped)
            else:
                # 其余情况，保持不变
                paf_subset.append(cols)

        if not paf_subset:
            return None
        records = []
        for i, line in enumerate(paf_subset):
            qstart, qend, tstart, tend = int(line[2]), int(line[3]), int(line[7]), int(line[8])
            records.append((qstart, tstart, i))

        def get_lis(points, reverse=False):
            from bisect import bisect_left
            points = [(x, -y if reverse else y, i) for x, y, i in points]
            points.sort()
            tails, idxs, parent = [], [], [-1] * len(points)
            for i, (x, y, _) in enumerate(points):
                pos = bisect_left(tails, y)
                if pos == len(tails):
                    tails.append(y)
                    idxs.append(i)
                else:
                    tails[pos] = y
                    idxs[pos] = i
                if pos > 0:
                    parent[i] = idxs[pos - 1]
            lis = []
            k = idxs[-1]
            while k != -1:
                lis.append(points[k][2])
                k = parent[k]
            lis.reverse()
            return [paf_subset[i] for i in lis]

        def calc_stats(lis):
            q, t = zip(*[(int(l[2]), int(l[7])) for l in lis])
            r2 = linregress(q, t).rvalue ** 2
            return r2, len(lis) / len(records)

        lis_inc = get_lis(records, reverse=False)
        lis_dec = get_lis(records, reverse=True)

        r2_inc, cov_inc = calc_stats(lis_inc)
        r2_dec, cov_dec = calc_stats(lis_dec)

        if len(lis_inc) >= len(lis_dec):
            best, r2, cov, trend = lis_inc, r2_inc, cov_inc, 'increasing'
        else:
            best, r2, cov, trend = lis_dec, r2_dec, cov_dec, 'decreasing'

        if cov < cov_thresh or r2 < r2_thresh:
            return None

        qstarts = [int(line[2]) for line in best]
        tends = [int(line[8]) for line in best]
        tstarts = [int(line[7]) for line in best]
        qends = [int(line[3]) for line in best]

        qlen = len(qseq)
        tlen = len(tseq)
        if trend == 'decreasing':
            tseq = ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}[base] for base in reversed(tseq))
            tstarts = [tlen - int(t) for t in tends]
            tends = [tlen - int(t) for t in tstarts]

        q_pos = 'head' if sum(qstarts) < sum(qlen - q for q in qstarts) else 'tail'
        t_pos = 'head' if sum(tstarts) < sum(tlen - t for t in tstarts) else 'tail'

        if q_pos == 'head' and t_pos == 'tail':
            merged = tseq[:min(tstarts)] + qseq[min(qstarts):]
            print('head -- tail')
        elif q_pos == 'tail' and t_pos == 'head':
            merged = qseq[:max(qends)] + tseq[max(tends):]
            print('tail -- head')
        elif q_pos == 'head' and t_pos == 'head':
            merged = qseq[:max(qends)] + tseq[max(tends):]
            print('head -- head')
        elif q_pos == 'tail' and t_pos == 'tail':
            merged = tseq[:min(tstarts)] + qseq[min(qstarts):]
            print('tail -- tail')
        else:
            return None

        return merged

    file_paf_filter = os.path.join(wd, 'merge.filter.paf')
    dic_paf_info, coord_dict, aln_pos = filter_paf(merge_paf, args.min_alignment_length, args.min_contig_length)
    r2_result = compute_r2(coord_dict, threshold=args.r2_threshold)

    with open(os.path.join(wd, 'removed_contigs_match.ratio.txt'), 'w') as f:
        for key, values in dic_paf_info.items():
            r2_status = r2_result.get(key, False)
            str_values = "\t".join(map(str, values))
            f.write(f'{key[0]}\t{key[1]}\t{str_values}\tR2_pass={r2_status}\n')

    remove_contig = set()
    merged_contig = {}

    with open(merge_paf) as f:
        paf_lines = f.readlines()

    for contig_pair, values in dic_paf_info.items():
        match_ratio_pass = values[4] > args.match_ratio_target and values[5] > args.match_ratio_query
        collinear_pass = r2_result.get(contig_pair, False)

        alns_full = aln_pos.get(contig_pair, [])
        if alns_full:
            tlen = alns_full[0][2]
            alns = [(s, e) for s, e, _ in alns_full]
            contained_pass = is_contained(contig_pair, alns, tlen, margin_ratio=args.internal_margin_ratio,
                                          ratio_threshold=args.internal_ratio_threshold)
        else:
            contained_pass = False
        if contig_pair[1] not in remove_contig:
            if (match_ratio_pass and contained_pass and collinear_pass):
                remove_contig.add(contig_pair[1])

    list_retain = []
    list_remove = []

    with open(f'{wd}/removed.fa', 'w') as out_remove, \
            open(f'{wd}/retained.fa', 'w') as out_retain:
            # open(f'{wd}/merged.fa', 'w') as out_merged:
        for id, seq in dic_contig.items():
            if len(seq) >= args.min_chr_length:
                if id in remove_contig:
                    out_remove.write(f'>{id}\n{seq}\n')
                    list_remove.append(id)
                else:
                    out_retain.write(f'>{id}\n{seq}\n')
                    list_retain.append(id)
    print(list_retain)

    with open(f'{wd}/mT2T.fa', 'w') as out_merged:
        for contig_pair, value in dic_paf_info.items():
            qname, tname = contig_pair[1], contig_pair[0]
            qseq, tseq = dic_contig[qname], dic_contig[tname]
            if qname in list_retain and tname in list_retain:
                print('lis lds')
                merged_seq = try_merge_lis_lds(paf_lines, qname, tname, qseq, tseq, args.r2_threshold, args.lis_cov, dic_contig)
                if merged_seq:
                    out_merged.write(f'>merged_{tname}_{qname}\n{merged_seq}\n')
                    list_retain.remove(qname)
                    list_retain.remove(tname)

        for id in list_retain:
            out_merged.write(f'>{id}\n{dic_contig[id]}\n')


    return None