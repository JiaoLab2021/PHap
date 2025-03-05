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
----
time: 2024-09-29
author: pxxiao
version: 3.0
description:
添加一个参数，控制 ONT reads 的长度，用于 hifiasm UL 组装
'''


import argparse
import multiprocessing
import os
import subprocess
import random
import pysam
from collections import defaultdict
import pickle
import logging

# 设置日志记录配置
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


def execute_command(command, env=None):
    '''
    执行命令并处理可能的错误
    :param command: 要执行的命令字符串
    :return: None
    '''
    try:
        process = subprocess.Popen(command, shell=True, env=env)
        process.wait()
        if process.returncode != 0:
            logging.error(f'Command failed with return code {process.returncode}: {command}')
    except Exception as e:
        logging.error(f'An error occurred while executing command: {command}\nError: {e}')
    return None


def run_in_parallel(commands, num_processes, env=None):
    '''
    多线程运行命令
    :param commands: 要执行的命令列表
    :param num_processes: 并行执行的进程数
    :param env: 环境变量字典
    :return: None
    '''
    # processes = []
    pool = multiprocessing.Pool(processes=num_processes)
    for command in commands:
        pool.apply_async(execute_command, args=(command, env, ))
    pool.close()
    pool.join()
    return None


def output_pickle(dict_, to):
    with open(to, 'wb') as fpkl:
        pickle.dump(dict_, fpkl)
    logging.info(f'Successfully saved dictionary to {to}')


def load_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    logging.info(f'Successfully loaded dictionary from {file_path}')
    return data


def parse_contig_type_based_on_dosage(file_contig_type):
    dic_dosage_contig_type = {}
    for lines in open(file_contig_type, 'r'):
        if lines.startswith('contig_ID'):
            continue
        else:
            line = lines.strip().split()
            contig_type = line[2]
            type_list = ['haplotig', 'diplotig', 'triplotig', 'tetraplotig']
            if contig_type not in type_list:
                contig_type = 'haplotig'
            dic_dosage_contig_type[line[0]] = contig_type
    logging.info('Parsed contig types successfully')
    return dic_dosage_contig_type


def parse_group_cluster(group_file):
    dic_group_ctg = {}
    for lines in open(group_file):
        line = lines.strip().split()
        group_name = line[0]
        dic_group_ctg[group_name] = line[1:]
    logging.info('Parsed group clusters successfully')
    return dic_group_ctg


def parse_bam_tgs(bam_file):
    # 打开 bam 文件
    bam = pysam.AlignmentFile(bam_file, 'rb')

    # 初始化字典存储结果
    contig_reads_dict = defaultdict(set)

    # 遍历 BAM 文件中的每个比对
    for read in bam.fetch():
        contig = bam.get_reference_name(read.reference_id)
        contig_reads_dict[contig].add(read.query_name)

    # 关闭 BAM 文件
    bam.close()
    logging.info(f'Parsed BAM {bam_file} file successfully')
    return contig_reads_dict


def parse_bam_ngs(bam_file):
    # 打开 bam 文件
    bam = pysam.AlignmentFile(bam_file, 'rb')

    # 初始化字典存储结果
    contig_reads_dict = defaultdict(set)

    # 遍历 BAM 文件中的每个比对
    for read in bam.fetch():
        # contig = bam.get_reference_name(read.reference_id)
        # 跳过未比对的 reads
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        # 仅处理第一端的 reads
        if read.is_read1:
            contig = bam.get_reference_name(read.reference_id)
            if contig:  # 确保 contig 不为 None
                contig_reads_dict[contig].add(read.query_name)

        # 仅处理第一端 或独立的单端
        # if read.is_read1 or read.is_unmapped or read.mate_is_unmapped:
        #     contig_reads_dict[contig].add(read.query_name)
    # 关闭 BAM 文件
    bam.close()
    logging.info(f'Parsed BAM {bam_file} file successfully')
    return contig_reads_dict


def split_reads(reads, fraction, seed=None):
    """将 HiFi reads 随机分割成指定的比例"""
    reads_list = list(reads)
    num_reads = len(reads_list)
    split_point = int(num_reads * fraction)
    # 设置种子，可以对结果进行重现
    if seed is not None:
        random.seed(seed)
    random.shuffle(reads_list)
    return set(reads_list[:split_point]), set(reads_list[split_point:])


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
    contig_type_file = 'contig_type.pickle'
    group_file = 'group_contig_type.pickle'
    contig_hifi_pickle = 'contig_hifi.pickle'
    contig_ont_pickle = 'contig_ont.pickle'
    contig_hic_pickle = 'contig_hic.pickle'

    ### step1: 解析 unitig type
    if not os.path.isfile(contig_type_file):
        dic_contig_type = parse_contig_type_based_on_dosage(args.contig_type)
        output_pickle(dic_contig_type, contig_type_file)
    dic_contig_type = load_pickle_file(contig_type_file)
    logging.debug(f'dic_contig_type: {dic_contig_type}')

    ### step2: 解析 group file
    if not os.path.exists(group_file):
        dic_group_ctg = parse_group_cluster(args.group)
        output_pickle(dic_group_ctg, group_file)
    dic_group_ctg = load_pickle_file(group_file)
    logging.debug(f'dic_group_ctg: {dic_group_ctg}')

    ### step3: 解析 bam file
    ## HiFi
    if not os.path.exists(contig_hifi_pickle):
        dic_contig_hifi = parse_bam_tgs(args.bam_hifi)
        output_pickle(dic_contig_hifi, contig_hifi_pickle)
    dic_contig_hifi = load_pickle_file(contig_hifi_pickle)
    logging.debug(f'dic_contig_hifi: {dic_contig_hifi}')
    ## ONT
    if not os.path.exists(contig_ont_pickle):
        dic_contig_ont = parse_bam_tgs(args.bam_ont)
        output_pickle(dic_contig_ont, contig_ont_pickle)
    dic_contig_ont = load_pickle_file(contig_ont_pickle)
    logging.debug(f'dic_contig_ont: {dic_contig_ont}')
    ## Hi-C
    if not os.path.exists(contig_hic_pickle):
        dic_contig_hic = parse_bam_ngs(args.bam_hic)
        output_pickle(dic_contig_hic, contig_hic_pickle)
    dic_contig_hic = load_pickle_file(contig_hic_pickle)
    logging.debug(f'dic_contig_hic: dic_contig_hic')    # 数据量太大了，注释掉

    ### step4: 输出每个 Group 的 HiFi / ONT / Hi-C reads FASTQ
    dic_group_hifi = defaultdict(set)
    dic_group_ont = defaultdict(set)
    dic_group_hic = defaultdict(set)
    # dic_remaining_reads = defaultdict(set)

    contig_split_fraction = {
        'haplotig': 1.0,
        'diplotig': 0.5,
        'triplotig': 1 / 3,
        'tetraplotig': 0.25
    }

    for group, unitigs in dic_group_ctg.items():
        for unitig in unitigs:
            contig_type = dic_contig_type.get(unitig, 'haplotig')
            fraction = contig_split_fraction.get(contig_type, 1.0)
            ## HiFi
            reads_hifi = dic_contig_hifi.get(unitig, set())
            part1_hifi, part2_hifi = split_reads(reads_hifi, fraction, seed=args.seed)
            dic_group_hifi[group].update(part1_hifi)
            # dic_remaining_reads[group].update(part2_hifi)
            ## ONT
            reads_ont = dic_contig_ont.get(unitig, set())
            part1_ont, part2_ont = split_reads(reads_ont, fraction, seed=args.seed+1)
            dic_group_ont[group].update(part1_ont)
            ## Hi-C
            reads_hic = dic_contig_hic.get(unitig, set())
            part1_hic, part2_hic = split_reads(reads_hic, fraction, seed=args.seed+2)
            dic_group_hic[group].update(part1_hic)

    ont_reads = args.ont
    ont_reads_bn = os.path.basename(ont_reads)
    ont_quality = args.ont_quality
    ont_length = args.ont_length
    if args.ont_length != 1 or args.ont_quality != 0:
        ont_reads_filter = ont_reads_bn + '.filter.ont.fq'
        ont_reads_filter_gzip = ont_reads_bn + '.filter.ont.fq.gz'
        cmd = f'~/tools/Anaconda3/envs/chopper/bin/chopper -l {ont_length} -q {ont_quality} -t {args.threads} -i <(zcat {ont_reads}) > {ont_reads_filter}; pigp -p 10 {ont_reads_filter}'
        subprocess.run(cmd, shell=True, executable='/bin/bash', close_fds=True)

    commands = []
    process_seqkit = args.process * args.threads // 1
    for group, reads in dic_group_hifi.items():
        with open(f'{group}.HiFi.txt', 'w') as f:
            f.write('\n'.join(reads) + '\n')
        cmd = f'seqkit grep -j {args.threads} -f {group}.HiFi.txt {args.hifi} > {group}.HiFi.fq; pigz -p {args.threads} {group}.HiFi.fq'
        commands.append(cmd)
        logging.info(f'Executing command: {cmd}')
    for group, reads in dic_group_ont.items():
        with open(f'{group}.ONT.txt', 'w') as f:
            f.write('\n'.join(reads) + '\n')
        if args.ont_length != 1 or args.ont_quality != 0:
            cmd = f'seqkit grep -j {args.threads} -f {group}.ONT.txt {ont_reads_filter_gzip} > {group}.ONT.fq; pigz -p {args.threads} {group}.ONT.fq'
        else:
            cmd = f'seqkit grep -j {args.threads} -f {group}.ONT.txt {args.ont} > {group}.ONT.fq; pigz -p {args.threads} {group}.ONT.fq'
        commands.append(cmd)
        logging.info(f'Executing command: {cmd}')
    for group, reads in dic_group_hic.items():
        with open(f'{group}.Hi-C.1.txt', 'w') as f1:
            f1.write('/1\n'.join(reads) + '/1' + '\n')
        with open(f'{group}.Hi-C.2.txt', 'w') as f2:
            f2.write('/2\n'.join(reads) + '/2' + '\n')
        cmd1 = f'seqkit grep -j {args.threads} -f {group}.Hi-C.1.txt {args.hic1} > {group}.Hi-C.1.fq; pigz -p 5 {group}.Hi-C.1.fq'
        cmd2 = f'seqkit grep -j {args.threads} -f {group}.Hi-C.2.txt {args.hic2} > {group}.Hi-C.2.fq; pigz -p 5 {group}.Hi-C.2.fq'
        commands.append(cmd1)
        commands.append(cmd2)
        logging.info(f'Executing command: {cmd1}')
        logging.info(f'Executing command: {cmd2}')
    run_in_parallel(commands, process_seqkit)

    ### step5: 对每个 group 进行组装（hifiasm）
    commands = []
    for group in dic_group_hifi.keys():
        os.makedirs(f'{group}.asm', exist_ok=True)
        cmd_hifiasm = f'nohup time -v hifiasm -t {args.threads} -o {group}.asm/{group}.asm --ul {group}.ONT.fq.gz {group}.HiFi.fq.gz > {group}.asm/log_hifiasm_out 2> {group}.asm/log_hifiasm_err'
        commands.append(cmd_hifiasm)
        logging.info(f'Executing command: {cmd_hifiasm}')
    run_in_parallel(commands, args.process)

    ### step6: 对每个 group p_ctg 进行挂载（haphic）
    commands = []
    for group in dic_group_hifi.keys():
        # 创建目录
        os.makedirs(f'{group}.asm/scaffolding/01_hic_mapping', exist_ok=True)
        os.makedirs(f'{group}.asm/scaffolding/02_haphic', exist_ok=True)
        # 构建命令序列
        cmd_scaffolding = [
            f'awk \'{{if($0~/^S/) print ">"$2"\\n"$3}}\' {group}.asm/{group}.asm.bp.p_ctg.gfa > {group}.asm/scaffolding/01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa',
            f'bwa index {group}.asm/scaffolding/01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa',
            f'bwa mem -5SP -t {args.threads} {group}.asm/scaffolding/01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa {group}.Hi-C.1.fq.gz {group}.Hi-C.2.fq.gz | /home/pxxiao/tools/Anaconda3/envs/haphic/bin/samblaster | samtools view - -@ {args.threads} -S -h -b -F 3340 -o {group}.asm/scaffolding/01_hic_mapping/HiC.bam',
            f'/home/pxxiao/tools/Assembly-tools/73_HapHiC/HapHiC/utils/filter_bam {group}.asm/scaffolding/01_hic_mapping/HiC.bam 1 --nm 3 --threads {args.threads} | samtools view - -b -@ {args.threads} -o {group}.asm/scaffolding/01_hic_mapping/HiC.filtered.bam',
            f'cd {group}.asm/scaffolding/02_haphic',
            f'nohup time -v /home/pxxiao/tools/Assembly-tools/73_HapHiC/1.0.5/HapHiC/haphic pipeline ../01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa ../01_hic_mapping/HiC.filtered.bam 1 --threads {args.threads} --processes {args.process} --Nx 100 > log_haphic_out 2> log_haphic_err',
            f'cd 04.build',
            f'nohup time -v bash juicebox.sh > log_juicebox_out 2> log_juicebox_err'
        ]
        # 将命令加入到 commands 列表中
        commands.append(" && ".join(cmd_scaffolding))
        logging.info(f'Executing command sequence for {cmd_scaffolding}')
    # 构建包含 haphic 环境的环境变量字典
    env_vars = os.environ.copy()
    env_vars['PATH'] = '/home/pxxiao/tools/Anaconda3/envs/haphic/bin:' + env_vars['PATH']
    env_vars['CONDA_DEFAULT_ENV'] = 'haphic'
    # 执行命令序列
    run_in_parallel(commands, args.process, env_vars)


if __name__ == '__main__':
    main()