#!/usr/bin/env python
'''
time: 2024-06-12
author: pxxiao
version: 1.0
description: 对 allelic table 重新洗牌
----\n
time: 2024-06-17
author: pxxiao
version: 2.0
description: 考虑collapsed unitig 信息，每一行 allelic table haplotype 加和不应该超过设置的 haplotype 数量。
                如果 unitig 为 diplotig，那么代表两个 haplotype，triplotig 代表三个，tetraplotig 代表四个，haplotig 代表一个。
                如果 allelic table 加起来超过 4，则需要对这一行的 allelic table 进行过滤。
                如果这个 unitig 在前面处理过，那么这个 unitig 可以过滤掉。过滤掉一个 unitig 之后，如果还是超过 4，那么再过滤一个，直到小于等于 4。
----\n
time: 2024-06-18
author: pxxiao
version: 3.0
description: 在过滤的之前，需要增加一个条件：
                如果这一行的 allelic table haplotype 加和大于4，这一行除了 new_unitig 已经在前面出现过，并且 new_unitig 在后面还会出现，
                那么这一行allelic table 的信息就是冗余的，为了避免错误，这一行可以直接过滤掉，标记为 None。
----
time: 2024-07-19
author: pxxiao
version: 3.2
description: debugs:
before:
   4731 chr08   21500000        21600000        utg000087l      utg000061l      utg000598l
   4732 chr08   21600000        21700000        utg000087l      utg000061l      utg000598l      utg002852l
   4733 chr08   21700000        21800000        utg000061l      utg000087l      utg002852l
   4734 chr08   21800000        21900000        utg000061l      utg000087l      utg001068l      utg000274l
   4735 chr08   21900000        22000000        utg000061l      utg000087l      utg001068l
   4736 chr08   22000000        22100000        utg000087l      utg001068l      utg000061l

after:
   4078 chr08   20000000        20100000        utg000061l      utg000087l      utg000274l
   4079 chr08   21700000        21800000        utg000061l      utg000274l      utg002852l
   4080 chr08   21900000        22000000        utg000061l      utg000087l      utg001068l
   4081 chr08   22000000        22100000        utg000061l      utg000087l      utg001068l
'''


import argparse


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


def parse_allelic_table(file_allelic_table):
    ''' 解析 allelic table 文件 '''
    allelic_table = []
    with open(file_allelic_table, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            unitigs = parts[3:]
            allelic_table.append((chromosome, start, end, unitigs))
    return allelic_table


## 解析 contig type
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
    return dic_dosage_contig_type


def get_unitig_lengths(unitig_ids, dic_fasta_len):
    ''' 获取 unitig 的长度信息，可以从外部文件或预定义的字典中读取 '''
    # 假设这里有一个预定义的字典存储每个 unitig 的长度
    return {unitig_id: dic_fasta_len.get(unitig_id, 0) for unitig_id in unitig_ids}


def filter_allelic_table(allelic_table, unitig_types, contig_type_haplotypes, dic_fasta_len):
    ''' 根据 contig type，对 allelic table 进行过滤 '''
    seen_unitigs = set()        # 前面 allelic table 中出现过的 unitig
    filtered_table = []
    total_bins = len(allelic_table)
    print(f'allelic table: {allelic_table}')

    for idx, row in enumerate(allelic_table):
        print(f'pxp {row}')
        chrom, start, end, unitigs = row
        haplotype_count = 0
        filtered_unitigs = []
        removed_unitig = []

        # 计算每个 unitig 的长度并排序 -- 短到长，后面 pop 的时候，就是从长到短进行过滤
        unitig_lengths = get_unitig_lengths(unitigs, dic_fasta_len)
        sorted_unitigs = sorted(unitigs, key=lambda x: unitig_lengths[x], reverse=True)

        for unitig in sorted_unitigs:
            unitig_type = unitig_types.get(unitig, 'haplotig')
            haplotype_count += contig_type_haplotypes[unitig_type]
            if unitig in seen_unitigs:
                filtered_unitigs.append(unitig)

        # 检查是否所有的 unitig 都在前面出现过并且将会在后面出现
        if haplotype_count > 4:
            if idx > 0 and idx < total_bins - 1:
                prev_unitigs = set(allelic_table[idx - 1][3])
                next_unitigs = set(allelic_table[idx + 1][3])
                new_unitigs = set(unitigs) - prev_unitigs

                if all(unitig in next_unitigs for unitig in new_unitigs):
                    filtered_table.append(None)
                    continue

        while haplotype_count > 4 and filtered_unitigs:
            removed_unitig.append(filtered_unitigs.pop())
            haplotype_count -= contig_type_haplotypes[unitig_types[removed_unitig[-1]]]
        print(f'pxp {unitigs} {removed_unitig}')
        filtered_unitigs_set = set(unitigs) - set(removed_unitig)
        filtered_table.append((chrom, start, end, sorted(filtered_unitigs_set)))
        seen_unitigs.update(filtered_unitigs_set)

    return filtered_table


def correct_allelic_table(allelic_table, dic_fasta_len, dic_unitig_types, contig_type_haplotypes, top_n, search_range):
    ''' 修正 allelic table '''
    corrected_table = []
    total_bins = len(allelic_table)

    for i in range(total_bins):
        chromosome, start, end, unitigs = allelic_table[i]
        print(f"Processing bin: {chromosome}:{start}-{end}")
        original_unitigs = set(unitigs)

        # 获取前 search_range 个和后 search_range 个 bin 的 unitig
        prev_unitigs = set()
        for j in range(max(0, i - search_range), i):
            prev_unitigs.update(allelic_table[j][3])

        next_unitigs = set()
        for j in range(i + 1, min(total_bins, i + 1 + search_range)):
            next_unitigs.update(allelic_table[j][3])

        # # 获取前一个和后一个 bin 的 unitig
        # prev_unitigs = set(allelic_table[i - 1][3]) if i > 0 else set()
        # next_unitigs = set(allelic_table[i + 1][3]) if i < len(allelic_table) - 1 else set()

        # 找出当前 bin 中缺失的 unitig
        missing_unitigs = (prev_unitigs & next_unitigs) - original_unitigs

        if missing_unitigs:
            print(f"Missing unitigs in bin {chromosome}:{start}-{end}: {', '.join(missing_unitigs)}")

            # 获取 unitig 的长度信息
            unitig_lengths = get_unitig_lengths(original_unitigs | missing_unitigs, dic_fasta_len)

            # 将缺失的 unitig 添加进来
            new_unitigs = original_unitigs | missing_unitigs

            print('Debugs: new_unitigs ', new_unitigs)

            # 如果新集合的大小大于 top_n，移除最短的 unitig
            while len(new_unitigs) > top_n:
                shortest_unitig = min(new_unitigs, key=lambda x: unitig_lengths[x])
                print(f"Removed unitig from bin {chromosome}:{start}-{end}: {shortest_unitig}")
                new_unitigs.remove(shortest_unitig)

            corrected_table.append((chromosome, start, end, sorted(new_unitigs)))
        else:
            corrected_table.append((chromosome, start, end, sorted(original_unitigs)))

    print(corrected_table)

    return filter_allelic_table(corrected_table, dic_unitig_types, contig_type_haplotypes, dic_fasta_len)


def write_corrected_table1(corrected_table, output_file):
    ''' 将修正后的 table 写入文件 '''
    with open(output_file, 'w') as f:
        for chromosome, start, end, unitigs in corrected_table:
            line = f"{chromosome}\t{start}\t{end}\t" + "\t".join(unitigs) + "\n"
            f.write(line)


def write_corrected_table(corrected_table, output_file):
    ''' 将修正后的 table 写入文件 '''
    with open(output_file, 'w') as f:
        for row in corrected_table:
            if row is None:
                continue
            chromosome, start, end, unitigs = row
            line = f"{chromosome}\t{start}\t{end}\t" + "\t".join(unitigs) + "\n"
            f.write(line)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Refresh allelic table',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--allelic_table', type=str, help='Allelic table', required=True)
    parser.add_argument('--fasta', type=str, help='Fasta file', required=True)
    parser.add_argument('--contig_type', type=str, help='Contig type file', required=True)
    parser.add_argument('--wd', type=str, help='working directory', required=True)
    parser.add_argument('--top_n', type=int, help='The number of haplotypes', required=False, default=4)
    parser.add_argument('--search_range', type=int, help='Search range', required=False, default=5)
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    wd = args.wd
    dic_fasta = fasta_read(args.fasta)
    dic_fasta_len = {unitig_id: len(seq) for unitig_id, seq in dic_fasta.items()}
    dic_unitig_types = parse_contig_type_based_on_dosage(args.contig_type)

    # top_n = 4           # 设定的 haplotype 值
    # search_range = 5    # 设置搜索范围

    contig_type_haplotypes = {
        "diplotig": 2,
        "triplotig": 3,
        "tetraplotig": 4,
        "haplotig": 1
    }

    output_file = f'{wd}/corrected_allelic_table.txt'  # 输出文件路径

    allelic_table = parse_allelic_table(args.allelic_table)
    corrected_table = correct_allelic_table(allelic_table, dic_fasta_len, dic_unitig_types, contig_type_haplotypes, args.top_n, args.search_range)
    write_corrected_table(corrected_table, output_file)


if __name__ == '__main__':
    main()