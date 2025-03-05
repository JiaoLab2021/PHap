#!/usr/bin/env python
'''
time: 2024-06-11
author: pxxiao
version: 1.0
description: 根据 PAF，找到最长上升子序列 Longest Increasing Subsequence

利用最长上升子序列（LIS，Longest Increasing Subsequence）的算法可以帮助我们从一组可能包含噪声的比对中找出正确的连续比对区域。
这个方法适用于以下情况：你有多个比对，它们在 reference 上的起点是非递减的，并且你希望找出一组按顺序排列的比对，以此推断出正确的比对区域。

根据 reference 上的坐标信息，我们将所有的 alignments 划分为多个 LIS，然后根据 LIS 所包含的 alignment 的个数以及 LIS 的总的长度来选取最佳的 LIS 作为主要的比对区域。
----
time: 2024-06-13
author: pxxiao
version: 2.0
description: 添加对 PAF 过滤的步骤：有的 unitig 比对到 mT2T 的 match ratio 很低，也就是说只有部分可以比对上，对这部分 unitig 进行统计，并过滤.
这部分的 unitig 感觉是组装错误引起的：看了一下这些 unitig 与其他 unitig 的 Hi-C 信号，信号比较强的那些 unitig 都不在附近。
---
time: 2024-06-13 night
author: pxxiao
version: 3.0
description: 添加 merge LIS 的处理：如果 unitig 与 reference 的差异比较大，那么 alignments 在这里的间隔就比较大，这就会造成LIS 在这里断开，形成两个小的 LIS。
遇到这种情况，我们需要添加一个处理，判断 LIS 之间的距离，如果两个 LIS 包含的 alignment 大于阈值，且 LIS 之间的距离小于阈值，这时候可以考虑把小的 LIS 合并成大的 LIS。
如果这两个 LIS 之间，还有一些较小的 LIS，其包含的 alignment 数量小于阈值，也应该与这两个大的 LIS 合并进去。
---
time: 2024-06-26
author: pxxiao
version: 4.0
description: 添加 LIS 合并的条件：一个 LIS 中，alignment 的长度之和大于阈值，则认为这个 LIS 是有效的，可以进行合并。与 LIS 中 alignment 的数量互补。
---
time: 2024-06-27
author: pxxiao
version: 5.0
description: 过滤掉 contained alignment -- filter_contained_alignments
'''

import argparse
import os
from collections import defaultdict
from intervaltree import Interval, IntervalTree


def parse_paf_line(line, min_align_length, min_unitig_length):
    ''' 解析 PAF 文件的行 '''
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

    align_type = parts[16][-1]
    if align_type == 'S':
        return None
    elif query_len < min_unitig_length:
        return None
    elif match_len <= min_align_length:
        return None
    else:
        return parts


def filter_contained_alignments(alignments):   # 时间消耗：1m20.517s
    ''' 过滤包含 alignment，只对同一个 query 且在同一条染色体上的 alignments 进行处理 '''
    n = len(alignments)
    to_remove = set()
    for i in range(n):
        for j in range(n):
            if i != j and alignments[i][0] == alignments[j][0] and alignments[i][5] == alignments[j][5]:
                # 检查对齐i是否被对齐j包含
                if (int(alignments[j][7]) <= int(alignments[i][7]) and int(alignments[j][8]) >= int(alignments[i][8])):
                    to_remove.add(i)
                # 检查对齐j是否被对齐i包含
                elif (int(alignments[i][7]) <= int(alignments[j][7]) and int(alignments[i][8]) >= int(alignments[j][8])):
                    to_remove.add(j)
    return [alignments[i] for i in range(n) if i not in to_remove]


def read_paf(file_paf, min_align_length, min_unitig_length):
    ''' 读取 PAF 文件，返回所有对齐信息 '''
    alignments = defaultdict(list)
    with open(file_paf, 'r') as f:
        for line in f:
            alignment = parse_paf_line(line, min_align_length, min_unitig_length)
            if alignment:
                query_id = alignment[0]
                alignments[query_id].append(alignment)

    # 对每个 query 的 alignments 进行过滤，移除 contained alignment
    for query_id in list(alignments.keys()):
        filtered = filter_contained_alignments(alignments[query_id])
        if filtered:
            alignments[query_id] = filtered
        else:
            del alignments[query_id]  # 如果过滤后为空，则从字典中移除

    return alignments


# def find_best_reference_for_unitig(unitig_alignments):
#     return max(unitig_alignments.keys(), key=lambda unitig_id: (len(unitig_alignments[unitig_id]), sum(int(aln[9]) for aln in unitig_alignments[unitig_id])))
def find_best_reference_for_unitig(unitig_alignments, min_match_ratio, min_align_length, min_unitig_length, unitig_id):
    '''
    找到 unitig 的最佳匹配参考序列，并计算出匹配比率
    如果匹配比率低于阈值，则返回 None
    '''
    reference_stats = defaultdict(lambda: [0, 0])   # 每个 reference 的对齐数量和总匹配长度
    unitig_len = int(unitig_alignments[0][1])       # unitig 的长度

    if unitig_len < min_unitig_length:
        return None, None

    for alignment in unitig_alignments:
        match_len = int(alignment[9])
        if match_len < min_align_length:
            continue
        reference_id = alignment[5]
        reference_stats[reference_id][0] += 1           # 记录对齐数量
        reference_stats[reference_id][1] += match_len   # 记录总匹配长度

    # 找到比对数量最多的参考序列，如果数量相同，则选择总比对长度最大的那个
    best_reference = max(reference_stats.keys(), key=lambda ref: (reference_stats[ref][0], reference_stats[ref][1]))
    match_ratio = reference_stats[best_reference][1] / unitig_len

    if unitig_id == 'utg000211l':
        print(f'Debugs: {match_ratio}')

    if match_ratio < min_match_ratio:
        return None, None

    return best_reference, match_ratio


def find_lis_with_threshold(arr, min_distance):
    ''' 根据距离阈值，找到最长递增子序列 （LIS） '''
    n = len(arr)
    if n == 0:
        return []

    lis_segments = []
    current_lis = [arr[0]]

    for i in range(1, n):
        current_alignment = arr[i]
        prev_alignment = current_lis[-1]

        if int(current_alignment[7]) - int(prev_alignment[8]) <= min_distance:      # 避免有 alignment 很长，下一个 alignment 包含在 long alignment 里面
            # if current_alignment[7] == 20276257:
            #     print('xxp')
            current_lis.append(current_alignment)
        else:
            print(f'Debugs: current_alignment {current_alignment}')
            lis_segments.append(current_lis)
            current_lis = [current_alignment]

    lis_segments.append(current_lis)

    #
    lis_segments_sorted = sorted(lis_segments,
                                   key=lambda segment: int(segment[-1][7]) if segment and segment[-1] else 0)


    return lis_segments_sorted


def merge_lis_segments(lis_segments, min_lis_size, max_lis_distance):
    ''' 合并满足条件的 LIS 片段 '''
    merged_segments = []
    current_segment = []
    large_segments = []

    for lis in lis_segments:
        if len(current_segment) == 0:
            current_segment = lis
            continue

        last_alignment = current_segment[-1]        # 前面 LIS 的最后一个 alignment
        first_alignment = lis[0]                    # 当前面 LIS 的第一个 alignment
        distance = int(first_alignment[7]) - int(last_alignment[8])

        if len(current_segment) > min_lis_size and len(lis) > min_lis_size and distance < max_lis_distance:
            current_segment.extend(lis)
        else:
            if len(current_segment) >= min_lis_size:
                large_segments.append(current_segment)
            else:
                merged_segments.append(current_segment)
            current_segment = lis

    if current_segment:
        if len(current_segment) >= min_lis_size:
            large_segments.append(current_segment)
        else:
            merged_segments.append(current_segment)

    # for i in merged_segments:
    #     print(f'Debugs: i {i}')

    # 将较小的 LIS 插入到合适的位置
    for small_lis in merged_segments:
        inserted = False
        for i in range(len(large_segments) - 1):
            distance_to_previous = int(small_lis[0][7]) - int(large_segments[i][-1][8])
            distance_to_next = int(large_segments[i + 1][0][7]) - int(small_lis[-1][8])
            if distance_to_previous <= max_lis_distance and distance_to_next <= max_lis_distance:
                large_segments[i].extend(small_lis)
                inserted = True
                break
        if not inserted:
            large_segments.append(small_lis)

    return large_segments


def merge_lis_segments1(lis_segments, min_lis_size, min_lis_length, max_lis_distance):
    ''' 合并满足条件的 LIS 片段 '''
    merged_segments = []
    current_segment = []

    for lis in lis_segments:
        # print(f'Debugs {len(lis)}')
        if not current_segment:
            current_segment = lis
            continue

        last_alignment = current_segment[-1]
        first_alignment = lis[0]

        distance = abs(int(first_alignment[7]) - int(last_alignment[8]))
        # print('xpx', last_alignment[8], first_alignment[7], distance)

        total_length = sum(int(align[9]) for align in current_segment)  # Calculate total length of alignments

        # if distance <= max_lis_distance and len(current_segment) >= min_lis_size:
        if distance <= max_lis_distance and (len(current_segment) >= min_lis_size or total_length >= min_lis_length):
            current_segment.extend(lis)
        else:
            merged_segments.append(current_segment)
            current_segment = lis

    if current_segment:
        merged_segments.append(current_segment)

    # 处理小的 LIS 并将它们合并到合适的较大段中
    final_segments = []
    for segment in merged_segments:
        if len(segment) >= min_lis_size:
            final_segments.append(segment)
        else:
            if not final_segments:
                final_segments.append(segment)
            else:
                last_large_segment = final_segments[-1]
                distance_to_last_large = int(segment[0][7]) - int(last_large_segment[-1][8])

                if distance_to_last_large <= max_lis_distance:
                    last_large_segment.extend(segment)
                else:
                    final_segments.append(segment)
    # 排序
    for segment in final_segments:
        print(f'debugs: xpx segment {segment}')

    # for segment in final_segments:
        # 对segment内的alignment按终止坐标排序
        # segment.sort(key=lambda alignment: int(alignment[8]))
    return final_segments


def filter_alignments_by_lis(alignments, min_distance, min_lis_size, min_lis_length, max_lis_distance):
    ''' 根据最小距离划分 alignment 为 LIS '''
    ## 划分 LIS，两个相邻的 LIS 的距离阈值设置为：min_distance
    lis_segments = find_lis_with_threshold(alignments, min_distance)
    merged_segments = merge_lis_segments1(lis_segments, min_lis_size, min_lis_length, max_lis_distance)
    return merged_segments


def find_best_lis(lis_segments):
    ''' 找到匹配长度之和最长的 LIS '''
    ## 匹配长度之和最长的 LIS 为最佳 LIS
    best_lis = max(lis_segments, key=lambda lis: (len(lis), sum(int(aln[9]) for aln in lis)))
    return best_lis


def write_lis_to_file(lis_segments, output_file, best_lis_output_file):
    ''' 将 LIS 写入文件 '''
    with open(output_file, 'a+') as f:
        for i, lis in enumerate(lis_segments):
            # print(f'Debugs: lis_segments {lis_segments}')
            lis_start = lis[0][7]
            lis_end = lis[-1][8]
            lis_length = len(lis)
            f.write(f"LIS segment {i + 1}:\n")
            f.write(f"Range: {lis_start} - {lis_end}\n")
            f.write(f"Number of alignments: {lis_length}\n")
            f.write("Alignments:\n")
            for aln in lis:
                f.write(f"{aln}\n")
            f.write("\n")

    best_lis = find_best_lis(lis_segments)
    with open(best_lis_output_file, 'a+') as f2:
        for aln in best_lis:
            f2.write("\t".join(aln) + "\n")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--paf', required=True, type=str, help='The paf file')
    parser.add_argument('--min_align_length', type=int, default=1000, help="Minimum alignment length [1000]")
    parser.add_argument("--min_unitig_length", type=int, default=100000, help="Minimum unitig length [100000]")
    parser.add_argument('--min_alignment_distance', required=False, type=int, default=500000, help='The minimum distance between two alignments [500000]')
    parser.add_argument('--min_match_ratio', type=float, default=0.05, help='Minimum match ratio [0.05]')
    parser.add_argument('--min_lis_size', type=int, default=5, help='Minimum number of alignments in a LIS to be considered for mergeing [5]')
    parser.add_argument('--min_lis_length', type=int, default=1000000, help='Minimum length of alignments in a LIS to be considered for mergeing [1000000]')
    parser.add_argument('--max_lis_distance', type=int, default=3000000, help='Maximum distance between LIS to be considered for merge [3000000]')
    parser.add_argument('--best_lis_output', required=True, type=str, help='Output file for the best LIS in PAF format')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    # 示例使用
    file_paf = args.paf  # 替换为实际的 PAF 文件路径
    min_distance = args.min_alignment_distance
    min_align_length = args.min_align_length        # 最小的 alignment 大小
    min_unitig_length = args.min_unitig_length      # 最小的 unitig 长度
    min_match_ratio = args.min_match_ratio          # 最小的 match ratio
    min_lis_size = args.min_lis_size                # LIS 包含最少得 alignments 个数 （或者）
    min_lis_length = args.min_lis_length            # LIS 包含最小得 alignments 长度 （或者）
    max_lis_distance = args.max_lis_distance        # 两个 LIS 之间合并的最大距离
    output_file = 'lis_segments.txt'  # 输出文件路径
    best_lis_output_file = args.best_lis_output
    wd = os.path.dirname(best_lis_output_file)
    output_file = os.path.join(wd, output_file)

    ### step1: 解析 PAF 文件
    unitig_alignments = read_paf(file_paf, min_align_length, min_unitig_length)
    # print(f'Debugs: unitig_alignments {unitig_alignments}')

    ### step2: 寻找 unitig 的最佳匹配 reference 染色体
    filtered_alignments = defaultdict(list)
    for unitig_id, alignments in unitig_alignments.items():
        best_reference, match_ratio = find_best_reference_for_unitig(alignments, min_match_ratio, min_align_length, min_unitig_length, unitig_id)
        # print(f'Debugs: best_reference match_ratio {best_reference} {match_ratio}')
        if best_reference:
            filtered_alignments[unitig_id] = [aln for aln in alignments if aln[5] == best_reference]
        if unitig_id == 'utg000211l':
            print(f'utg000211l {best_reference} {match_ratio}')
        if unitig_id == 'utg000613l':
            print(f'utg000613l {best_reference} {match_ratio}')
    if 'utg000211l' in filtered_alignments:
        print(f'nihc utg000211l')
    print(f'Debugs: filtered_alignments {filtered_alignments}')

    ### step3: 寻找 unitig 的最佳 LIS
    if os.path.exists(best_lis_output_file):
        os.remove(best_lis_output_file)
    if os.path.exists(output_file):
        os.remove(output_file)

    for unitig_id, alignments in filtered_alignments.items():
        # Filter out contained alignments before finding LIS
        lis_segments = filter_alignments_by_lis(alignments, min_distance, min_lis_size, min_lis_length, max_lis_distance)
        write_lis_to_file(lis_segments, output_file, best_lis_output_file)


if __name__ == '__main__':
    main()