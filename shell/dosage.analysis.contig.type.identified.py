#!/usr/bin/env python
'''
time: 2024-04-27
author: pxxiao
version: 1.0
description: 根据 reads 深度，将 contig 划分为不同类型的 contig

haplotig [15, 44]
diplotig [45, 73]
triplotig [74, 102]
tetraplotig [103, 131]
replotig > 132

输入文件：cnv_winsize10000_step10000_hq.txt
utg000001l      1       10000   1       6       1.6084  1301636
utg000001l      10001   20000   1       45      13.1233 1301636
utg000001l      20001   30000   1       106     30.727  1301636
utg000001l      30001   40000   1       87      25.3734 1301636
utg000001l      40001   50000   1       63      18.182  1301636
utg000001l      50001   60000   1       67      19.3853 1301636
----
time: 2024-06-24
author: pxxiao
version: 2.0
description: 在浮点数的比较那里出现了问题，将浮点数四舍五入，保留两位小数。
'''


import argparse
import pandas as pd
import subprocess


def parse_args():
    parser = argparse.ArgumentParser('')
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--pandepth', action='store_true', help='input from pandepth or not')
    args = parser.parse_args()
    return args


def classify_contig_type(average_depth):
    # 四舍五入到两位小数
    average_depth = round(average_depth, 2)
    if 15 <= average_depth < 45:
        return 'haplotig'
    elif 45 <= average_depth < 74:
        return 'diplotig'
    elif 74 <= average_depth < 103:
        return 'triplotig'
    elif 103 <= average_depth < 132:
        return 'tetraplotig'
    elif average_depth >= 132:
        return 'replotig'
    else:
        return 'other'


def classify_contig_type_SD(average_depth):
    # 四舍五入到两位小数
    average_depth = round(average_depth, 2)
    if 17 <= average_depth < 49:
        return 'haplotig'
    elif 49 <= average_depth < 81:
        return 'diplotig'
    elif 81 <= average_depth < 113:
        return 'triplotig'
    elif 113 <= average_depth < 145:
        return 'tetraplotig'
    elif 145 <= average_depth < 177:
        return 'pentaplotig'
    elif 177 <= average_depth < 209:
        return 'hexaplotig'
    elif average_depth >= 209:
        return 'replotig'
    else:
        return 'other'


def classify_contig_type_SD_base21(average_depth):
    # 四舍五入到两位小数
    average_depth = round(average_depth, 2)
    if 11 <= average_depth < 32:
        return 'haplotig'
    elif 32 <= average_depth < 53:
        return 'diplotig'
    elif 53 <= average_depth < 74:
        return 'triplotig'
    elif 74 <= average_depth < 95:
        return 'tetraplotig'
    elif 95 <= average_depth < 116:
        return 'pentaplotig'
    elif 116 <= average_depth < 137:
        return 'hexaplotig'
    elif average_depth >= 137:
        return 'replotig'
    else:
        return 'other'


def classify_contig_normal(average_depth, base_value=28):
    # 四舍五入到两位小数
    average_depth = round(average_depth, 2)

    if 0.5 * base_value <= average_depth < 1.5 * base_value:
        return 'haplotig'
    elif 1.5 * base_value <= average_depth < 2.5 * base_value:
        return 'diplotig'
    elif 2.5 * base_value <= average_depth < 3.5 * base_value:
        return 'triplotig'
    elif 3.5 * base_value <= average_depth < 4.5 * base_value:
        return 'tetraplotig'
    elif average_depth >= 4.5 * base_value:
        return 'replotig'
    else:
        return 'other'


def main():
    args = parse_args()

    ## 获得每个 contig 的平均深度
    # 读取数据文件
    data = pd.read_csv(args.input_file, sep='\t', header=None)

    # 按照 contig ID 分组计算平均深度
    if args.pandepth:
        contig_depth = data.groupby(0)[7].mean().reset_index()
    else:
        contig_depth = data.groupby(0)[5].mean().reset_index()

    # 为结果添加列名
    contig_depth.columns = ['contig_ID', 'average_depth']

    # 添加 contig 类型列
    contig_depth['contig_type'] = contig_depth['average_depth'].apply(classify_contig_normal)

    # 打印结果
    print(contig_depth)

    # 将 DataFrame 以文件的形式输出
    contig_depth.to_csv("contig_depth.csv", index=False)
    contig_depth.to_csv("contig_depth.txt", sep="\t", index=False)

    #

    return None


if __name__ == '__main__':
    main()