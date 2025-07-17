#!/bin/bash

# My Bash Script
# Author: Pei-xuan Xiao
# Mail: 1107931486@qq.com
# Date: 2023/07/23
# Version: 1.0

# 显示帮助文档
show_help() {
    echo "Usage: $0 -g genome_file -i hifi_file -t threads [-o output_directory]"
    echo ""
    echo "Options:"
    echo "  -g genome_file      Path to the genome file"
    echo "  -i hifi_file        Path to the HiFi fastq file"
    echo "  -t threads          Number of threads to use (default: 30)"
    echo "  -o output_directory Directory to save output files (default: current directory)"
    echo ""
    echo " dosage.analysis.contig.type.identified.py 用来识别 contig 类型 "
    echo " grep -v \"#\" aln.sort.clean.pandepth.win.stat > temp && mv temp aln.sort.clean.pandepth.win.stat"
    echo " python dosage.analysis.contig.type.identified.py --input_file aln.sort.clean.pandepth.win.stat.txt --pandepth"
}

# 默认设置
output_directory="."
threads=30

# 获取脚本自身目录路径
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rscript_path="$script_dir/4p_rscript_dosage_reads_depth_no_text.R"

# 解析命令行参数
while getopts "hg:i:t:o:" opt; do
    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        g)
            genome=$OPTARG
            ;;
        i)
            hifi=$OPTARG
            ;;
        t)
            threads=$OPTARG
            ;;
        o)
            output_directory=$OPTARG
            ;;
        \?)
            show_help
            exit 1
            ;;
    esac
done

# 检查必需的参数是否已提供
if [ -z "$genome" ] || [ -z "$hifi" ]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# 创建输出目录
mkdir -p "$output_directory"

# HiFi mapping to genome
minimap2 -ax map-hifi -t $threads -N 1 --secondary=no -o "$output_directory/aln.sam" $genome $hifi

samtools view -@ $threads -bS "$output_directory/aln.sam" > "$output_directory/aln.bam"
rm "$output_directory/aln.sam"

samtools sort -@ $threads -o "$output_directory/aln.sort.bam" "$output_directory/aln.bam"
rm "$output_directory/aln.bam"

samtools view -h -@ $threads -F 3840 -bS "$output_directory/aln.sort.bam" | samtools sort -@ $threads -o "$output_directory/aln.sort.clean.bam" -
rm "$output_directory/aln.sort.bam"

# 计算Depth，窗口设置为 10kb
nohup time -v pandepth -i "$output_directory/aln.sort.clean.bam" -w 10000 -a -o "$output_directory/aln.sort.clean.pandepth" -t $threads > "$output_directory/log_pandepth_out" 2> "$output_directory/log_pandepth_err"

zcat "$output_directory/aln.sort.clean.pandepth.win.stat.gz" | grep -v RegionLength | cut -f 8 | sed 's/MeanDepth/Depth/g' > "$output_directory/dosage.win10000.txt"

$rscript_path -i "$output_directory/dosage.win10000.txt" -o "$output_directory/dosage.win10000"

echo "Done!^-^"