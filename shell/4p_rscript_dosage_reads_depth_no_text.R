#!/usr/bin/env Rscript

# Author: Pei-xuan Xiao
# Mail: 1107931486@qq.com
# Date: 2024/04/26
# Version: 1.0

# 导入必要的包
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(egg))
suppressPackageStartupMessages(library(dplyr))

# 定义参数选项
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input file containing depth data [Required]",
              dest="input_file"),
  make_option(c("-o", "--output"), type="character", default="plot.pdf",
              help="Output filename for the plot [Default: plot.pdf]",
              dest="output_file")
)

# 解析命令行参数
options(error=traceback)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opts = parse_args(parser)

# 检查是否提供了输入文件
if(is.null(opts$input_file)){
  cat("Error: Missing input file containing depth data!\n")
  print_help(parser)
  quit()
}

# 读取深度数据
depth <- read.table(opts$input_file, header=TRUE)
depth$Depth <- round(depth$Depth)
# 获取输入文件的基本名称
input_basename <- tools::file_path_sans_ext(basename(opts$input_file))
# 构造新文件名
new_filename <- paste0(input_basename, ".round.txt")
# 将修改后的数据写入新文件
write.table(depth, file = new_filename, sep = "\t", quote = FALSE, row.names = FALSE)
# 打印新文件名
cat(paste("Modified data saved as", new_filename, "\n"))

# 过滤深度大于 1000 的数据
# depth_filtered <- subset(depth, Depth <= 1000)

# 生成深度和数量的统计表
depth_summary <- depth %>%
  group_by(Depth) %>%
  summarise(Count = n())

# 将统计表保存为文件
depth_summary_file <- paste0(opts$output_file, ".depth_summary.txt")
write.table(depth_summary, file = depth_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("Depth summary saved as", depth_summary_file, "\n"))

# 创建 ggplot2 图形对象并添加密度曲线和频率曲线
p <- ggplot(data = depth, aes(x = Depth)) +
  # stat 计算每个类别的观察次数，并将其显示为条形的高度
  geom_bar(fill="lightblue", stat = "count") + 
  # egg 的一个主题格式
  theme_article() + 
  # 显示 x 轴的范围
  xlim(0,150) + 
  # 横纵坐标标签
  labs(x="Sequencing depth", y="Count") + 
  # haplotig
#  geom_segment(aes(x=29, y=0, xend=29, yend=13000), linetype="dashed", color="red") + 
#   geom_segment(aes(x=29, y=1300, xend=29, yend=13200), color="red") +
#  geom_text(aes(x=29, y = 13200, label = "haplotig (29X)"), vjust = 0.5, size=4) +
  # diplotig
#  geom_segment(aes(x=58, y=0, xend=58, yend=2000), linetype="dashed", color="red") + 
#   geom_segment(aes(x=58, y=1300, xend=58, yend=1520), color="red") +
#  geom_text(aes(x=58, y = 2200, label = "diplotig (58X)"), vjust = 0.5, size=4) +
  # triplotig
#  geom_segment(aes(x=87, y=0, xend=87, yend=1300), linetype="dashed", color="red") + 
#   geom_segment(aes(x=87, y=1300, xend=87, yend=1320), color="red") +
#  geom_text(aes(x=87, y = 1500, label = "triplotig (87X)"), vjust = 0.5, size=4) +
  # tetraplotig
#  geom_segment(aes(x=116, y=0, xend=116, yend=800), linetype="dashed", color="red") + 
#   geom_segment(aes(x=116, y=1300, xend=116, yend=820), color="red") +
#  geom_text(aes(x=116, y = 1000, label = "tetraplotig (116X)"), vjust = 0.5, size=4) + 
  # 修改横纵坐标字体
  theme(axis.text.x = element_text(color = "black"),  # 设置横坐标标签颜色为黑色
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"))  # 设置纵坐标标签颜色为黑色

# 保存图形到文件，分别保存为 PDF 和 PNG 格式
ggsave(filename = paste0(opts$output_file, ".pdf"), plot = p, width = 6, height = 4, device = "pdf")
ggsave(filename = paste0(opts$output_file, ".png"), plot = p, width = 6, height = 4, device = "png")

cat(paste("Plot saved as", opts$output_file, "\n"))
