# 安装和加载必要的R包
install.packages("gdata")  # 如果尚未安装
install.packages("gplots")  # 如果尚未安装
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomeInfoDbData")
BiocManager::install("DESeq2")
library(gdata)
library(gplots)
library(DESeq2)

# 设置工作路径
setwd("D:\\Gavin's\\Research\\出生队列\\Rawdata\\GTEx_data")
# 读取基因计数数据
# 读取数据文件
data <- read.csv("Countdata.csv", header = TRUE)

# 将第一列设置为行索引
rownames(data) <- data[, 1]
data <- data[, -1]  # 去除原本的第一列
# 转换为矩阵格式
mycounts <- as.matrix(data)
# 显示 count data 的前几行
head(mycounts)

# 设置样本组别和重复数
condition <- factor(c(rep("u1", 71), rep("u2", 71)), levels = c("u1", "u2"))

# 创建 colData
colData <- data.frame(condition = condition)

# 为 colData 设置行名
rownames(colData) <- colnames(mycounts)

# 创建 DESeq2 数据对象
dds <- DESeqDataSetFromMatrix(countData = mycounts, colData = colData, design = ~ condition)
# 检查处理效果
dds
#对原始dds进行normalize
dds <- DESeq(dds)
#显示dds信息
dds
#使用DESeq2包中的results()函数，提取差异分析的结果
#Usage:results(object, contrast, name, .....）
#将提取的差异分析结果定义为变量"res" 
#contrast: 定义谁和谁比较
res = results(dds, contrast=c("condition", "u1", "u2"))
#对结果res利用order()函数按pvalue值进行排序
#创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列
#order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
res = res[order(res$pvalue),]
#显示res结果首信息
head(res)
#对res矩阵进行总结，利用summary命令统计显示一共多少个genes上调和下调
summary(res)
#将分析的所有结果进行输出保存
write.csv(res, file="All_results.csv")
#显示显著差异的数目
table(res$padj<0.05)


