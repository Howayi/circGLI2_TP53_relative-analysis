BiocManager::install("EnhancedVolcano")

library(data.table)  # 用于高效处理大数据集
library(dplyr)       # 用于数据操作和转换
library(ggplot2)     # 画图图
library(pheatmap)    # 绘制热图
library(DESeq2)      # 差异分析一号选手
library(edgeR)       # 差异分析二号选手
library(limma)       # 差异分析三号选手
library(tinyarray)
library(statmod)
getwd("D:/Bioinfo/MLY/RNA-seq/")
# 读取两个 count 数据文件

normal <- fread("D:/Bioinfo/MLY/RNA-seq/normal.txt",header = T, sep = '\t', data.table = F)
tumor <- fread("D:/Bioinfo/MLY/RNA-seq/tumor.txt", header = T, sep = '\t', data.table = F)
mrna <- fread("D:/Bioinfo/MLY/RNA-seq/mrna_counts.txt", header = T, sep = '\t', data.table = F)
lncrna <- fread("D:/Bioinfo/MLY/RNA-seq/lncrna-counts.txt", header = T, sep = '\t', data.table = F)
merged_DEG <- fread("D:/Bioinfo/MLY/RNA-seq/merged_DEG_edgeR.txt", header = T, sep = '\t', data.table = F)
merged_rna <- fread("D:/Bioinfo/MLY/RNA-seq/merged_counts.txt", header = T, sep = '\t', data.table = F)


head(normal)
# 假设 normal 和 tumor 数据框已经存在
colnames(normal) <- c("gene_symbol", "counts_normal")
colnames(tumor) <- c("gene_symbol", "counts_tumor")
common_genes <- intersect(normal$gene_symbol, tumor$gene_symbol)

# 从 normal 和 tumor 中只保留共同基因
normal_filtered <- normal %>% filter(gene_symbol %in% common_genes)
tumor_filtered <- tumor %>% filter(gene_symbol %in% common_genes)

# 按基因符号合并两个数据框
merged_data <- merge(normal_filtered, tumor_filtered, by = "gene_symbol")

# 查看合并后的数据
head(merged_data)

# 使用全外连接（full join）合并两个表达矩阵
merged_data_O <- full_join(normal, tumor, by = "gene symbol", suffix = c("_normal", "_tumor"))

# 如果有 NA，替换为0
merged_data_O[is.na(merged_data_O)] <- 0

# 查看合并后的数据
head(merged_data_O)


# 创建表达矩阵 (只保留 counts 列)
expr_matrix <- as.matrix(merged_data[, c("counts_normal", "counts_tumor")])
rownames(expr_matrix) <- merged_data$`gene symbol`

# 创建设计矩阵
group <- factor(c("normal", "tumor"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

y <- DGEList(counts = merged_rna, group = group)
merged_data_standard <- calcNormFactors(merged_data)


merged_rna <- distinct(merged_rna, GeneSymbol, .keep_all = T)


# 计算归一化因子 (TMM normalization)
exprSet <- calcNormFactors(y)



bcv = 0.4#设置bcv为0.1
et <- exactTest(exprSet, dispersion=bcv^2)
allrna_edgeR=as.data.frame(topTags(et, n = nrow(exprSet$counts)))
head(allrna_edgeR)

write.csv(allrna_edgeR, 'merged_DEG_edgeR_cut_lowexpr.csv', quote = FALSE)
merged_edger_mc <- rbind(mrnaDEG_edgeR, DEG_edgeR)
write.csv(DEG_edgeR, 'D:/Bioinfo/MLY/RNA-seq/circDEG_edgeR.csv', quote = FALSE)

logFC = 2
P.Value = 0.01
k1 <- (allrna_edgeR$PValue < P.Value) & (allrna_edgeR$logFC < -logFC)
k2 <- (allrna_edgeR$PValue < P.Value) & (allrna_edgeR$logFC > logFC)
allrna_edgeR <- mutate(allrna_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(allrna_edgeR$change)
# down stable     up
#  634  30531   1904
deg_opt <- DEG_edgeR %>% filter(DEG_edgeR$change != "stable")

merged_data_heatmap <- merged_data %>% filter(rownames(merged_data) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(merged_data_heatmap)

p1 <- pheatmap(merged_data_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100))
p1

ggsave(filename = "./figure/heatmap_plot_edgeR.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()





# 火山图
logFC = 2
P.Value = 0.01
p <- ggplot(data = merged_DEG,
            aes(x = logFC,
                y = -log10(PValue))) +
  geom_point(alpha = 0.4, size = 3.5,
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("#4DBBD5", "grey", "#F39B7F"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./figure/volcano_plot_edgeR.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()









# 差异基因热图


# 将第一列设为行名，并删除第一列
# 使用 make.unique() 去重基因名
merged_rna_heatmap[, 1] <- make.unique(as.character(merged_rna_heatmap[, 1]))

rownames(merged_rna_heatmap) <- merged_rna_heatmap[, 1]  # 基因名作为行名
merged_rna_heatmap <- merged_rna_heatmap[, -1]  # 删除第一列，只保留数值型数据

pheatmap(merged_rna_heatmap)



deg_opt <- allrna_edgeR %>% filter(allrna_edgeR$change != "stable")
merged_rna_heatmap <- merged_rna %>% filter(rownames(merged_rna) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(merged_rna_heatmap)

p1 <- pheatmap(merged_rna, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-5, 5, length.out = 100))
p1

library(dplyr)

# 筛选出 pvalue < 0.01 的基因，去除稳定的基因
deg_opt <- allrna_edgeR %>%
  filter(change != "stable" & pvalue < 0.01)

# 找出 logFC 最大的5个基因（上调基因）
top_upregulated <- allrna_edgeR %>%
  arrange(desc(logFC)) %>%
  head(5)

# 找出 logFC 最小的5个基因（下调基因）
top_downregulated <- allrna_edgeR %>%
  arrange(logFC) %>%
  head(5)

# 合并上下调基因
top_genes <- bind_rows(top_upregulated, top_downregulated)

# 从表达矩阵中筛选出这些基因
merged_rna_heatmap <- merged_rna %>%
  filter(rownames(merged_rna) %in% rownames(top_genes))

# 绘制热图

rownames(merged_rna_heatmap) <- merged_rna_heatmap[, 1]  # 基因名作为行名
merged_rna_heatmap <- merged_rna_heatmap[, -1]  # 删除第一列，只保留数值型数据

merged_rna_heatmap <- fread("D:/Bioinfo/MLY/RNA-seq/log_exp_brca_top20_heatmap.txt")


p1 <- pheatmap(
  merged_rna_heatmap,  # 使用筛选后的矩阵
  angle_col = 45,
  cellwidth = 50,
  cellheight = 20,
  color = colorRampPalette(c('blue','white','red'))(500),  # 色块颜色从蓝到红，分为100个等级
  border_color = "black",  # 色块的边框颜色
  scale = "row",  # 按行进行归一化
  cluster_rows = TRUE,  # 对行聚类
  cluster_cols = FALSE,  # 不对列进行聚类
  legend = TRUE,  # 显示图例
  legend_breaks = c(-1, 0, 1),  # 图例的断点
  legend_labels = c("low", "", "high"),  # 图例的标签
  show_rownames = TRUE,  # 显示行名
  show_colnames = TRUE,  # 显示列名
  fontsize = 12  # 字体大小
)














write.table(exp_brca_heatmap,"D:/Bioinfo/MLY/RNA-seq/exp_brca_top20_heatmap.txt", sep="\t",quote=F,col.names = NA)
write.table(log_data,"D:/Bioinfo/MLY/RNA-seq/log_exp_brca_top20_heatmap.txt", sep="\t",quote=F,col.names = NA)
write.table(merged_data,"D:/Bioinfo/MLY/RNA-seq/merged_data.csv", sep="\t",quote=F,col.names = NA)
write.table(merged_data,"D:/Bioinfo/MLY/RNA-seq/merged_data.csv", sep="\t",quote=F,col.names = NA)








#=================解决excel自动转换基因名的问题/热图绘制======

merged_data <- fread("D:/Bioinfo/MLY/RNA-seq/merged_counts.txt", header = T, sep = '\t', data.table = F)

# 查看合并后的数据
head(merged_data_O)


# 创建表达矩阵 (只保留 counts 列)
expr_matrix <- as.matrix(merged_data[, c("7N", "15T")])
rownames(expr_matrix) <- merged_data$`gene symbol`

# 创建设计矩阵
group <- factor(c("normal", "tumor"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

y <- DGEList(counts = merged_data, group = group)
merged_data_standard <- calcNormFactors(merged_data)


merged_rna <- distinct(merged_rna, GeneSymbol, .keep_all = T)


# 计算归一化因子 (TMM normalization)
exprSet <- calcNormFactors(y)



bcv = 0.4#设置bcv为0.1
et <- exactTest(exprSet, dispersion=bcv^2)
allrna_edgeR=as.data.frame(topTags(et, n = nrow(exprSet$counts)))
head(allrna_edgeR)

write.csv(allrna_edgeR, 'D:/Bioinfo/MLY/RNA-seq/merged_DEG_edgeR_cut_lowexpr.csv', quote = FALSE)
merged_edger_mc <- rbind(mrnaDEG_edgeR, DEG_edgeR)

logFC = 2
P.Value = 0.01
k1 <- (allrna_edgeR$PValue < P.Value) & (allrna_edgeR$logFC < -logFC)
k2 <- (allrna_edgeR$PValue < P.Value) & (allrna_edgeR$logFC > logFC)
allrna_edgeR <- mutate(allrna_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(allrna_edgeR$change)
# down stable     up
#  634  30531   1904


deg_opt <- allrna_edgeR %>% filter(allrna_edgeR$change != "stable")

write.csv(deg_opt, 'D:/Bioinfo/MLY/RNA-seq/DEGs_nostable.csv', quote = FALSE)




























merged_rna_heatmap <- fread("D:/Bioinfo/MLY/RNA-seq/log_exp_brca_top20_heatmap.txt", header = T, sep = '\t', data.table = F)
# 将第一列设为行名，并删除第一列
# 使用 make.unique() 去重基因名
merged_rna_heatmap[, 1] <- make.unique(as.character(merged_rna_heatmap[, 1]))

rownames(merged_rna_heatmap) <- merged_rna_heatmap[, 1]  # 基因名作为行名
merged_rna_heatmap <- merged_rna_heatmap[, -1]  # 删除第一列，只保留数值型数据

pheatmap(merged_rna_heatmap)



deg_opt <- allrna_edgeR %>% filter(allrna_edgeR$change != "stable")
merged_rna_heatmap <- merged_rna %>% filter(rownames(merged_rna) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(merged_rna_heatmap)

p1 <- pheatmap(merged_rna, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-5, 5, length.out = 100))
p1

library(dplyr)

# 筛选出 pvalue < 0.01 的基因，去除稳定的基因
deg_opt <- allrna_edgeR %>%
  filter(change != "stable" & pvalue < 0.01)

# 找出 logFC 最大的5个基因（上调基因）
top_upregulated <- allrna_edgeR %>%
  arrange(desc(logFC)) %>%
  head(5)

# 找出 logFC 最小的5个基因（下调基因）
top_downregulated <- allrna_edgeR %>%
  arrange(logFC) %>%
  head(5)

# 合并上下调基因
top_genes <- bind_rows(top_upregulated, top_downregulated)

# 从表达矩阵中筛选出这些基因
merged_rna_heatmap <- merged_rna %>%
  filter(rownames(merged_rna) %in% rownames(top_genes))

# 绘制热图

rownames(merged_rna_heatmap) <- merged_rna_heatmap[, 1]  # 基因名作为行名
merged_rna_heatmap <- merged_rna_heatmap[, -1]  # 删除第一列，只保留数值型数据

merged_rna_heatmap <- fread("D:/Bioinfo/MLY/RNA-seq/log_exp_brca_top20_heatmap.txt")

breaks <- seq(0, 17, length.out = 1000)  # 选择一个合理的范围并根据需要修改

p1 <- pheatmap(
  merged_rna_heatmap,  # 使用筛选后的矩阵
  angle_col = 45,
  cellwidth = 62,
  cellheight = 16,
  color = colorRampPalette(c('#63ADEE','white','#F5867F'))(1000),  # 色块颜色从蓝到红，分为500个等级
  border_color = "#999999",  # 色块的边框颜色

  cluster_rows = TRUE,  # 对行聚类
  cluster_cols = TRUE,  # 不对列进行聚类
  treeheight_col = 10,
  legend = TRUE,  # 显示图例
  legend_breaks = c(0, 8, 16),  # 图例的断点
  legend_labels = c("0", "8", "16"),  # 图例的标签
  show_rownames = TRUE,  # 显示行名
  show_colnames = TRUE,  # 显示列名
  fontsize = 15  # 字体大小
  #breaks = breaks  # 应用自定义断点，提升阈值范围
)

p1
