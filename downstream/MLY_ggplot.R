#MLY 火山图
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
# 标签添加

merged_DEG <- fread("D:/Bioinfo/MLY/RNA-seq/merged_DEG_edgeR_cut_lowexpr.txt", header = T, sep = '\t', data.table = F)
# 按需求添加标签

DEG_limma_voom <- merged_DEG[,-1]

DEG_limma_voom <- cbind(DEG_limma_voom, GeneSymbol = DEG_limma_voom[, 1])

# 创建一个新的列label，并初始化为NA
DEG_limma_voom$label <- NA

# 根据symbol的值，为特定基因添加标签信息
DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "GLI2:chr2|121684936|121713006|+")] <- "GLI2:chr2|121684936|121713006|+;circGLI2"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "PGA4")] <- "PGA4"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "GLI2:chr2|121696167|121713006|+")] <- "GLI2:chr2|121696167|121713006|+"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "PGA3")] <- "PGA3"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "GKN1")] <- "GKN1"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "PGA5")] <- "PGA5"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "ATP4A")] <- "ATP4A"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "GIF")] <- "GIF"
DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "ATP4B")] <- "ATP4B"
DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "GLI2")] <- "GLI2"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "CLDN18")] <- "CLDN18"
#DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "KRT20")] <- "KRT20"
DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "MUC4")] <- "MUC4"
DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "TFF1")] <- "TFF1"
DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "TFF2")] <- "TFF2"
DEG_limma_voom$label[which(DEG_limma_voom$GeneSymbol == "MUC2")] <- "MUC2"



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



# 在普通火山图p的基础上，添加标签，并使用geom_label_repel函数进行标签的绘制


p0 <- p +
  geom_label_repel(
    data = DEG_limma_voom,
    aes(label = label),
    size = 5,                           # 设置标签大小
    box.padding = unit(0.5, "lines"),   # 设置标签内边距
    point.padding = unit(0.8, "lines"), # 设置标签与点的距离
    segment.color = "black",            # 设置标签边界线颜色
    show.legend = FALSE,                # 不显示图例
    max.overlaps = 10000                # 设置标签重叠的最大次数
  ) +
  xlab(expression(paste(log[2],"FC", sep = "")))
  theme(
    legend.title = element_text(size = 0),     # 设置图例标题的字体大小
    legend.text = element_text(size = 16),      # 设置图例标签的字体大小
    axis.title.x = element_text(size = 16),     # 设置x轴标题的字体大小
    axis.title.y = element_text(size = 16)      # 设置y轴标题的字体大小
  )

p0

ggsave(filename = "volcano_plot_tag1.pdf", plot = p0, device = "pdf", width = 6, height = 5)





# 差异基因热图_top10 gene

# 筛选上调中显著性top10的基因
up_data <- filter(DEG_limma_voom, change == 'up') %>%  # 从DEG_limma_voom中筛选出上调的基因
  distinct(GeneSymbol, .keep_all = TRUE) %>%               # 去除重复的基因，保留第一个出现的行
  top_n(11, -log10(PValue))                           # 选择-P.Value值最大的前10个基因

up_data <- up_data[-2,]


# 筛选下调中显著性top10的基因
down_data <- filter(DEG_limma_voom, change == 'down') %>%  # 从DEG_limma_voom中筛选出下调的基因
  distinct(GeneSymbol, .keep_all = TRUE) %>%                   # 去除重复的基因，保留第一个出现的行
  top_n(10, -log10(PValue))                               # 选择-P.Value值最大的前10个基因

head(up_data); head(down_data)
rownames(up_data) <- up_data[, 1]  # 基因名作为行名
up_data <- up_data[, -1]  # 删除第一列，只保留数值型数据
rownames(down_data) <- down_data[, 1]  # 基因名作为行名
down_data <- down_data[, -1]  # 删除第一列，只保留数值型数据
merged_top_DEG <- rbind(up_data, down_data)


merged_rna[, 1] <- make.unique(as.character(merged_rna[, 1]))
rownames(merged_rna) <- merged_rna[, 1]  # 基因名作为行名
merged_rna <- merged_rna[, -1]
rownames(merged_top_DEG) <- merged_top_DEG[, 1]  # 基因名作为行名
merged_top_DEG <- merged_top_DEG[, -1]


merged_rna_heatmap <- merged_rna %>% filter(rownames(merged_top_DEG))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(merged_rna_heatmap)

p1 <- pheatmap(up_data, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100))
p1


# 提取只包含已选择基因的表达矩阵
group <- factor(c("normal", "tumor"))
#selected_genes <- colnames(merged_top_DEG)  # 获取在 merged_top_DEG 中选择的基因列名

exp_brca_heatmap <- merged_rna %>% filter(rownames(merged_rna) %in% rownames(merged_top_DEG))
log_data <- log2(exp_brca_heatmap + 1)


annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(log_data)


colnames(log_data)[1] <- "Normal"
colnames(log_data)[2] <- "Tumor"

rownames(log_data)[19] <- "GLI2:chr2|121684936|121713006|+;circGLI2"

# 定义颜色分布的断点
breaks <- seq(0, 17, length.out = 1000)  # 选择一个合理的范围并根据需要修改

p1 <- pheatmap(
  log_data,  # 使用筛选后的矩阵
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
  fontsize = 15,  # 字体大小
  breaks = breaks  # 应用自定义断点，提升阈值范围
)


ggsave(filename = "heatmap_plot.pdf", plot = p1, device = "pdf", width = 10, height = 12)
dev.off()
