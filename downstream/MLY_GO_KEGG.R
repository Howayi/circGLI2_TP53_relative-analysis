#MLY质谱分析
#go, kegg分析
BiocManager::install("org.Mm.eg.db")
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggridges)
library(data.table)



ms_result <- fread("D:/Bioinfo/MLY/202412联合分析/MS_result.txt", header = T, sep = '\t', data.table = F)

uniqe_peptide = 2
ms_H_result <- (ms_result$`Unique peptides` >= 2)
ms_H_result <- ms_result %>% filter(ms_result$`Unique peptides` >=3)


GO_database <- 'org.Mm.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'mmu' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
gene <- bitr(ms_H_result$`Gene names`,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T )

# 将 GO 富集分析结果转换为数据框格式
GO_msH_results <- as.data.frame(GO)

# 查看前几行，了解数据框结构
head(GO_results)

KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
KEGG_results <- as.data.frame(KEGG)
head(KEGG_results)
#将ENTREZ重转为symbol：
KEGG_diff <- setReadable(KEGG,
                         OrgDb = GO_database,
                         keyType = "ENTREZID")
View(KEGG_diff@result)
#计算Rich Factor（富集因子）：
KEGG_diff2 <- mutate(KEGG_diff,
                     RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#计算Fold Enrichment（富集倍数）：
KEGG_diff2 <- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]
head(KEGG_diff2)
# 导出，改名和排序（p值大小升序）
write.table(GO_results,"D:/Bioinfo/MLY/202412联合分析/MS_newGO.xls",sep="\t",quote=F,col.names = NA)
write.table(KEGG_diff@result,"D:/Bioinfo/MLY/202412联合分析/MS_newKEGG.xls",sep="\t",quote=F,col.names = NA)


#GSEA分析
names(deg_opt) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
deg_opt_merge <- merge(deg_opt,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- deg_opt_merge$Log2FoldChange
names(GSEA_input) = deg_opt_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.05)#GSEA富集分析
GSEA_results <- as.data.frame(GSEA_KEGG@result)
write.table(GSEA_results,"D:/Bioinfo/BRCA/202411/GSEA.xls",sep="\t",quote=F,col.names = NA)



GO_msH_results$Description <- factor(GO_msH_results$Description,
                                     levels = GO_msH_results$Description[order(GO_msH_results$p.adjust)])

# 绘制图形
p_GO_msH_results <- ggplot(GO_msH_results, aes(y = Description, x = Count)) +
  geom_point(aes(size = Count, color = -1 * p.adjust)) + # 使用 p.adjust 保持一致
  scale_color_gradient(high = "red", low = "green") +
  labs(
    color = expression(p.adjust),
    size = "Count",
    x = "Count",
    y = "Pathway",
    title = "MS_GO_result"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),            # 全局文字大小
    plot.title = element_text(size = 16),     # 标题文字大小
    axis.title = element_text(size = 14),     # 坐标轴标题文字大小
    axis.text = element_text(size = 12),      # 坐标轴刻度文字大小
    legend.text = element_text(size = 12),    # 图例文字大小
    legend.title = element_text(size = 14)    # 图例标题文字大小
  )

# 显示图形
p_GO_msH_results

ggsave(filename = "D:/Bioinfo/MLY/Fig/气泡图_GO_msH_results.pdf", plot = p_GO_msH_results , device = "pdf", width = 13, height = 28)






GO_RNA_down_results$Description <- factor(GO_RNA_down_results$Description,
                                     levels = GO_RNA_down_results$Description[order(GO_RNA_down_results$p.adjust)])

p_GO_RNA_down_results <- ggplot(GO_RNA_down_results, aes(y = Description, x = Count)) +
  geom_point(aes(size = Count, color = -1 * p.adjust)) + # 使用 p.adjust 保持一致
  scale_color_gradient(high = "red", low = "green") +
  labs(
    color = expression(p.adjust),
    size = "Count",
    x = "Count",
    y = "Pathway",
    title = "RNA_GO_result"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),            # 全局文字大小
    plot.title = element_text(size = 16),     # 标题文字大小
    axis.title = element_text(size = 14),     # 坐标轴标题文字大小
    axis.text = element_text(size = 12),      # 坐标轴刻度文字大小
    legend.text = element_text(size = 12),    # 图例文字大小
    legend.title = element_text(size = 14)    # 图例标题文字大小
  )

# 显示图形
p_GO_RNA_down_results

ggsave(filename = "D:/Bioinfo/MLY/Fig/气泡图_GO_RNA_down_results.pdf", plot = p_GO_RNA_down_results , device = "pdf", width = 13, height = 45)











































#RNA-DEG GO KEGG分析
ms_result <- fread("D:/Bioinfo/MLY/202412联合分析/MS_result.txt", header = T, sep = '\t', data.table = F)
DEG_limma_voom <- fread("D:/Bioinfo/MLY/202412联合分析/RNA_DEG_LvsC.xls", header = T, sep = '\t', data.table = F)

logFC = 1
P.Value = 0.01
k1 <- (DEG_limma_voom$pvalue < P.Value) & (DEG_limma_voom$log2FoldChange < -logFC)
k2 <- (DEG_limma_voom$pvalue < P.Value) & (DEG_limma_voom$log2FoldChange > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma_voom$change)


RNA_down <- DEG_limma_voom %>% filter(DEG_limma_voom$change != "up")


GO_h_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_h_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
RNA_up_gene <- bitr(RNA_down$gene_name,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_h_database)

GO<-enrichGO( RNA_up_gene$ENTREZID,#GO富集分析
              OrgDb = GO_h_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T )

# 将 GO 富集分析结果转换为数据框格式
GO_RNA_down_results <- as.data.frame(GO)

# 查看前几行，了解数据框结构
head(GO_RNA_down_results)

KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
KEGG_results <- as.data.frame(KEGG)
head(KEGG_results)
#将ENTREZ重转为symbol：
KEGG_diff <- setReadable(KEGG,
                         OrgDb = GO_database,
                         keyType = "ENTREZID")
View(KEGG_diff@result)
#计算Rich Factor（富集因子）：
KEGG_diff2 <- mutate(KEGG_diff,
                     RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#计算Fold Enrichment（富集倍数）：
KEGG_diff2 <- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]
head(KEGG_diff2)
# 导出，改名和排序（p值大小升序）
write.table(GO_results,"D:/Bioinfo/MLY/202412联合分析/MS_newGO.xls",sep="\t",quote=F,col.names = NA)
write.table(KEGG_diff@result,"D:/Bioinfo/MLY/202412联合分析/MS_newKEGG.xls",sep="\t",quote=F,col.names = NA)



write.table(GO_msH_results,"D:/Bioinfo/MLY/202412联合分析/MS_GO_result.xls",sep="\t",quote=F,col.names = NA)
write.table(GO_RNA_down_results,"D:/Bioinfo/MLY/202412联合分析/RNA_GO_result.xls",sep="\t",quote=F,col.names = NA)
