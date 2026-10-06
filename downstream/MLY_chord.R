library(data.table)
library(dplyr)
library(biomaRt)
library(GO.db)
library(circlize)
library(tidyr)

# 读取数据
rna_go <- fread("D:/Bioinfo/MLY/202412联合分析/RNA_GO.xls")    # RNA-seq GO分析结果
ms_go <- read.csv("D:/Bioinfo/MLY/202412联合分析/MS_GO.csv") # 质谱 GO分析结果

# 提取通路名
rna_go_down_terms <- GO_RNA_down_results$ID

ms_go_terms <- GO_msH_results$ID
GO_ms_results <- GO_results
# 找到overlap通路
common_go <- intersect(rna_go_down_terms, ms_go_terms)

# 筛选共同通路的数据
rna_common <- GO_RNA_down_results[GO_RNA_down_results$ID %in% common_go, ]
ms_common <- GO_msH_results[GO_msH_results$ID %in% common_go, ]

# 基因-通路关系

rna_gene_go <- subset(rna_common, select = c(geneID, Description))
ms_gene_go <- subset(ms_common, select = c(geneID, Description))

#rna_common$Description
#rna_gene_go <- rna_common %>%
  select(Description, geneID) %>%
  unique()

#ms_gene_go <- ms_common %>%
  select(geneID, Description) %>%
  unique()

colnames(ms_gene_go)[colnames(ms_gene_go) == "geneID"] <- "geneName"
colnames(ms_gene_go)[colnames(ms_gene_go) == "ID"] <- "GOID"

# 转换鼠基因名为人基因名
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "asia")
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")



# 安装 homologene 包
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("oganm/homologene")

BiocManager::install("homologene")

# 加载 homologene
library(homologene)

mouse2human(c("Trp53"))

# 拆分基因名字符串
ms_gene_go <- ms_gene_go %>%
  separate_rows(geneName, sep = "/")

# 查看结果
head(ms_gene_go)


# 映射鼠基因到人基因
mapping <- mouse2human(ms_gene_go$geneName)
ms_gene_go_mapped <- merge(ms_gene_go, mapping, by.x = "geneName", by.y = "mouseGene", all.x = TRUE)

# 筛选成功映射的基因
ms_gene_go_mapped <- ms_gene_go_mapped[!is.na(ms_gene_go_mapped$humanGene), ]

# 更新基因名为人基因名
ms_gene_go$geneName <- ms_gene_go_mapped$humanGene
















#拆分RNA结果
rna_gene_go <- rna_gene_go %>%
  separate_rows(geneName, sep = "/")




# 合并 RNA-seq 和质谱数据
combined_gene_go <- bind_rows(rna_gene_go, ms_gene_go)
relation_matrix <- table(combined_gene_go$geneName, combined_gene_go$GOID)
relation_matrix <- as.matrix(relation_matrix)

# 绘制和弦图
library(circlize)

# 生成和弦图
chordDiagram(relation_matrix,
             transparency = 0.5,
             annotationTrack = c("name", "grid"),
             preAllocateTracks = 1)

# 外环标注
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  # 获取当前的细分名称
  sector.name <- get.cell.meta.data("sector.index")
  # 设置文本标注位置
  circos.text(x = get.cell.meta.data("xcenter"),
              y = get.cell.meta.data("ylim")[2] + 1,
              labels = sector.name,
              facing = "clockwise",
              niceFacing = TRUE,
              adj = c(0, 0.5),
              cex = 0.6)  # 文本大小可调整
}, bg.border = NA)

# 标题
title("Common GO Terms and Related Genes in RNA-seq and Mass Spec")


























#==============================================只show质谱结果===============================================================


ms_result <- fread("D:/Bioinfo/MLY/202412联合分析/MS_result.txt", header = T, sep = '\t', data.table = F)


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







# 读取数据
rna_go <- fread("D:/Bioinfo/MLY/202412联合分析/RNA_GO.xls")    # RNA-seq GO分析结果
ms_go <- read.csv("D:/Bioinfo/MLY/202412联合分析/MS_GO.csv") # 质谱 GO分析结果

# 提取通路名
rna_go_down_terms <- GO_RNA_down_results$ID

ms_go_terms <- GO_msH_results$ID

# 找到overlap通路
common_go <- intersect(rna_go_down_terms, ms_go_terms)

# 筛选共同通路的数据
rna_common <- GO_RNA_down_results[GO_RNA_down_results$ID %in% common_go, ]
ms_common <- GO_msH_results[GO_msH_results$ID %in% common_go, ]

# 基因-通路关系

rna_gene_go <- subset(rna_common, select = c(geneID, Description))
ms_gene_go <- subset(ms_common, select = c(geneID, Description))


#ms_gene_go <- read.csv("D:/Bioinfo/MLY/202412联合分析/GO_overlapping_in MS partial.csv")
ms_gene_go <- ms_gene_go %>%
  separate_rows(geneID, sep = "/")

relation_matrix1 <- table(ms_gene_go$geneID, ms_gene_go$Description)
relation_matrix1 <- as.matrix(relation_matrix1)

# 绘制和弦图
library(circlize)



# 111设置 PDF 输出
circos.clear()
pdf("chord_diagram_output.pdf", width = 7, height = 7.2)


# 绘制和弦图
chordDiagram(relation_matrix1,
             transparency = 0.5,
             annotationTrack = "grid",  # 仅保留外圈注释
             preAllocateTracks = 1)

# 添加外圈标注
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.par(start.degree = -90)
  circos.text(x = get.cell.meta.data("xcenter"),
              y = get.cell.meta.data("ylim")[1] + 1.35,
              labels = sector.name,
              facing = "bending",
              niceFacing = TRUE,
              adj = c(0.5, 0.5),
              font = 2,
              cex = 1.2)  # 调整文字大小
}, bg.border = NA)

# 添加标题
#title("Common GO Terms and Related Genes in RNA-seq and Mass Spec")

# 关闭 PDF
dev.off()



#==========================韦恩图======================
library(VennDiagram)
common_go <- intersect(rna_go_down_terms, ms_go_terms)

rna_go_down_terms <- as.data.frame(rna_go_down_terms)
ms_go_terms <- as.data.frame(ms_go_terms)

p1 <- venn.diagram(x=list(rna_go_down_terms,ms_go_terms),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('black','black'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
              label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色

             cex = 4, # 数字大小

             fontface = "bold",  # 字体粗细；加粗bold

             fill=c('#E64B35','#4DBBD5'), # 填充色 配色https://www.58pic.com/

             category.names = c("RNA-seq_GO", "MS_GO"),

             cat.dist = 0.02, # 标签距离圆圈的远近

             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)

             cat.cex = 2.75, #标签字体大小

             cat.fontface = "bold",  # 标签字体加粗

             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色

             cat.default.pos = "outer",  # 标签位置, outer内;text 外

             #output=TRUE,

             filename= NULL,# 文件保存

             #imagetype="png",  # 类型（tiff png svg）

             #resolution = 400,  # 分辨率

            # compression = "lzw"# 压缩算法

)

p1

ggsave("D:/Bioinfo/MLY/202412联合分析/GO_OVERLAP.pdf", plot = p1, device = "pdf", width = 8, height = 8)


grid.draw(data)
