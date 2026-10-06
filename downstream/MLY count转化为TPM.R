####MLY count转化为TPM




if(!require(GenomicFeatures))BiocManager::install("GenomicFeatures")
if(!require(GenomicFeatures))BiocManager::install("DGEobj.utils")
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("D:/Bioinfo/data/gencode.v36.annotation.gtf",format="gtf")
# 获取每个基因id的外显子数据
exons.list.per.gene <- exonsBy(txdb,by="gene")
# 对于每个基因，将所有外显子减少成一组非重叠外显子，计算它们的长度(宽度)并求和
exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))
# 得到geneid和长度数据
gfe <- data.frame(gene_id=names(exonic.gene.sizes),
                  length=exonic.gene.sizes)
head(gfe)[1:5,1:2]
#                               gene_id length
# ENSG00000000003.15 ENSG00000000003.15   4536
# ENSG00000000005.6   ENSG00000000005.6   1476
# ENSG00000000419.13 ENSG00000000419.13   1207
# ENSG00000000457.14 ENSG00000000457.14   6883
# ENSG00000000460.17 ENSG00000000460.17   5970
save(gfe,file = "D:/Bioinfo/data/gfe.Rdata")




exp <- read.csv("E:/X101SC19123287-Z01-J001-B1-16/X101SC19123287-Z01-J001-B1-16-result/3.Quant/gene_count.csv", header = T)
str(exp)
rownames(exp) <- exp$gene_id
exp <- exp[, -1]
all_numeric <- all(sapply(exp, is.numeric))


effLen = gfe$length
effLen <- effLen[rownames(exp)]
BiocManager::install("DGEobj.utils")
library(DGEobj.utils)






CC_res <- convertCounts(
  countsMatrix = exp,  # 数值型矩阵
  unit = "TPM",
  geneLength = effLen,  # 基因长度向量
  log = FALSE,
  normalize = "none",
  prior.count = NULL
)
head(CC_res)[1:3,1:3]
