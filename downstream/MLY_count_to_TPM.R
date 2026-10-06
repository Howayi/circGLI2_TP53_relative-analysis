####MLY count转化为TPM





# 加载必要的R包
library(GenomicFeatures)
library(dplyr)
library(readr)
library(DGEobj.utils)

# 步骤1：生成基因长度数据
txdb <- makeTxDbFromGFF("D:/Bioinfo/data/gencode.v36.annotation.gtf", format = "gtf")
exons.list.per.gene <- exonsBy(txdb, by = "gene")
exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))
gfe <- data.frame(gene_id = names(exonic.gene.sizes), length = exonic.gene.sizes)
save(gfe, file = "D:/Bioinfo/data/gfe.Rdata")

# 步骤2：读取count数据
exp <- read.csv("E:/X101SC19123287-Z01-J001-B1-16/X101SC19123287-Z01-J001-B1-16-result/3.Quant/gene_count.csv", header = T)
str(exp)
rownames(exp) <- exp$gene_id
exp <- exp[, -1]


all_numeric <- all(sapply(exp, is.numeric))
if (!all_numeric) {
  stop("Error: Some columns in exp are not numeric!")
}

# 步骤3：处理基因ID匹配问题
# 移除gfe$gene_id中的版本号（如果存在）
gfe$gene_id_clean <- sub("\\.[0-9]+$", "", gfe$gene_id)  # 移除.15等版本号
exp$gene_id_clean <- sub("\\.[0-9]+$", "", rownames(exp))  # 对exp的gene_id也做同样处理

# 取交集，确保exp和gfe的基因ID匹配
common_genes <- intersect(gfe$gene_id_clean, exp$gene_id_clean)
if (length(common_genes) == 0) {
  stop("Error: No common genes found between exp and gfe!")
}

# 过滤exp和gfe，只保留共有基因
exp <- exp[exp$gene_id_clean %in% common_genes, , drop = FALSE]
gfe <- gfe[gfe$gene_id_clean %in% common_genes, ]

# 确保effLen与exp的行名对齐
effLen <- gfe$length
names(effLen) <- gfe$gene_id_clean
effLen <- effLen[exp$gene_id_clean]

# 验证长度是否匹配
if (length(effLen) != nrow(exp) || any(is.na(effLen))) {
  stop("Error: effLen length does not match exp rows or contains NA values!")
}

# 验证effLen是否为正数值
if (any(effLen <= 0)) {
  stop("Error: effLen contains non-positive values!")
}


exp1 <- exp[, -3]




effLen = gfe$length
#转化
Counts2TPM <- function(exp1, effLen){
  rate <- log(exp1) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

CC_res <- apply(exp1, 2, Counts2TPM, effLen = effLen)















# 步骤5：合并TPM到差异分析结果
diff_result <- read.csv("E:/X101SC19123287-Z01-J001-B1-16/X101SC19123287-Z01-J001-B1-16-result/4.Differential/1.deglist/AGS_LvsAGS/AGS_LvsAGS_deg.csv", header = T)  # 替换为你的差异分析结果文件路径
tpm_data <- as.data.frame(CC_res)
tpm_data$gene_id <- exp$gene_id_clean  # 使用去掉版本号的gene_id

merged_data <- diff_result %>%
  mutate(gene_id_clean = sub("\\.[0-9]+$", "", gene_id)) %>%  # 确保diff_result的gene_id也去掉版本号
  left_join(tpm_data, by = "gene_id") %>%
  select(-gene_id_clean)  # 移除临时列



# 保存结果
write_csv(merged_data, "E:/Supplementary_table_6.csv")

# 验证TPM总和
colSums(CC_res, na.rm = TRUE)  # 每列应接近1e6
