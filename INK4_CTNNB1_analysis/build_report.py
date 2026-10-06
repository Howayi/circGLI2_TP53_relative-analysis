"""Write the interpretation report from checked evidence tables."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"


def table(frame):
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join("—" if pd.isna(v) else str(v).replace("|", "/") for v in row) + " |")
    return "\n".join(lines)


def main():
    summary = json.loads((ROOT / "validation.json").read_text())
    ann = pd.read_csv(DATA / "variant_annotations.tsv", sep="\t")
    paired = pd.read_csv(DATA / "paired_evidence_annotated.tsv", sep="\t")
    candidates = pd.read_csv(DATA / "tumor_only_candidates_review.tsv", sep="\t")
    alleles = pd.read_csv(DATA / "allele_evidence_annotated.tsv", sep="\t")
    coverage = pd.read_csv(DATA / "coverage_recomputed.tsv", sep="\t")
    counts = pd.crosstab(paired.target, paired.label).reindex(["CDKN2A", "CDKN2B", "CTNNB1"])
    overview = pd.DataFrame({"基因": counts.index, "候选变异记录": [16, 22, 590],
                             "肿瘤侧候选": counts.candidate_tumor_only.values,
                             "两侧均有ALT": counts.alt_in_both.values,
                             "正常侧判断能力不足": counts.tumor_alt_normal_underpowered.values,
                             "正常侧无覆盖": counts.tumor_alt_normal_no_coverage.values})
    main_sites = candidates[(candidates.patient == 1) & candidates.pos.isin([41239336, 41239897, 41240109, 41239922, 41239918])].copy()
    main_sites = main_sites.sort_values("T_alt", ascending=False)
    main_table = pd.DataFrame({"位点（GRCh38）": main_sites.variant_id,
                                "1T ALT/总AD": main_sites.T_alt.astype(str) + "/" + main_sites.T_depth.astype(str),
                                "1T VAF": main_sites.T_vaf_exact.map(lambda v: f"{v:.2%}"),
                                "1N ALT/总AD": main_sites.N_alt.astype(str) + "/" + main_sites.N_depth.astype(str),
                                "主转录本后果": main_sites.canonical_consequence.map({"synonymous_variant": "同义 p.Asp780=", "3_prime_UTR": "3′UTR"})})
    changed = paired[paired.label_changed_by_exact_p0].sort_values(["patient", "pos"])
    label_names = {"candidate_tumor_only": "肿瘤侧候选", "tumor_alt_normal_underpowered": "正常侧判断能力不足"}
    changed_table = pd.DataFrame({"样本": changed.patient.astype(str) + "T", "位点": changed.variant_id,
                                  "原始P0": changed.p_zero_normal_alt_original.map(lambda x: f"{x:.4g}"),
                                  "整数计数重算P0": changed.p_zero_exact.map(lambda x: f"{x:.5f}"),
                                  "原始标签": changed.label_original.map(label_names),
                                  "更正标签": changed.label.map(label_names)})
    miss = candidates[candidates.canonical_consequence == "missense_variant"].sort_values(["patient", "pos"])
    miss_table = pd.DataFrame({"样本": miss.patient.astype(str) + "T", "蛋白改变": miss.protein_change,
                              "ALT/总AD": miss.T_alt.astype(str) + "/" + miss.T_depth.astype(str),
                              "VAF": miss.T_vaf_exact.map(lambda v: f"{v:.2%}"),
                              "ALT正/反向": miss.T_alt_fwd.astype(str) + "/" + miss.T_alt_rev.astype(str)})
    ink = alleles[alleles.target.isin(["CDKN2A", "CDKN2B"]) & alleles["sample"].str.endswith("T") & alleles.alt_reads.ge(3)].merge(
        ann[["variant_id", "ARF_protein_change"]], on="variant_id")
    ink = ink[ink.canonical_consequence.eq("missense_variant") | ink.ARF_protein_change.notna()]
    ink_table = pd.DataFrame({"样本": ink["sample"], "基因": ink.target, "位点": ink.variant_id,
                             "主转录本改变": ink.protein_change, "ARF改变": ink.ARF_protein_change,
                             "ALT/总AD": ink.alt_reads.astype(str) + "/" + ink.depth.astype(str),
                             "ALT正/反向": ink.alt_fwd.astype(str) + "/" + ink.alt_rev.astype(str)})
    coding_canonical = ann[(ann.target == "CTNNB1") & ann.protein_position.isin([33, 37, 41, 45])]
    hotspots = alleles.merge(coding_canonical[["variant_id", "protein_position"]], on="variant_id")
    hotspots = hotspots[hotspots.alt_reads >= 3]
    hotspot_table = pd.DataFrame({"样本": hotspots["sample"], "蛋白后果": hotspots.protein_change,
                                 "位点": hotspots.variant_id, "ALT/总AD": hotspots.alt_reads.astype(str) + "/" + hotspots.depth.astype(str),
                                 "ALT正/反向": hotspots.alt_fwd.astype(str) + "/" + hotspots.alt_rev.astype(str)})

    report = f"""# RNA-seq 候选变异解释与可视化：精确 P0 更正版

生成日期：2026-09-30。参考坐标：GRCh38 / hg38。图形由 R 制作，图内文字为英文。

本版更新第 2、4、5 张图、标签与相关说明；第 1、3、6 张图从上一版直接复制，四种格式的文件内容完全一致。原始测序结果未修改。

## 主要结论

这些数据支持若干 RNA 层面的候选等位信号，尚不足以确认 DNA 体细胞突变。最强的 CTNNB1 信号主要是主转录本同义变异或 3′UTR 信号；不能把它们直接称为蛋白激活突变。蛋白改变候选的计数普遍较弱，CDKN2A/CDKN2B 阴性结论还受到正常侧覆盖不足的限制。

- 共 19 个样本、8 对 N/T；41T、44N、56N 未配对。
- 候选 VCF 有 628 条位点/ALT 记录，其中 CTNNB1 590、CDKN2A 16、CDKN2B 22。这不是 628 个经过验证的突变。
- 配对表仍有 727 条“患者 × 位点 × ALT”比较记录。按整数计数重算并双向更正标签后，`candidate_tumor_only` 为 **34 条**，全部属于 CTNNB1，对应 **31 个位点/ALT**。
- 34 条中有 16 条 ALT 仅来自一个比对方向、25 条 ALT < 10、5 条正常深度 < 20、6 条肿瘤深度 < 20；这些提示可以重叠。
- 原先 39 条中有 6 条退出，但另有 1 条原 `underpowered` 记录进入，因此正确计算为 **39 − 6 + 1 = 34**。这些仍是筛查候选，不是已经确认的体细胞突变。
- 患者 1 的候选由 17 变为 16，患者 28 由 5 变为 3，患者 33 由 6 变为 4。其他患者数量不变。

{table(overview)}

表中标签计数是患者—位点—ALT 记录数；“候选变异记录”是整个候选 VCF 的独立位点—ALT 数。其他分类见完整明细。图例 `Normal ALT only` 对应 N_ALT ≥3、T_ALT <3，允许肿瘤侧有 1–2 条 ALT；它不代表肿瘤侧绝对没有 ALT。

## 怎样读突变字段

- `chrom` / `pos`：染色体及 1-based 基因组坐标。
- `ref` / `alt`：参考与非参考等位序列；单碱基替换通常是 SNV，长度变化通常是小型插入/缺失。
- `depth`：保存在拆分前的原始总 AD（VCF `ADT`），是本流程过滤后的等位计数分母，并非直接来自 DNA 的覆盖度。
- `ref_reads` / `alt_reads`：REF 和当前 ALT 的支持计数；多等位拆分后，两者之和可小于原始总 AD。
- `vaf`：ALT/原始总 AD；原文件保留三位小数，本报告与图使用整数计数得到的精确比例。没有覆盖时为 NA。
- `alt_fwd` / `alt_rev`：ALT 的比对正向/反向计数。RNA 文库本身可能有方向偏倚，不能仅凭一个方向就判定是假阳性。
- `p_zero_normal_alt` / `p_zero_exact`：本更正版均使用原始整数计数重算，避免提前舍入 VAF。它是假设正常与肿瘤具有相同 ALT 比例时，正常侧观察不到 ALT 的插值概率，不是校正后的显著性 P 值、致病概率或体细胞变异置信度。
- `label`：使用精确计数 P0 的探索性分类；`label_original` 与 `p_zero_normal_alt_original` 保留原始标签和概率，`label_changed_by_exact_p0` 指示是否更正。VCF 没有经过正式基因型/体细胞变异判定，不能把其中的 QUAL 或 FILTER 当成完整检测结果。

原 VCF 已保存拆分前总 AD；829 条样本—等位记录的 REF+当前 ALT 小于总 AD，这是其他 ALT 计数造成的正常现象。

## 最值得先复核的 CTNNB1 信号

{table(main_table)}

`chr3:41239336 C>T` 的 1T 支持为 602/1244（48.39%），1N 为 0/840，两个方向分别为 158/444。按 CTNNB1 主转录本 ENST00000349496.11 本地注释，它是 **p.Asp780=**：第 780 个氨基酸仍为 Asp。强等位信号不等于蛋白功能改变。

这个位点也见于 3N（58/121，47.93%）和 33N（298/554，53.79%）；6T 为 10/991（1.01%）。这些跨样本观察提示它不能统一解释成肿瘤特异突变，需考虑遗传多态、RNA 等位表达、低比例背景或其他技术/生物学解释。

患者 1 的 N/T 差异值得核查样本身份，但**单个位点的 RNA 等位比例不一致不能证明样本错配**。建议用全 BAM 中分散的多个常见 SNP 做样本指纹比较，并结合肿瘤拷贝数/杂合性变化及 RNA 等位表达解释。Picard CrosscheckFingerprints 使用多个基因组位置与单倍型块汇总身份匹配证据；本次没有进行该检查。[Picard 官方说明](https://gatk.broadinstitute.org/hc/en-us/articles/360037435811-CrosscheckFingerprints-Picard)

`41239897 T>G`、`41240109 C>CTAAT` 与 `41239922 T>A` 在主转录本上属于 **3′UTR**。`C>CTAAT` 是锚定在 C 后的 4 bp 插入，不能因插入长度不是 3 的倍数就称为主转录本编码区移码。多个位点在其他患者的正常样本中也有明显 ALT 支持；同一区域相邻 SNV/indel 还需要检查局部重复序列与比对情况。

功能位置具有转录本依赖性。`41239897` 在一个其他编码转录本的 CDS 范围内，所以本报告明确使用“主转录本 3′UTR”，不声称所有异构体都没有蛋白后果。各位点的其他编码转录本锚点重叠信息已保存在完整注释表中；未完成全转录本 VEP 注释或致病性分类。

## CTNNB1 蛋白改变候选

更正后的肿瘤侧候选中，主转录本后果为：11 条同义、9 条错义、8 条 3′UTR、6 条不位于主转录本外显子。这是 34 条患者—位点记录的分类。

{table(miss_table)}

9 条错义记录都只有 3–8 个 ALT 支持，VAF 范围为 0.28%–2.07%，其中 6 条只来自一个方向。新增的 1T p.Ala87Val 为 3/874（0.34%），方向 0/3，支持仍然很弱；33T p.His118Tyr 和 28T p.Lys666Arg 因精确 P0 >0.05 退出本候选集合。图例仅区分一个或两个 ALT 比对方向，不再把重叠的 P0 状态当作互斥类别。

42T 的 p.Ala13Thr 为 8/939，方向 6/2；33T 的 p.Ser179Pro 为 6/752，方向 5/1；57T 的 p.Met194Ile 为 3/145，方向 2/1。它们可进入人工复核名单，但不能直接作为致病或驱动突变报告。

若用本报告透明设定的复核条件（更正后的肿瘤侧候选、N/T 总 AD 均 ≥20、T_ALT ≥10、ALT 两个方向各 ≥2），仍有 6 条患者—位点记录通过，均为 p.Asp780= 或主转录本 3′UTR 信号。这是展示用的复核分层，未经验证，不能替代正式调用和过滤。

在 CTNNB1 第 33、37、41、45 位残基附近，达到 ALT ≥3 的观察如下：

{table(hotspot_table)}

41T 的 p.Ser33Pro 只有 3/1726（0.17%），均来自反向，且缺少配对正常样本；6N 的 p.Ser33Phe 是正常样本中的低计数观察。其余两条是同义变化。当前数据没有提供强支持的这些残基的肿瘤错义信号。

## CDKN2A / CDKN2B：不能把“0 条肿瘤侧候选”理解为“没有异常”

CDKN2A 的 8 个配对正常样本中，目标区域达到深度 ≥10 的比例范围为 **0%–13.27%**，中位数 **6.23%**；CDKN2B 为 **9.65%–37.50%**，中位数 **20.98%**。覆盖度分母是所有目标外显子的合并区间加两侧 10 bp，包含 UTR 与不同转录本，不等于主转录本 CDS 的覆盖比例。

`chr9:21968200 C>G` 在多份 N/T 中 ALT 比例为 100%，例如 1N 为 11/11、28N 为 9/9、42N 为 5/5，且主 p16 和选定 ARF 转录本均属于 3′UTR。因此它不是这些数据中的可靠肿瘤特异蛋白改变；100% RNA ALT 也不等于已经证明 DNA 纯合。

肿瘤样本中达到 ALT ≥3 的 p16、ARF 或 p15 蛋白改变信号为：

{table(ink_table)}

这 10 条均只有 3–4 个 ALT 支持，并且全部只来自一个比对方向。41T 无配对正常样本，其他配对的正常側深度偏低。特别是 `chr9:21971033 G>T` 在 p16 上预测为 p.Ala109Asp，在 ARF 上预测为 p.Cys123Ter；但 41T 只有 3/226，方向 3/0，不能据此断言 ARF 已失活。

CDKN2A 的两个主要蛋白产品使用不同阅读框，所以不能把一个转录本的蛋白后果直接套用到另一个。[NCBI Gene 对 CDKN2A 转录本的说明](https://www.ncbi.nlm.nih.gov/gene/1029)

本分析仅针对 RNA 层面的 SNV/小型 indel 信号；不用于判断 CDKN2A/CDKN2B 的 DNA 纯合缺失、拷贝数变化、启动子甲基化或蛋白缺失。

## P0 精度问题

按原始整数计数，正确的插值式为：

```text
P0 = (1 - T_alt / T_depth) ^ N_depth
```

原始表中的概率与使用三位小数 T_vaf 的计算一致（最大绝对差约 0.000499，为输出概率舍入误差）。在低 VAF、正常深度很高时，先舍入 VAF 再求幂会改变标签，且可以向两个方向跨越阈值：

{table(changed_table)}

本版三张相关图与导出表使用更正后的标签。新增位点 `chr3:41224972 C>T` 的 1T 为 3/874、1N 为 0/877：精确 P0=0.0490225，而用三位小数 VAF 得到约 0.0717，因此由判断能力不足转为候选。退出的六条与这一新增记录均保存在 `P0_label_changes.tsv`。

更正仅改变 7 个配对标签及相关概率，所有原始 REF/ALT、深度、支持计数、位点与注释不变。H 盘原文件未修改；更正表保留原字段用于追溯。建议后续筛查脚本始终使用整数计数计算概率，最后才格式化展示。即使没有此精度问题，这个 P0 仍忽略肿瘤 VAF 估计的不确定性、RNA 等位表达差异、相关错误及多重比较；低深度下 VAF=1 会得到 P0=0，尤其容易显得过于确定。

运行说明还引用了 BCFtools 1.19 的 indel 默认阈值，但 VCF 头部记录实际为 1.24，需以实际版本核对默认行为。[BCFtools 1.24 文档](https://www.htslib.org/doc/1.24/bcftools.html#mpileup)

## 图像与阅读方法

1. `01_coverage_heatmap`：先看正常侧是否具备观察到 ALT 的覆盖基础；灰度/蓝度和格内百分数为覆盖比例。
2. `02_candidate_status_overview`：更正后的分类按患者、基因汇总，纵轴是比较记录数量，不是确诊突变负荷。
3. `03_paired_VAF_evidence`：主要位点 N/T 比例对照，蓝色为 N、橙色为 T；空心点表示深度 <20，NA 表示无覆盖。
4. `04_all_tumor_only_candidate_heatmap`：更正后的 31 个候选位点/ALT ×16 个配对样本，共 496 格、34 个橙色候选边框；×为单方向 ALT、○为低深度、灰色为无覆盖。移除原来用于 P0 不过阈值的三角标记。
5. `05_missense_candidate_evidence`：更正后的 9 条 CTNNB1 肿瘤侧错义候选的比例与计数，提示蛋白改变证据强度。
6. `06_ink4_protein_candidate_evidence`：CDKN2A/CDKN2B 的 p16、ARF、p15 候选计数与正常侧覆盖说明。

每幅图提供 PNG（300 dpi）、PDF、SVG 与 TIFF（600 dpi）。图文与数据均为描述性展示，没有将 reads 当作生物学重复进行组间显著性检验。

## 核对与复现

已核对 VCF 的 19 个样本名与顺序、全部 11,932 条样本—等位记录的 AD/ADT、正反向计数和、配对计数、VAF 舍入、20,224 个目标碱基和全部 57 个覆盖汇总单元。BCF 可完整解压并解析 24,950 条记录，样本头部一致。487 条主转录本编码 SNV 的 REF 与公共 CDS 匹配；p16 CDS 另外与公共基因组序列按外显子重建结果一致。

还直接核对了原始 `tumor_normal_pairs.tsv` 的全部 727 条记录与旧版派生表，每个原始字段一致；原始 `depth_summary.tsv` 的 57 行与重算结果的目标碱基数、最大深度及 114 个覆盖比例一致（按原文件三位小数精度）。本版另验证精确 P0 的两向阈值变化、34 条候选、31 个位点和 9 条错义记录，并逐项比对图中数据。第 1、3、6 张图的 PNG/PDF/SVG/TIFF 与上一版逐字节一致。

公共参考来自 Ensembl 的 GRCh38 基因与转录本信息，已缓存完整 JSON；主转录本为 CTNNB1 ENST00000349496.11、CDKN2A ENST00000304494.10、CDKN2B ENST00000276925.7；附加 ARF 为 ENST00000579755.2。参考为下载时 Ensembl 注释，与输入 GENCODE v48 的全转录本集合可能存在差异，因此所有后果限定在指定参考转录本上。[Ensembl 基因/转录本接口](https://rest.ensembl.org/documentation/info/symbol_lookup)

完整注释、计数、覆盖比例、候选复核标记与输入 SHA-256 摘要位于 `source_data`。索引仅用于文件访问，不被当作独立变异证据。未获得原始 BAM/FASTQ，未检查单分子独立性、局部比对图或样本基因型一致性；因此仍需要人工 IGV 复核和相应 DNA 验证。

在同目录运行：

```text
python prepare_analysis.py --input <results_directory>
Rscript --vanilla run_figures.R
python build_report.py
```

Python 准备脚本需要 pandas；R 脚本需要 data.table、ggplot2、patchwork、ggrepel、jsonlite、ragg、svglite、systemfonts 与 scales。`run_figures.R` 显式处理 UTF-8，以避免 Windows 中文标签乱码。公共参考已缓存，重新分析无需上传样本变异数据。
"""
    (ROOT / "REPORT.zh-CN.md").write_text(report, encoding="utf-8")
    print("Report written.")


if __name__ == "__main__":
    main()
