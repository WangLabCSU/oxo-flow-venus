# 生产级独立审计:临床解读模块(2026-10-27)

独立评审(双背景评审员:肿瘤分子病理 + 生物信息学,基于 commit c372ff8 + 本轮修复),
按科学性 / 专业性 / 严谨准确性 / 工程可靠性四维度执行,联网核实了 FDA TMB-H(2026-10-16)、
T-DXd 中国获批(2026-01-22)、selpercatinib(2022-10)、pamiparib(2021-05 卵巢-only)、
zolbetuximab(2024-12)等关键事实。

## 执行摘要

分级纪律(Predictive-only Tier I/II、benign→IV、TMB 双阈值、tumor-only 披露、
tier_basis 可审计)达到准临床工具素养;但审计发现 1 处数据损坏缺陷(M1)与
3 条"表面在服务、实际永不触发"的死键(M2-M4,含胃癌最高价值的 HER2 行),
判定"未达生产级"。以下修复后达到"研究级、正交确认前置"的可运行状态。

## 修复清单(commit 对应)

| # | 级别 | 缺陷 | 修复 |
|---|---|---|---|
| M1 | major | `clinical_tier.py` pysam 对 Number=. String INFO 返回标量,`;".join` 逐字符损坏 CLNSIG,germline 筛查恒 0 | isinstance 分支 + 集成测试(带 CLNSIG=Pathogenic) |
| M2 | major | 代码键 `ERBB2 amp` vs 表键 `ERBB2 amp/IHC+` → trastuzumab/T-DXd 永不命中 | 表键改 `ERBB2 amp`,requires 列承载 IHC3+/ISH+ 强制确认 |
| M3 | major | 代码 `NTRK1/2/3 fusion resolved` vs 表键 `NTRK fusion resolved` → larotrectinib/entrectinib 永不命中 | `fusion_nmpa_key()` 折叠 NTRK 前缀 |
| M4 | major | `BRCA1/2 pathogenic (germline-grade)` 无生成路径,制造 PARP 假象 | 改键 `BRCA1/2 loss (somatic, context)`,由 BRCA1/2 缺失调用触发,context_only |
| M5 | major | CNV/fusion 行 tier 硬编码 "II",与 SNV 侧 Predictive-only 纪律矛盾 | `actionable_cnv_genes.tsv` 新增 `amp_tier` 列(ERBB2=II,其余 III);fusion tier 与 NMPA 命中挂钩(有获批行才 II) |
| M6 | major | AA-interval 匹配任意 1 残基重叠即命中 | 重叠 ≥50% max(interval 长度) + 负例测试(749-750 不吃 746-750 热点) |
| M7 | major | therapy-matching 层零测试 | `tests/test_clinical_interpret.py`:NMPA 可达性静态契约(新增 `keying_modality` 列区分 dna/protein/status;蛋白级键豁免并须可见)+ 端到端行为测试 |
| m1 | minor | 配对队列 SNV/TMB 仍标 tumor-only | `--cohort-mode`(paired/tumor_only),两条互补规则分别传参 |
| m2 | minor | 空结果无表头空 TSV | 固定 schema + fillna("") |
| m3 | minor | 报告仍用 ≥10 旧 cutoff 文案,与 FDA ≥13 内部矛盾 | 双阈值措辞(≥13 label / 10-<13 legacy-only) |
| m4 | minor | 文档旧键名 `oncokb_token` | 改 `oncokb_token_file`(OncoKB 装饰性列/静默降级遗留为已知项) |
| m8 | minor | md_to_html 表格单元格未转义 | 单元格 html.escape |

## 保留的已知取舍(已披露,不属缺陷)

- 蛋白级键(HER2-low RC48、CLDN18.2 zolbetuximab)不做 DNA 自动映射 —— WGS 不可判定,
  走正交检测推荐(报告 + docs/20261024-gc01-03-clinical-baseline.md)。
- OncoKB 富集列下游暂不消费(token 未注册时整链 skip,诚实 skip 规则)。
- CIViC 治疗行不带 per-evidence-item level(m5,biomarker 行有 levels 汇总暂代)。
- interval 匹配的残留风险:同 op + ≥50% 重叠仍可能跨变体继承映射,tier_basis 可审计,
  报告 RUO 定位。

## 与主流工具的差距(路线图,非阻断)

患者级签发报告要素(RUO 声明/LoD/局限段)、per-evidence-item 治疗行、
allele-specific CNV(接 ASCAT)、HRD、融合 RNA 正交、OncoKB CNV/fusion 覆盖、
新抗原模块(独立调研已立项)。
