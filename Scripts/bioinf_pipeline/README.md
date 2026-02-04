# bioinf_pipeline（R）

一个“可拼装”的差异表达（DE）+ 富集分析小骨架，适合从**表达矩阵 + 分组信息**快速跑出：
- DEG 表（UP/DOWN/NOT）
- 火山图（PDF）
- GO/KEGG 富集（RDS + dotplot PDF，失败不影响主流程）

当前支持两类输入：
- **Raw counts** → `DESeq2`
- **TPM/FPKM/其他已归一化表达** → `limma`（推荐做 `log2(x+1)`，见下方 `cfg$expr$log_transform`）

物种支持：
- 人：`hs`（`org.Hs.eg.db`，KEGG: `hsa`）
- 鼠：`mm`（`org.Mm.eg.db`，KEGG: `mmu`）

> 说明：`make_cfg()` 里虽然列了 `geo/seurat`，但当前 `run_pipeline()` 只实现了 `counts_table` 和 `expr_table` 两种输入。

## 目录结构

- `run_pipeline.R`：最小可运行入口示例（默认走 `limma` + TPM 示例）
- `R/config.R`：配置（物种、对比、阈值等）
- `R/io_table.R`：读入 `meta.csv` 与表达矩阵，并按 `sample_id` 对齐
- `R/preprocess.R`：预处理（重复基因合并、过滤、可选 log 转换）
- `R/de_deseq2.R`：DESeq2 差异分析
- `R/de_limma.R`：limma 差异分析
- `R/idmap.R`：SYMBOL/ALIAS → ENTREZID 映射（用于富集）
- `R/enrich.R`：GO/KEGG 富集与 rank 工具
- `R/viz.R`：火山图/富集 dotplot/热图函数
- `R/export.R`：写 TSV、保存 RDS、保存 PDF
- `R/pipeline.R`：`run_pipeline()` 主流程

## 依赖（R 包）

按功能分组（缺包时脚本会提示 `Please install ...`）：
- I/O：`data.table`
- DESeq2 路线：`DESeq2`, `SummarizedExperiment`
- limma 路线：`limma`
- 注释/映射：`AnnotationDbi`, `org.Hs.eg.db` 或 `org.Mm.eg.db`
- 富集：`clusterProfiler`
- 作图：`ggplot2`, `ggrepel`, `enrichplot`（dotplot）, `pheatmap`（可选热图）

## 输入文件规范

### 1) `meta.csv`

必须包含：
- `sample_id`：样本名（必须与表达矩阵列名一致）
- `group`：分组列（默认 `group_col="group"`）

可额外包含：`batch`、协变量等（用于 DESeq2 的 `design` 公式）。

### 2) 表达矩阵（TSV/CSV，推荐 TSV）

要求：
- 第一列是基因列（默认列名 `SYMBOL`，可通过 `gene_col` 修改）
- 其余列是样本表达值，列名必须能与 `meta.csv` 的 `sample_id` 对上

两种输入类型：
- `counts_table`：raw counts（应为非负、接近整数）
- `expr_table`：TPM/FPKM/归一化表达（可按需要做 `log2(x+1)`）

脚本会自动按 `meta$sample_id` 与矩阵列名取交集并排序对齐；若完全对不上会直接报错。

## 配置与对比（contrasts）

`make_cfg()` 的常用参数：
- `species`：`"hs"` / `"mm"`
- `input_type`：`"counts_table"` 或 `"expr_table"`
- `outdir`：输出目录
- `group_col`：分组列名（默认 `"group"`）
- `design`：DESeq2 的设计公式字符串（默认 `"~ group"`）
- `contrasts`：对比列表，形如 `list(c("group","A","B"), c("group","C","D"))`（第 1 个元素必须等于 `group_col`）
- `cutoff`：阈值（默认 `padj < 0.05` 且 `|logFC| >= 1`）

对比命名规则：
- 每个 contrast 会生成一个 tag：`A_vs_B`，并用于输出文件前缀。

## 快速运行

建议在 `Scripts/bioinf_pipeline` 目录下运行（保证 `source("R/xxx.R")` 相对路径正确）。

### 方式 A：直接用 `run_pipeline.R`

编辑 `run_pipeline.R` 里的 `meta_file / expr_file / counts_file`，然后在 R 里执行：

```r
source("run_pipeline.R")
```

### 方式 B：在 R 会话里显式调用 `run_pipeline()`

先加载所有模块：

```r
source("R/config.R")
source("R/io_table.R")
source("R/preprocess.R")
source("R/design.R")
source("R/de_deseq2.R")
source("R/de_limma.R")
source("R/idmap.R")
source("R/enrich.R")
source("R/viz.R")
source("R/export.R")
source("R/pipeline.R")
```

#### 1) Raw counts → DESeq2

```r
cfg <- make_cfg(
  species   = "hs",
  input_type= "counts_table",
  outdir    = "outputs_counts",
  contrasts = list(c("group","HL","L"))
)

run_pipeline(
  cfg,
  input  = list(file="raw_counts.tsv", meta="meta.csv", gene_col="SYMBOL"),
  engine = "deseq2"
)
```

预处理要点（counts）：
- 同名基因会按行求和合并（`sum`）
- 默认过滤：至少 `min_samples=2` 个样本中 counts ≥ `min_count=10`
- 若非整数比例很高会给 warning（可能不是 raw counts）

#### 2) TPM/FPKM → limma

```r
cfg <- make_cfg(
  species   = "hs",
  input_type= "expr_table",
  outdir    = "outputs_tpm",
  contrasts = list(c("group","HL","L"))
)

# expr 预处理：设置为 TRUE 才会做 log2(x+1)；不设置/设为 FALSE 则保持原值
cfg$expr <- list(log_transform = TRUE)

run_pipeline(
  cfg,
  input  = list(file="tpm_matrix.tsv", meta="meta.csv", gene_col="SYMBOL"),
  engine = "limma"
)
```

注意（limma）：
- 当前实现只使用 `group` 做设计矩阵（`~0+group`），不会自动纳入 batch/协变量。

## 输出说明（每个对比 A vs B）

输出目录为 `cfg$outdir`，每个对比会生成：
- `A_vs_B_DEG.tsv`：DEG 表（列：`SYMBOL`, `logFC`, `P.Value`, `padj`, `change`）
- `A_vs_B_volcano.pdf`：火山图

仅当满足依赖且基因数足够时（不足/缺包会被捕获并跳过）：
- `A_vs_B_GO_UP.rds` + `A_vs_B_GO_UP_dot.pdf`
- `A_vs_B_KEGG_UP.rds` + `A_vs_B_KEGG_UP_dot.pdf`

DESeq2 额外保存：
- `A_vs_B_expr_vst.rds`：VST 表达矩阵（用于后续热图等）

## 常见问题

- **报错：No overlap between meta$sample_id and matrix colnames**：检查 `meta.csv` 的 `sample_id` 是否与矩阵列名完全一致（大小写/空格/前后缀）。
- **富集没有输出**：可能是缺少 `clusterProfiler/enrichplot/org.*.eg.db`，或 UP 基因数太少（<5），或 KEGG 访问受限/失败（当前设计为非致命跳过）。
- **counts 非整数 warning**：确保输入确实是 raw counts（不是 TPM/FPKM 或已归一化矩阵）。
