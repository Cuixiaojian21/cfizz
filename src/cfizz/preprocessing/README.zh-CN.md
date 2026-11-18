# `src/cfizz/preprocessing/` — 上游 Hi-C 预处理参考脚本

> 3 个参数化脚本，完成 Hi-C 数据从 **FASTQ 原始 reads** 到 **多分辨率平衡 `.mcool`** 的全链路预处理。
> 本目录是**参考实现**，**不是强制流程**。已有 HiC-Pro / distiller / 4DN 流程的，可以直接把 `.mcool` 喂给下游 cfizz。

## 目录内容

| 脚本 | 形态 | 功能 | 包装的 CLI 工具 |
|---|---|---|---|
| `fastq2pairs.sh` | bash | FASTQ → 高质量 pairs | fastp + bwa-mem2 + pairtools（parse/sort/dedup/select） |
| `pairs2mcool.sh` | bash | pairs → 多分辨率平衡 `.mcool` | cooler（cload pairs + zoomify + balance） |
| `stats.py` | python | 解析 pairtools stats，生成处理得率表 | 纯 Python 文本处理 |

设计原则：

- **保持原形态**：`.sh` 保持 bash，`.py` 保持 Python；不做不必要的“统一 Python CLI”封装。
- **完全参数化**：无硬编码样本名、绝对路径或物种名；全部通过 CLI flag 传入。
- **不重新实现算法**：fastp / bwa-mem2 / pairtools / cooler 全部调用 `PATH` 中的现成工具；脚本只做参数组装、skip-if-exists 容错和显式日志。
- **可跳过**：已有 `.mcool` 文件的用户可完全忽略本目录。

## 典型用法（3 步流水线）

以下命令均在仓库根目录执行。

### Step 1 — FASTQ → filtered pairs

```bash
bash src/cfizz/preprocessing/fastq2pairs.sh \
    --samples sampleA sampleB \
    --data-dir /data/fastq \
    --output-root /out \
    --bwa-index /ref/hg38.fa \
    --chrom-sizes /ref/hg38.chrom.exceptYMetc.sorted.sizes \
    --assembly hg38 \
    --tech-type hic \
    --nproc 16
```

每个样本的输入约定：

```text
/data/fastq/sampleA/sampleA_R1.fq.gz
/data/fastq/sampleA/sampleA_R2.fq.gz
```

输出文件：

```text
/out/sampleA/rawdata/sampleA/sampleA_R{1,2}.fq.gz      # fastp 清洗后
/out/sampleA/pairtools/sampleA.bam                      # bwa-mem2 比对
/out/sampleA/pairtools/sampleA.pairs.gz                 # pairtools parse
/out/sampleA/pairtools/sampleA.sorted.pairs.gz          # pairtools sort
/out/sampleA/pairtools/sampleA.nodups.pairs.gz          # pairtools dedup
/out/sampleA/pairtools/sampleA.nodups.UU.pairs.gz       # pairtools select（最终）
/out/sampleA/pairtools/sampleA.{stats,dedup.stats,select.stats}  # 3 个统计文件
```

### Step 2 — pairs → balanced `.mcool`

```bash
bash src/cfizz/preprocessing/pairs2mcool.sh \
    --samples sampleA sampleB \
    --output-root /out \
    --chrom-sizes /ref/hg38.chrom.exceptYMetc.sorted.sizes \
    --base-resolution 1000 \
    --resolutions "5000,10000,25000,50000,100000,250000,500000,1000000" \
    --nproc 8
```

输入（由 Step 1 产出）：

```text
/out/sampleA/pairtools/sampleA.nodups.UU.pairs.gz
```

输出：

```text
/out/sampleA/sampleA_1000.cool    # 1 kb 基础分辨率 .cool
/out/sampleA/sampleA_1000.mcool   # 多分辨率平衡 .mcool（cfizz 标准输入）
```

### Step 3 — 处理得率（可选，QC）

```bash
# 默认：Tab 分隔输出到 stdout（可直接贴 PPT / Excel）
python3 src/cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO"

# 人类可读表格 + 汇总区
python3 src/cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO" \
    --no-ppt

# 同时写文件（.txt 含汇总区 + .tsv 可贴 Excel）
python3 src/cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO" \
    --output /out/processing_stats.txt
```

## 参数

### `fastq2pairs.sh`

| Flag | 必填 | 默认值 | 说明 |
|---|---|---|---|
| `--samples` | ✅ | - | 空格分隔的样本 ID 列表 |
| `--data-dir` | ✅ | - | 包含 `<sample_id>/<sample_id>_R{1,2}.fq.gz` 的根目录 |
| `--output-root` | ✅ | - | 输出根目录，会创建 `<output-root>/<sample_id>/` |
| `--bwa-index` | ✅ | - | BWA-MEM2 索引前缀，如 `/ref/hg38.fa` |
| `--chrom-sizes` | ✅ | - | 染色体大小文件 |
| `--assembly` | | `hg38` | 装配名 |
| `--tech-type` | | `hic` | `hic` 或 `microc`，决定 pairtools select 条件 |
| `--nproc` | | `8` | CPU 线程数 |
| `--mem-sort` | | `100G` | pairtools sort 内存 |
| `--fastp-extra` | | 空 | 附加 fastp 参数（多 token 用引号包起来） |
| `--max-parallel` | | `1` | 同时处理的样本数（谨慎使用） |
| `--overwrite` | | 否 | 已存在产物也覆盖 |
| `--no-skip` | | 否 | 同 `--overwrite`（语义别名） |

### `pairs2mcool.sh`

| Flag | 必填 | 默认值 | 说明 |
|---|---|---|---|
| `--samples` | ✅ | - | 空格分隔的样本 ID 列表 |
| `--output-root` | ✅ | - | 输出根目录，应与 `fastq2pairs.sh` 一致 |
| `--chrom-sizes` | ✅ | - | 染色体大小文件 |
| `--base-resolution` | | `1000` | 基础分辨率（1 kb） |
| `--resolutions` | | `5000,10000,25000,50000,100000,250000,500000,1000000` | 逗号分隔的多分辨率层级 |
| `--nproc` | | `8` | balance 线程数 |
| `--max-parallel` | | `1` | 同时处理的样本数 |
| `--overwrite` | | 否 | 已存在产物也覆盖 |
| `--no-skip` | | 否 | 同 `--overwrite` |

### `stats.py`

| Flag | 必填 | 默认值 | 说明 |
|---|---|---|---|
| `--pairtools-root` | ✅ | - | 包含 `<sample_id>/pairtools/` 的根目录 |
| `--samples` | ✅ | - | 样本 ID 列表（空格分隔） |
| `--display-names` | | 与 samples 同 | 与 `--samples` 一一对应的显示名 |
| `--output` | | 不写 | 输出 `.txt` 路径；不指定则只打印到 stdout |
| `--ppt` / `--no-ppt` | | `--ppt` | stdout 格式：`--ppt` Tab 分隔；`--no-ppt` 人类可读表格 + 汇总区 |

`--output` 行为细节：

- 指定时同时写 `.txt`（人类可读 + 汇总区）和 `.tsv`（Tab 分隔，同名衍生）。
- 加 `--no-ppt` 时，`.txt` 走人类可读；`.tsv` 仍是 Tab 分隔。

## 行为细节

### `fastq2pairs.sh` 内部流水线

```text
6 步：
  1. fastp             QC + adapter 修剪
  2. bwa-mem2 + samtools view   FASTQ → BAM
  3. pairtools parse   BAM → pairs.gz
  4. pairtools sort    按 (chrom1, chrom2, pos1, pos2) 排序
  5. pairtools dedup   去重（nodups / unmapped / dups）
  6. pairtools select  按 hic/microc 条件过滤 + pairtools stats
```

默认开启 skip-if-exists：产物存在则跳过。并发由 `--max-parallel` 控制。每次运行输出带时间戳的日志。`set -euo pipefail` 保证任一步失败立即停止。

### `pairs2mcool.sh` 内部流水线

```text
3 步：
  1. cooler cload pairs  filtered pairs.gz → 1 kb .cool
  2. cooler zoomify      .cool → 多分辨率 .mcool
  3. cooler balance      对每个分辨率做平衡（iterate cooler ls）

输入：<output-root>/<sample_id>/pairtools/<sample_id>.nodups.UU.pairs.gz
输出：<output-root>/<sample_id>/<sample_id>_<base-res>.mcool
```

### `stats.py` 3 步得率表

| 步骤 | 解析的 stats 文件 | 输入 pairs | 输出 pairs | 描述 |
|---|---|---|---|---|
| 1. `pairtools parse` | `<sample_id>.dedup.stats`（`total_mapped`） | 总 pairs | 两侧都比对成功的 pairs | mapping 率 + 单侧 + 未比对占比 |
| 2. `pairtools dedup` | `<sample_id>.dedup.stats` | total_mapped | total_nodups | dups 率 |
| 3. `pairtools select` | `<sample_id>.select.stats` | total_nodups | final | 最终保留率 + 累计得率 |

汇总区按样本打印：原始总 pairs / 最终高质量 pairs / 总体得率。

## 依赖

| 工具 | 渠道 | 用途 |
|---|---|---|
| fastp | bioconda | Step 1.1 质控 |
| bwa-mem2 | bioconda | Step 1.2 比对 |
| samtools | bioconda | Step 1.2 SAM → BAM |
| pairtools | bioconda | Step 1.3–1.6 |
| pairix | bioconda | pairs 索引（可选） |
| cooler | bioconda | Step 2.1–2.3 |
| ucsc-fetchchromsizes | bioconda | 染色体大小文件下载（可选） |
| Python ≥ 3.9 | - | Step 3 `stats.py` |

> 主 `environment.yml` 已覆盖前 7 项。已有 `.mcool` 文件的可跳过本目录，仅使用 cfizz 核心环境。

## 注意事项

- **bwa-mem2 和 pairtools sort 都很吃内存**：除非内存充足，`--max-parallel` 保持 1。
- **产物已存在会跳过**：用 `--overwrite` 或 `--no-skip` 强制重跑。
- **路径不要有空格**：目录名和文件名不要带空格；部分 CLI 对空格处理不友好。
- **`--bwa-index` 是索引前缀**：`--bwa-index /ref/hg38.fa` 对应 `/ref/hg38.fa.{amb,ann,bwt,pac}` 等；不要把 `.amb` 后缀写进去。
- **microc 需要更严的过滤**：`--tech-type microc` 时，建议在 `--fastp-extra` 中加 `--trim_poly_g`；pairtools select 条件会更严格。
