# `cfizz/preprocessing/` — 上游 Hi-C 预处理参考脚本

> 3 个参数化脚本,完成 Hi-C 数据从 **fastq 原始 reads** 到 **多分辨率平衡 .mcool** 的全链路预处理。
> 本目录的脚本是 **参考实现**,**不是强制** —— 已有 HiC-Pro / distiller / 4DN 流程的,可以直接拿 `.mcool` 喂给下游 cfizz。

## 这是什么

| 脚本 | 形态 | 功能 | 包装的 CLI 工具 |
|------|------|------|----------------|
| `fastq2pairs.sh` | bash | FASTQ → 高质量 pairs | fastp + bwa-mem2 + pairtools(parse/sort/dedup/select) |
| `pairs2mcool.sh` | bash | pairs → 多分辨率平衡 .mcool | cooler(cload pairs + zoomify + balance) |
| `stats.py` | python | 解析 pairtools stats,生成处理得率表 | (纯 Python 文本处理) |

3 个脚本的设计哲学:

- **形态继承**:原来是 .sh 就保持 .sh,原来是 .py 就保持 .py —— 不做"统一 Python CLI 包"那种过度封装
- **完全参数化**:无任何硬编码的样本名 / 绝对路径 / 物种名;一切通过 CLI flag 传入
- **不重新实现算法**:fastp / bwa-mem2 / pairtools / cooler 全部走系统 PATH 的现成 CLI 工具;本脚本只做参数组装 + skip-if-exists 容错 + 显式日志
- **可跳过**:已有 .mcool 文件的,可以完全忽略本目录

## 典型用法(3 步流水线)

### Step 1 — fastq → filtered pairs

```bash
bash cfizz/preprocessing/fastq2pairs.sh \
    --samples sampleA sampleB \
    --data-dir /data/fastq \
    --output-root /out \
    --bwa-index /ref/hg38.fa \
    --chrom-sizes /ref/hg38.chrom.exceptYMetc.sorted.sizes \
    --assembly hg38 \
    --tech-type hic \
    --nproc 16
```

输入文件约定(每个样本):

```
/data/fastq/sampleA/sampleA_R1.fq.gz
/data/fastq/sampleA/sampleA_R2.fq.gz
```

输出文件:

```
/out/sampleA/rawdata/sampleA/sampleA_R{1,2}.fq.gz      # fastp 清洗后
/out/sampleA/pairtools/sampleA.bam                      # bwa-mem2 比对
/out/sampleA/pairtools/sampleA.pairs.gz                 # pairtools parse
/out/sampleA/pairtools/sampleA.sorted.pairs.gz          # pairtools sort
/out/sampleA/pairtools/sampleA.nodups.pairs.gz          # pairtools dedup
/out/sampleA/pairtools/sampleA.nodups.UU.pairs.gz       # pairtools select(最终)
/out/sampleA/pairtools/sampleA.{stats,dedup.stats,select.stats}  # 3 个统计文件
```

### Step 2 — pairs → balanced .mcool

```bash
bash cfizz/preprocessing/pairs2mcool.sh \
    --samples sampleA sampleB \
    --output-root /out \
    --chrom-sizes /ref/hg38.chrom.exceptYMetc.sorted.sizes \
    --base-resolution 1000 \
    --resolutions "5000,10000,25000,50000,100000,250000,500000,1000000" \
    --nproc 8
```

输入文件(由 Step 1 产出):

```
/out/sampleA/pairtools/sampleA.nodups.UU.pairs.gz
```

输出文件:

```
/out/sampleA/sampleA_1000.cool    # 1kb 基础分辨率 .cool
/out/sampleA/sampleA_1000.mcool   # 多分辨率平衡 .mcool(下游 cfizz 标准输入)
```

### Step 3 — 统计处理得率(可选,QC 用)

```bash
# 默认:Tab 分隔输出到 stdout(可直接贴 PPT/Excel)
python3 cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO"

# 人类可读表格 + 汇总区
python3 cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO" \
    --no-ppt

# 写文件(同时生成 .txt 含汇总区 + .tsv 可贴 Excel)
python3 cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO" \
    --output /out/processing_stats.txt
```

## 全部参数

### `fastq2pairs.sh`

| flag | 必填 | 默认值 | 说明 |
|------|------|--------|------|
| `--samples` | ✅ | - | 空格分隔的样本 ID 列表 |
| `--data-dir` | ✅ | - | 包含 `<sample_id>/<sample_id>_R{1,2}.fq.gz` 的根目录 |
| `--output-root` | ✅ | - | 输出根目录,会创建 `<output-root>/<sample_id>/` 子目录 |
| `--bwa-index` | ✅ | - | BWA-MEM2 索引前缀(如 `/ref/hg38.fa`) |
| `--chrom-sizes` | ✅ | - | 染色体大小文件 |
| `--assembly` | | `hg38` | 装配名 |
| `--tech-type` | | `hic` | `hic` 或 `microc`,决定 pairtools select condition |
| `--nproc` | | `8` | CPU 线程数 |
| `--mem-sort` | | `100G` | pairtools sort 内存 |
| `--fastp-extra` | | 空 | 附加 fastp 参数(空格分隔的多 token 用引号包起来) |
| `--max-parallel` | | `1` | 同时处理的样本数(谨慎 >1,bwa-mem2 + sort 吃内存) |
| `--overwrite` | | 否 | 已存在产物也覆盖 |
| `--no-skip` | | 否 | 同 `--overwrite`(语义别名) |

### `pairs2mcool.sh`

| flag | 必填 | 默认值 | 说明 |
|------|------|--------|------|
| `--samples` | ✅ | - | 空格分隔的样本 ID 列表(与 `fastq2pairs.sh` 的 `--samples` 一致) |
| `--output-root` | ✅ | - | 输出根目录(应与 `fastq2pairs.sh` 共享同一个) |
| `--chrom-sizes` | ✅ | - | 染色体大小文件 |
| `--base-resolution` | | `1000` | 基础分辨率(1kb) |
| `--resolutions` | | `5000,10000,25000,50000,100000,250000,500000,1000000` | 多分辨率层级,逗号分隔 |
| `--nproc` | | `8` | balance 线程 |
| `--max-parallel` | | `1` | 同时处理的样本数 |
| `--overwrite` | | 否 | 已存在产物也覆盖 |
| `--no-skip` | | 否 | 同 `--overwrite` |

### `stats.py`

| flag | 必填 | 默认值 | 说明 |
|------|------|--------|------|
| `--pairtools-root` | ✅ | - | 包含 `<sample_id>/pairtools/` 子目录的根目录 |
| `--samples` | ✅ | - | 样本 ID 列表(空格分隔) |
| `--display-names` | | 与 samples 同 | 与 `--samples` 一一对应的显示名,用于表格"样本"列 |
| `--output` | | 不写 | 输出 `.txt` 路径;不指定 → 只打印到 stdout |
| `--ppt` / `--no-ppt` | | `--ppt` | stdout 输出格式:`--ppt` Tab 分隔(可贴 PPT/Excel);`--no-ppt` 人类可读表格 + 汇总区 |

`--output` 行为细节:

- 指定时同时写 `.txt`(人类可读 + 汇总区)+ `.tsv`(Tab 分隔,同名衍生)
- 同一文件路径 `--output` 加 `--no-ppt` 时,`.txt` 走人类可读;`.tsv` 仍是 Tab 分隔

## 行为细节

### `fastq2pairs.sh` 内部流水线

```
6 步流水线:
  1. fastp             QC + adapter 修剪
  2. bwa-mem2 + samtools view   fastq → BAM
  3. pairtools parse   BAM → pairs.gz
  4. pairtools sort    按 (chrom1, chrom2, pos1, pos2) 排序
  5. pairtools dedup   去重(nodups / unmapped / dups 三路)
  6. pairtools select  按 hic/microc 条件过滤 + pairtools stats

skip-if-exists 默认开启:产物存在则跳过该步
并发:用 --max-parallel 控制同时跑多少个样本
每次运行会输出时间戳日志
set -euo pipefail:任一步失败立即停
```

### `pairs2mcool.sh` 内部流水线

```
3 步流水线:
  1. cooler cload pairs  filtered pairs.gz → 1kb .cool
  2. cooler zoomify      .cool → 多分辨率 .mcool
  3. cooler balance      对每个分辨率做平衡(iterate cooler ls)

输入文件:<output-root>/<sample_id>/pairtools/<sample_id>.nodups.UU.pairs.gz
        (该文件应由 fastq2pairs.sh 跑出来)
输出文件:<output-root>/<sample_id>/<sample_id>_<base-res>.mcool
        (下游 cfizz 的标准输入)
```

### `stats.py` 3 步得率表

| 步骤 | 解析的 stats 文件 | 输入 pairs | 输出 pairs | 描述 |
|------|------------------|-----------|-----------|------|
| 1. `pairtools parse` | `<sample_id>.dedup.stats`(`total_mapped` 字段) | 总 pairs | 两侧都比对成功的 pairs | 算 mapping 率 + 单侧 + 未比对占比 |
| 2. `pairtools dedup` | `<sample_id>.dedup.stats` | total_mapped | total_nodups | 算 dups 率 |
| 3. `pairtools select` | `<sample_id>.select.stats` | total_nodups | final | 算最终保留率 + 累计得率 |

汇总区:每个样本打印 `原始总 pairs` / `最终高质量 pairs` / `总体得率`。

## 依赖

| 工具 | 渠道 | 用途 |
|------|------|------|
| fastp | bioconda | Step 1.1 质控 |
| bwa-mem2 | bioconda | Step 1.2 比对 |
| samtools | bioconda | Step 1.2 sam → bam |
| pairtools | bioconda | Step 1.3-1.6 |
| pairix | bioconda | pairs 文件索引(可选) |
| cooler | bioconda | Step 2.1-2.3 |
| ucsc-fetchchromsizes | bioconda | 染色体大小文件下载(可选) |
| Python ≥ 3.9 | - | Step 3 stats.py |

> cfizz 主仓库的 `environment.yml` 已经覆盖前 7 项,激活环境即可。
> 已有 `.mcool` 文件的可以完全跳过本目录,依赖 `environment.yml` 的 B 类(cfizz 本体)即可。

## 从原版迁移

本目录的 3 个脚本来自用户原项目的 3 个家什,经**统一参数化**收进 cfizz 仓库,对应关系:

| 本目录 | 原版(用户 `G:\2_0_demo\script\`) | 改了什么 |
|--------|----------------------------------|----------|
| `fastq2pairs.sh` | `a_1_run_50_1_50_2.sh` | 50_1/50_2 → `--samples`;<br/>/mnt/wsl/PHYSICALDRIVE6p1/... → `--data-dir` / `--output-root` / `--bwa-index` / `--chrom-sizes`;<br/>16 核写死 → `--nproc 8` 默认 + `--max-parallel 1` 并发可控 |
| `pairs2mcool.sh` | `a_2_pairs2mcool_50_1_50_2.sh` | 同上;**L46-52 if/else 拆 50_1/50_2 路径已改为循环**(派生自 `--output-root/<sample_id>/pairtools/`) |
| `stats.py` | `a_1_z_processing_stats.py` | 模块级 4 项硬编码(BASE_DIR/PAIRS_RESULT_ROOT/SAMPLES/SAMPLE_DISPLAY_NAMES)删除,改为 argparse 必填;<br/>`--ppt / --no-ppt` 替换原 `--ppt` 单独 flag |

行为完全保留:6 步 / 3 步 / 3 行统计 表格,skip-if-exists,`set -euo pipefail`,时间戳日志,`wait_for_slot` 并发控制,`hic` / `microc` 两套 pairtools select condition,等等。

## 注意事项

- **bwa-mem2 + pairtools sort 都很吃内存**:跑多样本时,`--max-parallel` 谨慎调到 >1(默认 1,稳妥)
- **产物已存在会跳过**:`--overwrite` 或 `--no-skip` 强制重跑
- **路径不要有空格**:3 个脚本的目录名 / 文件名都不要带空格(部分 CLI 工具对空格处理不友好)
- **bwa-mem2 索引是索引前缀**:`--bwa-index /ref/hg38.fa` 对应索引文件 `/ref/hg38.fa.{amb,ann,bwt,pac}` 等;不要把 `.amb` 后缀也写进去
- **microc 需要更严的过滤条件**:`--tech-type microc` 时,`--fastp-extra` 建议加 `--trim_poly_g`,pairtools select 条件会更严
