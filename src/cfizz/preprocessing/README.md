# `src/cfizz/preprocessing/` — upstream Hi-C preprocessing reference scripts

> Three parameterized scripts that take Hi-C data from **raw FASTQ reads** to a **multi-resolution balanced `.mcool`**.
> They are a **reference implementation**, not a requirement. If you already have a HiC-Pro / distiller / 4DN workflow, feed its `.mcool` directly to downstream cfizz.

## What is here

| Script | Form | Purpose | Wrapped CLI tools |
|---|---|---|---|
| `fastq2pairs.sh` | bash | FASTQ → high-quality pairs | fastp + bwa-mem2 + pairtools (parse/sort/dedup/select) |
| `pairs2mcool.sh` | bash | pairs → multi-resolution balanced `.mcool` | cooler (cload pairs + zoomify + balance) |
| `stats.py` | python | parse pairtools stats and build a processing-yield table | pure Python text processing |

Design principles:

- **Keep the original form**: `.sh` stays bash, `.py` stays Python; no unnecessary "unified Python CLI" wrapper.
- **Fully parameterized**: no hard-coded sample names, absolute paths or species; everything comes from CLI flags.
- **Do not reimplement algorithms**: fastp / bwa-mem2 / pairtools / cooler are called from `PATH`; these scripts only assemble parameters, provide skip-if-exists tolerance and explicit logging.
- **Skippable**: if you already have `.mcool` files, ignore this directory entirely.

## Typical usage (3-step pipeline)

All commands below assume you are in the repository root.

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

Expected input per sample:

```text
/data/fastq/sampleA/sampleA_R1.fq.gz
/data/fastq/sampleA/sampleA_R2.fq.gz
```

Outputs:

```text
/out/sampleA/rawdata/sampleA/sampleA_R{1,2}.fq.gz      # fastp-trimmed reads
/out/sampleA/pairtools/sampleA.bam                      # bwa-mem2 alignment
/out/sampleA/pairtools/sampleA.pairs.gz                 # pairtools parse
/out/sampleA/pairtools/sampleA.sorted.pairs.gz          # pairtools sort
/out/sampleA/pairtools/sampleA.nodups.pairs.gz          # pairtools dedup
/out/sampleA/pairtools/sampleA.nodups.UU.pairs.gz       # pairtools select (final)
/out/sampleA/pairtools/sampleA.{stats,dedup.stats,select.stats}  # three stats files
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

Input (produced by Step 1):

```text
/out/sampleA/pairtools/sampleA.nodups.UU.pairs.gz
```

Outputs:

```text
/out/sampleA/sampleA_1000.cool    # 1 kb base-resolution .cool
/out/sampleA/sampleA_1000.mcool   # multi-resolution balanced .mcool (standard cfizz input)
```

### Step 3 — processing yield (optional, QC)

```bash
# Default: tab-separated stdout (easy to paste into slides or Excel)
python3 src/cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO"

# Human-readable table + summary block
python3 src/cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO" \
    --no-ppt

# Also write files (.txt with summary + .tsv for Excel)
python3 src/cfizz/preprocessing/stats.py \
    --pairtools-root /out \
    --samples sampleA sampleB \
    --display-names "WT" "KO" \
    --output /out/processing_stats.txt
```

## Parameters

### `fastq2pairs.sh`

| Flag | Required | Default | Description |
|---|---|---|---|
| `--samples` | ✅ | - | space-separated sample IDs |
| `--data-dir` | ✅ | - | root containing `<sample_id>/<sample_id>_R{1,2}.fq.gz` |
| `--output-root` | ✅ | - | output root; creates `<output-root>/<sample_id>/` |
| `--bwa-index` | ✅ | - | BWA-MEM2 index prefix, e.g. `/ref/hg38.fa` |
| `--chrom-sizes` | ✅ | - | chromosome sizes file |
| `--assembly` | | `hg38` | assembly name |
| `--tech-type` | | `hic` | `hic` or `microc`; controls the pairtools select condition |
| `--nproc` | | `8` | CPU threads |
| `--mem-sort` | | `100G` | memory for pairtools sort |
| `--fastp-extra` | | empty | extra fastp arguments (quote multi-token values) |
| `--max-parallel` | | `1` | number of samples processed concurrently (use with care) |
| `--overwrite` | | no | overwrite existing products |
| `--no-skip` | | no | alias of `--overwrite` |

### `pairs2mcool.sh`

| Flag | Required | Default | Description |
|---|---|---|---|
| `--samples` | ✅ | - | space-separated sample IDs |
| `--output-root` | ✅ | - | output root; should match `fastq2pairs.sh` |
| `--chrom-sizes` | ✅ | - | chromosome sizes file |
| `--base-resolution` | | `1000` | base resolution (1 kb) |
| `--resolutions` | | `5000,10000,25000,50000,100000,250000,500000,1000000` | comma-separated resolution levels |
| `--nproc` | | `8` | threads for balancing |
| `--max-parallel` | | `1` | number of samples processed concurrently |
| `--overwrite` | | no | overwrite existing products |
| `--no-skip` | | no | alias of `--overwrite` |

### `stats.py`

| Flag | Required | Default | Description |
|---|---|---|---|
| `--pairtools-root` | ✅ | - | root containing `<sample_id>/pairtools/` |
| `--samples` | ✅ | - | sample IDs (space-separated) |
| `--display-names` | | same as samples | display names aligned with `--samples` |
| `--output` | | no file | output `.txt`; without it, print to stdout only |
| `--ppt` / `--no-ppt` | | `--ppt` | stdout format: `--ppt` tab-separated; `--no-ppt` human-readable table + summary |

`--output` details:

- When set, writes both `.txt` (human-readable + summary) and `.tsv` (tab-separated) with the same stem.
- With `--no-ppt`, the `.txt` is human-readable; the `.tsv` remains tab-separated.

## Behavior details

### `fastq2pairs.sh` pipeline

```text
6 steps:
  1. fastp             QC + adapter trimming
  2. bwa-mem2 + samtools view   FASTQ → BAM
  3. pairtools parse   BAM → pairs.gz
  4. pairtools sort    sort by (chrom1, chrom2, pos1, pos2)
  5. pairtools dedup   deduplicate (nodups / unmapped / dups)
  6. pairtools select  filter by hic/microc condition + pairtools stats
```

Skip-if-exists is enabled by default: existing products are skipped. Concurrency is controlled by `--max-parallel`. Each run prints timestamped logs. `set -euo pipefail` stops the pipeline on the first failure.

### `pairs2mcool.sh` pipeline

```text
3 steps:
  1. cooler cload pairs  filtered pairs.gz → 1 kb .cool
  2. cooler zoomify      .cool → multi-resolution .mcool
  3. cooler balance      balance each resolution (iterate cooler ls)

Input:  <output-root>/<sample_id>/pairtools/<sample_id>.nodups.UU.pairs.gz
Output: <output-root>/<sample_id>/<sample_id>_<base-res>.mcool
```

### `stats.py` 3-step yield table

| Step | Stats file | Input pairs | Output pairs | Description |
|---|---|---|---|---|
| 1. `pairtools parse` | `<sample_id>.dedup.stats` (`total_mapped`) | total pairs | pairs mapped on both sides | mapping rate + one-side + unmapped fractions |
| 2. `pairtools dedup` | `<sample_id>.dedup.stats` | total_mapped | total_nodups | duplicate rate |
| 3. `pairtools select` | `<sample_id>.select.stats` | total_nodups | final | final retention + cumulative yield |

The summary block reports, per sample: raw total pairs / final high-quality pairs / overall yield.

## Dependencies

| Tool | Channel | Used for |
|---|---|---|
| fastp | bioconda | Step 1.1 QC |
| bwa-mem2 | bioconda | Step 1.2 alignment |
| samtools | bioconda | Step 1.2 SAM → BAM |
| pairtools | bioconda | Steps 1.3–1.6 |
| pairix | bioconda | pairs index (optional) |
| cooler | bioconda | Steps 2.1–2.3 |
| ucsc-fetchchromsizes | bioconda | chromosome sizes download (optional) |
| Python ≥ 3.9 | - | Step 3 `stats.py` |

> The main `environment.yml` covers the first seven tools. If you already have `.mcool` files, skip this directory and use only the core cfizz environment.

## Notes

- **bwa-mem2 and pairtools sort are memory-hungry**: keep `--max-parallel` at 1 unless you have enough RAM.
- **Existing products are skipped**: use `--overwrite` or `--no-skip` to force a rerun.
- **Avoid spaces in paths**: directory and file names should not contain spaces; some wrapped CLIs handle them poorly.
- **`--bwa-index` is an index prefix**: `--bwa-index /ref/hg38.fa` expects `/ref/hg38.fa.{amb,ann,bwt,pac}` etc.; do not include the `.amb` suffix.
- **microc needs stricter filtering**: with `--tech-type microc`, consider adding `--trim_poly_g` to `--fastp-extra`; the pairtools select condition is stricter.
