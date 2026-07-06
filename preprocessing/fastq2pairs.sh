#!/bin/bash
#
# fastq2pairs.sh — 将 FASTQ 处理成高质量 pairs 的 6 步流水线
#
# 流水线:fastp -> bwa-mem2+samtools view -> pairtools parse
#         -> pairtools sort -> pairtools dedup -> pairtools select+stats
#
# 该脚本是 script/a_1_run_50_1_50_2.sh 的参数化版本:去除了 /mnt/wsl/...
# 路径与 50_1/50_2 硬编码,改为 CLI flag 驱动。
#
# 用法示例:
#   bash cfizz/preprocessing/fastq2pairs.sh \
#       --samples sampleA sampleB \
#       --data-dir /path/to/1_1_rawdata \
#       --output-root /path/to/1_2_pairs_result \
#       --bwa-index /path/to/hg38.fa \
#       --chrom-sizes /path/to/hg38.chrom.exceptYMetc.sorted.sizes \
#       --assembly hg38 --tech-type hic --nproc 16
#
# 输入文件约定:<data-dir>/<sample_id>/<sample_id>_R1.fq.gz
#             <data-dir>/<sample_id>/<sample_id>_R2.fq.gz
# 输出布局:  <output-root>/<sample_id>/rawdata/<sample_id>/...
#           <output-root>/<sample_id>/pairtools/<sample_id>.nodups.UU.pairs.gz
#
# skip-if-exists 默认开启:每一步前若产物已存在则跳过;--no-skip 强制重跑。
#

set -euo pipefail

# ==================== 默认值 ====================
ASSEMBLY="hg38"
TECH_TYPE="hic"
NPROC=8
MEM_SORT="100G"
FASTP_EXTRA=()
OVERWRITE=0
NO_SKIP=0
MAX_PARALLEL=1

SAMPLES=()
DATA_DIR=""
OUTPUT_ROOT=""
BWA_IDX_PATH=""
SPCS_PATH=""

# ==================== usage ====================
usage() {
    cat <<'EOF'
fastq2pairs.sh — FASTQ -> 高质量 pairs 的 6 步流水线(参数化版本)

用法:
    bash cfizz/preprocessing/fastq2pairs.sh --samples <S1> [S2 ...] \
        --data-dir <DIR> --output-root <DIR> \
        --bwa-index <PREFIX> --chrom-sizes <FILE> \
        [--assembly hg38] [--tech-type hic|microc] \
        [--nproc N] [--mem-sort 100G] [--fastp-extra "..."] \
        [--max-parallel N] [--overwrite] [--no-skip]

必填参数:
    --samples         空格分隔的样本 ID 列表
    --data-dir        包含 <sample_id>/<sample_id>_R{1,2}.fq.gz 的根目录
    --output-root     输出根目录,会创建 <output-root>/<sample_id>/ 子目录
    --bwa-index       BWA-MEM2 索引前缀(如 /ref/hg38.fa)
    --chrom-sizes     染色体大小文件

可选参数:
    --assembly        装配名(默认 hg38)
    --tech-type       hic 或 microc(默认 hic)
    --nproc           CPU 线程数(默认 8)
    --mem-sort        pairtools sort 内存(默认 100G)
    --fastp-extra     附加 fastp 参数(空格分隔的多 token 用引号包起来)
    --max-parallel    同时处理的样本数(默认 1,bwa-mem2+sort 吃内存)
    --overwrite       已存在产物也覆盖
    --no-skip         同 --overwrite(语义别名)
    -h, --help        显示此帮助

行为说明:
    1. 6 步流水线:fastp -> bwa-mem2+samtools view -> pairtools parse
                  -> pairtools sort -> pairtools dedup -> pairtools select+stats
    2. skip-if-exists 默认开启:产物存在则跳过该步
    3. 并发:用 --max-parallel 控制同时跑多少个样本(谨慎 >1)
    4. 每次运行会输出时间戳日志
EOF
}

# ==================== 长选项解析 ====================
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --samples)
                shift
                while [[ $# -gt 0 && "${1:0:2}" != "--" ]]; do
                    SAMPLES+=("$1")
                    shift
                done
                ;;
            --samples=*)
                SAMPLES+=("${1#*=}")
                shift
                ;;
            --data-dir)
                DATA_DIR="$2"; shift 2 ;;
            --data-dir=*)
                DATA_DIR="${1#*=}"; shift ;;
            --output-root)
                OUTPUT_ROOT="$2"; shift 2 ;;
            --output-root=*)
                OUTPUT_ROOT="${1#*=}"; shift ;;
            --bwa-index)
                BWA_IDX_PATH="$2"; shift 2 ;;
            --bwa-index=*)
                BWA_IDX_PATH="${1#*=}"; shift ;;
            --chrom-sizes)
                SPCS_PATH="$2"; shift 2 ;;
            --chrom-sizes=*)
                SPCS_PATH="${1#*=}"; shift ;;
            --assembly)
                ASSEMBLY="$2"; shift 2 ;;
            --assembly=*)
                ASSEMBLY="${1#*=}"; shift ;;
            --tech-type)
                TECH_TYPE="$2"; shift 2 ;;
            --tech-type=*)
                TECH_TYPE="${1#*=}"; shift ;;
            --nproc)
                NPROC="$2"; shift 2 ;;
            --nproc=*)
                NPROC="${1#*=}"; shift ;;
            --mem-sort)
                MEM_SORT="$2"; shift 2 ;;
            --mem-sort=*)
                MEM_SORT="${1#*=}"; shift ;;
            --fastp-extra)
                shift
                while [[ $# -gt 0 && "${1:0:2}" != "--" ]]; do
                    FASTP_EXTRA+=("$1")
                    shift
                done
                ;;
            --fastp-extra=*)
                # 单 token 的扩展参数(空格分隔多个用 --fastp-extra 多次给)
                FASTP_EXTRA+=("${1#*=}")
                shift
                ;;
            --max-parallel)
                MAX_PARALLEL="$2"; shift 2 ;;
            --max-parallel=*)
                MAX_PARALLEL="${1#*=}"; shift ;;
            --overwrite|--no-skip)
                OVERWRITE=1; NO_SKIP=1; shift ;;
            -h|--help)
                usage; exit 0 ;;
            *)
                echo "未知 flag: $1" >&2
                echo "使用 --help 查看用法" >&2
                exit 1 ;;
        esac
    done
}

parse_args "$@"

# ==================== 必填校验 ====================
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "错误: 必须通过 --samples 指定至少一个样本" >&2
    exit 1
fi
if [[ -z "$DATA_DIR" ]]; then
    echo "错误: 必须通过 --data-dir 指定数据目录" >&2
    exit 1
fi
if [[ -z "$OUTPUT_ROOT" ]]; then
    echo "错误: 必须通过 --output-root 指定输出根目录" >&2
    exit 1
fi
if [[ -z "$BWA_IDX_PATH" ]]; then
    echo "错误: 必须通过 --bwa-index 指定 BWA-MEM2 索引前缀" >&2
    exit 1
fi
if [[ -z "$SPCS_PATH" ]]; then
    echo "错误: 必须通过 --chrom-sizes 指定染色体大小文件" >&2
    exit 1
fi

# ==================== 静态检查 ====================
if [[ ! -f "$SPCS_PATH" ]]; then
    echo "错误: 找不到染色体大小文件: $SPCS_PATH" >&2
    exit 1
fi
if [[ ! -f "${BWA_IDX_PATH}.amb" ]]; then
    echo "错误: 找不到 BWA-MEM2 索引文件: ${BWA_IDX_PATH}.amb" >&2
    exit 1
fi

# ==================== 根据技术类型设置 select condition ====================
case "$TECH_TYPE" in
    hic)
        SELECT_CONDITION="mapq1 >= 10 and mapq2 >= 10 and abs(pos2 - pos1) > 2000 and not (strand1 == '-' and strand2 == '+' and abs(pos2 - pos1) <= 2000) and chrom1 != '!' and chrom2 != '!'"
        ;;
    microc)
        SELECT_CONDITION="mapq1 >= 20 and mapq2 >= 20 and abs(pos2 - pos1) > 500 and len1 >= 50 and len1 <= 600 and len2 >= 50 and len2 <= 600 and chrom1 != '!' and chrom2 != '!'"
        ;;
    *)
        echo "错误: 未知 --tech-type: $TECH_TYPE(支持 hic / microc)" >&2
        exit 1
        ;;
esac

# ==================== 横幅 ====================
echo "=========================================="
echo "Hi-C 数据处理:从 FASTQ 到高质量 pairs"
echo "物种:        $ASSEMBLY"
echo "技术类型:    $TECH_TYPE"
echo "数据目录:    $DATA_DIR"
echo "输出根目录:  $OUTPUT_ROOT"
echo "CPU 线程数:  $NPROC"
echo "样本列表:    ${SAMPLES[*]}"
echo "并发样本数:  $MAX_PARALLEL"
echo "=========================================="
echo ""

# ==================== skip 检查辅助 ====================
# 根据 OVERWRITE/NO_SKIP 决定是否强制重跑
should_run() {
    local outfile="$1"
    if [[ "$OVERWRITE" -eq 1 ]]; then
        return 0
    fi
    if [[ -f "$outfile" ]]; then
        return 1
    fi
    return 0
}

# ==================== 单样本处理函数 ====================
process_sample() {
    local SAMPLE_ID="$1"

    echo "=========================================="
    echo "开始处理: $SAMPLE_ID"
    echo "=========================================="

    # 检查 FASTQ 文件
    if [[ ! -f "${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_R1.fq.gz" || \
          ! -f "${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_R2.fq.gz" ]]; then
        echo "错误: 找不到 FASTQ 文件" >&2
        echo "  期望: ${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_R1.fq.gz" >&2
        echo "  期望: ${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_R2.fq.gz" >&2
        return 1
    fi

    # 创建输出目录
    local OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE_ID}"
    mkdir -p "${OUTPUT_DIR}/rawdata/${SAMPLE_ID}"
    mkdir -p "${OUTPUT_DIR}/pairtools"

    echo "[$(date +"%Y-%m-%d %H:%M:%S")] 开始处理样本: $SAMPLE_ID"

    # ---------- Step 1: 质量控制 (fastp) ----------
    echo "Step 1: 质量控制 (fastp)..."
    if should_run "${OUTPUT_DIR}/rawdata/${SAMPLE_ID}/${SAMPLE_ID}_R1.fq.gz"; then
        fastp \
            -g \
            --detect_adapter_for_pe \
            -w "${NPROC}" \
            -i "${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_R1.fq.gz" \
            -I "${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_R2.fq.gz" \
            -o "${OUTPUT_DIR}/rawdata/${SAMPLE_ID}/${SAMPLE_ID}_R1.fq.gz" \
            -O "${OUTPUT_DIR}/rawdata/${SAMPLE_ID}/${SAMPLE_ID}_R2.fq.gz" \
            "${FASTP_EXTRA[@]}"
        echo "Step 1 完成"
    else
        echo "质量控制文件已存在,跳过"
    fi
    echo ""

    # ---------- Step 2: BWA-MEM2 比对 ----------
    echo "Step 2: BWA-MEM2 比对..."
    if should_run "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.bam"; then
        bwa-mem2 mem \
            -t "${NPROC}" \
            -SP \
            "${BWA_IDX_PATH}" \
            "${OUTPUT_DIR}/rawdata/${SAMPLE_ID}/${SAMPLE_ID}_R1.fq.gz" \
            "${OUTPUT_DIR}/rawdata/${SAMPLE_ID}/${SAMPLE_ID}_R2.fq.gz" \
            | samtools view -@ "${NPROC}" -bS > "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.bam"
        echo "Step 2 完成"
    else
        echo "BAM 文件已存在,跳过"
    fi
    echo ""

    # ---------- Step 3: 转换为 pairs (pairtools parse) ----------
    echo "Step 3: 转换为 pairs 格式 (pairtools parse)..."
    if should_run "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.pairs.gz"; then
        pairtools parse \
            -o "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.pairs.gz" \
            -c "${SPCS_PATH}" \
            --drop-sam \
            --drop-seq \
            --output-stats "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.stats" \
            --assembly "${ASSEMBLY}" \
            --no-flip \
            --add-columns mapq \
            --walks-policy mask \
            "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.bam"
        echo "Step 3 完成"
    else
        echo "pairs 文件已存在,跳过"
    fi
    echo ""

    # ---------- Step 4: 排序 pairs (pairtools sort) ----------
    echo "Step 4: 排序 pairs 文件 (pairtools sort)..."
    if should_run "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.sorted.pairs.gz"; then
        pairtools sort \
            --memory "${MEM_SORT}" \
            --nproc "${NPROC}" \
            -o "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.sorted.pairs.gz" \
            "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.pairs.gz"
        echo "Step 4 完成"
    else
        echo "排序文件已存在,跳过"
    fi
    echo ""

    # ---------- Step 5: 去重 (pairtools dedup) ----------
    echo "Step 5: 去重 (pairtools dedup)..."
    if should_run "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.dedup.stats"; then
        pairtools dedup \
            -p "${NPROC}" \
            --max-mismatch 3 \
            --mark-dups \
            --output \
                >( pairtools split \
                    --output-pairs "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.nodups.pairs.gz" \
                    --output-sam "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.nodups.bam" \
                 ) \
            --output-unmapped \
                >( pairtools split \
                    --output-pairs "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.unmapped.pairs.gz" \
                    --output-sam "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.unmapped.bam" \
                 ) \
            --output-dups \
                >( pairtools split \
                    --output-pairs "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.dups.pairs.gz" \
                    --output-sam "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.dups.bam" \
                 ) \
            --output-stats "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.dedup.stats" \
            "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.sorted.pairs.gz"
        echo "Step 5 完成"
    else
        echo "去重文件已存在,跳过"
    fi
    echo ""

    # ---------- Step 6: 选择高质量 pairs (pairtools select) ----------
    echo "Step 6: 选择高质量 pairs (pairtools select)..."
    if should_run "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.nodups.UU.pairs.gz"; then
        pairtools select \
            "$SELECT_CONDITION" \
            "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.nodups.pairs.gz" \
            -o "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.nodups.UU.pairs.gz"

        pairtools stats \
            "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.nodups.UU.pairs.gz" \
            -o "${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.select.stats"
        echo "Step 6 完成"
    else
        echo "高质量 pairs 文件已存在,跳过"
    fi
    echo ""

    echo "=========================================="
    echo "样本 $SAMPLE_ID 处理完成!"
    echo "最终 pairs 文件: ${OUTPUT_DIR}/pairtools/${SAMPLE_ID}.nodups.UU.pairs.gz"
    echo "=========================================="
    echo ""
}

# ==================== 并发控制 ====================
wait_for_slot() {
    local max_jobs="$1"
    while :; do
        local running
        running=$(jobs -r -p 2>/dev/null | wc -l)
        if (( running < max_jobs )); then
            break
        fi
        sleep 1
    done
}

# ==================== 主循环 ====================
for SAMPLE_ID in "${SAMPLES[@]}"; do
    if (( MAX_PARALLEL > 1 )); then
        wait_for_slot "$MAX_PARALLEL"
        ( process_sample "$SAMPLE_ID" ) &
    else
        process_sample "$SAMPLE_ID"
    fi
done

# 等所有并行任务结束
if (( MAX_PARALLEL > 1 )); then
    wait
fi

echo "所有样本处理完成!"
