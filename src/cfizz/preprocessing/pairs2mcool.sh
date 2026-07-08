#!/bin/bash
#
# pairs2mcool.sh — 将 nodups.UU.pairs.gz 转为多分辨率 mcool + balance
#
# 流水线:cooler cload pairs -> cooler zoomify -> cooler balance(逐分辨率)
#
# 该脚本是 script/a_2_pairs2mcool_50_1_50_2.sh 的参数化版本:去除了
# /mnt/wsl/... 路径与 50_1/50_2 硬编码的 if/else 拆分,改为 CLI flag 驱动,
# 对任意样本列表通用。
#
# 用法示例:
#   bash cfizz/preprocessing/pairs2mcool.sh \
#       --samples sampleA sampleB \
#       --output-root /path/to/1_2_pairs_result \
#       --chrom-sizes /path/to/hg38.chrom.exceptYMetc.sorted.sizes \
#       --base-resolution 1000 \
#       --resolutions 5000,10000,25000,50000,100000,250000,500000,1000000 \
#       --nproc 8 --max-parallel 2
#
# 输入文件约定:<output-root>/<sample_id>/pairtools/<sample_id>.nodups.UU.pairs.gz
# 输出文件:      <output-root>/<sample_id>/<sample_id>_1000.cool
#              <output-root>/<sample_id>/<sample_id>_1000.mcool(内含多分辨率)
#

set -euo pipefail

# ==================== 默认值 ====================
BASE_RESOLUTION="1000"
RESOLUTIONS="5000,10000,25000,50000,100000,250000,500000,1000000"
NPROC=8
MAX_PARALLEL=1
OVERWRITE=0
NO_SKIP=0

SAMPLES=()
OUTPUT_ROOT=""
CHROM_SIZES=""

# ==================== usage ====================
usage() {
    cat <<'EOF'
pairs2mcool.sh — pairs.gz -> 多分辨率 mcool + balance(参数化版本)

用法:
    bash cfizz/preprocessing/pairs2mcool.sh --samples <S1> [S2 ...] \
        --output-root <DIR> --chrom-sizes <FILE> \
        [--base-resolution 1000] [--resolutions 5000,...,1000000] \
        [--nproc N] [--max-parallel N] [--overwrite] [--no-skip]

必填参数:
    --samples         空格分隔的样本 ID 列表
    --output-root     输出根目录(应与 fastq2pairs.sh 共享同一个)
    --chrom-sizes     染色体大小文件

可选参数:
    --base-resolution 基础分辨率(默认 1000)
    --resolutions     多分辨率层级,逗号分隔(默认 5000,10000,25000,50000,100000,250000,500000,1000000)
    --nproc           balance 线程(默认 8)
    --max-parallel    同时处理的样本数(默认 1;cooler cload 也吃内存)
    --overwrite       已存在产物也覆盖
    --no-skip         同 --overwrite
    -h, --help        显示此帮助

行为说明:
    1. 3 步流水线:cooler cload pairs -> cooler zoomify -> cooler balance
    2. skip-if-exists 默认开启:产物存在则跳过该步
    3. 输入文件:<output-root>/<sample_id>/pairtools/<sample_id>.nodups.UU.pairs.gz
       该文件应由 fastq2pairs.sh 跑出来
    4. 输出文件:<output-root>/<sample_id>/<sample_id>_<base-res>.cool
               <output-root>/<sample_id>/<sample_id>_<base-res>.mcool
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
            --output-root)
                OUTPUT_ROOT="$2"; shift 2 ;;
            --output-root=*)
                OUTPUT_ROOT="${1#*=}"; shift ;;
            --chrom-sizes)
                CHROM_SIZES="$2"; shift 2 ;;
            --chrom-sizes=*)
                CHROM_SIZES="${1#*=}"; shift ;;
            --base-resolution)
                BASE_RESOLUTION="$2"; shift 2 ;;
            --base-resolution=*)
                BASE_RESOLUTION="${1#*=}"; shift ;;
            --resolutions)
                RESOLUTIONS="$2"; shift 2 ;;
            --resolutions=*)
                RESOLUTIONS="${1#*=}"; shift ;;
            --nproc)
                NPROC="$2"; shift 2 ;;
            --nproc=*)
                NPROC="${1#*=}"; shift ;;
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
if [[ -z "$OUTPUT_ROOT" ]]; then
    echo "错误: 必须通过 --output-root 指定输出根目录" >&2
    exit 1
fi
if [[ -z "$CHROM_SIZES" ]]; then
    echo "错误: 必须通过 --chrom-sizes 指定染色体大小文件" >&2
    exit 1
fi

# ==================== 静态检查 ====================
if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "错误: 找不到染色体大小文件 $CHROM_SIZES" >&2
    exit 1
fi
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "错误: 未找到依赖命令: $1" >&2; exit 1; }; }
need_cmd pairtools
need_cmd cooler

# ==================== 横幅 ====================
echo "=== Hi-C 数据处理:pairs.gz -> mcool (balance) ==="
echo "输出根目录:    $OUTPUT_ROOT"
echo "染色体大小:    $CHROM_SIZES"
echo "基础分辨率:    ${BASE_RESOLUTION}bp"
echo "多分辨率:      $RESOLUTIONS"
echo "平衡线程数:    $NPROC"
echo "并发样本数:    $MAX_PARALLEL"
echo "样本列表:      ${SAMPLES[*]}"
echo ""

# ==================== skip 检查辅助 ====================
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
# 路径派生约定:pairtools 文件在 <output-root>/<sample_id>/pairtools/
#              cool / mcool 输出到 <output-root>/<sample_id>/(即同根的 sample 子目录)
process_sample() {
    local sample_id="$1"
    local input_dir="${OUTPUT_ROOT}/${sample_id}/pairtools"
    local output_dir="${OUTPUT_ROOT}/${sample_id}"
    local input_pairs="${input_dir}/${sample_id}.nodups.UU.pairs.gz"

    echo "=== 开始处理样本: $sample_id ==="

    # 检查输入文件
    if [[ ! -f "$input_pairs" ]]; then
        echo "错误: 找不到输入文件 $input_pairs" >&2
        return 1
    fi

    echo "[$sample_id] 输入文件: $input_pairs"

    # 确保 cool / mcool 输出目录存在
    mkdir -p "$output_dir"

    # ---------- 步骤 1: 生成 1000bp cool 文件 ----------
    local bins_spec="${CHROM_SIZES}:${BASE_RESOLUTION}"
    local out_cool="${output_dir}/${sample_id}_${BASE_RESOLUTION}.cool"

    if should_run "$out_cool"; then
        echo "[$sample_id] cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${bins_spec} -> ${out_cool}"
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 "$bins_spec" "$input_pairs" "$out_cool"
    else
        echo "[$sample_id] cool 文件已存在,跳过 cload"
    fi

    # ---------- 步骤 2: 生成 mcool(多分辨率) ----------
    local mcool_file="${output_dir}/${sample_id}_${BASE_RESOLUTION}.mcool"

    if should_run "$mcool_file"; then
        echo "[$sample_id] cooler zoomify --resolutions $RESOLUTIONS -o $mcool_file $out_cool"
        cooler zoomify --resolutions "$RESOLUTIONS" -o "$mcool_file" "$out_cool"
    else
        echo "[$sample_id] mcool 文件已存在,跳过 zoomify"
    fi

    # ---------- 步骤 3: balance(逐分辨率) ----------
    echo "[$sample_id] 开始 balance..."
    for path in $(cooler ls "$mcool_file"); do
        # balance 加 --force 才能在已 balance 上重做;但默认行为保留为:产物 <bin-size>.cool 上有
        # 一个 .bal 之类的副作用文件,cooler balance 检测已 balance 时会跳过。这里保持原版行为。
        echo "[$sample_id] 执行: cooler balance --force $path (已有则可能提示)"
        cooler balance --force "$path" 2>/dev/null || cooler balance "$path"
    done

    echo "[$sample_id] 完成!"
    echo "[$sample_id] 输出: $mcool_file"
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
echo "=== 阶段 1:生成 cool 和 mcool ==="
for sample_id in "${SAMPLES[@]}"; do
    if (( MAX_PARALLEL > 1 )); then
        wait_for_slot "$MAX_PARALLEL"
        ( process_sample "$sample_id" ) &
    else
        process_sample "$sample_id"
    fi
done

if (( MAX_PARALLEL > 1 )); then
    wait
fi

echo "=== 全部完成 ==="
echo "输出文件:"
for sample_id in "${SAMPLES[@]}"; do
    echo "  ${sample_id}: ${OUTPUT_ROOT}/${sample_id}/*.mcool"
done
