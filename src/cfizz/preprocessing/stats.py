#!/usr/bin/env python3
"""
Hi-C 数据处理流程统计脚本(参数化版本)

处理流程:
  1. pairtools parse  - BAM -> pairs.gz (生成 .stats)
  2. pairtools dedup  - 去重 (生成 .dedup.stats)
  3. pairtools select - 质量过滤 (生成 .select.stats)

该脚本是 script/a_1_z_processing_stats.py 的参数化版本:去除了
/mnt/wsl/... 路径与 50_1/50_2 硬编码,改为 argparse flag 驱动。

3 行汇总:parse / dedup / select,同时附原始总 pairs / 最终高质量 pairs / 总体得率汇总区。

可选 flag(默认行为 = --ppt: 同时写 .txt + .tsv 文件,stdout Tab 分隔):
    --output      指定时同时写 .txt + .tsv(.tsv 自动同名衍生)
                 不传 --output 时只打印到 stdout,不写文件
    --ppt / --no-ppt
                 --ppt 启用 PPT 友好的 Tab 分隔输出到 stdout(默认)
                 --no-ppt 切换为人类可读表格输出(stdout + 汇总区)

用法示例:
    # 默认(等价于 --ppt):指定 --output,写 .txt + .tsv,stdout 是 Tab 分隔
    python3 cfizz/preprocessing/stats.py \
        --pairtools-root /path/to/1_2_pairs_result \
        --samples 50_1 50_2 \
        --display-names 2118465_T 2118465_N \
        --output /path/to/processing_statistics.txt

    # 仅打印到 stdout,默认 --ppt 模式:Tab 分隔,直接可贴 PPT/Excel
    python3 cfizz/preprocessing/stats.py \
        --pairtools-root /path/to/1_2_pairs_result \
        --samples sampleA sampleB

    # 仅打印,人类可读格式(等同 --no-ppt,有汇总区)
    python3 cfizz/preprocessing/stats.py \
        --pairtools-root /path/to/1_2_pairs_result \
        --samples 50_1 50_2 --no-ppt
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional


# ==================== 参数解析 ====================
def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='stats.py',
        description='Hi-C 数据处理流程统计(parse / dedup / select 三步 + 汇总)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            '必填:\n'
            '  --pairtools-root + --samples\n'
            '可选:\n'
            '  --display-names (与 --samples 一一对应)\n'
            '  --output        (指定时写 .txt + .tsv;不指定则只 print 到 stdout)\n'
            '  --ppt / --no-ppt  (stdout 输出格式;默认 --ppt = Tab 分隔可贴 PPT)'
        ),
    )
    parser.add_argument(
        '--pairtools-root', type=str, required=True,
        help='包含 <sample_id>/pairtools/ 子目录的根目录',
    )
    parser.add_argument(
        '--samples', type=str, nargs='+', required=True,
        help='要统计的样本 ID 列表(nargs="+",空格分隔)',
    )
    parser.add_argument(
        '--display-names', type=str, nargs='+', default=None,
        help='与 --samples 一一对应的显示名,用于表格中的"样本"列',
    )
    parser.add_argument(
        '--output', type=str, default=None,
        help='输出 .txt 路径;不指定则只打印到 stdout,不写文件',
    )
    parser.add_argument(
        '--ppt', action=argparse.BooleanOptionalAction, default=True,
        help='启用 PPT 友好的 Tab 分隔输出到 stdout(默认开启;--no-ppt 关闭)',
    )
    return parser


# ==================== 工具函数 ====================
def read_pairtools_stats(stats_file: str) -> Dict[str, float]:
    """读取 pairtools 生成的统计文件。

    返回字段可能是 int 或 float,所以用 float 容器。
    """
    stats: Dict[str, float] = {}
    if not os.path.exists(stats_file):
        return stats

    with open(stats_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                key = parts[0]
                value: Optional[float] = None
                try:
                    value = int(parts[1])
                except ValueError:
                    try:
                        value = float(parts[1])
                    except ValueError:
                        value = None
                if value is not None:
                    stats[key] = value
    return stats


def get_sample_stats(pairtools_root: str, sample_id: str) -> Dict[str, Dict[str, float]]:
    """获取单个样本各处理步骤的统计数据。

    返回结构:
        {
            'parse': {...},   # pairtools parse 后的统计
            'dedup': {...},   # pairtools dedup 后的统计
            'select': {...},  # pairtools select 后的统计
        }
    """
    sample_dir = Path(pairtools_root) / sample_id / "pairtools"
    stats: Dict[str, Dict[str, float]] = {}

    parse_stats_file = sample_dir / f"{sample_id}.stats"
    if parse_stats_file.exists():
        stats['parse'] = read_pairtools_stats(str(parse_stats_file))

    dedup_stats_file = sample_dir / f"{sample_id}.dedup.stats"
    if dedup_stats_file.exists():
        stats['dedup'] = read_pairtools_stats(str(dedup_stats_file))

    select_stats_file = sample_dir / f"{sample_id}.select.stats"
    if select_stats_file.exists():
        stats['select'] = read_pairtools_stats(str(select_stats_file))

    return stats


def format_number(n: float) -> str:
    """格式化数字,添加千分位分隔符。"""
    try:
        return f"{int(n):,}"
    except (ValueError, TypeError):
        return f"{n:,}"


def calculate_yield(current: float, previous: float) -> float:
    """计算得率 (百分比)。"""
    if not previous:
        return 0.0
    return (current / previous) * 100


def resolve_display_name(sample_id: str, display_names: Optional[List[str]], samples: List[str]) -> str:
    """根据 --display-names 与 --samples 一一对应关系取显示名。

    若未给 --display-names,直接返回 sample_id。
    若给错长度,fallback 到 sample_id。
    """
    if display_names is None:
        return sample_id
    if len(display_names) != len(samples):
        return sample_id
    try:
        idx = samples.index(sample_id)
        return display_names[idx]
    except ValueError:
        return sample_id


# ==================== 数据生成 ====================
def generate_sample_summary(
    sample_id: str,
    stats: Dict[str, Dict[str, float]],
    display_name: str,
) -> List[Dict]:
    """生成单个样本的处理流程摘要(3 行:parse / dedup / select)。

    pairtools 统计文件的字段说明:
    - total: 总 pairs 数量(每个 read pair 算一个)
    - total_mapped: 成功比对到基因组两侧的 reads 数量
    - total_nodups: 去重后的 pairs 数量
    - total_dups: 被识别为重复的 pairs 数量
    """
    results: List[Dict] = []

    parse_stats = stats.get('parse', {})
    dedup_stats = stats.get('dedup', {})

    # ----- Step 1: pairtools parse (BAM -> pairs.gz) -----
    if dedup_stats:
        total_pairs = dedup_stats.get('total', 0)
        total_unmapped = dedup_stats.get('total_unmapped', 0)
        total_single = dedup_stats.get('total_single_sided_mapped', 0)
        total_mapped = dedup_stats.get('total_mapped', 0)

        mapping_rate = calculate_yield(total_mapped, total_pairs) if total_pairs else 0
        single_rate = calculate_yield(total_single, total_pairs) if total_pairs else 0
        unmapped_rate = calculate_yield(total_unmapped, total_pairs) if total_pairs else 0

        results.append({
            'sample': display_name,
            'step': '1_pairtools_parse',
            'step_name': 'pairs.gz (parse)',
            'input_count': total_pairs,
            'output_count': total_mapped,
            'yield': mapping_rate,
            'description': (
                f'两端比对: {format_number(total_mapped)} ({mapping_rate:.1f}%) | '
                f'单侧: {format_number(total_single)} ({single_rate:.1f}%) | '
                f'未比对: {format_number(total_unmapped)} ({unmapped_rate:.1f}%)'
            ),
        })

    # ----- Step 2: pairtools dedup (去重) -----
    if dedup_stats:
        total_mapped = dedup_stats.get('total_mapped', 0)
        total_nodups = dedup_stats.get('total_nodups', 0)
        total_dups = dedup_stats.get('total_dups', 0)

        results.append({
            'sample': display_name,
            'step': '2_pairtools_dedup',
            'step_name': '去重 (dedup)',
            'input_count': total_mapped,
            'output_count': total_nodups,
            'yield': calculate_yield(total_nodups, total_mapped) if total_mapped else 0,
            'description': (
                f'去除重复 pairs: {format_number(total_dups)} '
                f'({calculate_yield(total_dups, total_mapped):.1f}%)'
            ),
        })

    # ----- Step 3: pairtools select (质量过滤) -----
    select_stats = stats.get('select', {})
    if select_stats:
        dedup_nodups = dedup_stats.get('total_nodups', 0)
        final_pairs = select_stats.get('total', 0)

        results.append({
            'sample': display_name,
            'step': '3_pairtools_select',
            'step_name': '质量过滤 (select)',
            'input_count': dedup_nodups,
            'output_count': final_pairs,
            'yield': calculate_yield(final_pairs, dedup_nodups) if dedup_nodups else 0,
            'description': (
                f'筛选高质量 pairs (保留 {calculate_yield(final_pairs, dedup_nodups):.1f}%)'
            ),
        })

    return results


# ==================== 渲染 ====================
def get_sample_original_totals(all_results: List[Dict]) -> Dict[str, float]:
    """提取每个样本的第一行 input_count 作为原始总数。"""
    samples_data: Dict[str, float] = {}
    for r in all_results:
        sample = r['sample']
        if sample not in samples_data:
            samples_data[sample] = r['input_count']
    return samples_data


def print_table(all_results: List[Dict], for_ppt: bool = False) -> None:
    """打印格式化表格(for_ppt=True 时输出 Tab 分隔,可贴 Excel/PPT)。"""
    samples_data = get_sample_original_totals(all_results)

    if for_ppt:
        print("样本\t步骤\t输入 pairs\t输出 pairs\t得率(%)\t累计占比(%)\t描述")
        for r in all_results:
            original_total = samples_data.get(r['sample'], 0)
            cumulative_pct = calculate_yield(r['output_count'], original_total) if original_total else 0
            print(
                f"{r['sample']}\t{r['step_name']}\t{format_number(r['input_count'])}\t"
                f"{format_number(r['output_count'])}\t{r['yield']:.2f}%\t"
                f"{cumulative_pct:.2f}%\t{r['description']}"
            )
    else:
        header = (
            f"{'样本':<15} {'步骤':<20} {'输入 pairs':>15} {'输出 pairs':>15} "
            f"{'得率(%)':>10} {'累计占比(%)':>12} | 描述"
        )
        print("=" * 135)
        print(header)
        print("=" * 135)

        current_sample = None
        for r in all_results:
            sample = r['sample']
            original_total = samples_data.get(sample, 0)
            cumulative_pct = calculate_yield(r['output_count'], original_total) if original_total else 0

            if sample != current_sample:
                current_sample = sample
                print("-" * 135)

            row = (
                f"{r['sample']:<15} {r['step_name']:<20} "
                f"{format_number(r['input_count']):>15} {format_number(r['output_count']):>15} "
                f"{r['yield']:>10.2f}% {cumulative_pct:>12.2f}% | {r['description']}"
            )
            print(row)

        print("=" * 135)


def print_summary(all_results: List[Dict]) -> None:
    """打印汇总区:每个样本的原始总 pairs / 最终高质量 pairs / 总体得率。"""
    if not all_results:
        return

    print()
    print("=" * 120)
    print("汇总统计")
    print("=" * 120)
    print()

    seen_samples = []
    for sample in (r['sample'] for r in all_results):
        if sample not in seen_samples:
            seen_samples.append(sample)

    for sample in seen_samples:
        sample_results = [r for r in all_results if r['sample'] == sample]
        if not sample_results:
            continue
        original_total = sample_results[0]['input_count']
        final_result = sample_results[-1]
        overall_yield = calculate_yield(final_result['output_count'], original_total)

        print(f"样本 {sample}:")
        print(f"  - 原始总 pairs: {format_number(original_total)} pairs")
        print(f"  - 最终高质量 pairs: {format_number(final_result['output_count'])} pairs")
        print(f"  - 总体得率: {overall_yield:.2f}%")
        print()


# ==================== 文件输出 ====================
def save_table_to_file(all_results: List[Dict], output_file: str) -> None:
    """保存表格到文件(.txt 人类可读 + .tsv 可贴 Excel/PPT,两者一起写)。

    - .txt 包含表格 + 汇总区
    - .tsv 文件名 = output_file 后缀从 .txt 替换为 .tsv(否则追加 .tsv 后缀)
    """
    samples_data = get_sample_original_totals(all_results)

    # 1. 保存普通文本格式(.txt)— 总是写
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 135 + "\n")
        f.write("Hi-C 数据处理流程统计表\n")
        f.write("=" * 135 + "\n\n")

        header = (
            f"{'样本':<15} {'步骤':<20} {'输入 pairs':>15} {'输出 pairs':>15} "
            f"{'得率(%)':>10} {'累计占比(%)':>12} | 描述\n"
        )
        f.write(header)
        f.write("-" * 135 + "\n")

        current_sample = None
        for r in all_results:
            sample = r['sample']
            original_total = samples_data.get(sample, 0)
            cumulative_pct = calculate_yield(r['output_count'], original_total) if original_total else 0

            if sample != current_sample:
                current_sample = sample
                f.write("-" * 135 + "\n")

            row = (
                f"{r['sample']:<15} {r['step_name']:<20} "
                f"{format_number(r['input_count']):>15} {format_number(r['output_count']):>15} "
                f"{r['yield']:>10.2f}% {cumulative_pct:>12.2f}% | {r['description']}\n"
            )
            f.write(row)

        f.write("=" * 135 + "\n")

        # 汇总区
        f.write("\n" + "=" * 120 + "\n")
        f.write("汇总统计\n")
        f.write("=" * 120 + "\n\n")

        seen_samples = []
        for sample in (r['sample'] for r in all_results):
            if sample not in seen_samples:
                seen_samples.append(sample)

        for sample in seen_samples:
            sample_results = [r for r in all_results if r['sample'] == sample]
            if not sample_results:
                continue
            original_total = sample_results[0]['input_count']
            final_result = sample_results[-1]
            overall_yield = calculate_yield(final_result['output_count'], original_total)

            f.write(f"样本 {sample}:\n")
            f.write(f"  - 原始总 pairs: {format_number(original_total)} pairs\n")
            f.write(f"  - 最终高质量 pairs: {format_number(final_result['output_count'])} pairs\n")
            f.write(f"  - 总体得率: {overall_yield:.2f}%\n\n")

    print(f"表格已保存至: {output_file}")

    # 2. TSV 格式(可贴 Excel/PPT)— 总是写
    base, ext = os.path.splitext(output_file)
    if ext.lower() == '.txt':
        tsv_file = base + '.tsv'
    else:
        tsv_file = output_file + '.tsv'

    with open(tsv_file, 'w', encoding='utf-8') as f:
        f.write("样本\t步骤\t输入 pairs\t输出 pairs\t得率(%)\t累计占比(%)\t描述\n")
        for r in all_results:
            original_total = samples_data.get(r['sample'], 0)
            cumulative_pct = calculate_yield(r['output_count'], original_total) if original_total else 0
            f.write(
                f"{r['sample']}\t{r['step_name']}\t{format_number(r['input_count'])}\t"
                f"{format_number(r['output_count'])}\t{r['yield']:.2f}%\t"
                f"{cumulative_pct:.2f}%\t{r['description']}\n"
            )

    print(f"TSV 格式已保存至: {tsv_file}")
    print(f"  -> 可复制到 Excel/PPT,直接粘贴即可形成表格")


# ==================== main ====================
def main() -> int:
    parser = build_argparser()
    args = parser.parse_args()

    # 校验 --samples / --display-names 长度一致
    if args.display_names is not None and len(args.display_names) != len(args.samples):
        print(
            f"错误: --display-names 数量({len(args.display_names)})与 "
            f"--samples 数量({len(args.samples)})不一致",
            file=sys.stderr,
        )
        return 1

    print("=" * 70)
    print("Hi-C 数据处理流程统计")
    print("=" * 70)
    print(f"pairtools 根目录: {args.pairtools_root}")
    print(f"统计样本: {', '.join(args.samples)}")
    if args.display_names:
        print(f"显示名映射: {dict(zip(args.samples, args.display_names))}")
    print()

    all_results: List[Dict] = []

    for sample_id in args.samples:
        print(f"处理样本: {sample_id}")
        stats = get_sample_stats(args.pairtools_root, sample_id)

        if not stats:
            print(f"  ⚠ 未找到统计文件,跳过")
            continue

        display_name = resolve_display_name(sample_id, args.display_names, args.samples)
        sample_results = generate_sample_summary(sample_id, stats, display_name)
        all_results.extend(sample_results)

        print(f"  ✓ 找到 {len(stats)} 个统计文件")

    if not all_results:
        print("\n❌ 未找到任何统计数据")
        return 1

    # stdout 输出格式:默认 --ppt = Tab 分隔可贴 PPT;--no-ppt = 人类可读 + 汇总区
    print()
    print_table(all_results, for_ppt=args.ppt)
    if not args.ppt:
        print_summary(all_results)

    # 文件输出(只要指定了 --output 就会同时写 .txt + .tsv)
    if args.output:
        save_table_to_file(all_results, args.output)

    return 0


if __name__ == '__main__':
    sys.exit(main())
