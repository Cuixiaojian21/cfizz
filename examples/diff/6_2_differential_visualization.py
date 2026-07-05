#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hiPSC_var vs hiPSC_nor 差异分析可视化(Compartment + TAD + Loops)
============================================================

基于 6_1 计算结果,绘制差异特征区域热图:
1. Compartment 差异区域(A_to_B / B_to_A) — 用 plot_multi_compartment + 5_2 已存 O/E 矩阵
2. TAD 边界差异(Unique / Boundary_shift),3 个窗口 — 用 quick_plot_integrated
3. Loops 差异(gain / lost) — 用 plot_multi_heatmap_with_loops

范围限制:chr17 范围内(O/E 矩阵只存了 chr17)
样本:hiPSC_var(肿瘤),hiPSC_nor(正常)
比较组:hiPSC_var--hiPSC_nor(默认)
"""

import sys
import os
# 脚本位于 cfizz/examples/diff/<name>.py → 上 4 层到 cfizz/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import argparse
import logging
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

# cfizz imports
from cfizz.utils.range_utils import (
    calc_smart_union_range,
    calc_anchor_range,
    calc_symmetric_range,
)
from cfizz.utils.bedpe_utils import (
    make_loop_region_id,
    write_single_loop_bedpe,
)

# viz — 全部用现有的库函数
from cfizz.viz.compartment import plot_multi_compartment, slice_compartment_region
from cfizz.viz.loop import plot_multi_heatmap_with_loops
from cfizz.api.integrated.quick_plot import quick_plot_integrated
from cfizz.api.integrated.heatmap_tracks import GenomeRange

# ============================================================================
# 路径配置(相对路径)
# ============================================================================

# 脚本位于 cfizz/examples/diff/<name>.py → 上 3 层 (到 /mnt/g/2_0_demo/cfizz/) 拼 demo/
DEMO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "demo")
DATA_DIR = os.path.join(DEMO_DIR, "data")
OUTPUT_ROOT = os.path.join(DEMO_DIR, "output")
# 5_1/5_2/6_x 共用一个 output 目录(6_x 读 5_1/5_2 产物)
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "5_1_primary_analysis_template")

SAMPLES = {
    "hiPSC_var": os.path.join(DATA_DIR, "hiPSC_var_chr17.mcool"),  # 肿瘤
    "hiPSC_nor": os.path.join(DATA_DIR, "hiPSC_nor_chr17.mcool"),  # 正常
}

SAMPLE_DISPLAY_NAMES = None  # 默认 None = 用 sample ID

def get_display_name(sample: str) -> str:
    if SAMPLE_DISPLAY_NAMES is None:
        return sample
    return SAMPLE_DISPLAY_NAMES.get(sample, sample)

# 6_1 产物(只读) — comparison 决定文件名
DIFF_ROOT = os.path.join(OUTPUT_DIR, "4_differential")
COMPARISON = "hiPSC_var--hiPSC_nor"  # 默认比较组(双横线)

# 5_2 已存的 O/E 矩阵(不重算)
VIZ_5_2_ROOT = os.path.join(OUTPUT_DIR, "2_visualization/compartment")

# 6_2 产物
VIZ_ROOT = f"{DIFF_ROOT}/{COMPARISON}/viz"

# === 限制范围:chr17(demo mcool 只含 chr17) ===
ALLOWED_CHROMS = ['chr17']

# === Compartment 配置 ===
COMPARTMENT_RESOLUTION = 100000  # 100kb
COMPARTMENT_HALF_SPAN_BINS = 45  # 上下游各 45 bins(45 × 100kb = 4.5Mb @ 100kb)

# === TAD 配置 ===
TAD_RESOLUTION = 10000  # 10kb
TAD_DEFAULT_FLANK_BP = 20 * TAD_RESOLUTION  # 200kb(默认 fallback)
TAD_WINDOW_MULTS = [5, 10, 50]  # 50kb / 100kb / 500kb

# === Loops 配置 ===
LOOP_RESOLUTION = 10000  # 10kb

# === 限制每个类型画图数量 ===
MAX_PER_TYPE = 5


# ============================================================================
# 日志配置
# ============================================================================

def setup_logging(output_root: str):
    log_dir = Path(output_root) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / "6_2_viz.log"
    logger = logging.getLogger("6_2_viz")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ============================================================================
# 读数据函数(comparison 拼路径,不要用 treatment_control 下划线)
# ============================================================================

def read_compartment_diff_regions(
    diff_root: str, comparison: str, diff_types, logger
) -> pd.DataFrame:
    """
    读 compartment 差异区域,合并为连续块,过滤 chr17

    路径用 comparison 拼(=`hiPSC_var--hiPSC_nor`),不是 `{t}_{c}` 下划线格式。
    """
    diff_file = Path(diff_root) / comparison / 'compartment' / f'compartment_{comparison}_diff.tsv'
    if not diff_file.exists():
        logger.warning(f"compartment diff 文件不存在: {diff_file}")
        return pd.DataFrame()

    df = pd.read_csv(diff_file, sep='\t')

    # 列名: 'compartment_change_hiPSC_var'
    change_col = [c for c in df.columns if 'compartment_change' in c][0]

    # 筛选差异类型
    df = df[df[change_col].isin(diff_types)]

    # 限制 chr17
    df = df[df['chrom'] == 'chr17']

    if len(df) == 0:
        logger.warning(f"compartment {diff_types}: 无 chr17 区域")
        return pd.DataFrame()

    # 合并连续区域
    region_groups = []
    for chrom in df['chrom'].unique():
        chrom_data = df[df['chrom'] == chrom].sort_values('start').reset_index(drop=True)
        start_idx = 0
        while start_idx < len(chrom_data):
            current_type = chrom_data.iloc[start_idx][change_col]
            end_idx = start_idx + 1
            while end_idx < len(chrom_data):
                next_row = chrom_data.iloc[end_idx]
                if (next_row[change_col] == current_type and
                        next_row['start'] == chrom_data.iloc[end_idx - 1]['end']):
                    end_idx += 1
                else:
                    break
            region_start = chrom_data.iloc[start_idx]['start']
            region_end = chrom_data.iloc[end_idx - 1]['end']
            center = (region_start + region_end) // 2
            region_groups.append({
                "chrom": chrom,
                "region_start": region_start,
                "region_end": region_end,
                "center": int(center),
                "change_type": current_type,
                "bin_count": end_idx - start_idx,
            })
            start_idx = end_idx

    regions_df = pd.DataFrame(region_groups)
    logger.info(f"compartment {diff_types}: {len(regions_df)} 个差异区域块(chr17)")
    return regions_df


def read_tad_diff_regions(
    diff_root: str, comparison: str, window_mult: int, logger
) -> pd.DataFrame:
    """
    读 TAD 边界差异区域,合并双向,过滤 chr17

    路径用 comparison 拼文件名:
        {comparison.split("--")[0]}_vs_{comparison.split("--")[1]}_boundary_classification_final.tsv
    """
    t, c = comparison.split('--')  # t='hiPSC_var', c='hiPSC_nor'
    tad_dir = Path(diff_root) / comparison / 'tad_boundary' / f'{window_mult}b'
    f1 = tad_dir / f'{t}_vs_{c}_boundary_classification_final.tsv'
    f2 = tad_dir / f'{c}_vs_{t}_boundary_classification_final.tsv'

    all_records = []
    for f, source in [(f1, f'{t}_vs_{c}'), (f2, f'{c}_vs_{t}')]:
        if not f.exists():
            logger.warning(f"TAD 文件不存在: {f}")
            continue
        df = pd.read_csv(f, sep='\t')
        df = df[df['chrom'] == 'chr17']
        df['source'] = source
        all_records.append(df)

    if not all_records:
        logger.warning(f"tad_boundary/{window_mult}b: 无 chr17 数据")
        return pd.DataFrame()

    all_df = pd.concat(all_records, ignore_index=True)

    # 合并双向记录
    results = []
    for idx, row in all_df.iterrows():
        plot_start, plot_end = calc_smart_union_range(
            s1_start=int(row['sample1_start']),
            s1_end=int(row['sample1_end']),
            s2_start=int(row['sample2_start']),
            s2_end=int(row['sample2_end']),
            s1_dist_up=float(row['sample1_dist_upstream']),
            s1_dist_down=float(row['sample1_dist_downstream']),
            s2_dist_up=float(row['sample2_dist_upstream']),
            s2_dist_down=float(row['sample2_dist_downstream']),
            default_flank=TAD_DEFAULT_FLANK_BP,
            resolution=TAD_RESOLUTION,
        )

        # 标签: final_classification 或 primary_label
        label = row.get('final_classification', row.get('primary_label', 'unknown'))

        results.append({
            'chrom': row['chrom'],
            'sample1_start': int(row['sample1_start']),
            'sample1_end': int(row['sample1_end']),
            'sample2_start': int(row['sample2_start']),
            'sample2_end': int(row['sample2_end']),
            'plot_start': plot_start,
            'plot_end': plot_end,
            'primary_label': label,
            'direction': row.get('direction', 'N/A'),
            'distance': row.get('distance', 0),
            'source': row['source'],
        })

    result_df = pd.DataFrame(results)
    logger.info(f"tad_boundary/{window_mult}b: {len(result_df)} 条记录(chr17)")
    return result_df


def read_loop_diff_regions(
    diff_root: str, comparison: str, diff_types, logger
) -> pd.DataFrame:
    """
    读 loops 差异区域(gain/lost),过滤 chr17

    路径用 comparison 拼:
        gain:  1_1_{treatment}_unique_loops.tsv
        lost:  1_1_{control}_unique_loops_vs_{treatment}.tsv
    """
    t, c = comparison.split('--')  # t='hiPSC_var', c='hiPSC_nor'
    loops_dir = Path(diff_root) / comparison / 'loops'
    dfs = []
    for diff_type in diff_types:
        if diff_type == 'gain':
            file = loops_dir / f'1_1_{t}_unique_loops.tsv'
        elif diff_type == 'lost':
            file = loops_dir / f'1_1_{c}_unique_loops_vs_{t}.tsv'
        else:
            continue
        if not file.exists():
            logger.warning(f"loops 文件不存在: {file}")
            continue
        df = pd.read_csv(file, sep='\t', header=0)
        df['diff_type'] = diff_type
        dfs.append(df)

    if not dfs:
        logger.warning(f"loops {diff_types}: 无数据")
        return pd.DataFrame()

    all_df = pd.concat(dfs, ignore_index=True)
    all_df = all_df[all_df['chrom1'] == 'chr17']
    logger.info(f"loops {diff_types}: {len(all_df)} 个 loops(chr17)")
    return all_df


# ============================================================================
# 绘图函数(照搬 5_2 compartment 逻辑,其余用库函数)
# ============================================================================

def plot_compartment_diff_region(
    row: pd.Series,
    comparison: str, output_dir, logger
):
    """
    画单个 compartment 差异区域热图(O/E + E1 柱状图)

    参数:
        row: 来自 read_compartment_diff_regions 的 1 行
             至少含 chrom / center / region_start / region_end
    """
    # 1. 从 row 拿 center 并扩展上下游各 25 bins(修 9 格 bug)
    center_bp = int(row['center'])
    chrom = str(row['chrom'])
    change_type = str(row['change_type'])
    start, end = calc_symmetric_range(
        center_bp,
        half_span_bins=COMPARTMENT_HALF_SPAN_BINS,  # 25
        resolution=COMPARTMENT_RESOLUTION,           # 100000
    )

    region_safe = f"{chrom}_{start}_{end}".replace(":", "_")
    output_prefix = str(Path(output_dir) / f"compartment_{region_safe}")

    logger.info(f"  compartment: {chrom}:{start}-{end} (center={center_bp}, {change_type})")

    try:
        results = []
        for idx, (sample_internal, mcool_path) in enumerate(SAMPLES.items()):
            # 1. eigenvector 路径(5_1 产物)
            eig_path = f"{OUTPUT_DIR}/1_computation/compartment/{sample_internal}/eigenvector.100k.tsv"
            if not os.path.exists(eig_path):
                logger.warning(f"    eigenvector 不存在: {eig_path}")
                continue

            # 2. O/E npy 路径(5_2 产物,直接 np.load)
            oe_npy_path = f"{VIZ_5_2_ROOT}/{sample_internal}/{sample_internal}_chr17_100k.oe.npy"
            if not os.path.exists(oe_npy_path):
                logger.warning(f"    O/E 文件不存在: {oe_npy_path}")
                continue

            # 3. 照搬 slice_compartment_region 逻辑(compartment.py:421)
            #    eig_df: bp 过滤;oe_matrix: bin 切片
            eig_df = pd.read_csv(eig_path, sep='\t')
            oe_matrix = np.load(oe_npy_path)

            # 过滤 eig_df(用 bp 位置)
            region_eig_df = eig_df[
                (eig_df['chrom'] == chrom) &
                (eig_df['start'] >= start) &
                (eig_df['end'] <= end)
            ].copy()

            # 切片 oe_matrix(用 bin 索引)
            start_bin = start // COMPARTMENT_RESOLUTION
            end_bin = end // COMPARTMENT_RESOLUTION

            # 防越界
            if end_bin > oe_matrix.shape[0]:
                end_bin = oe_matrix.shape[0]
            if start_bin < 0:
                start_bin = 0
            if start_bin >= oe_matrix.shape[0]:
                logger.warning(f"    起始位置超过矩阵范围")
                continue

            region_oe_matrix = oe_matrix[start_bin:end_bin, start_bin:end_bin]

            # 对齐
            if len(region_eig_df) != region_oe_matrix.shape[0]:
                min_len = min(len(region_eig_df), region_oe_matrix.shape[0])
                region_eig_df = region_eig_df.iloc[:min_len]
                region_oe_matrix = region_oe_matrix[:min_len, :min_len]

            # NaN 处理
            if region_eig_df['E1'].isna().any():
                region_eig_df['E1'] = region_eig_df['E1'].fillna(0)

            results.append({
                'status': 'success',
                'eig_df': region_eig_df,
                'oe_matrix': region_oe_matrix,
                'sample_name': get_display_name(sample_internal),
                'idx': idx,
            })

        if not results:
            logger.warning(f"    无有效数据")
            return None

        # 4. 画图(plot_multi_compartment 内部用 calculate_heatmap_layout)
        plot_multi_compartment(
            results=results,
            output_prefix=output_prefix,
            vmin=-2,
            vmax=2,
            plot_size=3.0,
            bar_height_ratio=0.3,
            start_pos=start,
            end_pos=end,
            chrom=chrom,
            group_name=change_type,
        )

        # 5. 找输出文件
        for ext in ['.png', '.svg']:
            output_file = Path(f"{output_prefix}{ext}")
            if output_file.exists():
                return output_file
        return None
    except Exception as e:
        import traceback
        logger.warning(f"    绘图失败: {e}\n{traceback.format_exc()}")
        return None


def plot_tad_diff_region(
    row: pd.Series, output_dir,
    resolution: int, tad_window_bp: int, logger
):
    """
    画 TAD 边界差异区域热图(2 sample 对比 + insulation 标注)

    quick_plot_integrated 内部已用 calculate_heatmap_layout,不需要额外 import。
    """
    start = row['plot_start']
    end = row['plot_end']
    label = row['primary_label']
    chrom = row['chrom']

    if end - start < tad_window_bp:
        logger.warning(f"  范围过小,跳过: {chrom}:{start}-{end}")
        return None

    region_safe = f"{chrom}_{start}_{end}".replace(':', '_')
    output_prefix = str(Path(output_dir) / f"tad_{region_safe}_{label}")

    logger.info(f"  TAD: {chrom}:{start}-{end} ({label})")

    try:
        # 构建 hics 列表(照搬 5_2 L290-305 multi 代码)
        hics = []
        for idx, (sample_internal, mcool_path) in enumerate(SAMPLES.items()):
            basename = Path(mcool_path).stem  # 'hiPSC_var_chr17'
            insulation_path = (
                f"{OUTPUT_DIR}/1_computation/tad/{sample_internal}/"
                f"1_0.{basename}.{resolution}.insulation.tsv"
            )
            flip = (idx == 1)  # 第二个 sample 翻转
            hics.append({
                'file': mcool_path,
                'triangle_ratio': 0.5,
                'cmap': 'Reds',
                'balance': False,
                'name': get_display_name(sample_internal),
                'flip_vertical': flip,
                'insulation_path': insulation_path,
                'window_size': tad_window_bp,
                'boundary_cmap': 'Blues',
                'boundary_alpha': 0.9,
            })

        quick_plot_integrated(
            hics=hics,
            tracks=[],
            region=GenomeRange(chrom, start, end),
            output=output_prefix,
            width_cm=8,
            gap_cm=0.2,
            left_margin_cm=1.0,
            right_margin_cm=2.0,
            dpi=300,
        )

        for ext in ['.png', '.svg']:
            output_file = Path(f"{output_prefix}{ext}")
            if output_file.exists():
                return output_file
        return None
    except Exception as e:
        import traceback
        logger.warning(f"    绘图失败: {e}\n{traceback.format_exc()}")
        return None


def plot_loop_diff_region(row: pd.Series, output_dir, logger):
    """
    画 loop 差异区域热图(2 sample 对比 + 1 个 loop 标注)
    """
    start, end = calc_anchor_range(
        int(row['start1']), int(row['end1']),
        int(row['start2']), int(row['end2']),
    )

    loop_id = make_loop_region_id(row)
    diff_type = row['diff_type']
    loop_tsv = Path(output_dir) / f"{diff_type}_{loop_id}_loops.tsv"
    write_single_loop_bedpe(row, loop_tsv)

    output_prefix = str(Path(output_dir) / f"{diff_type}_{loop_id}")

    logger.info(f"  Loops: {row['chrom1']}:{start}-{end} ({diff_type})")

    try:
        plot_multi_heatmap_with_loops(
            mcool_paths=list(SAMPLES.values()),
            loops_paths=[str(loop_tsv), str(loop_tsv)],
            chrom=str(row['chrom1']),
            start=start,
            end=end,
            resolution=LOOP_RESOLUTION,
            output_path=output_prefix,
            sample_names=[get_display_name(s) for s in SAMPLES.keys()],
            cmap='Reds',
            color_scale='linear',
            loop_color='blue',
            loop_alpha=0.6,
            loop_size=10,
            balance=True,
            dpi=1000,
        )
        for ext in ['.png', '.svg']:
            output_file = Path(f"{output_prefix}{ext}")
            if output_file.exists():
                return output_file
        return None
    except Exception as e:
        import traceback
        logger.warning(f"    绘图失败: {e}\n{traceback.format_exc()}")
        return None


# ============================================================================
# 主函数
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="hiPSC_var vs hiPSC_nor 差异分析可视化")
    parser.add_argument('--diff_root', type=str, default=DIFF_ROOT)
    parser.add_argument('--comparison', type=str, default=COMPARISON)
    parser.add_argument('--output_root', type=str, default=VIZ_ROOT)
    parser.add_argument('--max_per_type', type=int, default=MAX_PER_TYPE)
    parser.add_argument(
        '--diff_types_compartment', type=str, nargs='+',
        default=['A_to_B', 'B_to_A']
    )
    parser.add_argument(
        '--diff_types_tad', type=str, nargs='+',
        default=['Unique_boundary', 'Boundary_shift']
    )
    parser.add_argument(
        '--diff_types_loop', type=str, nargs='+',
        default=['gain', 'lost']
    )
    parser.add_argument('--tad_window_mults', type=int, nargs='+', default=TAD_WINDOW_MULTS)
    args = parser.parse_args()

    Path(args.output_root).mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output_root)

    print(f"\n{'='*70}")
    print(f"差异分析可视化(Compartment + TAD + Loops)")
    print(f"比较组: {args.comparison}")
    print(f"输入: {args.diff_root}/{args.comparison}")
    print(f"输出: {args.output_root}")
    print(f"{'='*70}")

    # 1. Compartment
    logger.info("=" * 60)
    logger.info("1. Compartment 差异可视化")
    logger.info("=" * 60)
    for diff_type in args.diff_types_compartment:
        regions = read_compartment_diff_regions(
            args.diff_root, args.comparison, [diff_type], logger
        )
        if len(regions) == 0:
            continue
        type_dir = Path(args.output_root) / 'compartment' / diff_type
        type_dir.mkdir(parents=True, exist_ok=True)
        for _, row in regions.head(args.max_per_type).iterrows():
            plot_compartment_diff_region(
                row=row,
                comparison=args.comparison,
                output_dir=type_dir,
                logger=logger,
            )

    # 2. TAD(3 窗口)
    logger.info("=" * 60)
    logger.info("2. TAD 边界差异可视化")
    logger.info("=" * 60)
    for wm in args.tad_window_mults:
        tad_window_bp = wm * TAD_RESOLUTION
        regions = read_tad_diff_regions(
            args.diff_root, args.comparison, wm, logger
        )
        if len(regions) == 0:
            continue
        for diff_type in args.diff_types_tad:
            subset = regions[regions['primary_label'] == diff_type] \
                if 'primary_label' in regions.columns else regions
            subset = subset.head(args.max_per_type)
            type_dir = Path(args.output_root) / 'tad_boundary' / f'{wm}b' / diff_type
            type_dir.mkdir(parents=True, exist_ok=True)
            for _, row in subset.iterrows():
                plot_tad_diff_region(
                    row=row,
                    output_dir=type_dir,
                    resolution=TAD_RESOLUTION,
                    tad_window_bp=tad_window_bp,
                    logger=logger,
                )

    # 3. Loops
    logger.info("=" * 60)
    logger.info("3. Loops 差异可视化")
    logger.info("=" * 60)
    for diff_type in args.diff_types_loop:
        regions = read_loop_diff_regions(
            args.diff_root, args.comparison, [diff_type], logger
        )
        if len(regions) == 0:
            continue
        type_dir = Path(args.output_root) / 'loops' / diff_type
        type_dir.mkdir(parents=True, exist_ok=True)
        for _, row in regions.head(args.max_per_type).iterrows():
            plot_loop_diff_region(row=row, output_dir=type_dir, logger=logger)

    print(f"\n{'='*70}")
    print(f"✅ 差异分析可视化完成!")
    print(f"输出: {args.output_root}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()