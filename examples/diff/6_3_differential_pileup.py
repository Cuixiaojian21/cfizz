#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hiPSC_var vs hiPSC_nor 差异 Pileup(Loops APA + TAD boundary pileup)
================================================================

基于 6_1 差异分析结果,做 2 类差异 pileup:
1. Loops APA:差异 loops(hiPSC_var unique / hiPSC_nor unique)在 2 sample mcool 上的 APA
2. TAD boundary pileup:差异 boundaries(hiPSC_var unique / hiPSC_nor unique)在 2 sample mcool 上的 pileup

每类输出 2 张图(gain / lost),每张图 2 子图(hiPSC_var / hiPSC_nor)
预期:hiPSC_var unique loop 在 hiPSC_var 上富集(强中心),在 hiPSC_nor 上弱;反之亦然

样本:hiPSC_var(肿瘤), hiPSC_nor(正常)
限制:chr17
"""

import sys
import os
# 脚本位于 cfizz/examples/diff/<name>.py → 上 4 层到 cfizz/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import argparse
import logging
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

# cfizz imports
from cfizz.api.integrated.apa_pileup import plot_multi_apa_heatmap
from cfizz.viz.pileup import plot_multi_tad_boundary_pileup

# ============================================================================
# 路径配置(相对路径)
# ============================================================================

# 脚本位于 cfizz/examples/diff/<name>.py → 上 3 层 (到 /mnt/g/2_0_demo/cfizz/) 拼 demo/
DEMO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "demo")
DATA_DIR = os.path.join(DEMO_DIR, "data")
OUTPUT_ROOT = os.path.join(DEMO_DIR, "output")
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "5_1_primary_analysis_template")  # 5_1/5_2/6_x 共用

SAMPLES = {
    "hiPSC_var": os.path.join(DATA_DIR, "hiPSC_var_chr17.mcool"),  # 肿瘤
    "hiPSC_nor": os.path.join(DATA_DIR, "hiPSC_nor_chr17.mcool"),  # 正常
}

# 6_1 产物(只读)
DIFF_ROOT = os.path.join(OUTPUT_DIR, "4_differential", "hiPSC_var--hiPSC_nor")

# 6_3 产物(新增)
VIZ_ROOT = f"{DIFF_ROOT}/viz"

# === APA 参数(沿用 5_3) ===
APA_RESOLUTION = 10_000
APA_WINDOW = 7
APA_CORNER_SIZE = 5
APA_MIN_DISTANCE = 0  # 设为0,让所有6_1差异loops都参与APA(不过滤近距离loop)

# === TAD pileup 参数(沿用 5_3) ===
TAD_PILEUP_RESOLUTION = 10_000
TAD_PILEUP_FLANK = 300_000
TAD_PILEUP_TOP_N = 1000
TAD_PILEUP_WINDOWS = [5, 10, 50]  # 50kb/100kb/500kb


# ============================================================================
# 日志配置
# ============================================================================

def setup_logging(output_root: str):
    log_dir = Path(output_root) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / "6_3_pileup.log"
    logger = logging.getLogger("6_3_pileup")
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
# 数据读取函数
# ============================================================================

def read_diff_loops(diff_root, comparison, diff_type, logger) -> pd.DataFrame:
    """
    读差异 loops(gain = hiPSC_var unique, lost = hiPSC_nor unique)
    """
    t, c = comparison.split('--')  # t=hiPSC_var, c=hiPSC_nor
    loops_dir = Path(diff_root) / 'loops'
    if diff_type == 'gain':
        file = loops_dir / f'1_1_{t}_unique_loops.tsv'
    elif diff_type == 'lost':
        file = loops_dir / f'1_1_{c}_unique_loops_vs_{t}.tsv'
    else:
        raise ValueError(f"diff_type must be 'gain' or 'lost', got {diff_type}")
    if not file.exists():
        logger.warning(f"loops 文件不存在: {file}")
        return pd.DataFrame()
    df = pd.read_csv(file, sep=r'\s+', header=0)
    df = df[df['chrom1'] == 'chr17']
    logger.info(f"  {diff_type} loops: {len(df)} 个(chr17)")
    return df


def read_diff_tad_boundaries(diff_root, comparison, window_mult, diff_type, logger) -> pd.DataFrame:
    """
    读差异 TAD boundaries(Unique_boundary 类型,区分 hiPSC_var/hiPSC_nor unique)
    """
    t, c = comparison.split('--')  # t=hiPSC_var, c=hiPSC_nor
    tad_dir = Path(diff_root) / 'tad_boundary' / f'{window_mult}b'

    if diff_type == 'gain':  # hiPSC_var unique
        # 在 hiPSC_var_vs_hiPSC_nor 文件里,sample1=hiPSC_var 的 Unique_boundary = hiPSC_var 独有
        file = tad_dir / f'{t}_vs_{c}_boundary_classification_final.tsv'
    elif diff_type == 'lost':  # hiPSC_nor unique
        # 在 hiPSC_nor_vs_hiPSC_var 文件里,sample1=hiPSC_nor 的 Unique_boundary = hiPSC_nor 独有
        file = tad_dir / f'{c}_vs_{t}_boundary_classification_final.tsv'
    else:
        raise ValueError(f"diff_type must be 'gain' or 'lost', got {diff_type}")

    if not file.exists():
        logger.warning(f"TAD 文件不存在: {file}")
        return pd.DataFrame()

    df = pd.read_csv(file, sep='\t')
    df = df[df['chrom'] == 'chr17']
    df = df[df['final_classification'] == 'Unique_boundary']
    logger.info(f"  {diff_type} {window_mult}b boundaries: {len(df)} 个(chr17)")
    return df


# ============================================================================
# 辅助函数
# ============================================================================

def _rename_for_pileup(boundaries_df: pd.DataFrame) -> pd.DataFrame:
    """
    6_1 分类结果列名 → pileup 期望的 chrom/start/end
    内部 extract_pileup_snippets 期望这 3 列(viz/pileup.py:163-164)
    """
    return boundaries_df.rename(
        columns={'sample1_start': 'start', 'sample1_end': 'end'}
    )[['chrom', 'start', 'end']].copy()


# ============================================================================
# 绘图函数
# ============================================================================

def plot_diff_loops_apa(loops_df, diff_type, comparison, output_dir, logger):
    """
    差异 loops APA:
    - gain: hiPSC_var unique loops 在 hiPSC_var / hiPSC_nor 上的 APA
    - lost: hiPSC_nor unique loops 在 hiPSC_var / hiPSC_nor 上的 APA
    """
    if len(loops_df) == 0:
        logger.warning(f"  {diff_type} loops 为空,跳过")
        return None

    # 写单 sample loops 文件(给 plot_multi_apa_heatmap 喂)
    # 注意:read_loops_bedpe 用 header=None,所以不能写 header
    loop_tsv = Path(output_dir) / f"{diff_type}_loops.tsv"
    # 只写 6 列 BEDPE 格式,不含 header
    bedpe_cols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    loops_df[bedpe_cols].to_csv(loop_tsv, sep='\t', index=False, header=False)

    output_prefix = str(Path(output_dir) / f"loops_apa_{diff_type}")

    logger.info(f"  画 {diff_type} APA: 2 sample 各自的 mcool + 同一份 {diff_type} loops")

    try:
        plot_multi_apa_heatmap(
            mcool_paths=[SAMPLES['hiPSC_var'], SAMPLES['hiPSC_nor']],  # 2 sample
            loops_paths=[str(loop_tsv), str(loop_tsv)],                  # 同一份 diff loops
            output_path=output_prefix,
            sample_names=['hiPSC_var', 'hiPSC_nor'],
            resolution=APA_RESOLUTION,
            window=APA_WINDOW,
            corner_size=APA_CORNER_SIZE,
            min_distance=APA_MIN_DISTANCE,  # 0=不过滤
            vmin=None,  # 自动计算
            vmax=None,  # 自动计算
            balance=True,
            plot_size=4.0,  # 沿用 5_3 4cm
        )
        for ext in ['.png', '.svg']:
            output_file = Path(f"{output_prefix}{ext}")
            if output_file.exists():
                return output_file
        return None
    except Exception as e:
        import traceback
        logger.warning(f"  {diff_type} APA 失败: {e}\n{traceback.format_exc()}")
        return None


def plot_diff_tad_pileup(boundaries_df, diff_type, window_mult, comparison, output_dir, logger):
    """
    差异 TAD boundary pileup:
    - gain: hiPSC_var unique boundaries 在 hiPSC_var / hiPSC_nor mcool 上的 pileup
    - lost: hiPSC_nor unique boundaries 在 hiPSC_var / hiPSC_nor mcool 上的 pileup
    """
    if len(boundaries_df) == 0:
        logger.warning(f"  {diff_type} {window_mult}b boundaries 为空,跳过")
        return None

    # 列名映射:6_1 → pileup 期望
    boundaries_pileup = _rename_for_pileup(boundaries_df)
    # T-6.18: bioframe 要求 start/end 是 int,6_1 输出是 float → cast
    boundaries_pileup = boundaries_pileup.assign(
        start=boundaries_pileup['start'].astype(int),
        end=boundaries_pileup['end'].astype(int),
    )

    output_prefix = str(Path(output_dir) / f"tad_pileup_{diff_type}_{window_mult}b")

    logger.info(f"  画 {diff_type} {window_mult}b TAD pileup: 2 sample + 同一份 diff boundaries (OE+log2 风格)")

    try:
        # T-6.18: 算 per-sample P(s) 曲线(参 cooltools CTCF 教程 cell [18])
        import cooler as cooler_lib
        import cooltools
        from cfizz.analyze.compartment import get_view_df

        expected_dfs = []
        for sample_name in ['hiPSC_var', 'hiPSC_nor']:
            clr = cooler_lib.Cooler(f"{SAMPLES[sample_name]}::resolutions/{TAD_PILEUP_RESOLUTION}")
            view_df = get_view_df(clr)
            exp = cooltools.expected_cis(clr, view_df=view_df, nproc=8)
            expected_dfs.append(exp)

        # boundaries_list 是 DataFrame 列表(每 sample 一份)
        plot_multi_tad_boundary_pileup(
            mcool_paths=[SAMPLES['hiPSC_var'], SAMPLES['hiPSC_nor']],
            boundaries_list=[boundaries_pileup, boundaries_pileup],  # 列名映射后的 boundaries
            output_path=output_prefix,
            sample_names=['hiPSC_var', 'hiPSC_nor'],
            flank=TAD_PILEUP_FLANK,
            resolution=TAD_PILEUP_RESOLUTION,
            vmin=-1,                                  # T-6.18: OE 风格固定 vmin
            vmax=1,                                   # T-6.18: OE 风格固定 vmax
            balance=True,                             # T-6.18: OE 路径需要 balance 矩阵
            method='mean',
            color_scale='log2',                       # T-6.18: OE + log2
            color_scale_for_cbar='log2',              # T-6.18: cbar 刻度跟 heatmap 一致
            cmap='coolwarm',                          # T-6.18: OE 风格发散型 cmap
            plot_size=4.0,                            # 保持 4cm
            cbar_label='log2(obs/exp)',               # T-6.18: 跟 color_scale='log2' 一致
            expected_dfs=expected_dfs,                # T-6.18: per-sample P(s)
        )
        for ext in ['.png', '.svg']:
            output_file = Path(f"{output_prefix}{ext}")
            if output_file.exists():
                return output_file
        return None
    except Exception as e:
        import traceback
        logger.warning(f"  {diff_type} {window_mult}b TAD pileup 失败: {e}\n{traceback.format_exc()}")
        return None


# ============================================================================
# 主函数
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="hiPSC_var vs hiPSC_nor 差异 Pileup")
    parser.add_argument('--diff_root', type=str, default=DIFF_ROOT)
    parser.add_argument('--comparison', type=str, default='hiPSC_var--hiPSC_nor')
    parser.add_argument('--output_root', type=str, default=VIZ_ROOT)
    parser.add_argument('--diff_types', type=str, nargs='+', default=['gain', 'lost'])
    parser.add_argument('--tad_window_mults', type=int, nargs='+', default=TAD_PILEUP_WINDOWS)
    args = parser.parse_args()

    Path(args.output_root).mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output_root)

    print(f"\n{'='*70}")
    print(f"差异 Pileup(Loops APA + TAD boundary pileup)")
    print(f"比较组: {args.comparison}")
    print(f"输入: {args.diff_root}")
    print(f"输出: {args.output_root}")
    print(f"{'='*70}")

    # 1. Loops APA
    logger.info("=" * 60)
    logger.info("1. Loops APA 差异 Pileup")
    logger.info("=" * 60)
    for diff_type in args.diff_types:
        loops_df = read_diff_loops(args.diff_root, args.comparison, diff_type, logger)
        outdir = Path(args.output_root) / 'loops_apa' / diff_type
        outdir.mkdir(parents=True, exist_ok=True)
        plot_diff_loops_apa(loops_df, diff_type, args.comparison, outdir, logger)

    # 2. TAD boundary pileup
    logger.info("=" * 60)
    logger.info("2. TAD boundary 差异 Pileup")
    logger.info("=" * 60)
    for diff_type in args.diff_types:
        for window_mult in args.tad_window_mults:
            boundaries_df = read_diff_tad_boundaries(
                args.diff_root, args.comparison, window_mult, diff_type, logger
            )
            outdir = Path(args.output_root) / 'tad_pileup' / diff_type / f'{window_mult}b'
            outdir.mkdir(parents=True, exist_ok=True)
            plot_diff_tad_pileup(
                boundaries_df, diff_type, window_mult, args.comparison, outdir, logger
            )

    print(f"\n{'='*70}")
    print(f"✅ 差异 Pileup 完成!")
    print(f"输出: {args.output_root}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()