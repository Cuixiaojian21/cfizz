#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
差异分析(Compartment + TAD + Loops)
===================================

基于 5_1 主分析产物,做三类差异分析:
1. Compartment 差异(Stable_A/B/A_to_B/B_to_A)
2. TAD 边界差异(Unique/Boundary_shift/Stable),3 个窗口倍数
3. Loops 差异(gain/lost/common,HiCCUPS q-value 乘积)

每个类别输出:
- 主分析结果(diff_*.tsv)
- 统计图

样本:hiPSC_var(肿瘤),hiPSC_nor(正常)

Source: T-6.1 differential compute
"""

import sys
import os
# 脚本位于 cfizz/examples/diff/<name>.py → 上 4 层到 cfizz/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import argparse
import logging
from pathlib import Path
import pandas as pd
import numpy as np

# cfizz imports
from cfizz.analyze.compartment import (
    load_compartment_data,
    calculate_compartment_diff,
    save_plot_data as save_compartment_plot_data,
    load_plot_data as load_compartment_plot_data,
    plot_compartment_scatter,
    analyze_single_comparison as analyze_compartment,
    merge_diff_regions,
    COLORS_DICT,
)
from cfizz.analyze.tad import (
    calculate_boundary_context,
    generate_dual_direction_pairing,
    generate_final_classification,
    save_plot_data_tad,
    load_plot_data_tad,
    plot_tad_stacked_bar,
    analyze_single_comparison_window,
    TAD_COLORS,
    WINDOW_MULTIPLES,
    get_actual_window,
)
from cfizz.analyze.loop import (
    load_loops_data,
    find_matching_loops,
    calculate_sample_stats,
    save_plot_data_loops,
    load_plot_data_loops,
    plot_loops_stacked_bar,
    analyze_single_comparison_loops,
    LOOP_COLORS,
)

# ============================================================================
# 配置
# ============================================================================

# === demo 路径(相对路径) ===
# 脚本位于 cfizz/examples/diff/<name>.py → 上 3 层 (到 /mnt/g/2_0_demo/cfizz/) 拼 demo/
DEMO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "demo")
DATA_DIR = os.path.join(DEMO_DIR, "data")

# 样本 mcool 路径配置 (analyze 不感知具体路径，由 example 脚本管理)
# demo mcool: 50_1 → hiPSC_var (肿瘤), 50_2 → hiPSC_nor (正常)
SAMPLES = {
    "hiPSC_var": os.path.join(DATA_DIR, "hiPSC_var_chr17.mcool"),  # 肿瘤
    "hiPSC_nor": os.path.join(DATA_DIR, "hiPSC_nor_chr17.mcool"),  # 正常
}

# 可选: 显示名称映射。默认 None = 走最原汁原味(sample ID)
# 想要 2118465 风格就改成:
#   SAMPLE_DISPLAY_NAMES = {"hiPSC_var": "2118465_T", "hiPSC_nor": "2118465_N"}
SAMPLE_DISPLAY_NAMES = None


def get_display_name(sample: str) -> str:
    """获取样本的显示名称 — 默认用 sample ID,可选映射

    参数:
        sample: 样本 ID(如 "hiPSC_var")

    返回:
        显示名称 — 默认 = sample ID,配 SAMPLE_DISPLAY_NAMES 后 = 映射值
    """
    if SAMPLE_DISPLAY_NAMES is None:
        return sample
    return SAMPLE_DISPLAY_NAMES.get(sample, sample)

COMPARISONS = [("hiPSC_var", "hiPSC_nor")]  # (treatment, control) = (肿瘤, 正常)

# 5_1 产物根目录 (相对路径 → demo/output/5_1_primary_analysis_template/1_computation)
BASE_5_1 = Path(DEMO_DIR) / "output" / "5_1_primary_analysis_template" / "1_computation"

# TAD 窗口配置 - WINDOW_MULTIPLES 映射到实际文件名
# window_mult: (实际kb数, 文件中的kb字符串)
TAD_WINDOW_CONFIG = {
    50: ("500kb", "500kb"),
    10: ("100kb", "100kb"),
    5: ("50kb", "50kb"),
}

# ============================================================================
# 路径推导函数
# ============================================================================

def get_compartment_path(sample: str) -> Path:
    """获取 Compartment E1 文件路径"""
    return BASE_5_1 / "compartment" / sample / "eigenvector.100k.tsv"


# demo mcool 文件名是 hiPSC_var_chr17.mcool / hiPSC_nor_chr17.mcool
# process_tads 内部从 mcool 文件名推导 basename,5_1 实际产物用 {sample}_chr17 后缀
SAMPLE_BASENAMES = {
    "hiPSC_var": "hiPSC_var_chr17",
    "hiPSC_nor": "hiPSC_nor_chr17",
}


def get_tad_boundaries_path(sample: str, window_mult: int) -> Path:
    """获取 TAD boundaries 文件路径"""
    basename = SAMPLE_BASENAMES.get(sample, f"{sample}_1000")
    _, window_str = TAD_WINDOW_CONFIG[window_mult]
    return BASE_5_1 / "tad" / sample / f"2_0.{basename}.10000.{window_str}.boundaries.tsv"


def get_loop_path(sample: str) -> Path:
    """获取 Loops 文件路径"""
    basename = SAMPLE_BASENAMES.get(sample, f"{sample}_1000")
    return BASE_5_1 / "loop" / sample / f"{basename}.10k.loops.txt"


# ============================================================================
# 日志配置
# ============================================================================

def setup_logging(output_dir: Path, comparison: str):
    """配置日志"""
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"6_1_{comparison}.log"
    logger = logging.getLogger(f"6_1_{comparison}")
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
# Compartment 差异分析
# ============================================================================

def run_compartment_diff(treatment: str, control: str, comparison: str, output_dir: Path, logger: logging.Logger):
    """运行 Compartment 差异分析"""
    treatment_name = get_display_name(treatment)
    control_name = get_display_name(control)

    logger.info(f"\n{'='*60}")
    logger.info(f"1. Compartment 差异分析: {treatment_name} vs {control_name}")
    logger.info(f"{'='*60}")

    # 路径
    control_e1_path = get_compartment_path(control)
    treatment_e1_path = get_compartment_path(treatment)

    logger.info(f"  对照组: {control_e1_path}")
    logger.info(f"  处理组: {treatment_e1_path}")

    if not control_e1_path.exists():
        logger.error(f"  ❌ 对照组文件不存在: {control_e1_path}")
        return False
    if not treatment_e1_path.exists():
        logger.error(f"  ❌ 处理组文件不存在: {treatment_e1_path}")
        return False

    # 调用 analyze_compartment_diff
    compartment_output = output_dir / comparison / "compartment"
    compartment_output.mkdir(parents=True, exist_ok=True)

    analyze_compartment(
        comparison=comparison,
        output_root=str(output_dir / comparison),
        run_mode='all',
        control_e1_path=str(control_e1_path),
        treatment_e1_path=str(treatment_e1_path)
    )

    logger.info(f"  ✅ Compartment 差异分析完成")
    return True


# ============================================================================
# TAD 边界差异分析
# ============================================================================

def run_tad_diff(treatment: str, control: str, comparison: str, output_dir: Path, logger: logging.Logger):
    """运行 TAD 边界差异分析"""
    treatment_name = get_display_name(treatment)
    control_name = get_display_name(control)

    logger.info(f"\n{'='*60}")
    logger.info(f"2. TAD 边界差异分析: {treatment_name} vs {control_name}")
    logger.info(f"{'='*60}")

    for window_mult in WINDOW_MULTIPLES:
        actual_window = get_actual_window(window_mult)
        window_str = f"{window_mult}b"

        logger.info(f"\n  窗口 {window_str} ({actual_window//1000}kb)...")

        # 路径
        control_boundaries_path = get_tad_boundaries_path(control, window_mult)
        treatment_boundaries_path = get_tad_boundaries_path(treatment, window_mult)

        if not control_boundaries_path.exists():
            logger.warning(f"    ❌ 对照组文件不存在: {control_boundaries_path}")
            continue
        if not treatment_boundaries_path.exists():
            logger.warning(f"    ❌ 处理组文件不存在: {treatment_boundaries_path}")
            continue

        # 调用 analyze_tad_diff
        tad_output = output_dir / comparison / "tad_boundary"
        tad_output.mkdir(parents=True, exist_ok=True)

        analyze_single_comparison_window(
            comparison=comparison,
            window_mult=window_mult,
            output_root=str(tad_output),
            run_mode='all',
            treatment_boundaries_path=str(treatment_boundaries_path),
            control_boundaries_path=str(control_boundaries_path),
            treatment_name=treatment_name,
            control_name=control_name
        )

        logger.info(f"    ✅ 窗口 {window_str} 完成")

    return True


# ============================================================================
# Loops 差异分析
# ============================================================================

def run_loop_diff(treatment: str, control: str, comparison: str, output_dir: Path, logger: logging.Logger):
    """运行 Loops 差异分析"""
    treatment_name = get_display_name(treatment)
    control_name = get_display_name(control)

    logger.info(f"\n{'='*60}")
    logger.info(f"3. Loops 差异分析: {treatment_name} vs {control_name}")
    logger.info(f"{'='*60}")

    # 路径
    control_loop_path = get_loop_path(control)
    treatment_loop_path = get_loop_path(treatment)

    logger.info(f"  对照组: {control_loop_path}")
    logger.info(f"  处理组: {treatment_loop_path}")

    if not control_loop_path.exists():
        logger.error(f"  ❌ 对照组文件不存在: {control_loop_path}")
        return False
    if not treatment_loop_path.exists():
        logger.error(f"  ❌ 处理组文件不存在: {treatment_loop_path}")
        return False

    # 调用 analyze_loop_diff
    loop_output = output_dir / comparison / "loops"
    loop_output.mkdir(parents=True, exist_ok=True)

    analyze_single_comparison_loops(
        comparison=comparison,
        output_root=str(output_dir / comparison),
        run_mode='all',
        control_loops_path=str(control_loop_path),
        treatment_loops_path=str(treatment_loop_path)
    )

    logger.info(f"  ✅ Loops 差异分析完成")
    return True


# ============================================================================
# 主函数
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="差异分析(Compartment + TAD + Loops)")
    parser.add_argument('--output_dir', type=str, default=None,
                       help='输出目录,默认: demo/output/5_1_primary_analysis_template/4_differential')
    parser.add_argument('--run_mode', type=str, default='all', choices=['all', 'compute', 'plot'],
                       help='运行模式: all=计算+绘图, compute=只计算, plot=只绘图')
    args = parser.parse_args()

    # 默认输出路径(相对路径 → demo/output/...)
    if args.output_dir is None:
        output_dir = BASE_5_1.parent / "4_differential"
    else:
        output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}")
    print(f"差异分析(Compartment + TAD + Loops)")
    print(f"{'='*70}")
    print(f"输出目录: {output_dir}")
    print(f"运行模式: {args.run_mode}")
    print(f"比较组: {COMPARISONS}")
    print(f"{'='*70}")

    # 对每个比较组运行分析
    for treatment, control in COMPARISONS:
        treatment_name = get_display_name(treatment)
        control_name = get_display_name(control)
        comparison = f"{treatment}--{control}"  # 用 sample ID 拼,不受 SAMPLE_DISPLAY_NAMES 影响

        print(f"\n处理比较组: {comparison}")

        # 设置日志
        logger = setup_logging(output_dir, comparison)

        logger.info(f"=" * 60)
        logger.info(f"开始差异分析: {comparison}")
        logger.info(f"=" * 60)

        # 1. Compartment 差异
        try:
            run_compartment_diff(treatment, control, comparison, output_dir, logger)
        except Exception as e:
            logger.error(f"Compartment 差异分析失败: {e}")

        # 2. TAD 边界差异
        try:
            run_tad_diff(treatment, control, comparison, output_dir, logger)
        except Exception as e:
            logger.error(f"TAD 边界差异分析失败: {e}")

        # 3. Loops 差异
        try:
            run_loop_diff(treatment, control, comparison, output_dir, logger)
        except Exception as e:
            logger.error(f"Loops 差异分析失败: {e}")

        logger.info(f"\n{'='*60}")
        logger.info(f"差异分析完成!")
        logger.info(f"输出目录: {output_dir / comparison}")
        logger.info(f"{'='*60}")

    print(f"\n{'='*70}")
    print(f"✅ 全部差异分析完成!")
    print(f"输出目录: {output_dir}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()