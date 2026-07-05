#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example 5_2: Primary Analysis Visualization

基于 5_1 的计算结果，可视化 4 个特征：
  1. Heatmap (generate_heatmap / generate_multi_heatmap)
  2. Compartment (plot_compartment / generate_multi_compartment)
  3. TAD (plot_heatmap_with_tad_boundaries / plot_multi_heatmap_with_tad_boundaries)
  4. Loop (plot_heatmap_with_loops / plot_multi_heatmap_with_loops / quick_plot_integrated)

一句话解决问题(在 cfizz/ 根目录):
  python examples/5_2_primary_analysis_visualization.py

可视化区域(chr17 demo mcool):
  - Heatmap/TAD/Loop: chr17:10_000_000 - 12_000_000 (2Mb @ 10kb)
  - Compartment: chr17 全染色体 (0-83M, 100kb)

设计哲学:
  - 只用 viz / api 层，不动任何 analyze 层
  - 不导入 cfizz.io.paths (已放弃)
  - 不导入 hicviz 任何模块
  - 复用 5_1 产物，只算 O/E 矩阵自己算
"""

import sys
import os
# 脚本位于 cfizz/examples/<name>.py → 上 3 层 (到 /mnt/g/2_0_demo/)到 cfizz/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import time

import matplotlib
matplotlib.use('Agg')


# === 数据(相对路径) ===
DEMO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "demo")
DATA_DIR = os.path.join(DEMO_DIR, "data")
OUTPUT_ROOT = os.path.join(DEMO_DIR, "output")
# 5_1/5_2 共用一个 output 目录(因为 5_2 读 5_1 产物)
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "5_1_primary_analysis_template")

# FASTA: 直接用绝对路径 (5_2 暂不用,保留供未来扩展)
FASTA_PATH = "/mnt/g/2_0_demo/1_0_support/hg38.fa"  # ⚠️ 用户需自己准备 (3GB)

# === 双 sample (hiPSC demo) ===
SAMPLES = {
    "hiPSC_var": os.path.join(DATA_DIR, "hiPSC_var_chr17.mcool"),  # 肿瘤
    "hiPSC_nor": os.path.join(DATA_DIR, "hiPSC_nor_chr17.mcool"),  # 正常
}

# === 5_1 产物路径(5_1 已跑通,直接读) ===
COMP_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/compartment"
TAD_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/tad"
LOOP_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/loop"

# === 5_2 产物路径(本脚本新增) ===
VIZ_DIR = f"{OUTPUT_DIR}/2_visualization"


# === 辅助函数 ===

def viz_already_done(viz_dir: str, min_files: int = 1) -> bool:
    """
    检查 viz_dir 是否有 PNG/SVG 产物
    至少 min_files 个文件就认为已完成
    """
    if not os.path.isdir(viz_dir):
        return False
    files = [f for f in os.listdir(viz_dir)
             if f.endswith('.png') or f.endswith('.svg')]
    return len(files) >= min_files


# === 核心函数 ===

def visualize_heatmap(samples, viz_dir):
    """
    可视化基础热图(单 sample + 多 sample + quick_plot_integrated)
    范围:chr17:43.5M-43.7M, res=10kb
    """
    from cfizz.viz.heatmap import generate_heatmap, generate_multi_heatmap
    from cfizz.api.integrated.quick_plot import quick_plot_integrated
    from cfizz.api.integrated.heatmap_tracks import GenomeRange

    results = []
    chrom = "chr17"
    start_pos = 10_000_000   # 改为 chr17:10M-12M (2Mb)
    end_pos = 12_000_000
    resolution = 10_000

    # 单 sample
    for sample_name, mcool_path in samples.items():
        sample_viz_dir = f"{viz_dir}/heatmap/{sample_name}"
        os.makedirs(sample_viz_dir, exist_ok=True)
        start = time.time()
        try:
            generate_heatmap(
                file_path=mcool_path,
                chrom=chrom,
                resolution=resolution,
                start_pos=start_pos,
                end_pos=end_pos,
                output_dir=sample_viz_dir,
                balance=False,
                cmap='Reds',
                formats=('png', 'svg'),
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")

    # 多 sample
    multi_viz_dir = f"{viz_dir}/heatmap/multi"
    os.makedirs(multi_viz_dir, exist_ok=True)
    start = time.time()
    try:
        generate_multi_heatmap(
            file_paths=list(samples.values()),
            sample_names=list(samples.keys()),
            chrom=chrom,
            resolution=resolution,
            start_pos=start_pos,
            end_pos=end_pos,
            output_dir=multi_viz_dir,
            balance=False,
            cmap='Reds',
            plot_size=4, 
            formats=('png', 'svg'),
        )
        elapsed = time.time() - start
        results.append(("multi", 'success', elapsed))
        print(f"  ✅ multi: {elapsed:.1f}s")
    except Exception as e:
        elapsed = time.time() - start
        results.append(("multi", 'error', elapsed))
        print(f"  ❌ multi: {e}")

    # quick_plot_integrated (多 sample 整合图, n_tracks=0 只画热图)
    qp_viz_dir = f"{viz_dir}/heatmap/quick_plot_integrated"
    os.makedirs(qp_viz_dir, exist_ok=True)
    start = time.time()
    try:
        hics = []
        for i, (sample_name, mcool_path) in enumerate(samples.items()):
            flip = (i == 1)  # 第二个 sample 翻转对齐
            hics.append({
                'file': mcool_path,
                'name': sample_name,
                'cmap': 'Reds',
                'color_scale': 'linear',
                'balance': False,
                'resolution': resolution,
                'flip_vertical': flip,
                'triangle_ratio': 1,  # 全三角
            })
        output = f"{qp_viz_dir}/heatmap_quick_plot"
        quick_plot_integrated(
            hics=hics,
            region=GenomeRange(chrom=chrom, start=start_pos, end=end_pos),
            output=output,
            n_tracks=0,
            dpi=300,
            width_cm=8.0,
        )
        elapsed = time.time() - start
        results.append(("quick_plot_integrated", 'success', elapsed))
        print(f"  ✅ quick_plot_integrated: {elapsed:.1f}s")
    except Exception as e:
        elapsed = time.time() - start
        results.append(("quick_plot_integrated", 'error', elapsed))
        print(f"  ❌ quick_plot_integrated: {e}")

    return results


def visualize_compartment(samples, viz_dir, comp_dir):

    """
    可视化 Compartment(单 sample + 多 sample)
    范围:chr17 全染色体(0 - 83M), res=100kb
    流程:
      1. load_or_compute_oe_matrix(force_recompute=True) 算 O/E 矩阵(缓存 npy)
      2. plot_compartment 单 sample 画图(读路径)
      3. generate_multi_compartment 多 sample 画图(读路径列表)
    """
    from cfizz.viz.compartment import (
        plot_compartment,
        generate_multi_compartment,
    )
    from cfizz.analyze.oe import load_or_compute_oe_matrix

    results = []
    chrom = "chr17"
    start_pos = 0
    end_pos = 83_000_000  # chr17 全长(chr17 demo mcool 仅此 chr)
    resolution = 100_000

    # 准备多 sample 数据
    multi_eig_paths = []
    multi_oe_paths = []
    multi_sample_names = []

    # 单 sample
    for sample_name, mcool_path in samples.items():
        sample_viz_dir = f"{viz_dir}/compartment/{sample_name}"
        os.makedirs(sample_viz_dir, exist_ok=True)

        # Step A: 5_1 已算的 eigenvector 路径(只读不重算)
        eig_path = f"{comp_dir}/{sample_name}/eigenvector.100k.tsv"
        if not os.path.exists(eig_path):
            print(f"  ❌ {sample_name}: 5_1 eigenvector 缺失: {eig_path}")
            results.append((sample_name, 'error', 0))
            continue

        # Step B: 调 load_or_compute_oe_matrix(force_recompute=True)算 O/E 矩阵
        # 缓存到 sample_viz_dir(避免污染 5_1 的 1_computation 目录)
        start = time.time()
        try:
            oe_matrix, oe_metadata = load_or_compute_oe_matrix(
                mcool_path=mcool_path,
                chrom=chrom,
                resolution=resolution,
                output_dir=sample_viz_dir,
                sample_name=sample_name,
                balance=False,
                force_recompute=True,  # 5_2 自己算,5_1 没算,跟 2_1 无依赖
            )
            print(f"  [{sample_name}] O/E 矩阵 shape: {oe_matrix.shape}, "
                  f"elapsed: {time.time()-start:.1f}s")

            # Step C: oe.npy 路径(给 plot_compartment 用)
            oe_npy_path = os.path.join(
                sample_viz_dir,
                f"{sample_name}_{chrom}_{resolution//1000}k.oe.npy"
            )

            # Step D: 画图(plot_compartment 自己内部读路径 + 切片)
            plot_compartment(
                eig_tsv_path=eig_path,
                oe_npy_path=oe_npy_path,
                output_dir=sample_viz_dir,
                chrom=chrom,
                resolution=resolution,
                start_pos=start_pos,
                end_pos=end_pos,
                sample_name=sample_name,
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")

            # 累积到 multi
            multi_eig_paths.append(eig_path)
            multi_oe_paths.append(oe_npy_path)
            multi_sample_names.append(sample_name)
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")

    # 多 sample(generate_multi_compartment 读路径列表)
    if len(multi_eig_paths) >= 2:
        multi_viz_dir = f"{viz_dir}/compartment/multi"
        os.makedirs(multi_viz_dir, exist_ok=True)
        start = time.time()
        try:
            generate_multi_compartment(
                eig_tsv_paths=multi_eig_paths,
                oe_npy_paths=multi_oe_paths,
                output_dir=multi_viz_dir,
                sample_names=multi_sample_names,
                chrom=chrom,
                resolution=resolution,
                start_pos=start_pos,
                end_pos=end_pos,
                plot_size=4
            )
            elapsed = time.time() - start
            results.append(("multi", 'success', elapsed))
            print(f"  ✅ multi: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append(("multi", 'error', elapsed))
            print(f"  ❌ multi: {e}")

    return results


def visualize_tad(samples, viz_dir, tad_dir):
    """
    可视化 TAD(单 sample + 多 sample)
    范围:chr17:43.5M-44.5M (1Mb, 可看到 5-10 个 TAD), res=10kb, window=100kb

    使用 quick_plot_integrated + insulation_path (参考 3_2_tad_visualization.py)
    """
    from cfizz.api.integrated.quick_plot import quick_plot_integrated
    from cfizz.api.integrated.heatmap_tracks import GenomeRange

    results = []
    chrom = "chr17"
    start_pos = 10_000_000   # 改为 chr17:10M-12M (2Mb)
    end_pos = 12_000_000
    resolution = 10_000
    window_size = 100_000

    # 收集所有 sample 的 insulation paths
    valid_pairs = []
    for sample_name, mcool_path in samples.items():
        insulation_path = f"{tad_dir}/{sample_name}/1_0.{os.path.basename(mcool_path).split('.')[0]}.{resolution}.insulation.tsv"
        if os.path.exists(insulation_path):
            valid_pairs.append((sample_name, mcool_path, insulation_path))

    # 单 sample (triangle_ratio=1)
    for sample_name, mcool_path, insulation_path in valid_pairs:
        sample_viz_dir = f"{viz_dir}/tad/{sample_name}"
        os.makedirs(sample_viz_dir, exist_ok=True)

        start = time.time()
        try:
            output = f"{sample_viz_dir}/{sample_name}_tad"
            quick_plot_integrated(
                hics=[{
                    'file': mcool_path,
                    'triangle_ratio': 1,  # 全三角
                    'cmap': 'Reds',
                    'balance': False,
                    'name': sample_name,
                    'insulation_path': insulation_path,
                    'window_size': window_size,
                    'boundary_cmap': 'Blues_r',
                    'boundary_alpha': 0.9
                }],
                tracks=[],
                region=GenomeRange(chrom=chrom, start=start_pos, end=end_pos),
                output=output,
                width_cm=8,
                gap_cm=0.2,
                left_margin_cm=1.0,
                right_margin_cm=2.0,
                dpi=300
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")

    # 多 sample (各自用自己的 boundary, triangle_ratio=0.5 半三角)
    multi_viz_dir = f"{viz_dir}/tad/multi"
    os.makedirs(multi_viz_dir, exist_ok=True)
    
    if len(valid_pairs) >= 2:
        start = time.time()
        try:
            # hiPSC_var 不翻转, hiPSC_nor 翻转用于对齐比较, 配色统一用 Reds
            hics = []
            for i, (sample_name, mcool_path, insulation_path) in enumerate(valid_pairs):
                flip = (i == 1)  # 第二个 sample 翻转
                hics.append({
                    'file': mcool_path,
                    'triangle_ratio': 0.5,  # 半三角
                    'cmap': 'Reds',  # 配色统一
                    'balance': False,
                    'name': sample_name,  # 去掉 "(shared boundary)"
                    'flip_vertical': flip,
                    'insulation_path': insulation_path,  # 各自用自己的 boundary
                    'window_size': window_size,
                    'boundary_cmap': 'Blues_r',
                    'boundary_alpha': 0.9
                })
            
            output = f"{multi_viz_dir}/tad_multi"
            quick_plot_integrated(
                hics=hics,
                tracks=[],
                region=GenomeRange(chrom=chrom, start=start_pos, end=end_pos),
                output=output,
                width_cm=8,
                gap_cm=0.2,
                left_margin_cm=1.0,
                right_margin_cm=2.0,
                dpi=300
            )
            elapsed = time.time() - start
            results.append(("multi", 'success', elapsed))
            print(f"  ✅ multi: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append(("multi", 'error', elapsed))
            print(f"  ❌ multi: {e}")

    return results


def visualize_loop(samples, viz_dir, loop_dir):
    """
    可视化 Loop(单 sample + 多 sample + quick_plot_integrated)
    范围:chr17:43.5M-43.7M (200kb, 与 heatmap 同区域), res=10kb
    """
    from cfizz.viz.loop import (
        plot_heatmap_with_loops,
        plot_multi_heatmap_with_loops,
    )
    from cfizz.api.integrated.quick_plot import quick_plot_integrated
    from cfizz.api.integrated.heatmap_tracks import GenomeRange

    results = []
    chrom = "chr17"
    start_pos = 10_000_000   # 改为 chr17:10M-12M (2Mb)
    end_pos = 12_000_000
    resolution = 10_000

    # 单 sample
    for sample_name, mcool_path in samples.items():
        sample_viz_dir = f"{viz_dir}/loop/{sample_name}"
        os.makedirs(sample_viz_dir, exist_ok=True)

        # Loop 实际产物路径(固定用 10k 的 loops 文件)
        loops_path = f"{loop_dir}/{sample_name}/{os.path.basename(mcool_path).split('.')[0]}.10k.loops.txt"

        if not os.path.exists(loops_path):
            print(f"  ❌ {sample_name}: loops.txt 缺失: {loops_path}")
            results.append((sample_name, 'error', 0))
            continue

        start = time.time()
        try:
            output_path = f"{sample_viz_dir}/{sample_name}_loops"
            plot_heatmap_with_loops(
                mcool_path=mcool_path,
                loops_path=loops_path,
                output_path=output_path,
                chrom=chrom,
                start=start_pos,
                end=end_pos,
                resolution=resolution,
                loop_color='blue',
                loop_alpha=0.6,
                loop_size=2,
                balance=False,
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")

    # 多 sample
    multi_viz_dir = f"{viz_dir}/loop/multi"
    os.makedirs(multi_viz_dir, exist_ok=True)
    valid_pairs = []
    for sample_name, mcool_path in samples.items():
        loops_path = f"{loop_dir}/{sample_name}/{os.path.basename(mcool_path).split('.')[0]}.10k.loops.txt"
        if os.path.exists(loops_path) and os.path.getsize(loops_path) > 0:
            valid_pairs.append((sample_name, mcool_path, loops_path))

    if len(valid_pairs) >= 2:
        start = time.time()
        try:
            output_path = f"{multi_viz_dir}/loops_multi"
            plot_multi_heatmap_with_loops(
                mcool_paths=[p[1] for p in valid_pairs],
                loops_paths=[p[2] for p in valid_pairs],
                output_path=output_path,
                sample_names=[p[0] for p in valid_pairs],
                chrom=chrom,
                start=start_pos,
                end=end_pos,
                resolution=resolution,
                loop_color='blue',
                loop_alpha=0.6,
                loop_size=2,
                balance=False,
                plot_size=4,
            )
            elapsed = time.time() - start
            results.append(("multi", 'success', elapsed))
            print(f"  ✅ multi: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append(("multi", 'error', elapsed))
            print(f"  ❌ multi: {e}")

    # quick_plot_integrated(2 sample 整合图,n_tracks=0 只画热图+loop 标注)
    # hiPSC_var 不翻转, hiPSC_nor 翻转用于对齐比较
    qp_viz_dir = f"{viz_dir}/loop/quick_plot_integrated"
    os.makedirs(qp_viz_dir, exist_ok=True)
    if len(valid_pairs) >= 2:
        start = time.time()
        try:
            hics = []
            for i, (sample_name, mcool_path, loops_path) in enumerate(valid_pairs):
                flip = (i == 1)  # 第二个 sample 翻转对齐
                hics.append({
                    'file': mcool_path,
                    'name': sample_name,
                    'cmap': 'Reds',  # 配色统一
                    'color_scale': 'linear',
                    'balance': False,
                    'resolution': resolution,
                    'flip_vertical': flip,
                    'loops_path': loops_path,
                    'loop_color': 'blue',
                    'loop_alpha': 0.8,
                    'loop_size': 15,
                })
            output = f"{qp_viz_dir}/quick_plot_integrated"
            quick_plot_integrated(
                hics=hics,
                region=GenomeRange(chrom=chrom, start=start_pos, end=end_pos),
                output=output,
                n_tracks=0,
                dpi=3000,
                width_cm=8.0,
            )
            elapsed = time.time() - start
            results.append(("quick_plot_integrated", 'success', elapsed))
            print(f"  ✅ quick_plot_integrated: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append(("quick_plot_integrated", 'error', elapsed))
            print(f"  ❌ quick_plot_integrated: {e}")

    return results


# === 主函数 ===

def main():
    print("=" * 70)
    print("Example 5_2: Primary Analysis Visualization")
    print("=" * 70)
    print(f"  样本: {list(SAMPLES.keys())}")
    print(f"  输入: {OUTPUT_DIR}/1_computation/  (5_1 产物)")
    print(f"  输出: {VIZ_DIR}/")
    print()
    print("4 个块,13 个可视化(已有产物会跳过):")
    print("  1. Heatmap         (chr17:10-12M, 10kb)  3 个: hiPSC_var, hiPSC_nor, multi")
    print("  2. Compartment     (chr17 全染色体, 100kb)     3 个: hiPSC_var, hiPSC_nor, multi")
    print("  3. TAD             (chr17:43.5-44.5M, 10kb)  3 个: hiPSC_var, hiPSC_nor, multi")
    print("  4. Loop            (chr17:10-12M, 10kb)  4 个: hiPSC_var, hiPSC_nor, multi, quick_plot_integrated")
    print()

    os.makedirs(VIZ_DIR, exist_ok=True)
    overall_start = time.time()

    # Step 1: Heatmap
    print("=" * 70)
    print("Step 1: Heatmap 可视化")
    print("=" * 70)
    # 只有 multi 和 quick_plot_integrated 都存在时才跳过整个块
    if (viz_already_done(f"{VIZ_DIR}/heatmap/hiPSC_var") 
        and viz_already_done(f"{VIZ_DIR}/heatmap/hiPSC_nor")
        and viz_already_done(f"{VIZ_DIR}/heatmap/multi")
        and viz_already_done(f"{VIZ_DIR}/heatmap/quick_plot_integrated")):
        print("  ⏭️  heatmap 全部已有产物,跳过整块")
        heatmap_results = []
    else:
        heatmap_results = visualize_heatmap(SAMPLES, VIZ_DIR)

    # Step 2: Compartment
    print("\n" + "=" * 70)
    print("Step 2: Compartment 可视化")
    print("=" * 70)
    if viz_already_done(f"{VIZ_DIR}/compartment/multi"):
        print("  ⏭️  compartment/multi 已有产物,跳过整块")
        comp_results = []
    elif viz_already_done(f"{VIZ_DIR}/compartment/hiPSC_var") and viz_already_done(f"{VIZ_DIR}/compartment/hiPSC_nor"):
        print("  ⏭️  compartment/hiPSC_var + hiPSC_nor 已有产物,跳过单 sample")
        comp_results = []
    else:
        comp_results = visualize_compartment(SAMPLES, VIZ_DIR, COMP_COMPUTATION_DIR)

    # Step 3: TAD(库 bug 已修,必跑)
    print("\n" + "=" * 70)
    print("Step 3: TAD 可视化")
    print("=" * 70)
    if viz_already_done(f"{VIZ_DIR}/tad/multi") and viz_already_done(f"{VIZ_DIR}/tad/hiPSC_var") and viz_already_done(f"{VIZ_DIR}/tad/hiPSC_nor"):
        print("  ⏭️  tad 全部已有产物,跳过整块")
        tad_results = []
    else:
        tad_results = visualize_tad(SAMPLES, VIZ_DIR, TAD_COMPUTATION_DIR)

    # Step 4: Loop
    print("\n" + "=" * 70)
    print("Step 4: Loop 可视化")
    print("=" * 70)
    if (viz_already_done(f"{VIZ_DIR}/loop/quick_plot_integrated")
        and viz_already_done(f"{VIZ_DIR}/loop/multi")
        and viz_already_done(f"{VIZ_DIR}/loop/hiPSC_var")
        and viz_already_done(f"{VIZ_DIR}/loop/hiPSC_nor")):
        print("  ⏭️  loop 全部已有产物,跳过整块")
        loop_results = []
    else:
        loop_results = visualize_loop(SAMPLES, VIZ_DIR, LOOP_COMPUTATION_DIR)

    # 总结(显示 skipped vs success)
    total_elapsed = time.time() - overall_start
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    def count_success(results):
        return sum(1 for r in results if r[1] == 'success')

    def count_skipped(block_name, results):
        if block_name == "Heatmap":
            return "skipped" if not results else f"{count_success(results)}/{len(results)}"
        elif block_name == "Compartment":
            return "skipped" if not results else f"{count_success(results)}/{len(results)}"
        elif block_name == "TAD":
            return "skipped" if not results else f"{count_success(results)}/{len(results)}"
        elif block_name == "Loop":
            return "skipped" if not results else f"{count_success(results)}/{len(results)}"

    print(f"  Heatmap:     {count_skipped('Heatmap', heatmap_results)}")
    print(f"  Compartment: {count_skipped('Compartment', comp_results)}")
    print(f"  TAD:         {count_skipped('TAD', tad_results)}")
    print(f"  Loop:        {count_skipped('Loop', loop_results)}")
    print(f"\n  总耗时: {total_elapsed:.1f}s")
    print("\n  产物结构:")
    print(f"    {VIZ_DIR}/")
    print("      heatmap/{hiPSC_var, hiPSC_nor, multi}/")
    print("      compartment/{hiPSC_var, hiPSC_nor, multi}/")
    print("      tad/{hiPSC_var, hiPSC_nor, multi}/")
    print("      loop/{hiPSC_var, hiPSC_nor, multi, quick_plot_integrated}/")

    print("\n" + "=" * 70)
    print("✅ Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
