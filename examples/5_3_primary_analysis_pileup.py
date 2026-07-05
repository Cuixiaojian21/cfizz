#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example 5_3: Primary Analysis Pileup - 3 类累计分析

只做累计分析，3 类:
1. Saddle plot (A/B compartment 交互)
2. TAD pileup (TAD 边界附近信号堆叠)
3. APA (Aggregate Peak Analysis, loop 中心信号堆叠)

一句话解决问题(在 cfizz/ 根目录):
  python examples/5_3_primary_analysis_pileup.py

设计哲学:
  - 只用 viz / api 层，不动任何 analyze 层
  - 不导入 cfizz.io.paths (已放弃)
  - 不导入 hicviz 任何模块
  - Saddle: 读 5_1 已算的全基因组 eigenvector.tsv
  - TAD pileup: 读 5_1 已算的 insulation.tsv (不重算!)
  - APA: 读 5_1 已算的 loops.txt (不重算!)
  - 区域 chr17 demo: SADDLE 全染色体, TAD pileup +43.5M-44.5M, APA +43.5M-43.7M
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
# 5_1/5_2/5_3 共用一个 output 目录(因为 5_3 读 5_1/5_2 产物)
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "5_1_primary_analysis_template")

# FASTA: 直接用绝对路径 (5_3 saddle 算 oe 矩阵需要 fasta)
FASTA_PATH = "/mnt/g/2_0_demo/1_0_support/hg38.fa"  # ⚠️ 用户需自己准备 (3GB)

# === 双 sample (hiPSC demo) ===
SAMPLES = {
    "hiPSC_var": os.path.join(DATA_DIR, "hiPSC_var_chr17.mcool"),  # 肿瘤
    "hiPSC_nor": os.path.join(DATA_DIR, "hiPSC_nor_chr17.mcool"),  # 正常
}

# === 累计分析输出路径(本脚本新增) ===
PILEUP_DIR = f"{OUTPUT_DIR}/3_pileup"


# === Saddle plot 参数 ===
SADDLE_CHROM = "chr17"  # chr17 demo mcool 只有 chr17
SADDLE_RESOLUTION = 100_000  # 跟 5_1 eigenvector 同 res
SADDLE_N_BINS = 98
SADDLE_VMIN = -1
SADDLE_VMAX = 1
SADDLE_NPROC = 8


# === TAD pileup 参数 ===
TAD_PILEUP_CHROM = "chr17"
TAD_PILEUP_START = 43_500_000
TAD_PILEUP_END = 44_500_000  # 跟 5_2 TAD 可视化同区域 (1Mb, 可堆叠 5-10 个 TAD)
TAD_PILEUP_RESOLUTION = 10_000
TAD_PILEUP_WINDOW = 100_000  # 用 100kb window(跟 5_2 同)
TAD_PILEUP_FLANK = 300_000
TAD_PILEUP_TOP_N = 1000
TAD_PILEUP_DPI = 300


# === APA 参数 ===
APA_CHROM = "chr17"
APA_START = 43_500_000
APA_END = 43_700_000  # 跟 5_2 Loop 可视化同区域 (200kb)
APA_RESOLUTION = 10_000
APA_WINDOW = 7
APA_CORNER_SIZE = 5
APA_MIN_DISTANCE = 20
APA_VMIN = 0
APA_VMAX = 0.01


def compute_saddle(samples, output_dir, n_bins, vmin, vmax, nproc):
    """
    Saddle plot (A/B compartment 强度矩阵)
    流程:
      1. 读 5_1 已算的全基因组 eigenvector.tsv (必需!)
      2. 调 generate_single_saddle / generate_multi_saddle 算 saddle
    """
    from cfizz.api.integrated.saddle_plot import (
        generate_single_saddle,
        generate_multi_saddle,
    )

    # 5_1 Compartment 产物路径(读已算的全基因组 eigenvector)
    COMP_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/compartment"

    # 单 sample heatmap_size = 4cm
    # multi heatmap_size = 4cm (用户要求所有累计分析图统一 4cm)
    HEATMAP_SIZE_SINGLE = 4.0
    HEATMAP_SIZE_MULTI = 4.0

    results = []

    # 单 sample
    for sample_name, mcool_path in samples.items():
        sample_outdir = f"{output_dir}/saddle/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)
        start = time.time()
        try:
            # 读 5_1 已算的 eigenvector(全基因组,不是 chr1 npy!)
            eigenvector_path = f"{COMP_COMPUTATION_DIR}/{sample_name}/eigenvector.100k.tsv"
            if not os.path.exists(eigenvector_path):
                print(f"  ❌ {sample_name}: 5_1 eigenvector 缺失: {eigenvector_path}")
                results.append((sample_name, 'error', 0))
                continue

            # mcool resolution 必须跟 eigenvector 一致(都是 100kb)
            generate_single_saddle(
                cool_file=f"{mcool_path}::resolutions/100000",
                eigenvector_file=eigenvector_path,
                output_dir=sample_outdir,
                sample_name=sample_name,
                cache_dir=f"{sample_outdir}/cache",
                n_bins=n_bins,
                contact_type='cis',
                heatmap_size=HEATMAP_SIZE_SINGLE,
                vmin=vmin,
                vmax=vmax,
                nproc=nproc,
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")

    # 多 sample
    multi_outdir = f"{output_dir}/saddle/multi"
    os.makedirs(multi_outdir, exist_ok=True)
    start = time.time()
    try:
        cool_files = [f"{p}::resolutions/100000" for p in samples.values()]
        eigenvector_files = [
            f"{COMP_COMPUTATION_DIR}/{sn}/eigenvector.100k.tsv"
            for sn in samples.keys()
        ]
        # 验证 eigenvector 全部存在
        if not all(os.path.exists(ef) for ef in eigenvector_files):
            missing = [ef for ef in eigenvector_files if not os.path.exists(ef)]
            print(f"  ❌ multi: eigenvector 缺失 {missing}")
            results.append(("multi", 'error', 0))
            return results

        generate_multi_saddle(
            cool_files=cool_files,
            eigenvector_files=eigenvector_files,
            output_dir=multi_outdir,
            sample_names=list(samples.keys()),
            cache_dir=f"{multi_outdir}/cache",
            n_bins=n_bins,
            contact_type='cis',
            heatmap_size=HEATMAP_SIZE_MULTI,  # 多 sample 用 4cm
            vmin=vmin,
            vmax=vmax,
            n_cols=len(samples),
            n_rows=1,
            max_workers=2,
            nproc=nproc,
        )
        elapsed = time.time() - start
        results.append(("multi", 'success', elapsed))
        print(f"  ✅ multi: {elapsed:.1f}s")
    except Exception as e:
        elapsed = time.time() - start
        results.append(("multi", 'error', elapsed))
        print(f"  ❌ multi: {e}")

    return results


def compute_tad_pileup(samples, output_dir, flank, top_n, dpi):
    """
    TAD 边界附近信号堆叠
    流程:
      1. 读 5_1 已算的 boundaries.tsv (直接读,不是 insulation + 自己提取)
      2. plot_multi_tad_boundary_pileup(...) 画 pileup
    """
    from cfizz.viz.pileup import plot_multi_tad_boundary_pileup
    import pandas as pd

    # 5_1 TAD 产物路径(读已算的,不平白无故重算)
    TAD_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/tad"

    results = []

    # 单 sample
    for sample_name, mcool_path in samples.items():
        sample_outdir = f"{output_dir}/tad_pileup/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)
        start = time.time()
        try:
            # 读 5_1 已算的 boundaries.tsv (关键:不是 insulation + 自己提取!)
            # basename 从 mcool 文件名提取,如 "hiPSC_var_chr17"
            mcool_basename = os.path.basename(mcool_path)  # "hiPSC_var_chr17.mcool"
            basename = mcool_basename.split('.')[0]  # "hiPSC_var_chr17"
            # 文件名格式: 2_0.{basename}.{resolution}.{window}.boundaries.tsv
            boundaries_path = f"{TAD_COMPUTATION_DIR}/{sample_name}/2_0.{basename}.10000.100kb.boundaries.tsv"

            if not os.path.exists(boundaries_path):
                print(f"  ❌ {sample_name}: 5_1 boundaries 缺失: {boundaries_path}")
                results.append((sample_name, 'error', 0))
                continue

            boundaries = pd.read_csv(boundaries_path, sep='\t')

            # 画 pileup (关键参数 method='mean', color_scale='linear',跟需求一致)
            output_path = f"{sample_outdir}/{sample_name}_tad_pileup"
            plot_multi_tad_boundary_pileup(
                mcool_paths=[mcool_path],
                boundaries_list=[boundaries],
                output_path=output_path,
                sample_names=[sample_name],
                flank=flank,
                resolution=10000,
                dpi=dpi,
                balance=True,
                method='mean',  # ← 关键! 需求用 mean
                color_scale='linear',  # ← 关键! 需求用 linear
                top_n=top_n,
                cbar_label='mean normalized contacts',  # T-6.16: 跟 color_scale='linear' 一致
                plot_size=4.0,                       # T-6.17: 显式 4cm
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")

    # 多 sample (同样读 5_1 已算的 boundaries.tsv)
    multi_outdir = f"{output_dir}/tad_pileup/multi"
    os.makedirs(multi_outdir, exist_ok=True)
    start = time.time()
    try:
        valid_boundaries = []
        valid_mcool = []
        valid_names = []
        for sample_name, mcool_path in samples.items():
            mcool_basename = os.path.basename(mcool_path)  # "hiPSC_var_chr17.mcool"
            basename = mcool_basename.split('.')[0]  # "hiPSC_var_chr17"
            # 文件名格式: 2_0.{basename}.{resolution}.{window}.boundaries.tsv
            boundaries_path = f"{TAD_COMPUTATION_DIR}/{sample_name}/2_0.{basename}.10000.100kb.boundaries.tsv"
            if not os.path.exists(boundaries_path):
                continue
            boundaries = pd.read_csv(boundaries_path, sep='\t')
            valid_boundaries.append(boundaries)
            valid_mcool.append(mcool_path)
            valid_names.append(sample_name)

        if len(valid_boundaries) >= 2:
            output_path = f"{multi_outdir}/tad_pileup_multi"
            plot_multi_tad_boundary_pileup(
                mcool_paths=valid_mcool,
                boundaries_list=valid_boundaries,
                output_path=output_path,
                sample_names=valid_names,
                flank=flank,
                resolution=10000,
                dpi=dpi,
                balance=True,
                method='mean',  # ← 关键! 需求用 mean
                color_scale='linear',  # ← 关键! 需求用 linear
                top_n=top_n,
                cbar_label='mean normalized contacts',  # T-6.16: 跟 color_scale='linear' 一致
                plot_size=4.0,                       # T-6.17: 显式 4cm
            )
            elapsed = time.time() - start
            results.append(("multi", 'success', elapsed))
            print(f"  ✅ multi: {elapsed:.1f}s")
        else:
            results.append(("multi", 'error', 0))
            print(f"  ❌ multi: 少于 2 sample 有 boundaries")
    except Exception as e:
        elapsed = time.time() - start
        results.append(("multi", 'error', elapsed))
        print(f"  ❌ multi: {e}")

    return results


def compute_apa(samples, output_dir, chrom, start_pos, end_pos,
                resolution, window, corner_size, min_distance, vmin, vmax):
    """
    APA (Aggregate Peak Analysis)
    流程:
      1. 读 5_1 已算的 loops.txt (不重算!)
      2. plot_multi_apa_heatmap(...) 画 APA
    """
    from cfizz.api.integrated.apa_pileup import plot_multi_apa_heatmap

    # 5_1 Loop 产物路径(读已算的)
    LOOP_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/loop"

    results = []

    # 单 sample
    for sample_name, mcool_path in samples.items():
        sample_outdir = f"{output_dir}/apa/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)
        start = time.time()
        try:
            # 读 5_1 已算的 loops
            basename = os.path.basename(mcool_path).split('.')[0]
            loops_path = f"{LOOP_COMPUTATION_DIR}/{sample_name}/{basename}.{resolution//1000}k.loops.txt"

            if not os.path.exists(loops_path) or os.path.getsize(loops_path) == 0:
                print(f"  ❌ {sample_name}: 5_1 loops 缺失或空: {loops_path}")
                results.append((sample_name, 'error', 0))
                continue

            # 画 APA (vmin/vmax 传 None，让函数自动计算)
            output_path = f"{sample_outdir}/{sample_name}_apa"
            plot_multi_apa_heatmap(
                mcool_paths=[mcool_path],
                loops_paths=[loops_path],
                output_path=output_path,
                sample_names=[sample_name],
                resolution=resolution,
                window=window,
                corner_size=corner_size,
                min_distance=min_distance,
                vmin=None,  # 自动计算
                vmax=None,  # 自动计算
                balance=True,
                n_processes=None,
                plot_size=4.0,  # 统一 4cm
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")

    # 多 sample (同样读 5_1 已算的)
    multi_outdir = f"{output_dir}/apa/multi"
    os.makedirs(multi_outdir, exist_ok=True)
    start = time.time()
    try:
        valid_loops = []
        valid_mcool = []
        valid_names = []
        for sample_name, mcool_path in samples.items():
            basename = os.path.basename(mcool_path).split('.')[0]
            loops_path = f"{LOOP_COMPUTATION_DIR}/{sample_name}/{basename}.{resolution//1000}k.loops.txt"
            if os.path.exists(loops_path) and os.path.getsize(loops_path) > 0:
                valid_loops.append(loops_path)
                valid_mcool.append(mcool_path)
                valid_names.append(sample_name)

        if len(valid_loops) >= 2:
            output_path = f"{multi_outdir}/apa_multi"
            plot_multi_apa_heatmap(
                mcool_paths=valid_mcool,
                loops_paths=valid_loops,
                output_path=output_path,
                sample_names=valid_names,
                resolution=resolution,
                window=window, corner_size=corner_size,
                min_distance=min_distance, vmin=None, vmax=None,  # 自动计算
                balance=True, n_processes=None, plot_size=4.0,  # 统一 4cm
            )
            elapsed = time.time() - start
            results.append(("multi", 'success', elapsed))
            print(f"  ✅ multi: {elapsed:.1f}s")
        else:
            results.append(("multi", 'skipped', 0))
            print(f"  ⏭️ multi: 少于 2 sample 有 loops")
    except Exception as e:
        elapsed = time.time() - start
        results.append(("multi", 'error', elapsed))
        print(f"  ❌ multi: {e}")

    return results


def main():
    print("=" * 70)
    print("Example 5_3: Primary Analysis Pileup - 3 类累计分析")
    print("=" * 70)
    print(f"  样本: {list(SAMPLES.keys())}")
    print(f"  输出: {PILEUP_DIR}/")
    print()
    print("3 类累计分析, 9 个产物:")
    print("  1. Saddle plot     (chr17, 100kb, n_bins=98)    3 个: hiPSC_var, hiPSC_nor, multi")
    print("  2. TAD pileup      (chr17:43.5-44.5M, 10kb, 100kb) 3 个: hiPSC_var, hiPSC_nor, multi")
    print("  3. APA             (chr17:43.5-43.7M, 10kb)        3 个: hiPSC_var, hiPSC_nor, multi")
    print()
    print("关键决策:")
    print("  - Saddle: 读 5_1 全基因组 eigenvector.tsv")
    print("  - TAD: 读 5_1 已算的 insulation.tsv (不重算!)")
    print("  - APA: 读 5_1 已算的 loops.txt (不重算!)")
    print()

    os.makedirs(PILEUP_DIR, exist_ok=True)
    overall_start = time.time()

    # Step 1: Saddle
    print("=" * 70)
    print("Step 1: Saddle plot")
    print("=" * 70)
    saddle_results = compute_saddle(
        SAMPLES, PILEUP_DIR,
        SADDLE_N_BINS, SADDLE_VMIN, SADDLE_VMAX, SADDLE_NPROC,
    )

    # Step 2: TAD pileup
    print("\n" + "=" * 70)
    print("Step 2: TAD pileup")
    print("=" * 70)
    tad_results = compute_tad_pileup(
        SAMPLES, PILEUP_DIR,
        TAD_PILEUP_FLANK, TAD_PILEUP_TOP_N, TAD_PILEUP_DPI,
    )

    # Step 3: APA
    print("\n" + "=" * 70)
    print("Step 3: APA")
    print("=" * 70)
    apa_results = compute_apa(
        SAMPLES, PILEUP_DIR, APA_CHROM,
        APA_START, APA_END,
        APA_RESOLUTION, APA_WINDOW, APA_CORNER_SIZE,
        APA_MIN_DISTANCE, APA_VMIN, APA_VMAX,
    )

    # 总结
    total_elapsed = time.time() - overall_start
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    def count_success(results):
        return sum(1 for r in results if r[1] == 'success')

    print(f"  Saddle: {count_success(saddle_results)}/{len(saddle_results)} 成功")
    print(f"  TAD:    {count_success(tad_results)}/{len(tad_results)} 成功")
    print(f"  APA:    {count_success(apa_results)}/{len(apa_results)} 成功")
    print(f"\n  总耗时: {total_elapsed:.1f}s")
    print(f"\n  产物结构:")
    print(f"    {PILEUP_DIR}/")
    print(f"      saddle/{{hiPSC_var, hiPSC_nor, multi}}/")
    print(f"      tad_pileup/{{hiPSC_var, hiPSC_nor, multi}}/")
    print(f"      apa/{{hiPSC_var, hiPSC_nor, multi}}/")

    print("\n" + "=" * 70)
    print("✅ Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
