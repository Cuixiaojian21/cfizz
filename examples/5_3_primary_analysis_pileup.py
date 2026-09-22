#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example 5_3: Primary Analysis Pileup - 4 类累计分析

只做累计分析，4 类:
1. Saddle plot (A/B compartment 交互)
2. TAD linear pileup (TAD 边界附近信号堆叠, balance+mean+linear)
3. TAD O/E pileup (TAD 边界附近 obs/exp 偏离, OE + mean + log2)
4. APA O/E pileup (loop 中心 obs/exp 偏离, OE + median + linear, 用户推荐配置)

一句话解决问题(在 cfizz/ 根目录):
  python examples/5_3_primary_analysis_pileup.py

设计哲学:
  - 只用 viz / api 层，不动任何 analyze 层
  - 不导入 cfizz.io.paths (已放弃)
  - 不导入 hicviz 任何模块
  - Saddle: 读 5_1 已算的全基因组 eigenvector.tsv
  - TAD O/E: 读 5_1 已算的 boundaries.tsv + per-sample cooltools.expected_cis (OE + mean + log2)
  - APA O/E: 读 5_1 已算的 loops.txt + per-sample cooltools.expected_cis + cooltools.pileup (OE + median + linear)
  - 区域 chr17 demo: SADDLE 全染色体, TAD O/E 全基因组 top1000, APA O/E 全基因组 loops
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
FASTA_PATH = "/mnt/e/1_1_LiMin/1_0_supports/hg38/hg38.fa"  # ⚠️ 用户需自己准备 (3GB)

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
TAD_PILEUP_NPROC = 8  # cooltools.expected_cis 并行核数(O/E pileup 必需)


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

# === APA O/E pileup 参数 (用户指定: OE + median + linear) ===
APA_OE_FLANK = 70_000          # = 7 bins @ 10kb resolution, 对应 APA_WINDOW
APA_OE_NPROC = 8               # cooltools.expected_cis 并行核数
APA_OE_DPI = 300
APA_OE_VMIN = 0.5              # obs/exp 中心 1 (linear 直接显示)
APA_OE_VMAX = 2.0
APA_OE_CMAP = 'coolwarm'
APA_OE_COLOR_SCALE = 'linear'  # ← 用户指定 linear (无 log 变换)
APA_OE_CBAR_LABEL = 'median obs/exp'  # ← linear 时无括号


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


def compute_tad_oe_pileup(samples, output_dir, flank, top_n, dpi, nproc=8):
    """
    TAD boundary pileup (O/E log2 风格)

    与 compute_tad_pileup 的差别:
      - 输入矩阵: balance 矩阵逐元素除以 expected(distance) (O/E)
      - 配色: coolwarm (发散型), vmin/vmax = ±1 (log2 obs/exp 合理范围)
      - 颜色变换: log2(obs/exp)
      - colorbar 标签: 'log2(obs/exp)'

    流程:
      1. 读 5_1 已算的 boundaries.tsv (不重算!)
      2. 算 per-sample P(s) DataFrame (cooltools.expected_cis)
      3. plot_multi_tad_boundary_pileup 带 expected_dfs 画 O/E pileup
    """
    from cfizz.analyze.oe import compute_expected_cis_per_sample
    from cfizz.viz.pileup import plot_multi_tad_boundary_pileup
    import pandas as pd

    # 5_1 TAD 产物路径(读已算的, 不平白无故重算)
    TAD_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/tad"
    WINDOW_KB = f"{TAD_PILEUP_WINDOW // 1000}kb"  # 100kb (跟现有 TAD pileup 一致)

    results = []

    # 1) 收集各 sample 的 boundaries 与 mcool
    valid_mcools, valid_bd, valid_names = [], [], []
    for sample_name, mcool_path in samples.items():
        mcool_basename = os.path.basename(mcool_path)
        basename = mcool_basename.split('.')[0]
        boundaries_path = f"{TAD_COMPUTATION_DIR}/{sample_name}/2_0.{basename}.10000.{WINDOW_KB}.boundaries.tsv"
        if not os.path.exists(boundaries_path):
            print(f"  ❌ {sample_name}: 5_1 boundaries 缺失: {boundaries_path}")
            results.append((sample_name, 'error', 0))
            continue
        valid_mcools.append(mcool_path)
        valid_bd.append(pd.read_csv(boundaries_path, sep='\t'))
        valid_names.append(sample_name)

    if not valid_mcools:
        return results

    # 2) per-sample 算 expected_cis (一次算完, 单 sample 与 multi 复用)
    print(f"  → per-sample 算 P(s) DataFrame (nproc={nproc})...")
    p_start = time.time()
    samples_dict = dict(zip(valid_names, valid_mcools))
    expected_dfs = compute_expected_cis_per_sample(
        samples_dict, resolution=10000, nproc=nproc, balance=True,
    )
    print(f"  ✓ P(s) 算完 {len(expected_dfs)} 个 sample ({time.time() - p_start:.1f}s)")

    # 3) 单 sample
    for i, sample_name in enumerate(valid_names):
        sample_outdir = f"{output_dir}/tad_oe_pileup/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)
        start = time.time()
        try:
            output_path = f"{sample_outdir}/{sample_name}_tad_oe_pileup"
            plot_multi_tad_boundary_pileup(
                mcool_paths=[valid_mcools[i]],
                boundaries_list=[valid_bd[i]],
                expected_dfs=[expected_dfs[i]],
                output_path=output_path,
                sample_names=[sample_name],
                flank=flank,
                resolution=10000,
                vmin=-1,
                vmax=1,
                cmap='coolwarm',
                color_scale='log2',
                color_scale_for_cbar='log2',
                cbar_label='log2(mean obs/exp)',
                balance=True,
                method='mean',
                top_n=top_n,
                dpi=dpi,
                n_processes=nproc,
                plot_size=4.0,
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {type(e).__name__}: {e}")

    # 4) 多 sample (共用前面算好的 P(s))
    if len(valid_mcools) >= 2:
        multi_outdir = f"{output_dir}/tad_oe_pileup/multi"
        os.makedirs(multi_outdir, exist_ok=True)
        start = time.time()
        try:
            output_path = f"{multi_outdir}/tad_oe_pileup_multi"
            plot_multi_tad_boundary_pileup(
                mcool_paths=valid_mcools,
                boundaries_list=valid_bd,
                expected_dfs=expected_dfs,
                output_path=output_path,
                sample_names=valid_names,
                flank=flank,
                resolution=10000,
                vmin=-1,
                vmax=1,
                cmap='coolwarm',
                color_scale='log2',
                color_scale_for_cbar='log2',
                cbar_label='log2(mean obs/exp)',
                balance=True,
                method='mean',
                top_n=top_n,
                dpi=dpi,
                n_processes=nproc,
                plot_size=4.0,
            )
            elapsed = time.time() - start
            results.append(('multi', 'success', elapsed))
            print(f"  ✅ multi: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append(('multi', 'error', elapsed))
            print(f"  ❌ multi: {type(e).__name__}: {e}")
    else:
        results.append(('multi', 'skipped', 0))
        print(f"  ⏭️ multi: 少于 2 sample 有 boundaries")

    return results


def compute_apa_oe_pileup(samples, output_dir, dpi, nproc=8):
    """
    APA O/E pileup (用户推荐配置: OE + median + linear)

    配置:
      - 输入矩阵: balance 矩阵逐元素除以 expected(distance) (O/E)
      - method: median (抗 outlier, 推荐用于跨 sample 对比)
      - color_scale: linear (直接显示 obs/exp 倍数, 不 log 变换)
      - 配色: coolwarm (发散型)
      - vmin/vmax: 0.5 / 2.0 (obs/exp 中心 1 的合理区间)
      - colorbar 标签: 'median obs/exp' (无 log, 无括号)

    流程:
      1. 读 5_1 已算的 loops.txt (不重算!)
      2. 算 per-sample P(s) DataFrame (cooltools.expected_cis)
      3. compute_apa_oe_pileup_single/multi 画 O/E APA
    """
    from cfizz.analyze.oe import compute_expected_cis_per_sample
    from cfizz.api.integrated.apa_pileup import (
        compute_apa_oe_pileup_single,
        compute_apa_oe_pileup_multi,
    )

    # 5_1 Loop 产物路径(读已算的)
    LOOP_COMPUTATION_DIR = f"{OUTPUT_DIR}/1_computation/loop"

    results = []

    # 1) 收集各 sample 的 loops 路径
    valid_mcools, valid_loops_paths, valid_names = [], [], []
    for sample_name, mcool_path in samples.items():
        basename = os.path.basename(mcool_path).split('.')[0]
        loops_path = f"{LOOP_COMPUTATION_DIR}/{sample_name}/{basename}.{APA_RESOLUTION//1000}k.loops.txt"
        if not os.path.exists(loops_path) or os.path.getsize(loops_path) == 0:
            print(f"  ❌ {sample_name}: 5_1 loops 缺失或空: {loops_path}")
            results.append((sample_name, 'error', 0))
            continue
        valid_mcools.append(mcool_path)
        valid_loops_paths.append(loops_path)
        valid_names.append(sample_name)

    if not valid_mcools:
        return results

    # 2) per-sample 算 P(s) (一次算完, 单 sample 与 multi 复用)
    print(f"  → per-sample 算 P(s) DataFrame (nproc={nproc})...")
    p_start = time.time()
    samples_dict = dict(zip(valid_names, valid_mcools))
    expected_dfs = compute_expected_cis_per_sample(
        samples_dict, resolution=APA_RESOLUTION, nproc=nproc, balance=True,
    )
    print(f"  ✓ P(s) 算完 {len(expected_dfs)} 个 sample ({time.time() - p_start:.1f}s)")

    # 3) 单 sample
    for i, sample_name in enumerate(valid_names):
        sample_outdir = f"{output_dir}/apa_oe/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)
        start = time.time()
        try:
            output_path = f"{sample_outdir}/{sample_name}_apa_oe"
            compute_apa_oe_pileup_single(
                mcool_path=valid_mcools[i],
                loops_path=valid_loops_paths[i],
                expected_df=expected_dfs[i],
                output_path=output_path,
                sample_name=sample_name,
                flank=APA_OE_FLANK,
                resolution=APA_RESOLUTION,
                balance=True,
                method='median',  # ← 用户指定 median (替代 mean)
                vmin=APA_OE_VMIN, vmax=APA_OE_VMAX,
                cmap=APA_OE_CMAP,
                color_scale=APA_OE_COLOR_SCALE,
                cbar_label=APA_OE_CBAR_LABEL,
                dpi=dpi,
                plot_size=4.0,
            )
            elapsed = time.time() - start
            results.append((sample_name, 'success', elapsed))
            print(f"  ✅ {sample_name}: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {type(e).__name__}: {e}")

    # 4) 多 sample (共用前面算好的 P(s))
    if len(valid_mcools) >= 2:
        multi_outdir = f"{output_dir}/apa_oe/multi"
        os.makedirs(multi_outdir, exist_ok=True)
        start = time.time()
        try:
            output_path = f"{multi_outdir}/apa_oe_multi"
            compute_apa_oe_pileup_multi(
                mcool_paths=valid_mcools,
                loops_paths=valid_loops_paths,
                expected_dfs=expected_dfs,
                output_path=output_path,
                sample_names=valid_names,
                flank=APA_OE_FLANK,
                resolution=APA_RESOLUTION,
                balance=True,
                method='median',  # ← 用户指定 median (替代 mean)
                vmin=APA_OE_VMIN, vmax=APA_OE_VMAX,
                cmap=APA_OE_CMAP,
                color_scale=APA_OE_COLOR_SCALE,
                cbar_label=APA_OE_CBAR_LABEL,
                dpi=dpi,
                plot_size=4.0,
            )
            elapsed = time.time() - start
            results.append(('multi', 'success', elapsed))
            print(f"  ✅ multi: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            results.append(('multi', 'error', elapsed))
            print(f"  ❌ multi: {type(e).__name__}: {e}")
    else:
        results.append(('multi', 'skipped', 0))
        print(f"  ⏭️ multi: 少于 2 sample 有 loops")

    return results


def main():
    print("=" * 70)
    print("Example 5_3: Primary Analysis Pileup - 3 类累计分析")
    print("=" * 70)
    print(f"  样本: {list(SAMPLES.keys())}")
    print(f"  输出: {PILEUP_DIR}/")
    print()
    print("4 类累计分析, 12 个产物:")
    print("  1. Saddle plot     (chr17, 100kb, n_bins=98)       3 个: hiPSC_var, hiPSC_nor, multi")
    print("  2. TAD linear      (chr17:43.5-44.5M, 10kb, 100kb) 3 个: hiPSC_var, hiPSC_nor, multi")
    print("  3. TAD O/E pileup  (chr17 全基因组 top1000, log2(obs/exp)) 3 个: hiPSC_var, hiPSC_nor, multi")
    print("  4. APA O/E pileup  (OE + median + linear, 推荐配置) 3 个: hiPSC_var, hiPSC_nor, multi")
    print()
    print("关键决策:")
    print("  - Saddle: 读 5_1 全基因组 eigenvector.tsv")
    print("  - TAD O/E: 读 5_1 已算的 boundaries.tsv + per-sample cooltools.expected_cis (OE + mean + log2)")
    print("  - APA O/E: 读 5_1 已算的 loops.txt + per-sample cooltools.expected_cis + cooltools.pileup (OE + median + linear)")
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

    # Step 2: TAD pileup (linear balance 风格)
    print("\n" + "=" * 70)
    print("Step 2: TAD pileup (linear balance)")
    print("=" * 70)
    tad_results = compute_tad_pileup(
        SAMPLES, PILEUP_DIR,
        TAD_PILEUP_FLANK, TAD_PILEUP_TOP_N, TAD_PILEUP_DPI,
    )

    # Step 2.5: TAD O/E pileup (log2(obs/exp) + coolwarm)
    print("\n" + "=" * 70)
    print("Step 2.5: TAD O/E pileup (log2(obs/exp) + coolwarm)")
    print("=" * 70)
    tad_oe_results = compute_tad_oe_pileup(
        SAMPLES, PILEUP_DIR,
        TAD_PILEUP_FLANK, TAD_PILEUP_TOP_N, TAD_PILEUP_DPI,
        nproc=TAD_PILEUP_NPROC,
    )

    # Step 3: APA O/E pileup (OE + median + linear, 用户推荐配置)
    print("\n" + "=" * 70)
    print("Step 3: APA O/E pileup (OE + median + linear)")
    print("=" * 70)
    apa_results = compute_apa_oe_pileup(
        SAMPLES, PILEUP_DIR,
        dpi=APA_OE_DPI,
        nproc=APA_OE_NPROC,
    )

    # 总结
    total_elapsed = time.time() - overall_start
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    def count_success(results):
        return sum(1 for r in results if r[1] == 'success')

    print(f"  Saddle:      {count_success(saddle_results)}/{len(saddle_results)} 成功")
    print(f"  TAD linear:  {count_success(tad_results)}/{len(tad_results)} 成功")
    print(f"  TAD O/E:     {count_success(tad_oe_results)}/{len(tad_oe_results)} 成功")
    print(f"  APA O/E:     {count_success(apa_results)}/{len(apa_results)} 成功 (OE + median + linear)")
    print(f"\n  总耗时: {total_elapsed:.1f}s")
    print(f"\n  产物结构:")
    print(f"    {PILEUP_DIR}/")
    print(f"      saddle/{{hiPSC_var, hiPSC_nor, multi}}/")
    print(f"      tad_pileup/{{hiPSC_var, hiPSC_nor, multi}}/")
    print(f"      tad_oe_pileup/{{hiPSC_var, hiPSC_nor, multi}}/")
    print(f"      apa_oe/{{hiPSC_var, hiPSC_nor, multi}}/")

    print("\n" + "=" * 70)
    print("✅ Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
