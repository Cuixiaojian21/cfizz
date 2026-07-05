#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example 5_1: Primary Analysis Template - 计算 3 个特征

一句话解决问题(在 cfizz/ 根目录):
  python examples/5_1_primary_analysis_template.py

功能:
  1. Compartment (eigenvector + GC cov) - 100kb
  2. TAD (insulation + boundaries + tads) - 10kb + 3 windows
  3. Loop (loops.txt BEDPE 6列) - 10kb

输入(相对路径):
  demo/data/hiPSC_nor_chr17.mcool   # chr17 demo mcool (10kb + 100kb)
  demo/data/hiPSC_var_chr17.mcool   # 同上
  FASTA_PATH 用户需自己准备 (cfizz 默认值: /mnt/g/2_0_demo/1_0_support/hg38.fa)
  ⚠️ 克隆 cfizz 后必须改成自己的 hg38.fa 路径 (3GB)

产物结构:
  demo/output/5_1_primary_analysis_template/1_computation/
    compartment/{hiPSC_nor, hiPSC_var}/  # eigenvector.100k.tsv + gc_cov.100k.tsv
    tad/{hiPSC_nor, hiPSC_var}/          # insulation + boundaries + tads
    loop/{hiPSC_nor, hiPSC_var}/         # loops.txt

设计哲学:
  - 只用 analyze 层，不走 viz 层
  - 不导入 cfizz.io.paths (已放弃)
  - 不导入 hicviz 任何模块
  - 单进程顺序跑，简单够用
  - demo mcool 只含 chr17,跑 eigenvector 时其他 chr 空 bin 自动跳过
"""

import sys
import os
# 脚本位于 cfizz/examples/<name>.py → 上 3 层到 /mnt/g/2_0_demo/ (cfizz 的父目录)
# 需要在 PYTHONPATH 里能找到 cfizz 包, 所以加的是 cfizz 的父目录
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import time


# === 数据(相对路径) ===
DEMO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "demo")
DATA_DIR = os.path.join(DEMO_DIR, "data")
OUTPUT_ROOT = os.path.join(DEMO_DIR, "output")
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "5_1_primary_analysis_template")

# FASTA: 直接用绝对路径 (用户克隆 cfizz 后需改成自己的 hg38.fa 路径)
FASTA_PATH = "/mnt/g/2_0_demo/1_0_support/hg38.fa"  # ⚠️ 用户需自己准备 (3GB)

# === 双 sample (hiPSC demo) ===
SAMPLES = {
    "hiPSC_var": os.path.join(DATA_DIR, "hiPSC_var_chr17.mcool"),  # 肿瘤
    "hiPSC_nor": os.path.join(DATA_DIR, "hiPSC_nor_chr17.mcool"),  # 正常
}

# === 计算参数 ===
# Compartment
COMPARTMENT_RESOLUTIONS = [100_000]  # 只用 100kb

# TAD
TAD_RESOLUTION = 10_000
TAD_WINDOWS = [50_000, 100_000, 500_000]

# Loop
LOOP_RESOLUTION = 10_000
LOOP_PEAK_WIDTHS = [2]
LOOP_WINDOW_WIDTHS = [5]
NPROC = 26


# === 核心函数 ===

def compute_compartment_features(samples, output_dir, resolutions, fasta_path):
    """
    算 compartment 特征(eigenvector + GC cov)
    对每个 sample × 每个 resolution 算一次
    """
    from cfizz.analyze.compartment import process_compartment

    results = []
    for sample_name, mcool_path in samples.items():
        if not os.path.exists(mcool_path):
            print(f"  ❌ mcool 缺失: {mcool_path}")
            results.append((sample_name, None, None, 'error', 0))
            continue

        sample_outdir = f"{output_dir}/1_computation/compartment/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)

        for resolution in resolutions:
            start = time.time()
            try:
                eig_df = process_compartment(
                    mcool_path=mcool_path,
                    resolution=resolution,
                    fasta_path=fasta_path,
                    output_dir=sample_outdir,
                    n_eigs=3,
                )
                elapsed = time.time() - start
                # 验证产物
                res_str = f"{resolution // 1000}k"  # 100_000 → "100k"
                eig_path = os.path.join(sample_outdir, f"eigenvector.{res_str}.tsv")
                gc_path = os.path.join(sample_outdir, f"gc_cov.{res_str}.tsv")
                if os.path.exists(eig_path) and os.path.getsize(eig_path) > 0 \
                   and os.path.exists(gc_path) and os.path.getsize(gc_path) > 0:
                    results.append((sample_name, resolution, res_str, 'success', elapsed))
                    print(f"  ✅ {sample_name} @ {res_str}: {elapsed:.1f}s "
                          f"({len(eig_df)} bins)")
                else:
                    missing = []
                    if not os.path.exists(eig_path): missing.append(eig_path)
                    if not os.path.exists(gc_path): missing.append(gc_path)
                    results.append((sample_name, resolution, res_str, 'partial', elapsed))
                    print(f"  ⚠️ {sample_name} @ {res_str}: 缺失 {missing}")
            except Exception as e:
                elapsed = time.time() - start
                results.append((sample_name, resolution, None, 'error', elapsed))
                print(f"  ❌ {sample_name} @ {resolution}: {e}")
    return results


def compute_tad_features(samples, output_dir, resolution, windows):
    """
    算 TAD 特征(insulation + boundaries + tads)
    对每个 sample 算一次,产出 3 windows 的文件
    """
    from cfizz.analyze.tad import process_tads

    results = []
    for sample_name, mcool_path in samples.items():
        if not os.path.exists(mcool_path):
            print(f"  ❌ mcool 缺失: {mcool_path}")
            results.append((sample_name, 'error', 0))
            continue

        sample_outdir = f"{output_dir}/1_computation/tad/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)

        start = time.time()
        try:
            boundaries_dict, tads_dict = process_tads(
                cooler_path=mcool_path,
                resolution=resolution,
                windows=windows,
                output_dir=sample_outdir,
                max_tad_length=3_000_000,
                threshold_method='otsu',
                nproc=NPROC,
                verbose=False,
            )
            elapsed = time.time() - start

            # 验证产物(7 个文件:1 insulation + 3 boundaries + 3 tads)
            # process_tads 内部从 mcool 文件名推导 basename, demo mcool basename=hiPSC_var_chr17
            mcool_basename = os.path.basename(mcool_path).replace('.mcool', '')
            expected_files = [
                f"1_0.{mcool_basename}.{resolution}.insulation.tsv",
            ]
            for w in windows:
                w_str = f"{w // 1000}kb"
                expected_files.append(f"2_0.{mcool_basename}.{resolution}.{w_str}.boundaries.tsv")
                expected_files.append(f"2_1.{mcool_basename}.{resolution}.{w_str}.tads.tsv")

            missing = [f for f in expected_files
                       if not os.path.exists(os.path.join(sample_outdir, f))]

            if not missing:
                n_tads_total = sum(len(tads_dict[w]) for w in windows)
                results.append((sample_name, 'success', elapsed))
                print(f"  ✅ {sample_name}: {elapsed:.1f}s, "
                      f"共 {n_tads_total:,} TADs (3 windows)")
            else:
                results.append((sample_name, 'partial', elapsed))
                print(f"  ⚠️ {sample_name}: 缺失 {missing}")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed))
            print(f"  ❌ {sample_name}: {e}")
    return results


def compute_loop_features(samples, output_dir, resolution, peak_widths, window_widths):
    """
    算 Loop 特征(loops.txt BEDPE 6 列)
    对每个 sample 算一次
    """
    from cfizz.analyze.loop import process_loops

    results = []
    for sample_name, mcool_path in samples.items():
        if not os.path.exists(mcool_path):
            print(f"  ❌ mcool 缺失: {mcool_path}")
            results.append((sample_name, 'error', 0, None))
            continue

        sample_outdir = f"{output_dir}/1_computation/loop/{sample_name}"
        os.makedirs(sample_outdir, exist_ok=True)

        start = time.time()
        try:
            output_file = process_loops(
                mcool_path=mcool_path,
                resolution=resolution,
                output_dir=sample_outdir,
                peak_widths=peak_widths,
                window_widths=window_widths,
                only_anchors=True,
                nproc=NPROC,
                verbose=True,
            )
            elapsed = time.time() - start

            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                with open(output_file) as f:
                    n_loops = sum(1 for _ in f if _.strip())
                results.append((sample_name, 'success', elapsed, n_loops))
                print(f"  ✅ {sample_name}: {elapsed:.1f}s, {n_loops:,} loops")
            else:
                results.append((sample_name, 'error', elapsed, 0))
                print(f"  ❌ {sample_name}: loops.txt 未生成")
        except Exception as e:
            elapsed = time.time() - start
            results.append((sample_name, 'error', elapsed, 0))
            print(f"  ❌ {sample_name}: {e}")
    return results


# === 主函数 ===

def main():
    print("=" * 70)
    print("Example 5_1: Primary Analysis Template - 计算 3 个特征")
    print("=" * 70)
    print(f"  样本: {list(SAMPLES.keys())}")
    print(f"  输出: {OUTPUT_DIR}")
    print()
    print("流程:")
    print("  1. Compartment(eigenvector + GC cov)")
    print("  2. TAD(insulation + boundaries + tads)")
    print("  3. Loop(loops.txt)")
    print()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    overall_start = time.time()

    # Step 1: Compartment
    print("=" * 70)
    print("Step 1: Compartment 计算")
    print("=" * 70)
    comp_results = compute_compartment_features(
        SAMPLES, OUTPUT_DIR, COMPARTMENT_RESOLUTIONS, FASTA_PATH
    )

    # Step 2: TAD
    print("\n" + "=" * 70)
    print("Step 2: TAD 计算")
    print("=" * 70)
    tad_results = compute_tad_features(
        SAMPLES, OUTPUT_DIR, TAD_RESOLUTION, TAD_WINDOWS
    )

    # Step 3: Loop
    print("\n" + "=" * 70)
    print("Step 3: Loop 计算")
    print("=" * 70)
    loop_results = compute_loop_features(
        SAMPLES, OUTPUT_DIR, LOOP_RESOLUTION, LOOP_PEAK_WIDTHS, LOOP_WINDOW_WIDTHS
    )

    # 总结
    total_elapsed = time.time() - overall_start
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    comp_success = sum(1 for r in comp_results if r[3] == 'success')
    tad_success = sum(1 for r in tad_results if r[1] == 'success')
    loop_success = sum(1 for r in loop_results if r[1] == 'success')

    print(f"  Compartment: {comp_success}/{len(comp_results)} 成功")
    print(f"  TAD:         {tad_success}/{len(tad_results)} 成功")
    print(f"  Loop:        {loop_success}/{len(loop_results)} 成功")
    print(f"\n  总耗时: {total_elapsed:.1f}s")
    print(f"\n  产物结构:")
    print(f"    {OUTPUT_DIR}/1_computation/")
    print(f"      compartment/{{hiPSC_nor, hiPSC_var}}/  # eigenvector.{{res}}.tsv + gc_cov.{{res}}.tsv")
    print(f"      tad/{{hiPSC_nor, hiPSC_var}}/          # insulation + boundaries + tads")
    print(f"      loop/{{hiPSC_nor, hiPSC_var}}/         # loops.txt")
    print(f"\n💡 后续可视化需要 oe.npy 矩阵,目前未生成。")
    print(f"   如需 oe.npy,运行 examples/5_2_primary_analysis_visualization.py (chrom=chr17, region=43.5M-43.7M)")

    print("\n" + "=" * 70)
    print("✅ Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
