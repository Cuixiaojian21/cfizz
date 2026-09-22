#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example 5_0: Distance Decay (P(s)) 对比分析

基于 2 sample mcool, 算全基因组距离衰减曲线并对比.

参数 (demo 数据):
  - resolution = 10_000 bp (10kb, 与 5_3 一致)
  - max_distance = 50_000_000 (50Mb, demo 数据 chr17 短)
  - nproc = 4
  - balance = False (raw, clr_weight_name=None)
  - fit_window = (0.2Mb, 5Mb) (demo 数据较短, 适配窗口)
  - excluded_chroms = ('chrM', 'chrY')

产物:
  demo/output/5_0_distance_decay_comparison/7_distance_decay/
    distance_decay_loglog.png/svg     # log-log 对比曲线
    distance_decay_alpha.tsv         # 每 sample α + R² + 分段斜率
    distance_decay_comparison.tsv    # 两两 sample α +段斜率 diff
    distance_decay_data.tsv          # 完整 P(s) 数据 (~5000 行/sample)

跑法 (在 cfizz/ 根目录):
  PYTHONPATH=src python examples/5_0_distance_decay_comparison.py

设计哲学:
  - 用 cfizz.viz.run_distance_decay 一站式 (4 个产物)
  - 不引入额外分析, 不重复造 cooltools.expected_cis 包装
  - 与 g_6d_distance_decay.py reference 行为一致
  - 编号5_0 表示它在5_x 序列里靠前 (距离衰减 → 特征分析 的逻辑顺序)
"""

import sys
import os
# 脚本位于 cfizz_v208620_updata/examples/<name>.py → 上 2 层到 cfizz_v208620_updata/, 再加 src/
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

import time

import matplotlib
matplotlib.use('Agg')


# === Demo 数据 (与 5_3 一致) ===
DEMO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "demo")
DATA_DIR = os.path.join(DEMO_DIR, "data")
OUTPUT_ROOT = os.path.join(DEMO_DIR, "output")
# 5_0 独立输出 (不与 5_1/5_2/5_3 共享目录)
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "5_0_distance_decay_comparison")

# === 双 sample (hiPSC demo) ===
SAMPLES = {
    "hiPSC_var": os.path.join(DATA_DIR, "hiPSC_var_chr17.mcool"),  # 肿瘤
    "hiPSC_nor": os.path.join(DATA_DIR, "hiPSC_nor_chr17.mcool"),  # 正常
}

# === Distance decay 参数 (适合 demo 数据) ===
RESOLUTION = 10_000           # 10kb
MAX_DISTANCE = 50_000_000     # 50Mb (chr17 ~83Mb, 留 buffer)
NPROC = 4
FIT_WINDOW = (0.2e6, 5e6)    # demo 数据较短, 缩窗口

# === 产物路径 (与 5_3 风格一致, 7_xxx 子目录) ===
DECAY_DIR = f"{OUTPUT_DIR}/7_distance_decay"


def main():
    print("=" * 70)
    print("Example 5_0: Distance Decay (P(s)) 对比分析")
    print("=" * 70)
    print(f"  样本: {list(SAMPLES.keys())}")
    print(f"  输出: {DECAY_DIR}/")
    print()
    print("参数:")
    print(f"  resolution = {RESOLUTION:,} bp")
    print(f"  max_distance = {MAX_DISTANCE/1e6:.0f} Mb")
    print(f"  fit_window = {FIT_WINDOW[0]/1e6}-{FIT_WINDOW[1]/1e6} Mb")
    print(f"  nproc = {NPROC}")
    print(f"  excluded = ('chrM', 'chrY')")
    print()

    os.makedirs(DECAY_DIR, exist_ok=True)
    overall_start = time.time()

    # 一站式调用 (cfizz.viz.run_distance_decay 是 cfizz.analyze.distance + cfizz.viz.distance_decay_plot 的端到端 wrapper)
    from cfizz.viz import run_distance_decay

    run_distance_decay(
        samples=SAMPLES,
        output_dir=DECAY_DIR,
        resolution=RESOLUTION,
        max_distance=MAX_DISTANCE,
        nproc=NPROC,
        fit_window=FIT_WINDOW,
        balance=False,
        excluded_chroms=('chrM', 'chrY'),
    )

    total = time.time() - overall_start
    print()
    print("=" * 70)
    print(f"✅ Done! 总耗时: {total:.1f}s")
    print("=" * 70)
    print("产物:")
    print(f"  {DECAY_DIR}/distance_decay_loglog.png/svg")
    print(f"  {DECAY_DIR}/distance_decay_alpha.tsv")
    print(f"  {DECAY_DIR}/distance_decay_comparison.tsv")
    print(f"  {DECAY_DIR}/distance_decay_data.tsv")


if __name__ == "__main__":
    main()