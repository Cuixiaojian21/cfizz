"""
Distance decay log-log 曲线图与 TSV 输出模块.

整合自 g_6d_distance_decay.py 参考脚本.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional


# cfizz 通用 7 色调色板 (用于 distance_decay 曲线 fallback 配色, 跨 sample 名稳定)
COLOR_PALETTE: list = [
    '#264653',  # 深青蓝
    '#27736F',  # 青绿
    '#299D92',  # 浅青绿
    '#8AB17C',  # 草绿
    '#E8C56B',  # 金黄
    '#F2A361',  # 橙黄
    '#E66F51',  # 橙红
]


# cfizz 标准 3 sample 配色 (与 g_6d_distance_decay.py 一致)
DEFAULT_SAMPLE_COLORS: Dict[str, str] = {
    'NC':  '#264653',
    'HH':  '#E66F51',
    'PGS': '#299D92',
}


def get_sample_color(sample: str, sample_colors: Optional[Dict[str, str]] = None) -> str:
    """
    解析 sample 的颜色.

    优先级:
      1. ``sample_colors`` dict 中显式映射 (用户指定)
      2. ``DEFAULT_SAMPLE_COLORS`` dict (NC/HH/PGS)
      3. ``COLOR_PALETTE`` 按 sample 名 hash 索引 (任意 sample 名都拿不同色)

    Returns
    -------
    color : str
        '#hex' 颜色字符串.
    """
    if sample_colors and sample in sample_colors:
        return sample_colors[sample]
    if sample in DEFAULT_SAMPLE_COLORS:
        return DEFAULT_SAMPLE_COLORS[sample]
    # fallback: 按 hash 索引 COLOR_PALETTE, 同一样本名跨运行稳定
    idx = hash(sample) % len(COLOR_PALETTE)
    return COLOR_PALETTE[idx]


def plot_distance_decay_curves(
    agg_dict: Dict[str, pd.DataFrame],
    output_path: str,
    sample_colors: Optional[Dict[str, str]] = None,
    dpi: int = 2000,
) -> None:
    """
    画 log-log 距离衰减对比曲线, 每 sample 1 条线.

    Parameters
    ----------
    agg_dict : dict
        ``{sample: groupby('dist_bp').mean() DataFrame}`` 输出, 含
        ``dist_bp``, ``normalized_count``, ``normalized_count_smoothed`` 列.
    output_path : str
        输出前缀 (无扩展名), 会存 ``.png`` + ``.svg``.
    sample_colors : dict, optional
        ``{sample: '#hex'}``, 默认 ``DEFAULT_SAMPLE_COLORS``.
    dpi : int
        图片 dpi, 默认 2000.
    """
    # cfizz 标准尺寸: 6.4cm × 5cm
    fig, ax = plt.subplots(figsize=(6.4 / 2.54, 5.0 / 2.54), dpi=dpi)

    for sample, df in agg_dict.items():
        mask = (df['dist_bp'] > 0) & (df['normalized_count'] > 0)
        x_log = np.log10(df.loc[mask, 'dist_bp'].values)
        y_raw_log = np.log10(df.loc[mask, 'normalized_count'].clip(lower=1e-10).values)
        # 每 sample 不同色 (user 传 sample_colors 优先, 否则按 COLOR_PALETTE 索引)
        sample_color = get_sample_color(sample, sample_colors)

        # log10 域 savgol 平滑 (>20 点启用, 与 g_6d 一致)
        if len(x_log) > 20:
            from scipy.signal import savgol_filter
            w = min(21, max(3, (len(x_log) // 2) * 2 + 1))
            y_plot = savgol_filter(y_raw_log, w, 3)
        else:
            y_plot = y_raw_log

        ax.plot(x_log, y_plot, label=sample, color=sample_color,
                linewidth=0.8, zorder=5)

        # scatter overlay (<=500 点叠加)
        if len(x_log) <= 500:
            ax.scatter(x_log, y_raw_log, s=0.8,
                       color=sample_color,
                       edgecolor='white', linewidth=0.1, alpha=0.7, zorder=6)

    ax.set_xlabel('Log10(Genomic Distance) (bp)', fontsize=6, labelpad=0.01)
    ax.set_ylabel('Log10(Interaction Frequency)', fontsize=6, labelpad=0.01)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.spines['left'].set_linewidth(0.4)
    ax.tick_params(axis='both', labelsize=5, length=1, pad=1)
    ax.legend(loc=(1.01, 0.6), fontsize=5, frameon=False)

    plt.subplots_adjust(top=0.96, bottom=0.16, right=0.75, left=0.125)

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    fig.savefig(f"{output_path}.png", dpi=dpi, bbox_inches='tight')
    fig.savefig(f"{output_path}.svg", bbox_inches='tight')
    plt.close(fig)


def save_decay_alpha_table(
    alpha_dict: Dict[str, Dict],
    output_path: str,
    fit_window: Tuple[float, float] = (0.2e6, 10e6),
    segments: Optional[List[Tuple[float, float]]] = None,
) -> None:
    """
    保存每 sample α + R² + 分段斜率 TSV.

    Parameters
    ----------
    alpha_dict : dict
        ``{sample: fit_decay_alpha_segments() 输出}``.
    output_path : str
        输出 .tsv 路径.
    fit_window : tuple of (lo, hi) in bp
        拟合窗口 (用于 TSV 注释).
    segments : list of (lo_mb, hi_mb), optional
        分段区间 (用于 TSV 注释).
    """
    segs = segments if segments is not None else [
        (0.02, 0.2), (0.2, 1), (1, 5), (5, 20), (20, 50),
    ]
    alpha_df = pd.DataFrame.from_dict(alpha_dict, orient='index')
    alpha_df.index.name = 'sample'
    alpha_df = alpha_df.reset_index()

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w') as f:
        f.write("# Distance decay power-law fit (cfizz distance_decay)\n")
        f.write(f"# normalized_count = count.avg / total_contacts per sample (contact probability proportion)\n")
        f.write(f"# fit window: {fit_window[0]/1e6}-{fit_window[1]/1e6} Mb\n")
        f.write(f"# segments (Mb): {segs}\n")
        f.write("# columns: sample, A_raw/alpha_raw/R2_raw (on normalized_count), slope_<lo>-<hi>Mb\n")
        alpha_df.to_csv(f, sep='\t', index=False)


def save_decay_comparison_table(
    alpha_dict: Dict[str, Dict],
    output_path: str,
    segments: Optional[List[Tuple[float, float]]] = None,
) -> None:
    """
    保存两两 sample α + 分段斜率 diff TSV.

    Parameters
    ----------
    alpha_dict : dict
        ``{sample: fit_decay_alpha_segments() 输出}``.
    output_path : str
        输出 .tsv 路径.
    segments : list of (lo_mb, hi_mb), optional
        分段区间 (用于读取 alpha_dict 中的 slope key).
    """
    segs = segments if segments is not None else [
        (0.02, 0.2), (0.2, 1), (1, 5), (5, 20), (20, 50),
    ]
    comp_rows = []
    samples_list = list(alpha_dict.keys())
    for i in range(len(samples_list)):
        for j in range(i + 1, len(samples_list)):
            s1, s2 = samples_list[i], samples_list[j]
            row = {
                'sample1': s1,
                'sample2': s2,
                'alpha1_raw': alpha_dict[s1]['alpha_raw'],
                'alpha2_raw': alpha_dict[s2]['alpha_raw'],
                'alpha_diff_raw': alpha_dict[s2]['alpha_raw'] - alpha_dict[s1]['alpha_raw'],
            }
            for lo_mb, hi_mb in segs:
                key = f'slope_{lo_mb}-{hi_mb}Mb'
                row[f'{key}_1'] = alpha_dict[s1].get(key, np.nan)
                row[f'{key}_2'] = alpha_dict[s2].get(key, np.nan)
                row[f'{key}_diff'] = row[f'{key}_2'] - row[f'{key}_1']
            comp_rows.append(row)
    comp_df = pd.DataFrame(comp_rows)

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    comp_df.to_csv(output_path, sep='\t', index=False)


def save_decay_data_table(
    agg_dict: Dict[str, pd.DataFrame],
    output_path: str,
) -> None:
    """
    保存完整 P(s) 数据 TSV (每 sample ~5000 行, 跨 region 平均后).

    Parameters
    ----------
    agg_dict : dict
        ``{sample: groupby('dist_bp').mean() DataFrame}`` 输出.
    output_path : str
        输出 .tsv 路径.
    """
    data_rows = []
    for sample, df in agg_dict.items():
        sub = df[['dist_bp', 'normalized_count', 'normalized_count_smoothed', 'total_contacts']].copy()
        sub.insert(0, 'sample', sample)
        data_rows.append(sub)
    data_df = pd.concat(data_rows, ignore_index=True)
    data_df.columns = ['sample', 'dist_bp', 'normalized_count', 'normalized_count_smoothed', 'total_contacts']

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    data_df.to_csv(output_path, sep='\t', index=False)


def run_distance_decay(
    samples: Dict[str, str],
    output_dir: str,
    resolution: int = 1_000,
    max_distance: int = 100_000_000,
    nproc: int = 8,
    fit_window: Tuple[float, float] = (0.2e6, 10e6),
    balance: bool = False,
    excluded_chroms: Optional[Tuple[str, ...]] = ('chrM', 'chrY'),
    sample_colors: Optional[Dict[str, str]] = None,
) -> None:
    """
    端到端 distance decay pipeline: per-sample P(s) + 拟合 + 画图 + 写 TSV.

    一站式调用, 直接得到 4 类产物:
        ``distance_decay_loglog.png/svg``
        ``distance_decay_alpha.tsv``
        ``distance_decay_comparison.tsv``
        ``distance_decay_data.tsv``

    Parameters
    ----------
    samples : dict
        ``{sample_name: mcool_path}``.
    output_dir : str
        输出目录.
    其余参数: 同 compute_distance_decay_per_sample() / plot_distance_decay_curves().
    """
    from cfizz.analyze.distance import (
        compute_distance_decay_per_sample,
        DEFAULT_SEGMENTS_MB,
    )

    os.makedirs(output_dir, exist_ok=True)

    print(f"[distance_decay] per-sample 计算 ({len(samples)} 个, nproc={nproc})...")
    data_dict, agg_dict, alpha_dict = compute_distance_decay_per_sample(
        samples,
        resolution=resolution,
        nproc=nproc,
        max_distance=max_distance,
        fit_window=fit_window,
        balance=balance,
        excluded_chroms=excluded_chroms,
    )
    for s, df in agg_dict.items():
        print(f"  ✓ {s}: {len(df)} 个距离点 (跨 region 平均后)")

    print(f"[distance_decay] 画 log-log 图...")
    plot_distance_decay_curves(
        agg_dict,
        output_path=f"{output_dir}/distance_decay_loglog",
        sample_colors=sample_colors,
    )

    print(f"[distance_decay] 保存 3 类 TSV...")
    save_decay_alpha_table(
        alpha_dict,
        output_path=f"{output_dir}/distance_decay_alpha.tsv",
        fit_window=fit_window,
    )
    save_decay_comparison_table(
        alpha_dict,
        output_path=f"{output_dir}/distance_decay_comparison.tsv",
    )
    save_decay_data_table(
        agg_dict,
        output_path=f"{output_dir}/distance_decay_data.tsv",
    )
    print(f"[distance_decay] 产物: {output_dir}/")