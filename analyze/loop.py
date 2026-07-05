"""
Chromatin loop analysis module.

This module provides functions for chromatin loop detection
and analysis using peak calling methods.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)


def call_loops_hiccups(
    matrix: np.ndarray,
    resolution: int,
    window_size: int = 5,          # 1_7 标准
    peak_width: int = 2,           # 1_7 标准
    only_anchors: bool = True,     # 新加(1_7 有)
    nproc: int = 1,                # 新加(多进程)
    use_cache: bool = True,        # 新加(缓存)
    cache_path: str = None,        # 新加(缓存路径)
    FDR_threshold: float = 0.1,
    max_distance: int = 20_000_000  # 1_7 标准(20Mb)
) -> pd.DataFrame:
    """
    Call chromatin loops using a simplified HiCCUPS-style algorithm.
    (1_7 对齐版: 支持多进程并行 + 灵活分层缓存)
    
    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    resolution : int
        Bin resolution (bp)
    window_size : int
        Window size for peak detection (bins) [1_7 标准: 5]
    peak_width : int
        Width for peak area [1_7 标准: 2]
    only_anchors : bool
        Only report anchor regions [1_7 有]
    nproc : int
        Number of parallel processes [1_7 有]
    use_cache : bool
        Use cached results if available
    cache_path : str
        Path to cache file
    FDR_threshold : float
        False discovery rate threshold
    max_distance : int
        Maximum loop distance (bp) [1_7 标准: 20Mb]
        
    Returns
    -------
    loops : pd.DataFrame
        DataFrame with loop coordinates
    """
    import os
    import multiprocessing
    
    # 缓存检查(灵活分层)
    if use_cache and cache_path and os.path.exists(cache_path):
        try:
            cached_df = pd.read_csv(cache_path)
            if len(cached_df) > 0:
                return cached_df
        except Exception:
            pass  # 缓存读取失败,继续计算
    
    n = matrix.shape[0]
    max_distance_bins = max_distance // resolution
    
    # Normalize matrix
    normalized = _normalize_for_loop_calling(matrix)
    
    # Find peaks
    loops = []
    
    # 多进程并行(外层 for i 分配到多进程)
    if nproc > 1:
        def process_row(args):
            i, n, window_size, max_distance_bins, normalized, peak_width, resolution = args
            row_loops = []
            for j in range(i + window_size, min(n, i + max_distance_bins)):
                loop = _find_peak_at_position(
                    i, j, n, normalized, peak_width, resolution, window_size
                )
                if loop:
                    row_loops.append(loop)
            return row_loops
        
        # 准备任务参数
        tasks = [(i, n, window_size, max_distance_bins, normalized, peak_width, resolution) 
                 for i in range(n)]
        
        # 多进程执行
        with multiprocessing.Pool(nproc) as pool:
            results = pool.map(process_row, tasks)
        
        # 合并结果
        for row_loops in results:
            loops.extend(row_loops)
    else:
        # 单进程模式
        for i in range(n):
            for j in range(i + window_size, min(n, i + max_distance_bins)):
                loop = _find_peak_at_position(
                    i, j, n, normalized, peak_width, resolution, window_size
                )
                if loop:
                    loops.append(loop)
    
    if not loops:
        result_df = pd.DataFrame(columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'score'])
    else:
        loops_df = pd.DataFrame(loops)
        loops_df = loops_df.sort_values('score', ascending=False).reset_index(drop=True)
        
        # Remove nearby duplicates
        loops_df = _remove_nearby_loops(loops_df, window_size)
        result_df = loops_df
    
    # 保存缓存
    if use_cache and cache_path:
        try:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            result_df.to_csv(cache_path, index=False)
        except Exception:
            pass  # 缓存保存失败不影响返回
    
    return result_df


def _find_peak_at_position(i, j, n, normalized, peak_width, resolution, window_size):
    """检查指定位置是否是峰值"""
    center_value = normalized[i, j]
    
    # Check if local maximum
    is_peak = True
    local_mean = 0
    local_count = 0
    
    for di in range(-peak_width, peak_width + 1):
        for dj in range(-peak_width, peak_width + 1):
            if di == 0 and dj == 0:
                continue
            ni, nj = i + di, j + dj
            if 0 <= ni < n and 0 <= nj < n:
                local_mean += normalized[ni, nj]
                local_count += 1
                if normalized[ni, nj] > center_value:
                    is_peak = False
    
    if is_peak and local_count > 0:
        local_mean /= local_count
        
        # Calculate z-score-like statistic
        if local_mean > 0:
            score = (center_value - local_mean) / np.sqrt(local_mean)
        else:
            score = 0
        
        if score > 1.5:  # Lower threshold for sensitivity
            return {
                'bin1': i,
                'bin2': j,
                'pos1': i * resolution,
                'pos2': j * resolution,
                'score': score,
                'value': normalized[i, j]
            }
    return None


def _normalize_for_loop_calling(matrix: np.ndarray) -> np.ndarray:
    """
    Normalize matrix for loop calling.
    
    Simple row-column normalization.
    """
    n = matrix.shape[0]
    normalized = np.zeros_like(matrix, dtype=float)
    
    row_sums = np.sum(matrix, axis=1)
    col_sums = np.sum(matrix, axis=0)
    total = np.sum(matrix)
    
    for i in range(n):
        for j in range(n):
            if row_sums[i] > 0 and col_sums[j] > 0:
                expected = (row_sums[i] * col_sums[j]) / total
                if expected > 0:
                    normalized[i, j] = matrix[i, j] / expected
    
    return normalized


def _remove_nearby_loops(loops_df: pd.DataFrame, min_distance: int) -> pd.DataFrame:
    """Remove loops that are too close to each other."""
    if len(loops_df) == 0:
        return loops_df
    
    keep = [0]
    last_pos1 = loops_df.iloc[0]['pos1']
    last_pos2 = loops_df.iloc[0]['pos2']
    
    for i in range(1, len(loops_df)):
        pos1 = loops_df.iloc[i]['pos1']
        pos2 = loops_df.iloc[i]['pos2']
        
        if pos1 - last_pos1 > min_distance or pos2 - last_pos2 > min_distance:
            keep.append(i)
            last_pos1 = pos1
            last_pos2 = pos2
    
    return loops_df.iloc[keep].reset_index(drop=True)


def calculate_loop_score(
    matrix: np.ndarray,
    anchor1_pos: int,
    anchor2_pos: int,
    resolution: int,
    window_size: int = 5
) -> float:
    """
    Calculate loop strength score at a specific position.
    
    Parameters
    ----------
    matrix : np.ndarray
        Contact matrix
    anchor1_pos : int
        First anchor position (bp)
    anchor2_pos : int
        Second anchor position (bp)
    resolution : int
        Bin resolution (bp)
    window_size : int
        Window size for local background
        
    Returns
    -------
    score : float
        Loop score
    """
    bin1 = anchor1_pos // resolution
    bin2 = anchor2_pos // resolution
    
    n = matrix.shape[0]
    
    if bin1 < 0 or bin1 >= n or bin2 < 0 or bin2 >= n:
        return 0
    
    # Get center value
    center = matrix[bin1, bin2]
    
    # Calculate local background
    background = 0
    count = 0
    
    for di in range(-window_size, window_size + 1):
        for dj in range(-window_size, window_size + 1):
            if di == 0 and dj == 0:
                continue
            ni, nj = bin1 + di, bin2 + dj
            if 0 <= ni < n and 0 <= nj < n:
                background += matrix[ni, nj]
                count += 1
    
    if count > 0:
        background /= count
    
    if background > 0:
        return center / background
    return 0


def filter_loops_by_distance(
    loops: pd.DataFrame,
    min_distance: int = 0,
    max_distance: int = 20000000
) -> pd.DataFrame:
    """
    Filter loops by genomic distance.
    
    Parameters
    ----------
    loops : pd.DataFrame
        DataFrame with loop coordinates
    min_distance : int
        Minimum loop distance (bp)
    max_distance : int
        Maximum loop distance (bp)
        
    Returns
    -------
    filtered : pd.DataFrame
        Filtered loops
    """
    if len(loops) == 0:
        return loops
    
    if 'start1' in loops.columns and 'end1' in loops.columns:
        distance = loops['start2'] - loops['end1']
    elif 'pos1' in loops.columns and 'pos2' in loops.columns:
        distance = loops['pos2'] - loops['pos1']
    else:
        return loops
    
    return loops[(distance >= min_distance) & (distance <= max_distance)]

# =============================================================================
# =============================================================================

def call_loops(
    input_cool: str,
    output_file: str,
    peak_widths: list = None,
    window_widths: list = None,
    only_anchors: bool = True,
    resolution: int = 10000,
    nproc: int = 1,
    verbose: bool = True
) -> str:
    """
    使用pyHICCUPS调用染色质环
    
    
    参数:
        input_cool (str): 输入的cool文件路径
        output_file (str): 输出文件路径
        peak_widths (List[int]): 峰值宽度列表，默认[2]（10kb分辨率推荐参数）
        window_widths (List[int]): 窗口宽度列表，默认[5]（10kb分辨率推荐参数）
        only_anchors (bool): 是否只输出锚点，默认True
        resolution (int): 分辨率，默认10000（10kb）
        nproc (int): 并行处理的核心数，默认1
        verbose (bool): 是否显示进度信息，默认True
    
    返回:
        str: 输出文件路径
    """
    import subprocess
    from pathlib import Path
    from cfizz.utils.coordinates import ensure_dir
    
    if peak_widths is None:
        peak_widths = [2]
    if window_widths is None:
        window_widths = [5]
    
    # 确保输出目录存在
    output_file = Path(output_file)
    ensure_dir(output_file.parent)
    
    # 构建pyHICCUPS命令
    cmd = [
        "pyHICCUPS",
        "-O", str(output_file),
        "-p", f"{input_cool}::/resolutions/{resolution}",
        "--pw", *map(str, peak_widths),
        "--ww", *map(str, window_widths),
        "--nproc", str(nproc)
    ]
    
    if only_anchors:
        cmd.append("--only-anchors")
    
    if verbose:
        print(f"运行命令: {' '.join(cmd)}")
    
    # 执行命令
    try:
        subprocess.run(cmd, check=True)
        if verbose:
            print(f"Loop检测完成，结果保存在: {output_file}")
        return str(output_file)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"pyHICCUPS运行失败: {str(e)}")


def process_loops(
    mcool_path: str,
    resolution: int,
    output_dir: str,
    peak_widths: list = None,
    window_widths: list = None,
    only_anchors: bool = True,
    nproc: int = 1,
    verbose: bool = True
) -> str:
    """
    处理Hi-C数据并调用染色质环
    
    
    参数:
        mcool_path (str): mcool文件路径
        resolution (int): 分辨率
        output_dir (str): 输出目录
        peak_widths (List[int]): 峰值宽度列表，默认[2]（10kb分辨率推荐参数）
        window_widths (List[int]): 窗口宽度列表，默认[5]（10kb分辨率推荐参数）
        only_anchors (bool): 是否只输出锚点，默认True
        nproc (int): 并行处理的核心数，默认1
        verbose (bool): 是否显示进度信息，默认True
    
    返回:
        str: 输出文件路径
    """
    import os
    from cfizz.utils.coordinates import ensure_dir
    
    if peak_widths is None:
        peak_widths = [2]
    if window_widths is None:
        window_widths = [5]
    
    # 确保输出目录存在
    ensure_dir(output_dir)
    
    # 构建输出文件路径
    basename = os.path.basename(mcool_path).split('.')[0]
    output_file = os.path.join(
        output_dir,
        f"{basename}.{resolution//1000}k.loops.txt"
    )
    
    if verbose:
        print(f"\n开始检测染色质环...")
    
    # 调用loops
    output_file = call_loops(
        input_cool=mcool_path,
        output_file=output_file,
        peak_widths=peak_widths,
        window_widths=window_widths,
        only_anchors=only_anchors,
        resolution=resolution,
        nproc=nproc,
        verbose=verbose
    )
    
    if verbose:
        print(f"  - 完成！")
    
    return output_file


# =============================================================================
# Loop Differential Analysis Functions (移植自 script/a_6_loops_diff.py)
# Source: a_6 L59, L90, L167, L188, L214
# =============================================================================

from pathlib import Path
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging


# 配色方案 (a_6 配色)
LOOP_COLORS = {
    'Stable': '#264653',      # 深青色
    'Low_FE': '#299D92',      # 青色
    'Medium_FE': '#F2A361',   # 橙色
    'High_FE': '#E66F51'      # 红橙色
}


def load_loops_data(file_path: str) -> pd.DataFrame:
    """
    加载loops数据文件。
    Source: a_6 L59
    
    Parameters
    ----------
    file_path : str
        HiCCUPS loops输出文件路径
        
    Returns
    -------
    pd.DataFrame
        过滤后的显著loops DataFrame (q_product < 0.01)
        
    行为特征 6 问:
    1. 输入: HiCCUPS原始输出文件
    2. 输出: 过滤后DataFrame [chrom1, start1, end1, chrom2, start2, end2, ...]
    3. 跟5_1产物关系: 读4_1 calls_loops输出的loops.txt，不重算
    4. 跟a_6算法关系: 完全一致，L59-87
    5. 调谁: pd.read_csv, dropna
    6. 错误处理: numeric_cols dropna
    """
    print(f"  正在读取: {file_path}")
    column_names = [
        'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
        'col7', 'col8', 'col9', 'col10',
        'donut_fe', 'pvalue1', 'donut_q',
        'll_fe', 'pvalue2', 'll_q'
    ]
    
    df = pd.read_csv(file_path, sep='\t', header=None, names=column_names)
    print(f"    -> Loops数量: {len(df)}")
    
    numeric_cols = ['start1', 'end1', 'start2', 'end2', 'donut_fe', 'donut_q', 'll_fe', 'll_q']
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    df = df.dropna(subset=numeric_cols)
    
    df['center1'] = (df['start1'] + df['end1']) / 2
    df['center2'] = (df['start2'] + df['end2']) / 2
    df['loop_distance'] = abs(df['center2'] - df['center1'])
    df['q_product'] = df['donut_q'] * df['ll_q']
    
    significant = df['q_product'] < 0.01
    print(f"    -> 显著loops数 (q_product < 0.01): {significant.sum()}")
    df_significant = df[significant].copy()
    
    return df_significant


def find_matching_loops(
    loops1: pd.DataFrame, 
    loops2: pd.DataFrame, 
    sample1_name: str, 
    sample2_name: str, 
    distance_threshold: int = 5
) -> dict:
    """
    基于KDTree最近距离匹配loops。
    Source: a_6 L90
    
    核心函数。使用cKDTree在bin空间最近邻匹配。
    
    Parameters
    ----------
    loops1, loops2 : pd.DataFrame
        两个样本的loops数据
    sample1_name, sample2_name : str
        样本名称
    distance_threshold : int
        距离阈值（bin数），默认5
        
    Returns
    -------
    dict
        包含配对结果和统计信息的字典
    """
    print(f"    开始基于KDTree最近距离匹配loops... (距离阈值: {distance_threshold}个bin)")
    
    loops1 = loops1.copy()
    loops2 = loops2.copy()
    
    loops1['bin1'] = (loops1['start1'] // 10000).astype(int)
    loops1['bin2'] = (loops1['start2'] // 10000).astype(int)
    loops2['bin1'] = (loops2['start1'] // 10000).astype(int)
    loops2['bin2'] = (loops2['start2'] // 10000).astype(int)
    
    for df in [loops1, loops2]:
        swap_mask = df['bin1'] > df['bin2']
        df.loc[swap_mask, ['bin1', 'bin2']] = df.loc[swap_mask, ['bin2', 'bin1']].values
    
    loops1['loop_key'] = loops1['chrom1'] + '_' + loops1['bin1'].astype(str) + '_' + loops1['bin2'].astype(str)
    loops2['loop_key'] = loops2['chrom1'] + '_' + loops2['bin1'].astype(str) + '_' + loops2['bin2'].astype(str)
    
    loops1['min_distance'] = float('inf')
    loops1['matched_idx2'] = None
    loops1['is_stable'] = False
    
    chroms1 = loops1['chrom1'].unique()
    matched_indices_2 = set()
    
    for chrom in chroms1:
        mask1 = loops1['chrom1'] == chrom
        mask2 = loops2['chrom1'] == chrom
        
        loops1_chrom = loops1[mask1]
        loops2_chrom = loops2[mask2]
        
        if len(loops2_chrom) == 0:
            continue
        
        coords2 = loops2_chrom[['bin1', 'bin2']].values.astype(float)
        tree = cKDTree(coords2)
        
        coords1 = loops1_chrom[['bin1', 'bin2']].values.astype(float)
        distances, indices = tree.query(coords1, k=1)
        
        for i, (dist, idx) in enumerate(zip(distances, indices)):
            loop1_idx = loops1_chrom.index[i]
            loop2_idx = loops2_chrom.index[idx]
            
            loops1.loc[loop1_idx, 'min_distance'] = dist
            loops1.loc[loop1_idx, 'matched_idx2'] = loop2_idx
            
            if dist <= distance_threshold:
                loops1.loc[loop1_idx, 'is_stable'] = True
                matched_indices_2.add(loop2_idx)
    
    stable_count = int(loops1['is_stable'].sum())
    changed_count = len(loops1) - stable_count
    
    common_loops = stable_count
    unique_to_1 = changed_count
    unique_to_2 = len(loops2) - len(matched_indices_2)
    union_loops = len(loops1) + len(loops2) - common_loops
    
    loops2['is_stable'] = loops2.index.isin(matched_indices_2)
    
    result = {
        f'{sample1_name}_total': len(loops1),
        f'{sample2_name}_total': len(loops2),
        'common_loops': common_loops,
        'unique_to_sample1': unique_to_1,
        'unique_to_sample2': unique_to_2,
        'union_loops': union_loops,
        'distance_threshold': distance_threshold,
        'loops1_df': loops1,
        'loops2_df': loops2
    }
    
    return result


def calculate_sample_stats(pairing_df: pd.DataFrame, sample_name: str) -> dict:
    """
    计算样本统计信息。
    Source: a_6 L167
    """
    sample_data = pairing_df[pairing_df['loop_source'] == sample_name]
    total = len(sample_data)
    
    stats = {
        'sample': sample_name,
        'total': total,
        'Stable': len(sample_data[sample_data['quality_group'] == 'Stable']),
        'Low_FE': len(sample_data[sample_data['quality_group'] == 'Low_FE']),
        'Medium_FE': len(sample_data[sample_data['quality_group'] == 'Medium_FE']),
        'High_FE': len(sample_data[sample_data['quality_group'] == 'High_FE'])
    }
    
    stats['Stable_pct'] = stats['Stable'] / total * 100 if total > 0 else 0
    stats['Low_FE_pct'] = stats['Low_FE'] / total * 100 if total > 0 else 0
    stats['Medium_FE_pct'] = stats['Medium_FE'] / total * 100 if total > 0 else 0
    stats['High_FE_pct'] = stats['High_FE'] / total * 100 if total > 0 else 0
    
    return stats


def save_plot_data_loops(stats_list: list, comparison_dir: Path) -> Path:
    """保存Loops绘图数据"""
    plot_data = pd.DataFrame(stats_list)
    plot_data_file = comparison_dir / "loops_plot_data.tsv"
    plot_data.to_csv(plot_data_file, sep='\t', index=False)
    return plot_data_file


def load_plot_data_loops(output_path: Path) -> dict:
    """加载Loops绘图数据"""
    plot_data_file = output_path / "loops_plot_data.tsv"
    
    if not plot_data_file.exists():
        raise FileNotFoundError(f"绘图数据文件不存在: {plot_data_file}")
    
    plot_data = pd.read_csv(plot_data_file, sep='\t')
    
    data_dict = {}
    for _, row in plot_data.iterrows():
        sample = row['sample']
        data_dict[sample] = {
            'total': int(row['total']),
            'Stable': int(row['Stable']),
            'Low_FE': int(row['Low_FE']),
            'Medium_FE': int(row['Medium_FE']),
            'High_FE': int(row['High_FE']),
            'Stable_pct': row['Stable_pct'],
            'Low_FE_pct': row['Low_FE_pct'],
            'Medium_FE_pct': row['Medium_FE_pct'],
            'High_FE_pct': row['High_FE_pct']
        }
    
    return data_dict


def plot_loops_stacked_bar(data_dict: dict, comparison: str, output_dir: Path) -> None:
    """
    绘制Loops分类堆积柱状图。
    Source: a_6 (参考a_5的绘图逻辑)
    """
    print(f"\n生成 {comparison} 堆积柱状图...")
    
    width_cm, height_cm = 6.4, 5
    fig, ax = plt.subplots(figsize=(width_cm / 2.54, height_cm / 2.54))
    
    samples = list(data_dict.keys())
    n_samples = len(samples)
    x_positions = np.arange(n_samples)
    width = 0.8
    
    categories = ['Stable', 'Low_FE', 'Medium_FE', 'High_FE']
    category_labels = ['stable (shared)', 'unique (low FE)', 'unique (medium FE)', 'unique (high FE)']
    
    bottom = np.zeros(n_samples)
    
    for i, cat in enumerate(categories):
        pct_values = [data_dict[sample].get(f'{cat}_pct', 0) for sample in samples]
        count_values = [data_dict[sample].get(cat, 0) for sample in samples]
        
        bars = ax.bar(x_positions, pct_values, width, bottom=bottom, label=category_labels[i], color=LOOP_COLORS[cat])
        
        for j, (bar, pct_val, count_val) in enumerate(zip(bars, pct_values, count_values)):
            if pct_val > 1:
                height = bar.get_height()
                y_pos = bottom[j] + height / 2
                ax.text(bar.get_x() + bar.get_width()/2, y_pos, f'{pct_val:.1f}%({int(count_val)})',
                       ha='center', va='center', fontsize=5, color='black')
        
        bottom = [b + p for b, p in zip(bottom, pct_values)]
    
    ax.set_xlabel('Sample')
    ax.set_ylabel('Percentage (%)')
    ax.set_xticks(x_positions)
    ax.set_xticklabels(samples)
    ax.set_ylim(0, 100)
    ax.legend(title='loop change', loc=(1.01, 0.6), frameon=False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    plt.subplots_adjust(top=0.96, bottom=0.16, right=0.75, left=0.125)
    
    comparison_clean = comparison.replace('--', '_')
    output_png = output_dir / f"loops_{comparison_clean}_stacked_bar.png"
    output_svg = output_dir / f"loops_{comparison_clean}_stacked_bar.svg"
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_svg, format='svg', bbox_inches='tight')
    plt.close()
    
    print(f"  -> PNG: {output_png}")
    print(f"  -> SVG: {output_svg}")


def analyze_single_comparison_loops(
    comparison: str,
    output_root: str,
    run_mode: str,
    control_loops_path: str = None,
    treatment_loops_path: str = None
) -> None:
    """
    分析单个比较组的Loops差异。
    Source: a_6 L214
    
    主流程函数，整合loops配对和质量分类。
    """
    def setup_logger(output_path, comparison):
        log_dir = output_path / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"loops_{comparison}.log"
        logger = logging.getLogger(f"loops_diff_{comparison}")
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
    
    def parse_comparison_func(comparison_str):
        if '--' in comparison_str:
            parts = comparison_str.split('--')
            return parts[0].strip(), parts[1].strip()
        return comparison_str, comparison_str
    
    treatment, control = parse_comparison_func(comparison)
    
    comparison_clean = comparison.replace('--', '_')
    comparison_dir = Path(output_root) / "loops"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    
    logger = setup_logger(comparison_dir, comparison)
    
    logger.info("=" * 70)
    logger.info(f"Loops差异分析: {comparison}")
    logger.info(f"处理组: {treatment}")
    logger.info(f"对照组: {control}")
    logger.info(f"运行模式: {run_mode}")
    logger.info("=" * 70)
    
    try:
        if run_mode in ['all', 'compute']:
            if control_loops_path is None or treatment_loops_path is None:
                logger.warning("请提供 control_loops_path 和 treatment_loops_path 参数")
                return
            
            logger.info("\n1. 加载loops数据")
            loops_control = load_loops_data(control_loops_path)
            loops_treatment = load_loops_data(treatment_loops_path)
            
            logger.info(f"\n2. 执行 {treatment} vs {control} 匹配分析...")
            match_result = find_matching_loops(loops_control, loops_treatment, control, treatment, distance_threshold=5)
            
            logger.info(f"    并集总数: {match_result['union_loops']}")
            logger.info(f"    共有loops: {match_result['common_loops']} ({match_result['common_loops']/match_result['union_loops']*100:.1f}%)")
            logger.info(f"    {treatment}独有: {match_result['unique_to_sample2']} ({match_result['unique_to_sample2']/match_result['union_loops']*100:.1f}%)")
            logger.info(f"    {control}独有: {match_result['unique_to_sample1']} ({match_result['unique_to_sample1']/match_result['union_loops']*100:.1f}%)")
            
            loops1_with_key = match_result['loops1_df']
            loops2_with_key = match_result['loops2_df']
            
            logger.info(f"\n3. 保存独有loops")
            unique_treatment = loops2_with_key[~loops2_with_key['is_stable']].copy()
            unique_control = loops1_with_key[~loops1_with_key['is_stable']].copy()
            
            unique_treatment_file = comparison_dir / f"1_1_{treatment}_unique_loops.tsv"
            unique_control_file = comparison_dir / f"1_1_{control}_unique_loops_vs_{treatment}.tsv"
            
            unique_treatment.to_csv(unique_treatment_file, sep='\t', index=False)
            unique_control.to_csv(unique_control_file, sep='\t', index=False)
            
            logger.info(f"    {treatment}独有loops保存至: {unique_treatment_file}")
            logger.info(f"    {control}独有loops保存至: {unique_control_file}")
            
            logger.info(f"\n4. 生成分类文件")
            loops1_with_key['fe_product'] = loops1_with_key['donut_fe'] * loops1_with_key['ll_fe']
            loops2_with_key['fe_product'] = loops2_with_key['donut_fe'] * loops2_with_key['ll_fe']
            
            unique_loops1 = loops1_with_key[~loops1_with_key['is_stable']].copy()
            unique_loops2 = loops2_with_key[~loops2_with_key['is_stable']].copy()
            
            unique1_fe = unique_loops1['fe_product']
            unique2_fe = unique_loops2['fe_product']
            
            if len(unique1_fe) > 0:
                fe_q1_1, fe_q2_1 = unique1_fe.quantile(0.33), unique1_fe.quantile(0.67)
            else:
                fe_q1_1 = fe_q2_1 = 0
            
            if len(unique2_fe) > 0:
                fe_q1_2, fe_q2_2 = unique2_fe.quantile(0.33), unique2_fe.quantile(0.67)
            else:
                fe_q1_2 = fe_q2_2 = 0
            
            def classify_unique_loop_fe(fe_product):
                if fe_product <= min(fe_q1_1, fe_q1_2):
                    return "Low_FE"
                elif fe_product <= min(fe_q2_1, fe_q2_2):
                    return "Medium_FE"
                else:
                    return "High_FE"
            
            loops1_with_key['quality_group'] = loops1_with_key.apply(
                lambda row: 'Stable' if row['is_stable'] else classify_unique_loop_fe(row['fe_product']),
                axis=1
            )
            loops2_with_key['quality_group'] = loops2_with_key.apply(
                lambda row: 'Stable' if row['is_stable'] else classify_unique_loop_fe(row['fe_product']),
                axis=1
            )
            
            logger.info(f"\n5. 保存绘图数据")
            stats_list = [
                calculate_sample_stats(loops1_with_key.assign(loop_source=control), control),
                calculate_sample_stats(loops2_with_key.assign(loop_source=treatment), treatment)
            ]
            save_plot_data_loops(stats_list, comparison_dir)
        
        if run_mode in ['all', 'plot']:
            logger.info(f"\n6. 生成统计图")
            plot_data_file = comparison_dir / "loops_plot_data.tsv"
            
            if not plot_data_file.exists():
                logger.error(f"绘图数据文件不存在: {plot_data_file}")
                return
            
            data_dict = load_plot_data_loops(comparison_dir)
            plot_loops_stacked_bar(data_dict, comparison_clean, comparison_dir)
        
        logger.info("\n" + "=" * 70)
        logger.info(f"分析完成！")
        logger.info(f"输出目录: {comparison_dir}")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error(f"分析失败: {str(e)}")
        raise
