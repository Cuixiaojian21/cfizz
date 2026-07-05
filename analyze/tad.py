"""
TAD (Topologically Associating Domain) analysis module.

This module provides functions for TAD identification and analysis
using insulation score and boundary detection methods.

Key functions:
  - identify_boundaries(): Identify TAD boundaries from insulation scores
  - extract_tads(): Extract TAD regions from boundaries
  - process_tads(): Complete TAD analysis pipeline

"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
import warnings
import os
from pathlib import Path

warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# Boundary identification
# =============================================================================

def identify_boundaries(
    insulation_df: pd.DataFrame,
    window_size: int,
    threshold_method: str = 'otsu'
) -> pd.DataFrame:
    """
    Identify TAD boundaries from insulation scores using threshold method.

    Uses Otsu's method or Li's method to automatically determine
    the threshold for boundary identification.

    Parameters
    ----------
    insulation_df : pd.DataFrame
        Insulation scores DataFrame from cooltools
    window_size : int
        Window size in bp (e.g., 100000 for 100kb)
    threshold_method : str
        Threshold method: 'otsu' or 'li'

    Returns
    -------
    pd.DataFrame
        Boundaries DataFrame with columns [chrom, start, end, ...]

    Notes
    -----
    The insulation_df must contain:
    - 'chrom', 'start', 'end' columns
    - 'is_boundary_{window_size}' column (from cooltools)

    Correspondence with hicviz:
    >>> from hicviz.core.tads import identify_boundaries
    >>> # Equivalent to hicviz's identify_boundaries()
    """
    # Use is_boundary column from cooltools (most reliable)
    is_boundary_col = f'is_boundary_{window_size}'

    if is_boundary_col not in insulation_df.columns:
        # Try boundary_strength column
        strength_col = f'boundary_strength_{window_size}'
        if strength_col in insulation_df.columns:
            # Use threshold method on boundary strength
            strength_values = insulation_df[strength_col].dropna().values
            if len(strength_values) > 0:
                from skimage.filters import threshold_li, threshold_otsu
                if threshold_method == 'otsu':
                    threshold = threshold_otsu(strength_values)
                else:
                    threshold = threshold_li(strength_values)
                boundaries = insulation_df[insulation_df[strength_col] >= threshold].copy()
                return boundaries
        raise ValueError(f"Neither '{is_boundary_col}' nor 'boundary_strength_{window_size}' found in insulation_df. "
                        f"Available columns: {insulation_df.columns.tolist()}")

    # Filter by is_boundary column
    boundaries = insulation_df[insulation_df[is_boundary_col]].copy()

    return boundaries


def identify_boundaries_simple(
    insulation_df: pd.DataFrame,
    window_size: int,
    threshold: float = 0.1
) -> pd.DataFrame:
    """
    Identify TAD boundaries using a simple threshold.

    This is a simpler alternative to identify_boundaries() that uses
    a fixed threshold value instead of automatic threshold detection.

    Parameters
    ----------
    insulation_df : pd.DataFrame
        Insulation scores DataFrame
    window_size : int
        Window size in bp
    threshold : float
        Minimum boundary strength threshold

    Returns
    -------
    pd.DataFrame
        Boundaries DataFrame
    """
    # Try is_boundary first
    is_boundary_col = f'is_boundary_{window_size}'
    if is_boundary_col in insulation_df.columns:
        return insulation_df[insulation_df[is_boundary_col]].copy()

    # Fall back to boundary_strength
    strength_col = f'boundary_strength_{window_size}'
    if strength_col in insulation_df.columns:
        boundaries = insulation_df[insulation_df[strength_col] >= threshold].copy()
        return boundaries

    raise ValueError(f"Neither '{is_boundary_col}' nor '{strength_col}' found")


# =============================================================================
# TAD extraction
# =============================================================================

def extract_tads(
    insulation_df: pd.DataFrame,
    window_size: int,
    max_tad_length: int = 3_000_000
) -> pd.DataFrame:
    """
    Extract TAD regions from insulation scores using boundary information.
    
    Uses bioframe.merge() to merge non-boundary regions into TADs.
    
    Parameters
    ----------
    insulation_df : pd.DataFrame
        Insulation scores DataFrame with is_boundary column
    window_size : int
        Window size in bp
    max_tad_length : int
        Maximum TAD length in bp (default 3Mb)
        
    Returns
    -------
    pd.DataFrame
        TADs DataFrame with columns [chrom, start, end]
        
    Notes
    -----
    The insulation_df must contain:
    - 'chrom', 'start', 'end' columns
    - 'is_boundary_{window_size}' column

    Correspondence with hicviz:
    >>> from hicviz.core.tads import extract_tads
    >>> # Equivalent to hicviz's extract_tads()
    """
    import bioframe

    # Get the is_boundary column name
    is_boundary_col = f'is_boundary_{window_size}'

    if is_boundary_col not in insulation_df.columns:
        raise ValueError(f"Column '{is_boundary_col}' not found in insulation_df. "
                        f"Available columns: {insulation_df.columns.tolist()}")

    # Filter for non-boundary regions
    non_boundary = insulation_df[insulation_df[is_boundary_col] == False].copy()

    # Use bioframe.merge to merge adjacent non-boundary regions
    if len(non_boundary) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])

    # Merge regions
    tads = bioframe.merge(non_boundary[['chrom', 'start', 'end']])

    # Filter by maximum length
    tad_lengths = tads['end'] - tads['start']
    tads = tads[tad_lengths <= max_tad_length].reset_index(drop=True)

    return tads[['chrom', 'start', 'end']]


def calculate_tad_statistics(
    tads: pd.DataFrame
) -> Dict[str, Any]:
    """
    Calculate statistics for a TAD set.

    Parameters
    ----------
    tads : pd.DataFrame
        DataFrame with TAD coordinates

    Returns
    -------
    dict
        Dictionary with statistics
    """
    if len(tads) == 0:
        return {
            'n_tads': 0,
            'mean_size': 0,
            'median_size': 0,
            'min_size': 0,
            'max_size': 0
        }

    sizes = tads['end'] - tads['start']

    return {
        'n_tads': len(tads),
        'mean_size': float(np.mean(sizes)),
        'median_size': float(np.median(sizes)),
        'min_size': float(np.min(sizes)),
        'max_size': float(np.max(sizes)),
        'std_size': float(np.std(sizes))
    }



# =============================================================================
# Wrapper functions (for compatibility)
# =============================================================================

def calculate_insulation(
    matrix: np.ndarray,
    window_sizes: List[int] = [5, 10, 20]
) -> pd.DataFrame:
    """
    Calculate insulation scores for TAD boundary identification.
    
    This is a wrapper around cfizz.io.insulation functions for convenience.
    
    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    window_sizes : list
        List of window sizes (in bins)
        
    Returns
    -------
    pd.DataFrame
        DataFrame with insulation scores and boundary calls
    """
    from cfizz.io.insulation import calculate_insulation_score
    return calculate_insulation_score(matrix, window_sizes)


# Alias for backwards compatibility
def call_tads(boundaries, matrix, resolution, min_tad_size=5, max_tad_size=100):
    """
    Call TADs based on boundaries (deprecated, use extract_tads instead).
    """
    raise DeprecationWarning(
        "call_tads() is deprecated, use extract_tads() instead"
    )


# =============================================================================
# TAD Differential Analysis Functions (移植自 script/a_5_tad_diff.py)
# Source: a_5 L42, L119, L226, L340, L401, L430, L454, L467, L523, L559
# =============================================================================

from pathlib import Path
from scipy.spatial import KDTree
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 配色方案 (a_5 L76-79)
TAD_COLORS = {
    'Stable_boundary': '#264653',
    'Boundary_shift': '#299D92',
    'Unique_boundary': '#E66F51'
}

# 窗口倍数 (a_5 L42)
WINDOW_MULTIPLES = [50, 10, 5]

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 5


def get_actual_window(window_mult: int, resolution: int = 10000) -> int:
    """将倍数转换为实际的窗口大小(bp)"""
    return window_mult * resolution


def create_tad_annotation_fast(tad_file: str, chrom_sizes_file: str, output_file: str, bin_size: int = 10000):
    """
    创建全基因组指定大小bin注释文件。
    Source: a_5 L119
    
    Parameters
    ----------
    tad_file : str
        TAD文件路径
    chrom_sizes_file : str
        染色体大小文件路径
    output_file : str
        输出文件路径
    bin_size : int
        bin大小，默认10000
    """
    import logging
    logger = logging.getLogger(__name__)
    logger.info(f"正在处理: {tad_file}")
    
    tad_df = pd.read_csv(tad_file, sep='\t', header=0)
    tad_df['start'] = pd.to_numeric(tad_df['start'])
    tad_df['end'] = pd.to_numeric(tad_df['end'])
    
    chrom_sizes = {}
    if Path(chrom_sizes_file).exists():
        with open(chrom_sizes_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    chrom_sizes[parts[0]] = int(parts[1])
    
    chroms = sorted(tad_df['chrom'].unique())
    all_bins = []
    
    for chrom in chroms:
        if chrom not in chrom_sizes:
            continue
        chrom_size = chrom_sizes[chrom]
        tad_chrom = tad_df[tad_df['chrom'] == chrom].copy().sort_values('start').reset_index(drop=True)
        num_bins = int(np.ceil(chrom_size / bin_size))
        
        tad_bin_ranges = [(int(tad['start'] // bin_size), int((tad['end'] - 1) // bin_size)) for _, tad in tad_chrom.iterrows()]
        
        bins_with_annotations = []
        for i in range(num_bins):
            bin_start = i * bin_size
            bin_end = min((i + 1) * bin_size, chrom_size)
            is_tad = any(sb <= i <= eb for sb, eb in tad_bin_ranges)
            bins_with_annotations.append({'bin_idx': i, 'start': bin_start, 'end': bin_end, 'annotation': 'TAD' if is_tad else 'Boundary'})
        
        filtered_bins = []
        i = 0
        while i < len(bins_with_annotations):
            if bins_with_annotations[i]['annotation'] == 'TAD':
                filtered_bins.append(bins_with_annotations[i])
                i += 1
            else:
                start_idx = i
                while i < len(bins_with_annotations) and bins_with_annotations[i]['annotation'] == 'Boundary':
                    i += 1
                end_idx = i
                boundary_length = end_idx - start_idx
                
                if boundary_length == 1:
                    filtered_bins.append(bins_with_annotations[start_idx])
                elif boundary_length == 2:
                    filtered_bins.append(bins_with_annotations[start_idx])
                    filtered_bins.append(bins_with_annotations[start_idx + 1])
                else:
                    first_bin = bins_with_annotations[start_idx].copy()
                    filtered_bins.append(first_bin)
                    for j in range(start_idx + 1, end_idx - 1):
                        na_bin = bins_with_annotations[j].copy()
                        na_bin['annotation'] = 'NA'
                        filtered_bins.append(na_bin)
                    filtered_bins.append(bins_with_annotations[end_idx - 1])
        
        for bin_info in filtered_bins:
            all_bins.append({'chrom': chrom, 'start': bin_info['start'], 'end': bin_info['end'], 'annotation': bin_info['annotation']})
    
    result_df = pd.DataFrame(all_bins)
    result_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"  -> 生成 {len(result_df)} 个bins")


def calculate_boundary_context(boundaries: pd.DataFrame) -> pd.DataFrame:
    """
    计算边界上下文信息。
    Source: a_5 L226
    
    Parameters
    ----------
    boundaries : pd.DataFrame
        边界DataFrame
        
    Returns
    -------
    pd.DataFrame
        带上下游距离的边界DataFrame
    """
    boundaries = boundaries.copy().reset_index(drop=True)
    boundaries['start'] = boundaries['start'].astype(float)
    boundaries['end'] = boundaries['end'].astype(float)
    
    boundaries['dist_to_upstream'] = np.nan
    boundaries['dist_to_downstream'] = np.nan
    boundaries['adaptive_threshold'] = np.nan
    
    for chrom in boundaries['chrom'].unique():
        chrom_mask = boundaries['chrom'] == chrom
        chrom_boundaries = boundaries[chrom_mask].copy().reset_index(drop=True)
        chrom_boundaries = chrom_boundaries.sort_values('start').reset_index(drop=True)
        
        for i in range(len(chrom_boundaries)):
            current = chrom_boundaries.iloc[i]
            
            if i == 0:
                dist_up = np.inf
                dist_down = (chrom_boundaries.iloc[1]['start'] - current['end']) if len(chrom_boundaries) > 1 else np.inf
            elif i == len(chrom_boundaries) - 1:
                dist_down = np.inf
                dist_up = current['start'] - chrom_boundaries.iloc[i-1]['end']
            else:
                dist_up = current['start'] - chrom_boundaries.iloc[i-1]['end']
                dist_down = chrom_boundaries.iloc[i+1]['start'] - current['end']
            
            adaptive_threshold = min(abs(dist_up), abs(dist_down)) / 2 if (np.isfinite(dist_up) and np.isfinite(dist_down)) else np.inf
            
            idx = boundaries[chrom_mask].index[i]
            boundaries.loc[idx, 'dist_to_upstream'] = dist_up
            boundaries.loc[idx, 'dist_to_downstream'] = dist_down
            boundaries.loc[idx, 'adaptive_threshold'] = adaptive_threshold
    
    return boundaries


def generate_dual_direction_pairing(boundaries1: pd.DataFrame, boundaries2: pd.DataFrame,
                                     sample1_name: str, sample2_name: str,
                                     output_file1: str, output_file2: str,
                                     logger: logging.Logger) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    生成双向配对文件。
    Source: a_5 L340
    
    核心函数，使用KDTree双向配对。
    
    Parameters
    ----------
    boundaries1, boundaries2 : pd.DataFrame
        两个样本的边界数据
    sample1_name, sample2_name : str
        样本名称
    output_file1, output_file2 : str
        输出文件路径
    logger : logging.Logger
        日志记录器
        
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        双向配对结果
    """
    def pair_one_direction(b_from, b_to, from_name, to_name):
        paired_details = []
        for chrom in b_from['chrom'].unique():
            if chrom not in b_to['chrom'].unique():
                continue
            
            b_from_chrom = b_from[b_from['chrom'] == chrom].copy().sort_values('start').reset_index(drop=True)
            b_to_chrom = b_to[b_to['chrom'] == chrom].copy().sort_values('start').reset_index(drop=True)
            
            positions_to = b_to_chrom['start'].values.reshape(-1, 1)
            tree = KDTree(positions_to)
            
            for _, row_from in b_from_chrom.iterrows():
                pos = row_from['start']
                dist, idx = tree.query([pos])
                row_to = b_to_chrom.iloc[idx]
                
                distance = row_to['start'] - row_from['start']
                
                paired_details.append({
                    'chrom': chrom,
                    'sample1': from_name,
                    'sample1_start': row_from['start'],
                    'sample1_end': row_from['end'],
                    'sample1_dist_upstream': row_from.get('dist_to_upstream', np.nan),
                    'sample1_dist_downstream': row_from.get('dist_to_downstream', np.nan),
                    'sample1_adaptive_threshold': row_from.get('adaptive_threshold', np.nan),
                    'sample2': to_name,
                    'sample2_start': row_to['start'],
                    'sample2_end': row_to['end'],
                    'sample2_dist_upstream': row_to.get('dist_to_upstream', np.nan),
                    'sample2_dist_downstream': row_to.get('dist_to_downstream', np.nan),
                    'sample2_adaptive_threshold': row_to.get('adaptive_threshold', np.nan),
                    'distance': distance,
                    'direction': 'Downstream' if distance > 0 else ('Upstream' if distance < 0 else 'No_shift'),
                    'abs_distance': abs(distance),
                    'min_threshold': min(row_from.get('adaptive_threshold', np.inf), row_to.get('adaptive_threshold', np.inf))
                })
        
        return pd.DataFrame(paired_details)
    
    logger.info(f"    方向1: {sample1_name} -> {sample2_name}")
    paired_df1 = pair_one_direction(boundaries1, boundaries2, sample1_name, sample2_name)
    paired_df1.to_csv(output_file1, sep='\t', index=False)
    logger.info(f"      {sample1_name}视角: {len(paired_df1)}个配对")
    
    logger.info(f"    方向2: {sample2_name} -> {sample1_name}")
    paired_df2 = pair_one_direction(boundaries2, boundaries1, sample2_name, sample1_name)
    paired_df2.to_csv(output_file2, sep='\t', index=False)
    logger.info(f"      {sample2_name}视角: {len(paired_df2)}个配对")
    
    return paired_df1, paired_df2


def generate_final_classification(paired_df: pd.DataFrame, sample1_name: str, sample2_name: str,
                                   output_file: str, logger: logging.Logger) -> Tuple[pd.DataFrame, pd.Series]:
    """
    生成最终分类文件。
    Source: a_5 L401
    
    分类: Stable_boundary / Boundary_shift / Unique_boundary
    """
    final_classification = []
    for _, row in paired_df.iterrows():
        abs_distance = row['abs_distance']
        min_threshold = row['min_threshold']
        
        if abs_distance == 0:
            classification = 'Stable_boundary'
        elif abs_distance < min_threshold:
            classification = 'Boundary_shift'
        else:
            classification = 'Unique_boundary'
        
        final_classification.append(classification)
    
    paired_df_final = paired_df.copy()
    paired_df_final['final_classification'] = final_classification
    paired_df_final.to_csv(output_file, sep='\t', index=False)
    
    classification_counts = paired_df_final['final_classification'].value_counts()
    logger.info(f"      最终分类统计: {dict(classification_counts)}")
    
    return paired_df_final, classification_counts


def save_plot_data_tad(paired_df1_final: pd.DataFrame, paired_df2_final: pd.DataFrame,
                       treatment: str, control: str, output_path: Path, logger: logging.Logger) -> Path:
    """保存TAD绘图专用数据文件。Source: a_5 L430"""
    def get_stats(df, view_name):
        total = len(df)
        stats = {'view': view_name, 'total': total}
        for cat in ['Stable_boundary', 'Boundary_shift', 'Unique_boundary']:
            count = len(df[df['final_classification'] == cat])
            stats[cat] = count
            stats[f'{cat}_pct'] = (count / total * 100) if total > 0 else 0
        return stats
    
    stats1 = get_stats(paired_df1_final, treatment)
    stats2 = get_stats(paired_df2_final, control)
    
    plot_data = pd.DataFrame([stats1, stats2])
    plot_data_file = output_path / f"tad_plot_data.tsv"
    plot_data.to_csv(plot_data_file, sep='\t', index=False)
    logger.info(f"  -> 绘图数据: {plot_data_file}")
    
    return plot_data_file


def load_plot_data_tad(plot_data_file: Path, logger: logging.Logger) -> pd.DataFrame:
    """加载TAD绘图数据。Source: a_5 L454"""
    logger.info(f"读取绘图数据: {plot_data_file}")
    
    if not plot_data_file.exists():
        raise FileNotFoundError(f"绘图数据文件不存在: {plot_data_file}")
    
    df = pd.read_csv(plot_data_file, sep='\t')
    logger.info(f"  -> 数据行数: {len(df)}")
    
    return df


def plot_tad_stacked_bar(plot_data: pd.DataFrame, comparison: str, window_mult: int,
                          output_dir: Path, logger: logging.Logger) -> None:
    """
    绘制TAD边界分类堆积柱状图。
    Source: a_5 L467
    """
    logger.info(f"\n生成 {comparison} ({window_mult}b) 堆积柱状图...")
    
    width_cm, height_cm = 6.4, 5
    fig, ax = plt.subplots(figsize=(width_cm / 2.54, height_cm / 2.54))
    
    views = plot_data['view'].tolist()
    n_views = len(views)
    x_positions = np.arange(n_views)
    width = 0.8
    
    categories = ['Stable_boundary', 'Boundary_shift', 'Unique_boundary']
    category_labels = ['stable (unchanged)', 'unique (boundary shift)', 'unique (lost and gain)']
    
    bottom = np.zeros(n_views)
    
    for i, cat in enumerate(categories):
        pct_values = plot_data[f'{cat}_pct'].tolist()
        count_values = plot_data[cat].tolist()
        
        bars = ax.bar(x_positions, pct_values, width, bottom=bottom, label=category_labels[i], color=TAD_COLORS[cat])
        
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
    ax.set_xticklabels(views)
    ax.set_ylim(0, 100)
    ax.legend(title='boundary change', loc=(1.01, 0.6), frameon=False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    plt.subplots_adjust(top=0.96, bottom=0.16, right=0.75, left=0.125)
    
    comparison_clean = comparison.replace('--', '_')
    output_png = output_dir / f"tad_{comparison_clean}_{window_mult}b_stacked_bar.png"
    output_svg = output_dir / f"tad_{comparison_clean}_{window_mult}b_stacked_bar.svg"
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_svg, format='svg', bbox_inches='tight')
    plt.close()
    
    logger.info(f"  -> PNG: {output_png}")
    logger.info(f"  -> SVG: {output_svg}")


def generate_analysis_report_tad(paired_df1_final: pd.DataFrame, paired_df2_final: pd.DataFrame,
                                  treatment: str, control: str, comparison: str,
                                  output_path: Path, logger: logging.Logger) -> None:
    """
    生成TAD分析报告。
    Source: a_5 L523
    
    注意: 用户明确不写.md报告，此函数保留但可能不调用。
    """
    from datetime import datetime
    report_file = output_path / f"tad_{comparison}_analysis_report.txt"
    
    class_counts1 = paired_df1_final['final_classification'].value_counts()
    class_counts2 = paired_df2_final['final_classification'].value_counts()
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(f"{comparison} - TAD边界双向配对分析报告\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"方向1 ({treatment} -> {control}): {len(paired_df1_final)}个配对\n")
        f.write(f"方向2 ({control} -> {treatment}): {len(paired_df2_final)}个配对\n\n")
        f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    logger.info(f"  -> 分析报告: {report_file}")


def analyze_single_comparison_window(comparison: str, window_mult: int, output_root: str, run_mode: str,
                                     treatment_boundaries_path: str = None, control_boundaries_path: str = None,
                                     treatment_name: str = None, control_name: str = None) -> None:
    """
    分析单个比较组的单个窗口。
    Source: a_5 L559
    
    主流程函数。
    """
    def setup_logger(output_path, comparison):
        log_dir = output_path / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"tad_{comparison}_{window_mult}b.log"
        logger = logging.getLogger(f"tad_diff_{comparison}_{window_mult}b")
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
            return parts[0], parts[1]
        return comparison_str, comparison_str
    
    treatment, control = parse_comparison_func(comparison)
    actual_window = get_actual_window(window_mult)
    
    comparison_clean = comparison.replace('--', '_')
    output_path = Path(output_root) / f"{window_mult}b"
    output_path.mkdir(parents=True, exist_ok=True)
    
    logger = setup_logger(output_path, comparison)
    
    logger.info("=" * 70)
    logger.info(f"TAD边界差异分析: {comparison} ({window_mult}b={actual_window//1000}kb)")
    logger.info(f"运行模式: {run_mode}")
    logger.info("=" * 70)
    
    try:
        if run_mode in ['all', 'compute']:
            if treatment_boundaries_path is None or control_boundaries_path is None:
                logger.warning("请提供 treatment_boundaries_path 和 control_boundaries_path 参数")
                return
            
            logger.info(f"\n1. 加载边界数据")
            treatment_boundaries = pd.read_csv(treatment_boundaries_path, sep='\t')
            control_boundaries = pd.read_csv(control_boundaries_path, sep='\t')
            
            logger.info(f"  {treatment}: {len(treatment_boundaries)} 个边界")
            logger.info(f"  {control}: {len(control_boundaries)} 个边界")
            
            logger.info(f"\n2. 计算边界上下文")
            treatment_with_ctx = calculate_boundary_context(treatment_boundaries)
            control_with_ctx = calculate_boundary_context(control_boundaries)
            
            logger.info(f"\n3. 双向配对")
            output_file1 = output_path / f"{treatment}_vs_{control}_{window_mult}b_pairing1.tsv"
            output_file2 = output_path / f"{treatment}_vs_{control}_{window_mult}b_pairing2.tsv"
            
            paired_df1, paired_df2 = generate_dual_direction_pairing(
                treatment_with_ctx, control_with_ctx,
                treatment, control,
                str(output_file1), str(output_file2),
                logger
            )
            
            logger.info(f"\n4. 最终分类")
            final_file1 = output_path / f"{treatment}_vs_{control}_boundary_classification_final.tsv"
            paired_df1_final, _ = generate_final_classification(paired_df1, treatment, control, str(final_file1), logger)
            
            final_file2 = output_path / f"{control}_vs_{treatment}_boundary_classification_final.tsv"
            paired_df2_final, _ = generate_final_classification(paired_df2, control, treatment, str(final_file2), logger)
            
            logger.info(f"\n5. 保存绘图数据")
            save_plot_data_tad(paired_df1_final, paired_df2_final, treatment, control, output_path, logger)
        
        if run_mode in ['all', 'plot']:
            logger.info(f"\n6. 生成统计图")
            plot_data_file = output_path / "tad_plot_data.tsv"
            
            if not plot_data_file.exists():
                logger.error(f"绘图数据文件不存在: {plot_data_file}")
                return
            
            plot_data = load_plot_data_tad(plot_data_file, logger)
            plot_tad_stacked_bar(plot_data, comparison_clean, window_mult, output_path, logger)
        
        logger.info("\n" + "=" * 70)
        logger.info(f"分析完成！")
        logger.info(f"输出目录: {output_path}")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error(f"分析失败: {str(e)}")
        raise


# =============================================================================
# Complete pipeline
# =============================================================================

def process_tads(
    cooler_path: str,
    resolution: int,
    windows: List[int],
    output_dir: str,
    max_tad_length: int = 3_000_000,
    threshold_method: str = 'otsu',
    nproc: int = 1,
    verbose: bool = False
) -> Tuple[Dict[int, pd.DataFrame], Dict[int, pd.DataFrame]]:
    """
    Complete TAD analysis pipeline.

    This function combines all TAD analysis steps:
    1. Compute insulation scores
    2. Identify boundaries for each window
    3. Extract TADs for each window
    4. Save results to files

    Parameters
    ----------
    cooler_path : str
        Path to mcool file
    resolution : int
        Resolution in bp
    windows : list
        List of window sizes in bp (e.g., [50000, 100000, 500000])
    output_dir : str
        Output directory for results
    max_tad_length : int
        Maximum TAD length in bp
    threshold_method : str
        Threshold method for boundary identification: 'otsu' or 'li'
    nproc : int
        Number of parallel processes
    verbose : bool
        Print progress information

    Returns
    -------
    Tuple[Dict, Dict]
        - boundaries_dict: {window_size: boundaries_df}
        - tads_dict: {window_size: tads_df}

    Output files
    ------------
    The following files are saved to output_dir:
    - 1_0.{sample}.{resolution}.insulation.tsv           # Insulation scores
    - 2_0.{sample}.{resolution}.{window}kb.boundaries.tsv  # Boundaries
    - 2_1.{sample}.{resolution}.{window}kb.tads.tsv       # TADs

    Correspondence with hicviz:
    >>> from hicviz.core.tads import process_tads
    >>> boundaries, tads = process_tads(
    ...     mcool_path="sample.mcool",
    ...     resolution=10000,
    ...     windows=[50000, 100000, 500000],
    ...     output_dir="./output",
    ...     nproc=20
    ... )
    """
    from cfizz.io.insulation import compute_insulation_from_cooler

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Get sample name from file
    basename = os.path.basename(cooler_path).split('.')[0]

    if verbose:
        print(f"Processing TADs for {basename}")
        print(f"  Resolution: {resolution} bp")
        print(f"  Windows: {[w//1000 for w in windows]}kb")

    # Step 1: Compute insulation scores
    if verbose:
        print("  Computing insulation scores...")

    insulation_df = compute_insulation_from_cooler(
        cooler_path=cooler_path,
        resolution=resolution,
        windows=windows,
        nproc=nproc,
        verbose=verbose
    )

    # Save insulation scores
    insulation_file = os.path.join(
        output_dir,
        f"1_0.{basename}.{resolution}.insulation.tsv"
    )
    insulation_df.to_csv(insulation_file, index=False, sep='\t')

    if verbose:
        print(f"    Saved: {insulation_file}")

    # Step 2 & 3: Process each window
    boundaries_dict = {}
    tads_dict = {}

    for window_size in windows:
        window_str = f"{window_size // 1000}kb"

        if verbose:
            print(f"  Processing {window_str}...")

        # Identify boundaries
        boundaries = identify_boundaries(
            insulation_df,
            window_size,
            threshold_method
        )
        boundaries_dict[window_size] = boundaries

        # Extract TADs
        tads = extract_tads(
            insulation_df,
            window_size,
            max_tad_length
        )
        tads_dict[window_size] = tads

        # Save boundaries
        boundaries_file = os.path.join(
            output_dir,
            f"2_0.{basename}.{resolution}.{window_str}.boundaries.tsv"
        )
        boundaries.to_csv(boundaries_file, index=False, sep='\t')

        # Save TADs
        tads_file = os.path.join(
            output_dir,
            f"2_1.{basename}.{resolution}.{window_str}.tads.tsv"
        )
        tads.to_csv(tads_file, index=False, sep='\t')

        if verbose:
            print(f"    Boundaries: {len(boundaries)}, TADs: {len(tads)}")
            print(f"    Saved: {boundaries_file}")
            print(f"    Saved: {tads_file}")

    return boundaries_dict, tads_dict
