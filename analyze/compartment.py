"""
Compartment analysis module.

This module provides functions for A/B compartment analysis
using eigenvector decomposition of Hi-C matrices.

"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Any
from scipy.linalg import eigh
import warnings
import os

warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# =============================================================================

def calc_gc_cov(
    bins: pd.DataFrame,
    fasta_path: str,
    out_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Calculate GC content for genomic bins using bioframe.

    This function uses bioframe.frac_gc() to compute the GC content
    for each genomic bin in the provided bins DataFrame.

    Parameters
    ----------
    bins : pd.DataFrame
        DataFrame with columns: chrom, start, end (bin coordinates)
    fasta_path : str
        Path to reference genome FASTA file
    out_path : str, optional
        Path to save GC content results

    Returns
    -------
    gc_cov : pd.DataFrame
        DataFrame with columns: chrom, start, end, gc (GC content)

    Examples
    --------
    >>> clr = cooler.Cooler("sample.mcool::resolutions/100000")
    >>> bins = clr.bins()[:]
    >>> gc_cov = calc_gc_cov(bins, "hg38.fa")
    >>> print(f"GC range: {gc_cov['gc'].min():.3f} to {gc_cov['gc'].max():.3f}")
    """
    try:
        import bioframe
    except ImportError:
        raise ImportError("bioframe required: pip install bioframe")

    # Load genome FASTA
    hg_genome = bioframe.load_fasta(fasta_path)

    # Calculate GC content
    gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg_genome)

    # Save if path provided
    if out_path:
        gc_cov.to_csv(out_path, index=False, sep='\t')

    return gc_cov


def get_view_df(clr) -> pd.DataFrame:
    """
    Get chromosome view DataFrame from a cooler object.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object

    Returns
    -------
    view_df : pd.DataFrame
        DataFrame with columns: chrom, start, end, name
    """
    return pd.DataFrame({
        'chrom': clr.chromnames,
        'start': 0,
        'end': clr.chromsizes.values,
        'name': clr.chromnames
    })


# =============================================================================
# =============================================================================

def compute_eigenvector_cooltools(
    clr,
    gc_cov: pd.DataFrame,
    view_df: pd.DataFrame,
    n_eigs: int = 3
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Compute eigenvector using cooltools.eigs_cis.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object
    gc_cov : pd.DataFrame
        GC content DataFrame
    view_df : pd.DataFrame
        Chromosome view DataFrame
    n_eigs : int
        Number of eigenvectors to compute

    Returns
    -------
    eigenvalues : np.ndarray
        Eigenvalues for each chromosome
    eigenvector_df : pd.DataFrame
        DataFrame with columns: chrom, start, end, weight, E1, E2, ...
    """
    import cooltools

    # Compute eigenvectors using cooltools
    eigenvalues, eigenvector_df = cooltools.eigs_cis(
        clr,
        gc_cov,
        view_df=view_df,
        n_eigs=n_eigs
    )

    return eigenvalues, eigenvector_df


# =============================================================================
# Complete Process Functions
# =============================================================================

def process_compartment(
    mcool_path: str,
    resolution: int,
    fasta_path: Optional[str] = None,
    output_dir: str = ".",
    n_eigs: int = 3
) -> pd.DataFrame:
    """
    Process compartment calculation for Hi-C data.

    This is the main entry point for compartment analysis, corresponding
    to cfizz.analyze.compartment.process_compartment().

    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    resolution : int
        Resolution in bp (e.g., 1000000 for 1Mb)
    fasta_path : str, optional
        Path to reference genome FASTA file
    output_dir : str
        Output directory for results
    n_eigs : int
        Number of eigenvectors to compute

    Returns
    -------
    eig_df : pd.DataFrame
        Eigenvector results DataFrame with columns: chrom, start, end, E1, E2, ...

    Examples
    --------
    >>> eig_df = process_compartment(
    ...     "sample.mcool",
    ...     resolution=1000000,
    ...     fasta_path="hg38.fa",
    ...     output_dir="./output"
    ... )
    >>> print(f"E1 range: {eig_df['E1'].min():.3f} to {eig_df['E1'].max():.3f}")
    """
    import cooler

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Extract basename for file naming
    basename = os.path.basename(mcool_path).split('.')[0]

    print(f"\n[process_compartment] Computing compartment...")
    print(f"  mcool: {mcool_path}")
    print(f"  resolution: {resolution}")

    # Open cooler
    clr = cooler.Cooler(f'{mcool_path}::resolutions/{resolution}')
    bins = clr.bins()[:]

    # Calculate resolution string
    res_str = f"{resolution // 1000}k" if resolution < 1e6 else f"{resolution // 1e6}M"
    gc_cov_path = os.path.join(output_dir, f'gc_cov.{res_str}.tsv')

    # Calculate or load GC content
    if fasta_path is None:
        # Use default GC content (0.5)
        print(f"  No fasta provided, using default GC=0.5")
        gc_cov = pd.DataFrame({
            'chrom': bins['chrom'],
            'start': bins['start'],
            'end': bins['end'],
            'gc': 0.5
        })
    elif os.path.exists(gc_cov_path):
        print(f"  Loading existing GC content: {gc_cov_path}")
        gc_cov = pd.read_csv(gc_cov_path, sep='\t')
    else:
        print(f"  Calculating GC content with bioframe.frac_gc...")
        gc_cov = calc_gc_cov(bins, fasta_path, gc_cov_path)
        print(f"  GC content saved to: {gc_cov_path}")

    # Get chromosome view
    view_df = get_view_df(clr)

    # Compute eigenvectors
    print(f"  Computing eigenvectors (n_eigs={n_eigs})...")
    eigenvalues, eigenvector_df = compute_eigenvector_cooltools(
        clr, gc_cov, view_df, n_eigs=n_eigs
    )

    # Save results
    eig_path = os.path.join(output_dir, f'eigenvector.{res_str}.tsv')
    eigenvector_df.to_csv(eig_path, sep='\t', index=False)
    print(f"  Eigenvector saved to: {eig_path}")

    return eigenvector_df


# =============================================================================
# Original functions (kept for reference/backward compatibility)
# =============================================================================

def compute_compartment(
    oe_matrix: np.ndarray,
    n_components: int = 1
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute A/B compartment using eigenvector decomposition (custom method).

    The first eigenvector (E1) reflects the chessboard pattern of Hi-C data,
    where positive values indicate A compartments (open chromatin) and
    negative values indicate B compartments (closed chromatin).

    Note: This is a custom implementation. For cooler data, use
    process_compartment() which uses cooltools.eigs_cis for better results.
    """
    # Replace NaN values with 1 (neutral O/E value)
    oe_matrix = np.nan_to_num(oe_matrix, nan=1.0)

    # Ensure symmetric matrix
    if not np.allclose(oe_matrix, oe_matrix.T):
        oe_matrix = (oe_matrix + oe_matrix.T) / 2

    # Calculate correlation matrix
    n = oe_matrix.shape[0]
    corr_matrix = np.ones((n, n))

    for i in range(n):
        for j in range(n):
            if i != j:
                corr_matrix[i, j] = oe_matrix[i, j] * oe_matrix[j, i]

    # Replace diagonal with average of off-diagonal
    mask = ~np.eye(n, dtype=bool)
    valid_values = oe_matrix[mask]
    diag_value = np.nanmean(valid_values) if len(valid_values) > 0 else 1.0
    np.fill_diagonal(corr_matrix, diag_value)

    # Handle any remaining NaN/Inf values
    corr_matrix = np.nan_to_num(corr_matrix, nan=1.0, posinf=1.0, neginf=-1.0)

    # Eigenvalue decomposition
    eigenvalues, eigenvectors = eigh(corr_matrix)

    # Sort by eigenvalue (descending)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return eigenvectors[:, :n_components], eigenvalues[:n_components]


def compute_eigenvector(
    matrix: np.ndarray,
    n_components: int = 1,
    method: str = "classic"
) -> pd.DataFrame:
    """
    Compute eigenvectors for compartment analysis.
    """
    if method == "classic":
        eigenvectors, eigenvalues = compute_compartment(matrix, n_components)

    n = matrix.shape[0]
    data = {
        'start': np.arange(n),
        'end': np.arange(n) + 1
    }

    for i in range(n_components):
        data[f'E{i+1}'] = eigenvectors[:, i]


def assign_compartments(
    eigenvector: np.ndarray,
    gc_content: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, str]:
    """
    Assign A/B compartment labels based on eigenvector.
    """
    if gc_content is not None:
        correlation = np.corrcoef(eigenvector, gc_content)[0, 1]
        if correlation < 0:
            eigenvector = -eigenvector
            orientation = "inverted"
        else:
            orientation = "direct"
    else:
        orientation = "unknown"

    compartment_labels = np.where(eigenvector > 0, 1, -1)
    return compartment_labels, orientation


def calculate_compartment_strength(
    matrix: np.ndarray,
    compartment_labels: np.ndarray,
    n_bins: int = 10
) -> Dict[str, float]:
    """
    Calculate compartment strength metrics.
    """
    n = matrix.shape[0]

    within_a = within_b = between = 0
    count_a = count_b = count_between = 0

    for i in range(n):
        for j in range(i + 1, n):
            if compartment_labels[i] == 1 and compartment_labels[j] == 1:
                within_a += matrix[i, j]
                count_a += 1
            elif compartment_labels[i] == -1 and compartment_labels[j] == -1:
                within_b += matrix[i, j]
                count_b += 1
            else:
                between += matrix[i, j]
                count_between += 1

    avg_within_a = within_a / count_a if count_a > 0 else 0
    avg_within_b = within_b / count_b if count_b > 0 else 0
    avg_between = between / count_between if count_between > 0 else 0

    return {
        'avg_within_A': avg_within_a,
        'avg_within_B': avg_within_b,
        'avg_between': avg_between,
        'strength_A': avg_within_a / avg_between if avg_between > 0 else 0,
        'strength_B': avg_within_b / avg_between if avg_between > 0 else 0
    }


# =============================================================================
# =============================================================================

def calculate_saddle_matrix(
    matrix: np.ndarray,
    track: pd.DataFrame,
    n_bins: int = 98,
    qrange: Tuple[float, float] = (0.05, 0.95)
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate saddle score matrix using digitization and aggregation.

    Parameters
    ----------
    matrix : np.ndarray
        O/E normalized contact matrix
    track : pd.DataFrame
        DataFrame with compartment scores. Must have 'chrom', 'start', 'end', and E1 column.
    n_bins : int
        Number of bins for saddle
    qrange : tuple
        Quantile range for digitization

    Returns
    -------
    saddledata : np.ndarray
        Saddle matrix (n_bins x n_bins)
    track_values : np.ndarray
        Original track values (filtered)
    binedges : np.ndarray
        Bin edges for the track values
    """
    track_values = track['E1'].values
    n = min(len(track_values), matrix.shape[0])

    track_values = track_values[:n]
    matrix_cropped = matrix[:n, :n]

    binedges = np.percentile(track_values, np.linspace(0, 100, n_bins + 1))
    digitized = np.digitize(track_values, binedges[1:-1])

    saddledata = np.zeros((n_bins, n_bins))
    counts = np.zeros((n_bins, n_bins))

    for i in range(n):
        for j in range(n):
            di = digitized[i]
            dj = digitized[j]
            if 0 <= di < n_bins and 0 <= dj < n_bins:
                saddledata[di, dj] += matrix_cropped[i, j]
                counts[di, dj] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        saddledata = np.where(counts > 0, saddledata / counts, np.nan)

    return saddledata, track_values, binedges


def saddle_score(
    matrix: np.ndarray,
    track: pd.DataFrame,
    n_bins: int = 98,
    qrange: Tuple[float, float] = (0.075, 0.975)
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate saddle score matrix (legacy function).
    """
    return calculate_saddle_matrix(matrix, track, n_bins, qrange)


# =============================================================================
# Complete Analysis Pipeline
# =============================================================================

def calculate_compartment_analysis(
    mcool_path: str,
    chrom: str,
    resolution: int,
    n_bins_saddle: int = 10,
    fasta_path: Optional[str] = None,
    output_dir: str = "."
) -> Dict[str, Any]:
    """
    Perform complete compartment analysis pipeline.

    This function orchestrates the full compartment analysis workflow:
    1. Compute E1 eigenvector using cooltools
    2. Calculate saddle matrix
    3. Calculate compartment statistics

    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    chrom : str
        Chromosome name (e.g., "chr1")
    resolution : int
        Resolution in bp
    n_bins_saddle : int
        Number of bins for saddle plot
    fasta_path : str, optional
        Path to reference genome FASTA file
    output_dir : str
        Output directory

    Returns
    -------
    results : dict
        Dictionary containing:
        - 'eigenvector_df': E1 eigenvector DataFrame
        - 'saddle_matrix': Saddle matrix
        - 'saddle_stats': Saddle statistics
        - 'compartment_stats': A/B compartment metrics

    Examples
    --------
    >>> results = calculate_compartment_analysis(
    ...     "sample.mcool", chrom="chr1", resolution=1000000
    ... )
    >>> print(f"A compartment bins: {(results['eigenvector_df']['E1'] > 0).sum()}")
    """
    from cfizz.analyze.oe import calculate_oe_matrix_cooltools

    # Step 1: Calculate O/E matrix
    print(f"[calculate_compartment_analysis] Step 1: O/E calculation...")
    OE, _ = calculate_oe_matrix_cooltools(
        mcool_path=mcool_path,
        chrom=chrom,
        resolution=resolution,
        balance=True,
        smooth=True
    )

    # Step 2: Compute eigenvector using process_compartment
    print(f"[calculate_compartment_analysis] Step 2: Eigenvector computation...")
    eig_df = process_compartment(
        mcool_path=mcool_path,
        resolution=resolution,
        fasta_path=fasta_path,
        output_dir=output_dir,
        n_eigs=2
    )

    # Filter for our chromosome
    eig_df = eig_df[eig_df['chrom'] == chrom].copy()

    # Step 3: Calculate saddle matrix
    print(f"[calculate_compartment_analysis] Step 3: Saddle matrix calculation...")
    saddle_matrix, _, _ = calculate_saddle_matrix(
        matrix=OE,
        track=eig_df,
        n_bins=n_bins_saddle
    )

    # Step 4: Calculate statistics
    E1 = eig_df['E1'].values[:OE.shape[0]]
    n_a = np.sum(E1 > 0)
    n_b = np.sum(E1 <= 0)

    E1_pos = E1 > 0
    E1_neg = E1 <= 0

    strength_a = np.nanmean(OE[np.ix_(E1_pos, E1_pos)]) if np.sum(E1_pos) > 0 else np.nan
    strength_b = np.nanmean(OE[np.ix_(E1_neg, E1_neg)]) if np.sum(E1_neg) > 0 else np.nan

    results = {
        'eigenvector_df': eig_df,
        'saddle_matrix': saddle_matrix,
        'saddle_stats': {
            'n_bins': n_bins_saddle,
            'mean_aa': np.nanmean(saddle_matrix[:n_bins_saddle//2, :n_bins_saddle//2]),
            'mean_bb': np.nanmean(saddle_matrix[n_bins_saddle//2:, n_bins_saddle//2:])
        },
        'compartment_stats': {
            'n_A_bins': n_a,
            'n_B_bins': n_b,
            'pct_A': 100 * n_a / len(E1) if len(E1) > 0 else 0,
            'pct_B': 100 * n_b / len(E1) if len(E1) > 0 else 0,
            'strength_A': strength_a,
            'strength_B': strength_b
        }
    }

    return results


# =============================================================================
# Differential Analysis Functions (移植自 script/a_4_compartment_diff.py)
# Source: a_4 L48, L143, L188, L240, L258, L271, L368
# =============================================================================

import logging
from pathlib import Path
from typing import List, Optional
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math

# 配置字体和绘图参数 (a_4 L61-65)
FONT_CONFIG = {
    'family': 'Arial',
    'weight': 'regular',
    'size': 5,
}

CM = 1 / 2.54  # cm转inch (a_4 L48)

# 配色方案 (a_4 L69)
COLOR_PALETTE = ['#264653', '#27736F', '#299D92', '#8AB17C', '#E8C56B', '#F2A361', '#E66F51']

# Compartment状态颜色映射 (a_4 L72-77)
COMPARTMENT_COLORS = {
    'Stable_A': COLOR_PALETTE[0],   # 深青色 - A compartment稳定
    'Stable_B': COLOR_PALETTE[1],   # 青色 - B compartment稳定
    'A_to_B': COLOR_PALETTE[5],     # 橙色 - A转B
    'B_to_A': COLOR_PALETTE[6],     # 红橙色 - B转A
}

# 兼容性别名 (旧代码用 COLORS_DICT)
COLORS_DICT = COMPARTMENT_COLORS

def max_to_top(max_value: float) -> float:
    """计算绘图时的最大值上限 (a_4 L113-121)"""
    all_weishu = len(str(int(max_value)))
    top_value = math.ceil(max_value / (10 ** (all_weishu - 1))) * (10 ** (all_weishu - 1))
    if top_value < max_value * 1.2:
        return top_value
    else:
        top_value = math.ceil(max_value / (10 ** (all_weishu - 2))) * (10 ** (all_weishu - 2))
        return top_value


def parse_comparison_func(comparison_str: str) -> Tuple[str, str]:
    """解析比较组字符串，返回 (treatment, control)"""
    if '--' in comparison_str:
        parts = comparison_str.split('--')
        return parts[0], parts[1]
    return comparison_str, comparison_str


def load_compartment_data(file_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    加载compartment数据文件。
    Source: a_4 L143
    
    Parameters
    ----------
    file_path : str
        E1 TSV文件路径，列: chrom, start, end, E1
    logger : logging.Logger
        日志记录器
        
    Returns
    -------
    pd.DataFrame
        清理后的DataFrame，只包含 chrom, start, end, E1 列
        
    行为特征 6 问:
    1. 输入: TSV文件，列名可能有大小写差异
    2. 输出: DataFrame [chrom, start, end, E1]
    3. 跟5_1产物关系: 读eigenvector.100k.tsv，不重算
    4. 跟a_4算法关系: 完全一致，L143-185
    5. 调谁: pd.read_csv, dropna
    6. 错误处理: 列名映射，dropna处理
    """
    logger.info(f"正在读取文件: {file_path}")
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        expected_cols = {'chrom', 'start', 'end', 'E1'}
        actual_cols = set(df.columns)
        
        if not expected_cols.issubset(actual_cols):
            logger.warning(f"列名不匹配。期望: {expected_cols}, 实际: {actual_cols}")
            col_mapping = {}
            for col in actual_cols:
                col_lower = col.lower()
                if 'chr' in col_lower or 'chrom' in col_lower:
                    col_mapping['chrom'] = col
                elif 'start' in col_lower:
                    col_mapping['start'] = col
                elif 'end' in col_lower:
                    col_mapping['end'] = col
                elif 'e1' in col_lower:
                    col_mapping['E1'] = col
            
            if col_mapping:
                df = df.rename(columns=col_mapping)
        
        df_clean = df.dropna(subset=['E1']).copy()
        
        for col in ['start', 'end']:
            if col in df_clean.columns:
                df_clean.loc[:, col] = pd.to_numeric(df_clean[col], errors='coerce')
        
        df_clean = df_clean.dropna(subset=['start', 'end']).copy()
        
        logger.info(f"  -> 有效数据行数: {len(df_clean)}")
        logger.info(f"  -> 染色体数量: {df_clean['chrom'].nunique()}")
        logger.info(f"  -> E1值范围: [{df_clean['E1'].min():.3f}, {df_clean['E1'].max():.3f}]")
        
        return df_clean[['chrom', 'start', 'end', 'E1']].copy()
        
    except Exception as e:
        logger.error(f"读取文件失败 {file_path}: {str(e)}")
        raise


def calculate_compartment_diff(
    df_control: pd.DataFrame, 
    df_treatment: pd.DataFrame,
    treatment_name: str, 
    logger: logging.Logger
) -> pd.DataFrame:
    """
    计算两个样本间的E1差异，逐bin AB状态对比。
    Source: a_4 L188
    
    核心函数，内含 classify_compartment_change 内部函数。
    分类4类: Stable_A / Stable_B / A_to_B / B_to_A
    
    Parameters
    ----------
    df_control : pd.DataFrame
        对照组E1数据
    df_treatment : pd.DataFrame
        处理组E1数据
    treatment_name : str
        处理组名称
    logger : logging.Logger
        日志记录器
        
    Returns
    -------
    pd.DataFrame
        差异数据，包含 E1_control, E1_{treatment_name}, compartment_change_{treatment_name}
        
    行为特征 6 问:
    1. 输入: 两个E1 DataFrame，通过chrom/start/end合并
    2. 输出: 合并后的DataFrame，新增dE1和compartment_change列
    3. 跟5_1产物关系: 读两个eigenvector.tsv，不重算
    4. 跟a_4算法关系: 完全一致，L188-237
    5. 调谁: pd.merge, numpy操作
    6. 错误处理: 空DataFrame时返回空DataFrame
    """
    logger.info(f"\n计算 {treatment_name} 与 control 的差异...")
    
    merged = pd.merge(df_control, df_treatment, on=['chrom', 'start', 'end'],
                      suffixes=('_control', f'_{treatment_name}'),
                      how='inner')
    
    if len(merged) == 0:
        logger.warning(f"  -> 警告: 没有找到匹配的坐标区域")
        return pd.DataFrame()
    
    logger.info(f"  -> 匹配区域数: {len(merged)}")
    
    e1_control = merged['E1_control']
    e1_sample = merged[f'E1_{treatment_name}']
    
    dE1_col = f'dE1_{treatment_name}_vs_control'
    merged[dE1_col] = e1_sample - e1_control
    
    abs_dE1_col = f'abs_dE1_{treatment_name}_vs_control'
    rel_dE1_col = f'rel_dE1_{treatment_name}_vs_control'
    
    merged[abs_dE1_col] = np.abs(merged[dE1_col])
    
    merged[rel_dE1_col] = np.where(
        np.abs(e1_control) > 1e-10,
        np.abs(merged[dE1_col]) / np.abs(e1_control),
        0
    )
    
    def classify_compartment_change(e1_ctrl, e1_sample):
        """内部函数: 分类E1变化"""
        if e1_ctrl > 0 and e1_sample > 0:
            return 'Stable_A'
        elif e1_ctrl < 0 and e1_sample < 0:
            return 'Stable_B'
        elif e1_ctrl > 0 and e1_sample < 0:
            return 'A_to_B'
        elif e1_ctrl < 0 and e1_sample > 0:
            return 'B_to_A'
        else:
            return 'Boundary'
    
    merged[f'compartment_change_{treatment_name}'] = merged.apply(
        lambda row: classify_compartment_change(row['E1_control'], row[f'E1_{treatment_name}']),
        axis=1
    )
    
    return merged


def save_plot_data(
    df_diff: pd.DataFrame, 
    treatment_name: str, 
    output_path: Path, 
    logger: logging.Logger
) -> Path:
    """
    保存绘图专用数据文件。
    Source: a_4 L240
    
    Parameters
    ----------
    df_diff : pd.DataFrame
        差异数据
    treatment_name : str
        处理组名称
    output_path : Path
        输出目录
    logger : logging.Logger
        日志记录器
        
    Returns
    -------
    Path
        保存的绘图数据文件路径
    """
    plot_columns = [
        'chrom', 'start', 'end',
        'E1_control',
        f'E1_{treatment_name}',
        f'compartment_change_{treatment_name}'
    ]
    
    df_plot = df_diff[plot_columns].copy()
    
    plot_data_file = output_path / f"compartment_plot_data.tsv"
    df_plot.to_csv(plot_data_file, sep='\t', index=False)
    logger.info(f"  -> 绘图数据: {plot_data_file}")
    
    return plot_data_file


def load_plot_data(plot_data_file: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    加载绘图专用数据文件。
    Source: a_4 L258
    
    Parameters
    ----------
    plot_data_file : Path
        绘图数据文件路径
    logger : logging.Logger
        日志记录器
        
    Returns
    -------
    pd.DataFrame
        绘图数据
    """
    logger.info(f"读取绘图数据: {plot_data_file}")
    
    if not plot_data_file.exists():
        raise FileNotFoundError(f"绘图数据文件不存在: {plot_data_file}\n请先运行 'compute' 或 'all' 模式生成数据")
    
    df = pd.read_csv(plot_data_file, sep='\t')
    logger.info(f"  -> 数据行数: {len(df)}")
    
    return df


def plot_compartment_scatter(df: pd.DataFrame, treatment_name: str, comparison: str,
                             output_dir: Path, logger: logging.Logger) -> None:
    """绘制compartment差异散点图（带边缘分布） (a_4 L271-359)"""
    logger.info(f"\n生成 {comparison} 散点图...")

    figure = plt.figure(figsize=(8 * CM, 6.4 * CM), dpi=600)
    figure.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.95, wspace=0.1, hspace=0.1)

    grid = plt.GridSpec(6, 6, wspace=0.05, hspace=0.05, figure=figure)
    ax1 = figure.add_subplot(grid[0, 0:5])
    ax2 = figure.add_subplot(grid[1:6, 5])
    ax3 = figure.add_subplot(grid[1:6, 0:5])

    e1_control = df['E1_control']
    e1_sample = df[f'E1_{treatment_name}']
    max_val = max(e1_control.max(), e1_sample.max())
    top_value = max_to_top(max_val)
    bottom_value = top_value * 0.1 * -1
    bottom_value = bottom_value if bottom_value < -1.5 else -1.5

    sns.kdeplot(data=df, x='E1_control', color=COLOR_PALETTE[2], fill=True,
                common_norm=False, legend=False, alpha=0.5, linewidth=0.5, cut=0, ax=ax1)
    sns.kdeplot(data=df, y=f'E1_{treatment_name}', color=COLOR_PALETTE[5], fill=True,
                common_norm=False, legend=False, alpha=0.5, linewidth=0.5, cut=0, ax=ax2)

    x_line = np.linspace(bottom_value, top_value, 1000)
    plt.plot(x_line, x_line, color="#bbbbbb", linewidth=0.1)

    change_col = f'compartment_change_{treatment_name}'
    for status, color in COMPARTMENT_COLORS.items():
        subset = df[df[change_col] == status]
        if len(subset) > 0:
            ax3.scatter(x=subset['E1_control'], y=subset[f'E1_{treatment_name}'],
                       s=1, alpha=0.5, color=color, edgecolors=color, linewidths=0.5)

    total = len(df)

    legend_elements = []
    for status, color in COMPARTMENT_COLORS.items():
        count = len(df[df[change_col] == status])
        if count > 0:
            pct = count / total * 100
            legend_elements.append(mlines.Line2D([0], [0], marker='o', color='w',
                                              markerfacecolor=color, markersize=4,
                                              label=f'{status}\nn={count:,}\n({pct:.1f}%)\n'))

    if legend_elements:
        ax2.legend(handles=legend_elements, prop=FONT_CONFIG, labelspacing=0.4,
                  handleheight=1.5, handletextpad=0.2, loc=(1.01,0.3),
                  frameon=False, title='')

    ax1.spines[:].set_linewidth(0)
    ax1.tick_params(width=0.6, length=2.5, labelsize=6)
    ax1.set_xticks([])
    ax1.set_xlabel("")
    ax1.set_yticks([])
    ax1.set_ylabel("")

    ax2.spines[:].set_linewidth(0)
    ax2.tick_params(width=0.6, length=2.5, labelsize=6)
    ax2.set_xticks([])
    ax2.set_xlabel("")
    ax2.set_yticks([])
    ax2.set_ylabel("")

    ax3.spines[:].set_linewidth(0.4)
    ax3.tick_params(width=0.6, length=2.5, labelsize=7)

    _, control_name = parse_comparison_func(comparison)
    ax3.set_xlabel(f'Compartment of {control_name}', fontsize=7, x=0.55)
    ax3.set_ylabel(f'Compartment of {treatment_name}', fontsize=7, y=0.55)

    ax1.set_xlim(bottom_value, top_value)
    ax2.set_ylim(bottom_value, top_value)
    ax3.set_xlim(bottom_value, top_value)
    ax3.set_ylim(bottom_value, top_value)

    ax3.axhline(y=0, color='gray', linestyle='--', linewidth=0.3, alpha=0.5)
    ax3.axvline(x=0, color='gray', linestyle='--', linewidth=0.3, alpha=0.5)

    output_svg = output_dir / f"compartment_{comparison}_scatter.svg"
    output_png = output_dir / f"compartment_{comparison}_scatter.png"

    figure.savefig(output_svg, format='svg')
    figure.savefig(output_png, format='png', dpi=600)
    plt.close()

    logger.info(f"  -> 散点图已保存: {output_svg}")
    logger.info(f"  -> 散点图已保存: {output_png}")


def merge_diff_regions(
    diff_df: pd.DataFrame,
    treatment_name: str,
    change_type: str = 'A_to_B',
    min_region_length: int = 3,
    resolution: int = 100000
) -> pd.DataFrame:
    """
    合并连续差异区域。
    Source: a_8 L182-225 (注释说明)
    
    用于 A_to_B / B_to_A 找连续bin合并成region块。
    
    Parameters
    ----------
    diff_df : pd.DataFrame
        差异数据DataFrame
    treatment_name : str
        处理组名称
    change_type : str
        变化类型 ('A_to_B' 或 'B_to_A')
    min_region_length : int
        最小区域长度(bin数)
    resolution : int
        分辨率(bp)，默认100kb
        
    Returns
    -------
    pd.DataFrame
        合并后的差异区域
    """
    change_col = f'compartment_change_{treatment_name}'
    
    filtered = diff_df[diff_df[change_col] == change_type].copy()
    
    if len(filtered) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'n_bins', 'change_type'])
    
    filtered = filtered.sort_values(['chrom', 'start']).reset_index(drop=True)
    
    regions = []
    current_chrom = None
    current_start = None
    current_end = None
    current_count = 0
    
    for _, row in filtered.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        
        if current_chrom is None:
            current_chrom = chrom
            current_start = start
            current_end = end
            current_count = 1
        elif chrom != current_chrom or start > current_end + resolution:
            if current_count >= min_region_length:
                regions.append({
                    'chrom': current_chrom,
                    'start': current_start,
                    'end': current_end,
                    'n_bins': current_count,
                    'change_type': change_type
                })
            current_chrom = chrom
            current_start = start
            current_end = end
            current_count = 1
        else:
            current_end = end
            current_count += 1
    
    if current_count >= min_region_length:
        regions.append({
            'chrom': current_chrom,
            'start': current_start,
            'end': current_end,
            'n_bins': current_count,
            'change_type': change_type
        })
    
    return pd.DataFrame(regions)


def analyze_single_comparison(
    comparison: str, 
    output_root: str, 
    run_mode: str,
    control_e1_path: str = None,
    treatment_e1_path: str = None
):
    """
    分析单个比较组的Compartment差异。
    Source: a_4 L368
    
    主流程函数，整合上面所有函数，支持 run_mode='all'/'compute'/'plot'。
    
    Parameters
    ----------
    comparison : str
        比较组字符串，格式 "处理组--对照组"
    output_root : str
        输出根目录
    run_mode : str
        运行模式: 'all', 'compute', 'plot'
    control_e1_path : str, optional
        对照组E1文件路径，如果不提供则使用默认推导
    treatment_e1_path : str, optional
        处理组E1文件路径，如果不提供则使用默认推导
        
    行为特征 6 问:
    1. 输入: comparison字符串，E1文件路径
    2. 输出: diff.tsv, plot_data.tsv, scatter图
    3. 跟5_1产物关系: 读eigenvector.tsv，不重算
    4. 跟a_4算法关系: 完全一致，L368-420
    5. 调谁: load_compartment_data, calculate_compartment_diff等
    6. 错误处理: 空数据时打印错误并返回
    """
    def setup_logger(output_path, comparison):
        """设置日志记录器"""
        log_dir = output_path / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = log_dir / f"{comparison}.log"
        logger = logging.getLogger(f"compartment_diff_{comparison}")
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
    
    treatment, control = parse_comparison_func(comparison)
    
    output_path = Path(output_root) / "compartment"
    output_path.mkdir(parents=True, exist_ok=True)
    
    logger = setup_logger(output_path, comparison)
    
    logger.info("=" * 70)
    logger.info(f"Compartment差异分析: {comparison}")
    logger.info(f"处理组: {treatment}")
    logger.info(f"对照组: {control}")
    logger.info(f"运行模式: {run_mode}")
    logger.info("=" * 70)
    
    try:
        if run_mode in ['all', 'compute']:
            logger.info("\n1. 加载compartment数据")
            logger.info("-" * 70)
            
            if control_e1_path is None or treatment_e1_path is None:
                logger.warning("请提供 control_e1_path 和 treatment_e1_path 参数")
                return
            
            logger.info(f"  对照组({control})文件: {control_e1_path}")
            logger.info(f"  处理组({treatment})文件: {treatment_e1_path}")
            
            df_control = load_compartment_data(control_e1_path, logger)
            df_treatment = load_compartment_data(treatment_e1_path, logger)
            
            logger.info("\n2. 计算compartment差异")
            logger.info("-" * 70)
            
            df_diff = calculate_compartment_diff(df_control, df_treatment, treatment, logger)
            
            if len(df_diff) == 0:
                logger.error("没有有效的差异数据，跳后续分析")
                return
            
            logger.info("\n3. 保存差异数据")
            logger.info("-" * 70)
            
            diff_output = output_path / f"compartment_{comparison}_diff.tsv"
            df_diff.to_csv(diff_output, sep='\t', index=False)
            logger.info(f"  -> 完整差异数据: {diff_output}")
            
            save_plot_data(df_diff, treatment, output_path, logger)
        
        if run_mode in ['all', 'plot']:
            logger.info("\n4. 生成统计图")
            logger.info("-" * 70)
            
            plot_data_file = output_path / "compartment_plot_data.tsv"
            
            if not plot_data_file.exists():
                logger.error(f"绘图数据文件不存在: {plot_data_file}")
                return
            
            df_plot = load_plot_data(plot_data_file, logger)
            plot_compartment_scatter(df_plot, treatment, comparison, output_path, logger)
        
        logger.info("\n" + "=" * 70)
        logger.info(f"分析完成！")
        logger.info(f"输出目录: {output_path}")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error(f"分析失败: {str(e)}")
        raise
