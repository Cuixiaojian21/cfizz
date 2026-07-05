"""
Aggregate/Pileup visualization module.

Functions for creating aggregate plots:
- TAD boundary pileup
- Loop anchor pileup
- Generic pileup analysis

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from typing import List, Tuple, Dict, Optional, Any
from multiprocessing import Pool, cpu_count

from .layout import (
    setup_plot_style,
    setup_ratio_colormap, 
    calculate_heatmap_layout,
    save_figure_multi_format
)
from ..io.cooler import read_cooler


# =============================================================================
# =============================================================================

def setup_axes(ax):
    """
    设置坐标轴和边框 - 隐藏所有刻度和边框
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        matplotlib 的 axes 对象
        
    Returns
    -------
    tuple : (xmin, xmax, ymin, ymax) 坐标轴范围
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.tick_params(
        axis='both',
        bottom=False, top=False, left=False, right=False,
        labelbottom=False, labeltop=False, labelleft=False, labelright=False
    )
    # 无边框版本
    for spine in ['right', 'top', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(0)
    
    return xmin, xmax, ymin, ymax


def get_matrix_range(matrix, vmin=None, vmax=None):
    """
    获取矩阵的值范围，确保使用非零值
    
    本地实现(基于通用 numpy 数值范围算法)
    
    Parameters
    ----------
    matrix : np.ndarray
        输入矩阵
    vmin : float, optional
        指定的最小值
    vmax : float, optional
        指定的最大值
        
    Returns
    -------
    tuple : (vmin, vmax)
    """
    if vmin is None or vmax is None:
        # 获取所有非零值,然后排除 NaN(np.nonzero 不排除 NaN,会污染 min/percentile)
        nonzero = matrix[np.nonzero(matrix)]
        valid = nonzero[~np.isnan(nonzero)]
        
        if vmin is None:
            # 使用有效值的最小值
            vmin = np.min(valid) if valid.size > 0 else 0
            
        if vmax is None:
            if valid.size > 0:
                # 使用93百分位数
                perc = np.percentile(valid, 93)
                # 找到最接近perc的实际值
                vmax = valid[np.abs(valid - perc).argmin()]
            else:
                vmax = 1
                
    return vmin, vmax


def extract_snippets_from_npy(
    npy_path: str,
    boundaries: pd.DataFrame,
    flank: int,
    resolution: int,
    start_pos: int = 0
) -> np.ndarray:
    """
    从 npy 矩阵文件提取子矩阵用于 pileup 分析。
    
    Parameters
    ----------
    npy_path : str
        npy 文件路径（OE 矩阵）
    boundaries : pd.DataFrame
        DataFrame with columns: chrom, start, end
    flank : int
        侧翼大小（bp），如 300000 表示 300kb
    resolution : int
        分辨率（bp），如 10000 表示 10kb
    start_pos : int
        矩阵起始位置（bp），默认为 0
        
    Returns
    -------
    stack : np.ndarray
        子矩阵堆叠，形状为 (2*window, 2*window, n_boundaries)
        
    Notes
    -----
    window 大小（bin 数）：window = flank // resolution
    期望矩阵大小：(2*window, 2*window)
        
    Example
    -------
    >>> boundaries = pd.read_csv('boundaries.tsv', sep='\\t')
    >>> stack = extract_snippets_from_npy(
    ...     '50_1_chr2_10k.oe.npy', 
    ...     boundaries, 
    ...     300000, 
    ...     10000
    ... )
    >>> print(stack.shape)  # (60, 60, n_boundaries)
    """
    # 加载矩阵
    oe_matrix = np.load(npy_path)
    window = flank // resolution
    expected_size = 2 * window
    
    print(f"从 npy 提取子矩阵...")
    print(f"  npy_path: {npy_path}")
    print(f"  flank: {flank} bp ({flank//1000}kb)")
    print(f"  resolution: {resolution} bp ({resolution//1000}kb)")
    print(f"  window: {window} bins")
    print(f"  矩阵形状: {oe_matrix.shape}")
    print(f"  期望子矩阵大小: ({expected_size}, {expected_size})")
    print(f"  总边界数: {len(boundaries)}")
    
    submatrices = []
    valid_count = 0
    invalid_count = 0
    
    for _, row in boundaries.iterrows():
        # 计算边界中心位置
        boundary_center = (row['start'] + row['end']) // 2
        
        # 计算在矩阵中的 bin 索引
        center_bin = (boundary_center - start_pos) // resolution
        
        # 计算子矩阵的起止 bin
        start_bin = center_bin - window
        end_bin = center_bin + window
        
        # 检查边界
        if start_bin < 0 or end_bin > oe_matrix.shape[0]:
            invalid_count += 1
            continue
        
        # 提取子矩阵
        try:
            submatrix = oe_matrix[start_bin:end_bin, start_bin:end_bin]
            
            # 检查大小是否正确
            if submatrix.shape != (expected_size, expected_size):
                invalid_count += 1
                continue
            
            submatrices.append(submatrix)
            valid_count += 1
            
        except Exception as e:
            invalid_count += 1
            continue
    
    print(f"  有效子矩阵: {valid_count}")
    print(f"  无效子矩阵: {invalid_count}")
    
    if not submatrices:
        raise ValueError("没有提取到有效的子矩阵")
    
    stack = np.stack(submatrices, axis=2)
    return stack


def analyze_oe_boundary_pileup(
    npy_path: str,
    boundaries: pd.DataFrame,
    output_path: str,
    flank: int = 300000,
    resolution: int = 10000,
    start_pos: int = 0,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'RdBu_r',
    dpi: int = 300,
    method: str = 'mean',
    mask_diagonal: bool = True,
    mask_width: int = 1,
    title: str = None,
    # ↓↓↓ T-6.5 新加 4 参数 ↓↓↓
    cbar_label: Optional[str] = None,
    label_orientation: str = 'vertical',
    color_scale_for_cbar: str = 'linear',
    cbar_height: float = 0.15,
) -> None:
    """
    在 OE 矩阵上提取 TAD 边界并绘制累计热图。
    
    Parameters
    ----------
    npy_path : str
        npy 文件路径（OE 矩阵）
    boundaries : pd.DataFrame
        TAD 边界数据，包含列: chrom, start, end
    output_path : str
        输出路径（不含扩展名）
    flank : int
        侧翼大小（bp），默认 300000 = 300kb
    resolution : int
        分辨率（bp），默认 10000 = 10kb
    start_pos : int
        矩阵起始位置（bp），默认 0
    vmin : float, optional
        颜色范围最小值
    vmax : float, optional
        颜色范围最大值
    cmap : str
        颜色映射名称，默认 'RdBu_r'（适合 O/E 数据）
    dpi : int
        输出图像 DPI，默认 300
    method : str
        累计方法：'mean', 'median', 'sum'，默认 'mean'
    mask_diagonal : bool
        是否遮盖对角线，默认 True
    mask_width : int
        对角线遮盖宽度（bin），默认 1
    title : str, optional
        图表标题
        
    Returns
    -------
    None
        
    Example
    -------
    >>> boundaries = pd.read_csv('tad_boundaries.tsv', sep='\\t')
    >>> analyze_oe_boundary_pileup(
    ...     npy_path='50_1_chr2_10k.oe.npy',
    ...     boundaries=boundaries,
    ...     output_path='output/oe_boundary_pileup',
    ...     flank=300000,
    ...     resolution=10000
    ... )
    """
    from .layout import (
        setup_plot_style,
        calculate_heatmap_layout,
        save_figure_multi_format
    )
    
    # 设置绘图样式
    setup_plot_style()
    
    # 处理 cmap
    if isinstance(cmap, str):
        cmap = plt.colormaps[cmap]
    
    # 从 npy 提取子矩阵
    stack = extract_snippets_from_npy(
        npy_path=npy_path,
        boundaries=boundaries,
        flank=flank,
        resolution=resolution,
        start_pos=start_pos
    )
    
    # 计算累计矩阵
    pileup_matrix = calculate_pileup_matrix(
        stack, 
        method=method, 
        mask_diagonal=mask_diagonal, 
        mask_width=mask_width
    )
    
    # 计算 vmin/vmax
    if vmin is None or vmax is None:
        vmin, vmax = get_matrix_range(pileup_matrix, vmin, vmax)
    
    # 创建图表
    layout = calculate_heatmap_layout(n_plots=1, plot_size=4.0)
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    ax = fig.add_axes([
        layout['margin_left'] / layout['total_width'],
        layout['margin_bottom'] / layout['total_height'],
        layout['plot_width'],
        layout['plot_height']
    ])
    
    # 绘制热图
    pileup_matrix_masked = np.ma.masked_invalid(pileup_matrix)
    im = ax.imshow(
        pileup_matrix_masked, 
        cmap=cmap, 
        aspect='auto', 
        interpolation='none',
        vmin=vmin, 
        vmax=vmax
    )
    
    # 设置坐标轴和边框
    setup_axes(ax)
    
    # 设置标签
    window = flank // resolution
    ax.set_xlabel('boundary position')
    ax.set_ylabel('boundary position')
    
    # 添加 colorbar
    from .colorbar import setup_colorbar
    setup_colorbar(
        fig=fig,
        sc=im,
        vmin=vmin,
        vmax=vmax,
        balance=False,
        color_scale=color_scale_for_cbar,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.75,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=cbar_height,
        label_fontsize=6,
        tick_fontsize=5,
        label=cbar_label if cbar_label is not None else label,
        label_orientation=label_orientation,
    )
    
    # 保存图表
    save_figure_multi_format(fig, output_path, dpi=dpi)
    print(f"✅ OE TAD 边界累计图已保存: {output_path}")


def cooltools_tad_pileup(
    clr,
    boundaries: pd.DataFrame,
    flank: int = 300000,
    expected_df: Optional[pd.DataFrame] = None,
    view_df: Optional[pd.DataFrame] = None,
    method: str = 'mean',
    nproc: int = 1,
) -> np.ndarray:
    """统一接口的 TAD boundary pileup，内部调 cooltools.pileup（已 published）。

    跟 cooltools 官方教程对齐（min_diag='auto' + np.nanmean）。

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler 对象
    boundaries : pd.DataFrame
        TAD boundaries，列: chrom, start, end
    flank : int
        Flank size in bp（默认 300000 = 300kb）
    expected_df : pd.DataFrame, optional
        None = balance 子矩阵，非 None = obs/exp 子矩阵
    view_df : pd.DataFrame, optional
        染色体臂 view。None = 从 clr 自动读 chromsizes
    method : str
        累计方法: 'mean' / 'median' / 'sum'
    nproc : int
        cooltools.pileup 内部并行进程数

    Returns
    -------
    mtx : np.ndarray
        累计 2D 矩阵，shape (D, D)，D = 2*flank//binsize
        对角线 ± 2 bin = NaN（cooltools 内部 min_diag='auto' 行为）
    """
    import cooltools
    from cfizz.analyze.compartment import get_view_df

    if view_df is None:
        view_df = get_view_df(clr)

    stack = cooltools.pileup(
        clr=clr,
        features_df=boundaries,
        view_df=view_df,
        expected_df=expected_df,
        flank=flank,
        min_diag='auto',
        nproc=nproc,
    )
    print(f"stack.shape: {stack.shape}, n_features: {stack.shape[0]}")
    print(f"stack 第一条 snippet[0]: non-NaN count = {np.sum(~np.isnan(stack[0]))}, max = {np.nanmax(stack[0])}, min = {np.nanmin(stack[0])}")
    print(f"stack 全部: nan_count = {np.isnan(stack).sum()}, total = {stack.size}")
    expected_size = 2 * (flank // clr.binsize)
    if stack.shape[1] != expected_size or stack.shape[2] != expected_size:
        stack = stack[:, :expected_size, :expected_size]

    if method == 'mean':
        mtx = np.nanmean(stack, axis=0)
    elif method == 'median':
        mtx = np.nanmedian(stack, axis=0)
    elif method == 'sum':
        mtx = np.nansum(stack, axis=0)
    else:
        raise ValueError(f"Unsupported method: {method}")

    print(f"mtx: min = {np.nanmin(mtx)}, max = {np.nanmax(mtx)}, nonzero = {np.sum(mtx > 0)}, all_zero = {np.all(mtx == 0)}")
    return mtx


def calculate_pileup_matrix(
    stack: np.ndarray,
    method: str = 'mean',
    mask_diagonal: bool = True,
    mask_width: int = 1
) -> np.ndarray:
    """Calculate pileup matrix from stack of submatrices.
    
    Parameters
    ----------
    stack : np.ndarray
        Stack of submatrices with shape (n_bins, n_bins, n_samples)
    method : str
        Calculation method: 'mean', 'median', or 'sum'
    mask_diagonal : bool
        Whether to mask the diagonal
    mask_width : int
        Width of diagonal mask (in bins)
        
    Returns
    -------
    mtx : np.ndarray
        Aggregated pileup matrix
    """
    if method == 'mean':
        mtx = np.nanmean(stack, axis=2)
    elif method == 'median':
        mtx = np.nanmedian(stack, axis=2)
    elif method == 'sum':
        mtx = np.nansum(stack, axis=2)
    else:
        raise ValueError(f"Unsupported method: {method}")
    
    # Mask diagonal
    if mask_diagonal:
        n = mtx.shape[0]
        for i in range(n):
            for j in range(max(0, i - mask_width), min(n, i + mask_width + 1)):
                mtx[i, j] = 0
                mtx[j, i] = 0
    
    return mtx


def plot_pileup_heatmap(
    matrix: np.ndarray,
    flank: int,
    resolution: int,
    output_path: str,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    dpi: int = 3000,
    balance: bool = False,
    color_scale: str = 'linear',
    label: str = '',
    # ↓↓↓ T-6.5 新加 4 参数 ↓↓↓
    cbar_label: Optional[str] = None,
    label_orientation: str = 'vertical',
    color_scale_for_cbar: Optional[str] = None,
    cbar_height: float = 0.15,
) -> None:
    """Plot pileup heatmap.
    
    Parameters
    ----------
    matrix : np.ndarray
        Pileup matrix
    flank : int
        Flank size in bp
    resolution : int
        Resolution in bp
    output_path : str
        Output path (without extension)
    vmin : float, optional
        Minimum value for color scale
    vmax : float, optional
        Maximum value for color scale
    cmap : str
        Colormap name
    dpi : int
        DPI for output image
    balance : bool
        Whether data is balanced
    color_scale : str
        Color scale: 'linear' or 'log'
    label : str
        Colorbar label
    """
    # 设置绘图样式
    setup_plot_style()
    
    # 处理 cmap
    if isinstance(cmap, str):
        cmap = plt.colormaps[cmap]
    
    # 计算 vmin/vmax 使用非零值
    vmin, vmax = get_matrix_range(matrix, vmin, vmax)

    # Create figure with proper layout
    layout = calculate_heatmap_layout(n_plots=1, plot_size=4.0)
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))

    ax = fig.add_axes([
        layout['margin_left'] / layout['total_width'],
        layout['margin_bottom'] / layout['total_height'],
        layout['plot_width'],
        layout['plot_height']
    ])

    # Handle special values
    if color_scale == 'log':
        matrix = np.ma.masked_less_equal(matrix, 0)
    else:
        matrix = np.ma.masked_invalid(matrix)
    
    # 绘制热图
    if color_scale == 'log':
        im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='none',
                       norm=LogNorm(vmin=vmin, vmax=vmax))
    else:
        im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='none',
                       vmin=vmin, vmax=vmax)
    
    # 设置坐标轴和边框
    setup_axes(ax)
    
    # 设置标签
    window = flank // resolution
    ax.set_xlabel('boundary position')
    ax.set_ylabel('boundary position')
    
    # 添加 colorbar
    from .colorbar import setup_colorbar
    setup_colorbar(
        fig=fig,
        sc=im,
        vmin=vmin,
        vmax=vmax,
        balance=False,
        color_scale=color_scale_for_cbar,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.75,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=cbar_height,
        label_fontsize=6,
        tick_fontsize=5,
        label=cbar_label if cbar_label is not None else label,
        label_orientation=label_orientation,
    )
    
    # 保存图表
    save_figure_multi_format(fig, output_path, dpi=dpi)
    print(f"✅ OE TAD 边界累计图已保存: {output_path}")


def analyze_tad_boundary_pileup(
    mcool_path: str,
    tad_boundaries: pd.DataFrame,
    output_path: str,
    flank: int = 300000,
    resolution: int = 10000,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    dpi: int = 3000,
    balance: bool = True,
    method: str = 'mean',
    n_processes: Optional[int] = None,
    color_scale: str = 'linear',
    label: str = '',
    top_n: Optional[int] = None
) -> None:
    """Analyze TAD boundary pileup for a single sample.
    
    Parameters
    ----------
    mcool_path : str
        Path to mcool file
    tad_boundaries : pd.DataFrame
        TAD boundary data with columns: chrom, start, end
        Optionally: insulation_score for top_n selection
    output_path : str
        Output path (without extension)
    flank : int
        Flank size in bp (default: 300000 = 300kb)
    resolution : int
        Resolution in bp (default: 10000 = 10kb)
    vmin : float, optional
        Minimum value for color scale
    vmax : float, optional
        Maximum value for color scale
    cmap : str
        Colormap name (default: 'Reds')
    dpi : int
        DPI for output image (default: 3000)
    balance : bool
        Whether to use balanced matrix (default: True)
    method : str
        Calculation method: 'mean', 'median', 'sum' (default: 'mean')
    n_processes : int, optional
        Number of parallel processes
    color_scale : str
        Color scale: 'linear' or 'log' (default: 'linear')
    label : str
        Colorbar label text
    top_n : int, optional
        Only use top N boundaries with lowest insulation score
    
    Returns
    -------
    None
    """
    # Select top N boundaries if specified
    if top_n is not None and 'insulation_score' in tad_boundaries.columns:
        tad_boundaries = tad_boundaries.nsmallest(top_n, 'insulation_score')
        print(f"Using top {top_n} boundaries with lowest insulation score")
    
    # T-6.9: 调统一 cooltools_tad_pileup（风格 1:balance）
    import cooler as cooler_lib
    clr_obj = cooler_lib.Cooler(f"{mcool_path}::resolutions/{resolution}")
    mtx = cooltools_tad_pileup(
        clr=clr_obj,
        boundaries=tad_boundaries,
        flank=flank,
        expected_df=None,
        method=method,
        nproc=n_processes if n_processes is not None else 1,
    )
    
    # Apply log10 transformation if specified
    if color_scale == 'log10_linear':
        mtx = np.log10(np.ma.masked_less_equal(mtx, 0).filled(1e-10))
        colorbar_label = 'log10(Hi-C)'
        plot_color_scale = 'linear'
    else:
        colorbar_label = label if label else ''
        plot_color_scale = color_scale
    
    # Plot heatmap
    plot_pileup_heatmap(
        matrix=mtx,
        flank=flank,
        resolution=resolution,
        output_path=output_path,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        dpi=dpi,
        balance=balance,
        color_scale=plot_color_scale,
        label=colorbar_label
    )


def plot_multi_tad_boundary_pileup(
    mcool_paths: List[str],
    boundaries_list: List[pd.DataFrame],
    output_path: str,
    sample_names: List[str] = None,
    flank: int = 300000,
    resolution: int = 10000,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    dpi: int = 3000,
    balance: bool = True,
    method: str = 'sum',
    n_processes: int = None,
    color_scale: str = 'log10_linear',
    top_n: Optional[int] = None,
    plot_size: float = 4.0,             # T-6.17: 默认 4cm(跟其它 viz 统一)
    # ↓↓↓ T-6.5 新加 4 参数 ↓↓↓
    cbar_label: Optional[str] = None,
    label_orientation: str = 'vertical',
    color_scale_for_cbar: Optional[str] = None,
    cbar_height: float = 0.15,
    # ↓↓↓ 已有的 2 个参数 ↓↓↓
    expected_df: Optional[pd.DataFrame] = None,
    view_df: Optional[pd.DataFrame] = None,
    # ↓↓↓ T-6.12 新加 ↓↓↓
    expected_dfs: Optional[List[Optional[pd.DataFrame]]] = None,  # 每个 sample 自己的 expected_df
) -> None:
    """Plot multi-sample TAD boundary pileup comparison.
    
    Parameters
    ----------
    mcool_paths : list
        List of mcool file paths
    boundaries_list : list
        List of TAD boundary DataFrames
    output_path : str
        Output path (without extension)
    sample_names : list, optional
        List of sample names
    flank : int
        Flank size in bp (default: 300000)
    resolution : int
        Resolution in bp (default: 10000)
    vmin : float, optional
        Minimum value for color scale
    vmax : float, optional
        Maximum value for color scale
    cmap : str
        Colormap name (default: 'Reds')
    dpi : int
        DPI for output image (default: 3000)
    balance : bool
        Whether to use balanced matrix (default: True)
    method : str
        Calculation method (default: 'sum')
    n_processes : int, optional
        Number of parallel processes
    color_scale : str
        Color scale: 'linear', 'log', 'log10_linear' (default: 'log10_linear')
    top_n : int, optional
        Only use top N boundaries per sample
    plot_size : float
        Size of each subplot in inches (default: 1.8)
    expected_df : pd.DataFrame, optional
        Obs/exp expected DataFrame (single, shared by all samples for backward compat)
    expected_dfs : list of pd.DataFrame, optional
        Each sample's own expected_df. None = use expected_df (backward compat).
        e.g. expected_dfs=[df_50_1, df_50_2] means sample 0 uses df_50_1, sample 1 uses df_50_2.
        Element can be None = that sample uses balance path (no P(s) division).
    """
    # 设置绘图样式
    setup_plot_style()
    
    # 处理 cmap
    if isinstance(cmap, str):
        cmap = plt.colormaps[cmap]
    
    n_samples = len(mcool_paths)
    
    if sample_names is None:
        sample_names = [f"Sample{i+1}" for i in range(n_samples)]
    
    assert len(boundaries_list) == n_samples
    assert len(sample_names) == n_samples
    
    # Calculate all pileup matrices
    all_matrices = []
    for i in range(n_samples):
        boundaries = boundaries_list[i]
        
        # Select top N if specified
        if top_n is not None and 'insulation_score' in boundaries.columns:
            boundaries = boundaries.nsmallest(top_n, 'insulation_score')
            print(f"Sample {sample_names[i]}: using top {top_n} boundaries")
        
        # 决定该 sample 走 obs/exp 还是 balance
        # T-6.12: expected_dfs 优先,其次用 expected_df(backward compat)
        if expected_dfs is not None:
            sample_expected = expected_dfs[i]
        elif expected_df is not None:
            sample_expected = expected_df
        else:
            sample_expected = None
        if sample_expected is not None:
            # ====== 风格 2:obs/exp ======
            import cooler as cooler_lib
            import cooltools

            # view_df 处理(默认从 mcool 读 chromsizes)
            if view_df is None:
                from cfizz.analyze.compartment import get_view_df
                clr_obj = cooler_lib.Cooler(f"{mcool_paths[i]}::resolutions/{resolution}")
                view_df_local = get_view_df(clr_obj)
            else:
                clr_obj = cooler_lib.Cooler(f"{mcool_paths[i]}::resolutions/{resolution}")
                view_df_local = view_df

            # cooltools.pileup 算 obs/exp 子矩阵(返回 shape (n, D, D))
            stack = cooltools.pileup(
                clr=clr_obj,
                features_df=boundaries,
                view_df=view_df_local,
                expected_df=sample_expected,
                flank=flank,
            )
            # 截断到 DxD(cooltools 返回 61x61,OBO 1)
            expected_size = 2 * (flank // resolution)
            if stack.shape[1] != expected_size or stack.shape[2] != expected_size:
                stack = stack[:, :expected_size, :expected_size]
            # 关键修复:stack 是 (n, D, D),手动 mean 在 axis=0
            # 不自己写 mask 循环 — cooltools 内部 min_diag='auto' 已经把对角 ± 2 bin 设 NaN
            # 自己写 mask 会把对角设 0 → np.log2(0) = -inf → 1e-10 → -33 → clamp 到 vmin=-1 全深红(BUG)
            mtx = np.nanmean(stack, axis=0)
        else:
            # ====== T-6.9: 新路径(风格 1:balance),调统一 cooltools_tad_pileup ======
            import cooler as cooler_lib
            clr_obj = cooler_lib.Cooler(f"{mcool_paths[i]}::resolutions/{resolution}")
            mtx = cooltools_tad_pileup(
                clr=clr_obj,
                boundaries=boundaries,
                flank=flank,
                expected_df=None,
                method=method,
                nproc=n_processes if n_processes is not None else 1,
            )
        
        # 颜色缩放
        if color_scale == 'log10_linear':
            mtx = np.log10(np.ma.masked_less_equal(mtx, 0).filled(1e-10))
        elif color_scale == 'log2':  # ← 新加分支
            mtx = np.log2(np.ma.masked_less_equal(mtx, 0).filled(1e-10))
        
        all_matrices.append(mtx)
    
    # Get unified vmin/vmax from non-zero non-NaN values across all matrices
    # 排除 NaN:np.nonzero 不排除 NaN(NaN != 0 是 True),会污染 min/percentile
    # 跟 T-6.10 viz/pileup.py:60 get_matrix_range 修法一致
    all_values_list = []
    for m in all_matrices:
        if not np.any(m):
            continue
        nonzero = m[np.nonzero(m)]
        valid = nonzero[~np.isnan(nonzero)]
        if valid.size > 0:
            all_values_list.append(valid)
    all_values = np.concatenate(all_values_list) if all_values_list else np.array([])

    if all_values.size > 0:
        if vmin is None:
            vmin = np.min(all_values)
        if vmax is None:
            perc = np.percentile(all_values, 93)
            vmax = all_values[np.abs(all_values - perc).argmin()]
    else:
        vmin = 0
        vmax = 1
    
    if color_scale == 'log10_linear':
        auto_label = 'log10(Hi-C)'
        plot_color_scale = 'linear'
    elif color_scale == 'log2':
        auto_label = 'log2(obs/exp)'
        plot_color_scale = 'linear'
    else:
        auto_label = ''
        plot_color_scale = color_scale

    # ↓↓↓ T-6.5: 优先用入参 cbar_label,否则用自动算的 auto_label ↓↓↓
    final_label = cbar_label if cbar_label is not None else auto_label
    final_cbar_scale = color_scale_for_cbar if color_scale_for_cbar is not None else plot_color_scale
    
    # Create figure
    n_cols = min(8, n_samples)
    layout = calculate_heatmap_layout(n_plots=n_cols, plot_size=plot_size)
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # 设置字体
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5
    
    for i, (matrix, name) in enumerate(zip(all_matrices, sample_names)):
        col = i % n_cols
        row = i // n_cols
        left = layout['margin_left'] / layout['total_width'] + col * (layout['plot_width'] + layout['hspace'])
        bottom = layout['margin_bottom'] / layout['total_height'] + (layout['n_rows'] - 1 - row) * (layout['plot_height'] + layout['vspace'])
        ax = fig.add_axes([left, bottom, layout['plot_width'], layout['plot_height']])
        
        # Handle special values
        if plot_color_scale == 'log':
            matrix = np.ma.masked_less_equal(matrix, 0)
            im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='none',
                           norm=LogNorm(vmin=vmin, vmax=vmax))
        else:
            matrix = np.ma.masked_invalid(matrix)
            im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='none',
                           vmin=vmin, vmax=vmax)
        
        # 设置坐标轴和边框
        setup_axes(ax)
        
        # Set title (使用 Arial 字体)
        ax.text(0.5, 1.02, name, 
                transform=ax.transAxes, ha='center', va='bottom', 
                fontsize=5, fontname='Arial')
        
        if i == 0:
            ax.set_xlabel('boundary position')
            ax.set_ylabel('boundary position')
    
    # Add colorbar
    from .colorbar import setup_colorbar
    setup_colorbar(
        fig=fig,
        sc=im,
        vmin=vmin,
        vmax=vmax,
        balance=balance,
        color_scale=final_cbar_scale,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.75,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=cbar_height,
        label_fontsize=6,
        tick_fontsize=5,
        label=final_label,
        label_orientation=label_orientation,
    )
    
    # Save figure
    save_figure_multi_format(fig, output_path, dpi=dpi)


# =============================================================================
# Legacy functions for backward compatibility
# =============================================================================

def plot_pileup(
    matrices: List[np.ndarray],
    center_positions: List[int],
    resolution: int,
    window_size: int = 20,
    vmin: float = 0,
    vmax: float = 2,
    cmap: Optional[str] = None,
    figsize: Tuple[float, float] = (6, 5),
    title: Optional[str] = None,
    show_center: bool = True,
    show_diag: bool = True,
    colorbar_label: str = "Average contact frequency",
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Create aggregate pileup plot around specified positions.
    
    NOTE: This is a legacy function. For TAD boundary pileup, use
    analyze_tad_boundary_pileup() instead.
    
    Parameters
    ----------
    matrices : list
        List of contact matrices
    center_positions : list
        List of center positions (in bp) for each matrix
    resolution : int
        Bin resolution (bp)
    window_size : int
        Window size around center (in bins)
    vmin : float
        Minimum color value
    vmax : float
        Maximum color value
    cmap : str, optional
        Colormap name (default: blue-white-red for ratio data)
    figsize : tuple
        Figure size in inches
    title : str, optional
        Plot title
    show_center : bool
        Whether to show center lines
    show_diag : bool
        Whether to show diagonal line
    colorbar_label : str
        Colorbar label text
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    """
    # Aggregate matrices
    aggregate_matrix = np.zeros((2 * window_size, 2 * window_size))
    count = 0
    
    half_window = window_size // 2
    
    for matrix, pos in zip(matrices, center_positions):
        start_bin = pos // resolution
        
        # Calculate extraction bounds
        start_extract = max(0, start_bin - half_window)
        end_extract = min(matrix.shape[0], start_bin + half_window)
        
        extract_size = end_extract - start_extract
        
        if extract_size > 0:
            region = matrix[start_extract:end_extract, start_extract:end_extract]
            
            # Calculate padding
            pad_before = half_window - (start_bin - start_extract)
            pad_after = half_window - (end_extract - start_bin)
            
            # Create padded region
            padded = np.full((window_size, window_size), np.nan)
            actual_start = max(0, pad_before)
            actual_end = window_size - max(0, pad_after)
            padded[actual_start:actual_end, actual_start:actual_end] = region
            
            # Accumulate
            aggregate_matrix = np.nansum(np.dstack([aggregate_matrix, padded]), axis=2)
            count += 1
    
    if count == 0:
        raise ValueError("No valid positions found")
    
    aggregate_matrix /= count
    
    # Set colormap
    if cmap is None:
        cmap = setup_ratio_colormap()
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot aggregate matrix
    norm = Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(aggregate_matrix, cmap=cmap, norm=norm, origin='lower')
    
    # Mark center
    if show_center:
        center = window_size
        ax.axhline(y=center, color='white', linewidth=1, linestyle='--', alpha=0.7)
        ax.axvline(x=center, color='white', linewidth=1, linestyle='--', alpha=0.7)
    
    # Set axis labels
    ax.set_xlabel('Position relative to center', fontsize=10)
    ax.set_ylabel('Position relative to center', fontsize=10)
    
    # Set ticks
    tick_positions = [0, window_size // 2, window_size, window_size + window_size // 2, 2 * window_size]
    tick_labels = ['-{:.0f}'.format(window_size), '-{:.0f}'.format(window_size // 2), '0', 
                   '+{:.0f}'.format(window_size // 2), '+{:.0f}'.format(window_size)]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels, fontsize=8)
    
    # Add title
    if title:
        ax.set_title(f'{title} (n={count})', fontsize=10)
    else:
        ax.set_title(f'Aggregate Pileup (n={count})', fontsize=10)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(colorbar_label, fontsize=10)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


def plot_multi_pileup(
    pileup_results: List[Dict[str, Any]],
    n_cols: int = 3,
    plot_size: float = 2.5,
    vmin: float = 0,
    vmax: float = 2,
    cmap: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Plot multiple pileup results in a grid.
    
    Parameters
    ----------
    pileup_results : list
        List of dictionaries with 'matrix' and 'title'
    n_cols : int
        Number of columns in subplot grid
    plot_size : float
        Size of each subplot (inches)
    vmin : float
        Minimum color value
    vmax : float
        Maximum color value
    cmap : str, optional
        Colormap name
    figsize : tuple, optional
        Overall figure size
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    """
    n_results = len(pileup_results)
    n_rows = (n_results + n_cols - 1) // n_cols
    
    if figsize is None:
        fig_width = n_cols * (plot_size + 0.5)
        fig_height = n_rows * (plot_size + 0.5)
        figsize = (fig_width, fig_height)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    
    if n_results == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Set colormap
    if cmap is None:
        cmap = setup_ratio_colormap()
    
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    for i, result in enumerate(pileup_results):
        ax = axes[i]
        
        matrix = result.get('matrix')
        title = result.get('title', f'Pileup {i+1}')
        
        if matrix is None:
            ax.axis('off')
            continue
        
        im = ax.imshow(matrix, cmap=cmap, norm=norm, origin='lower')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=8)
    
    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    
    # Add shared colorbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Contact frequency', fontsize=10)
    
    plt.tight_layout(rect=[0, 0, 0.93, 1])
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig
