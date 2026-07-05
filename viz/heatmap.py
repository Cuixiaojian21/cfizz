"""
Heatmap visualization module.

Functions for plotting Hi-C contact matrices.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.collections import QuadMesh
import itertools
from typing import Optional, List, Tuple

from cfizz.viz.layout import (
    setup_plot_style,
    calculate_heatmap_layout,
    setup_axes,
    add_coordinate_labels,
    setup_colorbar,
    log2_and_mask,
)

from cfizz.utils.coordinates import (
    get_matrix_range,
)


# =============================================================================
# =============================================================================

def get_bin_index(query_start: int, query_end: int, matrix_start: int, resolution: int) -> int:
    """
    Calculate bin index for a given genomic position.
    
    
    Parameters
    ----------
    query_start : int
        Query start position
    query_end : int
        Query end position
    matrix_start : int
        Matrix start position
    resolution : int
        Resolution
        
    Returns
    -------
    int
        Bin index
    """
    if (query_end - query_start) != resolution:
        raise ValueError("Interval length must match resolution")
    idx = (query_start - matrix_start) // resolution
    if idx < 0:
        raise ValueError("Query start outside matrix range")
    return int(idx)


def mark_boundaries_from_insulation(
    insulation_data,
    matrix_data: np.ndarray,
    matrix_start: int,
    resolution: int
) -> np.ndarray:
    """
    Mark TAD boundaries from insulation data.
    
    
    Parameters
    ----------
    insulation_data : pd.DataFrame
        Insulation score data with 'start', 'end', 'is_boundary' columns
    matrix_data : np.ndarray
        Contact matrix
    matrix_start : int
        Matrix start position
    resolution : int
        Resolution
        
    Returns
    -------
    np.ndarray
        Marked matrix with boundaries
    """
    import pandas as pd
    marked_matrix = np.full_like(matrix_data, np.nan, dtype=float)
    boundaries = insulation_data[insulation_data['is_boundary']]
    if len(boundaries) == 0:
        return marked_matrix
    boundary_indices = []
    for _, boundary in boundaries.iterrows():
        idx = get_bin_index(boundary['start'], boundary['end'], matrix_start, resolution)
        boundary_indices.append(idx)
    if len(boundary_indices) == 1:
        i = boundary_indices[0]
        marked_matrix[i, :] = 1
        marked_matrix[:, i] = 1
    else:
        for k in range(len(boundary_indices)-1):
            i = boundary_indices[k]
            j = boundary_indices[k+1]
            marked_matrix[i, i:j+1] = 1
            marked_matrix[j, i:j+1] = 1
            marked_matrix[i:j+1, i] = 1
            marked_matrix[i:j+1, j] = 1
    return marked_matrix


# =============================================================================
# MAIN HEATMAP FUNCTIONS - Standard Hi-C visualizations
# =============================================================================

def plot_single_heatmap(
    matrix, 
    vmin=None, 
    vmax=None, 
    color_scale='linear', 
    cmap='Reds',
    start_pos=None, 
    end_pos=None, 
    chrom=None, 
    balance=False
):
    """
    Plot a single Hi-C matrix heatmap.
    
    This is a wrapper around plot_multi_heatmap for convenience.
    For single matrix visualization, this provides a simpler API.
    
    Parameters:
        matrix: np.ndarray, Hi-C contact matrix
        vmin, vmax: float, colorbar range
        color_scale: str, 'linear' or 'log'
        cmap: str or matplotlib.colors.Colormap, colormap name or object
        start_pos, end_pos: int, start and end positions (bp)
        chrom: str, chromosome name
        balance: bool, whether to use balanced matrix
        
    Returns:
        matplotlib.figure.Figure: The figure object
    """
    # Delegate to plot_multi_heatmap with single matrix
    return plot_multi_heatmap(
        matrices=[matrix],
        sample_names=None,
        ncols=1,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        color_scale=color_scale,
        plot_size=4.0,
        start_pos=start_pos,
        end_pos=end_pos,
        chrom=chrom,
        balance=balance
    )


def plot_heatmap(
    matrix: np.ndarray,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: str = "Reds",
    color_scale: str = "linear",
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300,
    figsize: Tuple[float, float] = None,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    chrom: Optional[str] = None,
    balance: bool = False
) -> plt.Figure:
    """
    Plot a single Hi-C heatmap (convenience wrapper).
    
    Parameters
    ----------
    matrix : np.ndarray
        Contact matrix
    vmin : float, optional
        Minimum value for color scale
    vmax : float, optional
        Maximum value for color scale
    cmap : str
        Colormap name
    color_scale : str
        Color scale ('linear' or 'log')
    title : str, optional
        Plot title
    xlabel, ylabel : str, optional
        Axis labels
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
    figsize : tuple, optional
        Figure size in inches
    start_pos, end_pos : int, optional
        Genomic positions for labels
    chrom : str, optional
        Chromosome name
    balance : bool
        Whether using balanced matrix
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    """
    fig = plot_single_heatmap(
        matrix=matrix,
        vmin=vmin,
        vmax=vmax,
        color_scale=color_scale,
        cmap=cmap,
        start_pos=start_pos,
        end_pos=end_pos,
        chrom=chrom,
        balance=balance
    )
    
    # Add labels
    ax = fig.axes[0] if fig.axes else None
    if ax is not None:
        if title:
            ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
    
    # Save if path provided
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


# =============================================================================
# Convenience wrappers: generate_* (read + plot + save)
# =============================================================================

def generate_multi_heatmap(
    file_paths: List[str],
    output_dir: str,
    sample_names: List[str],
    chrom: str,
    resolution: int,
    start_pos: int,
    end_pos: int,
    balance: bool = False,
    color_scale: str = 'linear',
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    sample_name: str = 'multi',
    plot_size: float = 4,
    output_prefix: str = None,
    formats: Tuple[str, ...] = ('png',),
    dpi: int = 300
) -> dict:
    """
    从多个 cooler 文件生成多子图热图（便利包装 = read + plot + save）。
    
    宪法对齐: §3 P3.1 单/多样本统一。内部调本地 read_matrix_from_cooler +
    plot_multi_heatmap + save_figure / save_figure_multi_format。
    
    Parameters
    ----------
    file_paths : List[str]
        mcool 文件路径列表
    output_dir : str
        输出目录（自动 mkdir）
    sample_names : List[str]
        样品名称列表（长度 = len(file_paths)）
    chrom : str
        染色体号（如 'chr1'）
    resolution : int
        分辨率（bp，如 1000000 = 1Mb）
    start_pos, end_pos : int
        起始/结束位置（bp）
    balance : bool
        是否使用 balance 矩阵
    color_scale : str
        'linear' / 'log'
    vmin, vmax : float, optional
        颜色条范围（None = 自动）
    cmap : str
        matplotlib colormap 名
    sample_name : str
        用于文件名（默认 'multi'）
    plot_size : float
        单子图宽度（英寸）
    output_prefix : str, optional
        手动指定输出文件前缀（默认 = generate_output_filename 自动）
    formats : Tuple[str]
        输出格式，默认 ('png',)；双格式用 ('png', 'svg')
    dpi : int
        输出 DPI
    
    Returns
    -------
    dict
        {format: output_path} 映射（如 {'png': '/path/to/multi.chr1.0-50M.1Mb.png'}）
    
    Examples
    --------
    >>> result = generate_multi_heatmap(
    ...     file_paths=['sample1.mcool', 'sample2.mcool'],
    ...     output_dir='./output',
    ...     sample_names=['WT', 'KO'],
    ...     chrom='chr1', resolution=1000000,
    ...     start_pos=0, end_pos=50000000,
    ...     balance=False, color_scale='log',
    ...     formats=('png', 'svg')
    ... )
    >>> print(result['png'], result['svg'])
    """
    from cfizz.viz.layout import generate_output_filename, save_figure, save_figure_multi_format
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 读所有矩阵
    matrices = []
    for fp in file_paths:
        M = read_matrix_from_cooler(fp, resolution, chrom, start_pos, end_pos, balance)
        if M is None:
            raise ValueError(f"无法从文件 {fp} 读取矩阵")
        matrices.append(M)
    
    # 2. 生成输出文件名前缀
    # 如果 end_pos 为 None，从第一个 mcool 文件获取染色体长度
    _end_pos = end_pos
    if _end_pos is None and file_paths:
        import cooler
        fp0 = file_paths[0]
        res_uri = f"{fp0}::/resolutions/{resolution}" if fp0.endswith('.mcool') else fp0
        clr = cooler.Cooler(res_uri)
        _end_pos = clr.chromsizes[chrom]
    
    if output_prefix is None:
        # region 是 tuple (start_pos, end_pos)，跟 layout.py generate_output_filename 签名对齐
        # 注意: sample_name='multi' 时直接用 'multi'，避免生成 'multi.multi'
        if sample_name and sample_name != 'multi':
            prefix_name = f"multi.{sample_name}"
        else:
            prefix_name = "multi"
        prefix = generate_output_filename(
            output_dir=output_dir,
            prefix=prefix_name,
            chrom=chrom,
            region=(start_pos, _end_pos),
            resolution=resolution,
            balance=balance,
            suffix=color_scale
        )
    else:
        prefix = os.path.join(output_dir, output_prefix)
    
    # 3. 绘图（本地 plot_multi_heatmap，跟宪法 P2.4-5 多样本对比对齐）
    fig = plot_multi_heatmap(
        matrices=matrices,
        sample_names=sample_names,
        vmin=vmin,
        vmax=vmax,
        color_scale=color_scale,
        cmap=cmap,
        start_pos=start_pos,
        end_pos=end_pos,
        chrom=chrom,
        balance=balance,
        plot_size=plot_size
    )
    
    # 4. 保存（支持单格式 / 双格式）
    # 注意: generate_output_filename 已经返回带扩展名的路径，不再重复添加
    if len(formats) == 1:
        output_path = prefix  # prefix 已含扩展名
        save_figure(fig, output_path, dpi=dpi)
        return {formats[0]: output_path}
    else:
        # 多格式 = 用 save_figure_multi_format，返回 list，转成 dict
        saved_list = save_figure_multi_format(fig, prefix, dpi=dpi, formats=formats)
        return {fmt: path for fmt, path in zip(formats, saved_list)}


def generate_heatmap(
    file_path: str,
    output_dir: str,
    chrom: str,
    resolution: int,
    start_pos: int,
    end_pos: int,
    sample_name: str = "",
    balance: bool = False,
    color_scale: str = 'linear',
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    plot_size: float = 4,
    output_prefix: str = None,
    formats: Tuple[str, ...] = ('png',),
    dpi: int = 300
) -> dict:
    """
    从单个 cooler 文件生成热图（便利包装）。
    
    宪法对齐: §3 P3.1 单/多样本统一。**此函数 = generate_multi_heatmap(file_paths=[file_path]) 的特殊情况**。
    内部直接调 generate_multi_heatmap(matrices=[M], sample_names=[sample_name])。
    
    Parameters
    ----------
    file_path : str
        单个 mcool 文件路径
    sample_name : str
        单个样品名称
    (其余参数同 generate_multi_heatmap)
    
    Returns
    -------
    dict
        {format: output_path} 映射
    
    Examples
    --------
    >>> result = generate_heatmap(
    ...     file_path='sample1.mcool',
    ...     output_dir='./output',
    ...     sample_name='WT',
    ...     chrom='chr1', resolution=1000000,
    ...     start_pos=0, end_pos=50000000,
    ...     formats=('png', 'svg')
    ... )
    """
    return generate_multi_heatmap(
        file_paths=[file_path],
        output_dir=output_dir,
        sample_names=[sample_name],
        chrom=chrom,
        resolution=resolution,
        start_pos=start_pos,
        end_pos=end_pos,
        balance=balance,
        color_scale=color_scale,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        sample_name=sample_name,  # 单 sample 用 sample_name 作文件名
        plot_size=plot_size,
        output_prefix=output_prefix,
        formats=formats,
        dpi=dpi
    )


# =============================================================================
# =============================================================================

def read_matrix_from_cooler(file_path, resolution, chrom, start_pos, end_pos, balance=False, suppress_stderr=False):
    """
    从cooler文件（支持.cool和.mcool）中读取指定区域的矩阵
    
    
    Parameters:
        file_path: str, cooler文件路径
        resolution: int, 分辨率（bp）
        chrom: str, 染色体号
        start_pos: int, 起始位置（bp）
        end_pos: int, 结束位置（bp）
        balance: bool, 是否使用balance矩阵
        suppress_stderr: 是否屏蔽cooler的报错输出
    
    Returns:
        np.ndarray or None: 读取的矩阵，如果读取失败则返回None
    """
    import contextlib, io
    import cooler
    try:
        if file_path.endswith('.mcool'):
            res_uri = f"{file_path}::/resolutions/{resolution}"
            clr = cooler.Cooler(res_uri)
        else:
            clr = cooler.Cooler(file_path)
        chrom_length = clr.chromsizes[chrom]
        start_pos = max(0, start_pos)
        end_pos = min(chrom_length, end_pos) if end_pos is not None else chrom_length
        start_bin = start_pos // clr.binsize
        end_bin = end_pos // clr.binsize
        if suppress_stderr:
            with contextlib.redirect_stderr(io.StringIO()):
                M = clr.matrix(balance=balance, sparse=False).fetch((chrom, start_bin*clr.binsize, end_bin*clr.binsize))
        else:
            M = clr.matrix(balance=balance, sparse=False).fetch((chrom, start_bin*clr.binsize, end_bin*clr.binsize))
        M[np.isnan(M)] = 0
        return M
    except Exception as e:
        return None


def plot_multi_heatmap(
    matrices: List[np.ndarray],
    sample_names: Optional[List[str]] = None,
    ncols: int = 2,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: str = "Reds",
    color_scale: str = "linear",
    plot_size: float = 4.0,
    save_path: Optional[str] = None,
    dpi: int = 300,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    chrom: Optional[str] = None,
    balance: bool = False
) -> plt.Figure:
    """
    Plot multiple Hi-C heatmaps.
    
    Parameters
    ----------
    matrices : list
        List of contact matrices
    sample_names : list, optional
        Names for each sample
    ncols : int
        Number of columns in subplot grid
    vmin, vmax : float, optional
        Color scale range
    cmap : str
        Colormap name
    color_scale : str
        Color scale ('linear' or 'log')
    plot_size : float
        Size of each subplot in inches
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    """
    setup_plot_style()
    if isinstance(cmap, str):
        cmap = plt.colormaps[cmap]
    
    n_samples = len(matrices)
    n_cols = min(8, n_samples)
    
    # plot_size is in cm (cfizz uses cm as the design unit)
    # layout function expects cm, matplotlib expects inches
    layout = calculate_heatmap_layout(n_plots=n_cols, plot_size=plot_size)
    
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # Merge all non-zero values from all matrices to calculate global range
    all_values = np.concatenate([M[np.nonzero(M)] for M in matrices if M is not None])
    vmin, vmax = get_matrix_range(all_values, vmin, vmax)
    
    # Plot subplots
    for i, M in enumerate(matrices):
        if M is None:
            continue
            
        col = i % n_cols
        left = layout['margin_left_relative'] + col * (layout['plot_width'] + layout['gap_relative'])
        bottom = layout['margin_bottom_relative']
        
        ax = fig.add_axes([left, bottom, layout['plot_width'], layout['plot_height']])
        
        # Handle special values
        if color_scale == 'log':
            M = np.ma.masked_less_equal(M, 0)
        else:
            M = np.ma.masked_invalid(M)
        
        # Plot heatmap
        if color_scale == 'log':
            sc = ax.imshow(M, cmap=cmap, aspect='auto', interpolation='none',
                          norm=LogNorm(vmin=vmin, vmax=vmax))
        else:
            sc = ax.imshow(M, cmap=cmap, aspect='auto', interpolation='none',
                          vmin=vmin, vmax=vmax)
        
        # Set axes
        xmin, xmax, ymin, ymax = setup_axes(ax)
        
        # Add sample name (only if non-empty string)
        if sample_names and i < len(sample_names) and sample_names[i]:
            ax.text(0.5, 1.02, sample_names[i], 
                    transform=ax.transAxes, ha='center', va='bottom', fontsize=5)
        
        # Add coordinate labels only on the first subplot
        if i == 0 and all(v is not None for v in [start_pos, end_pos, chrom]):
            add_coordinate_labels(ax, xmin, xmax, ymin, ymax, start_pos, end_pos, chrom, fontsize=5)
    
    # Add shared colorbar
    setup_colorbar(
        fig=fig,
        sc=sc,
        vmin=vmin,
        vmax=vmax,
        balance=balance,
        color_scale=color_scale,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.75,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=0.15,
        label_fontsize=6,
        tick_fontsize=5
    )
    
    # Save if path provided
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


def plot_oe_heatmap(
    matrix: np.ndarray,
    start: int = 0,
    resolution: int = 10000,
    vmin: float = 0.25,
    vmax: float = 4,
    cmap: str = "RdBu_r",
    title: Optional[str] = None,
    chrom: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Plot O/E (Observed/Expected) normalized heatmap.
    
    This function removes the distance decay effect to highlight local chromatin
    structure such as TADs and loops.
    
    Parameters
    ----------
    matrix : np.ndarray
        O/E normalized contact matrix
    start : int
        Starting genomic position in bp
    resolution : int
        Bin resolution in bp
    vmin : float
        Minimum value for color scale (default: 0.25)
    vmax : float
        Maximum value for color scale (default: 4)
    cmap : str
        Colormap name (default: "RdBu_r")
    title : str, optional
        Plot title
    chrom : str, optional
        Chromosome name
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    """
    setup_plot_style()
    
    # Calculate layout
    layout = calculate_heatmap_layout(n_plots=1, plot_size=4.0)
    
    # Create figure
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # Add heatmap area
    ax = fig.add_axes([layout['margin_left_relative'], 
                      layout['margin_bottom_relative'], 
                      layout['plot_width'], 
                      layout['plot_height']])
    
    # Mask invalid values
    masked_matrix = np.ma.masked_invalid(matrix)
    
    # Plot heatmap
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    im = ax.imshow(masked_matrix, cmap=cmap, aspect='auto', interpolation='none', norm=norm)
    
    # Set axes and borders
    xmin, xmax, ymin, ymax = setup_axes(ax)
    
    # Add coordinate labels
    if chrom:
        end = start + matrix.shape[0] * resolution
        add_coordinate_labels(ax, xmin, xmax, ymin, ymax, start, end, chrom)
    
    # Add title
    if title:
        ax.set_title(title, fontsize=5)
    
    # Set colorbar
    setup_colorbar(
        fig=fig,
        sc=im,
        vmin=vmin,
        vmax=vmax,
        balance=False,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.75,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=0.15,
        label_fontsize=6,
        tick_fontsize=5
    )
    
    # Save if path provided
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig
