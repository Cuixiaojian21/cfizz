"""
Chromatin loops visualization module.

Functions for plotting chromatin loop results:
- Loop annotation on heatmaps
- APA (Aggregate Peak Analysis)
- Loop comparison between samples

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from matplotlib.colors import Normalize
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
import warnings

from .heatmap import plot_heatmap

# Suppress warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


def plot_loops(
    matrix: np.ndarray,
    loops: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    resolution: int,
    vmin: float = 0,
    vmax: float = 1,
    cmap: str = "Reds",
    loop_color: str = "#d62728",
    loop_marker_size: float = 50,
    figsize: Tuple[float, float] = (10, 10),
    show_saddle: bool = False,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Plot heatmap with chromatin loops overlaid.
    
    This function creates a Hi-C heatmap with chromatin loop coordinates
    marked as points on the heatmap.
    
    Parameters
    ----------
    matrix : np.ndarray
        Contact matrix
    loops : pd.DataFrame
        DataFrame with loop coordinates, must contain columns:
        [chrom1, start1, end1, chrom2, start2, end2] or similar
    chrom : str
        Chromosome name
    start : int
        Start position (bp)
    end : int
        End position (bp)
    resolution : int
        Bin resolution (bp)
    vmin : float
        Minimum color value
    vmax : float
        Maximum color value
    cmap : str
        Colormap name
    loop_color : str
        Color for loop markers
    loop_marker_size : float
        Size of loop markers
    figsize : tuple
        Figure size in inches
    show_saddle : bool
        Whether to show saddle plot around loops (requires more computation)
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
        
    Examples
    --------
    >>> matrix = reader.fetch("chr1", start=0, end=10000000)
    >>> loops = pd.read_csv("loops.bedpe", sep="\t")
    >>> fig = plot_loops(matrix, loops, "chr1", 0, 10000000, 100000)
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap using standard imshow
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='none',
                   vmin=vmin, vmax=vmax)
    
    # Filter loops for region
    region_loops = loops[
        (loops['chrom1'] == chrom) &
        (loops['start1'] >= start) &
        (loops['end1'] <= end) &
        (loops['chrom2'] == chrom) &
        (loops['start2'] >= start) &
        (loops['end2'] <= end)
    ].copy()
    
    # Plot loops as points
    for _, loop in region_loops.iterrows():
        # Calculate bin indices
        x = (loop['start2'] - start) // resolution
        y = (loop['start1'] - start) // resolution
        
        if 0 <= x < matrix.shape[0] and 0 <= y < matrix.shape[0]:
            ax.scatter(x, y, color=loop_color, s=loop_marker_size, 
                      zorder=10, marker='o', edgecolors='white', linewidths=0.5)
    
    # Set axis labels
    ax.set_xlabel(f'{chrom} position', fontsize=10)
    ax.set_ylabel(f'{chrom} position', fontsize=10)
    
    # Add title
    if title:
        ax.set_title(f'{title} ({len(region_loops)} loops)')
    else:
        ax.set_title(f'Chromatin Loops ({len(region_loops)} detected)', fontsize=10)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Contact frequency', fontsize=10)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


# =============================================================================
# add_loops_to_heatmap_45deg - for rotated 45-degree heatmaps
# =============================================================================

def add_loops_to_heatmap_45deg(
    ax: plt.Axes,
    loops_data: pd.DataFrame,
    matrix: np.ndarray,
    start: int,
    resolution: int,
    loop_color: str = 'blue',
    loop_alpha: float = 0.6,
    loop_size: float = 50
):
    """
    Add loops annotations to a 45-degree rotated Hi-C heatmap.

    This function plots loops as pairs of circular markers at anchor positions,
    using the same 45-degree coordinate transformation as the rotated heatmap.

    Parameters
    ----------
    ax : plt.Axes
        Matplotlib axes object for the heatmap
    loops_data : pd.DataFrame
        Loops data with columns: chrom1, start1, end1, chrom2, start2, end2
    matrix : np.ndarray
        Hi-C contact matrix
    start : int
        Region start position (genomic coordinate in bp)
    resolution : int
        Matrix resolution (bp per bin)
    loop_color : str, default 'blue'
        Color for loop markers
    loop_alpha : float, default 0.6
        Transparency of loop markers
    loop_size : float, default 50
        Size of loop markers (scatter plot size parameter)
    """
    print(f"\n🔍 Adding {len(loops_data)} loops to 45-deg heatmap...")
    if len(loops_data) == 0:
        print("⚠️  No loops data to plot")
        return

    # Rotation matrix for 45-degree transformation (same as heatmap)
    rotation_matrix = np.array([[1, 0.5], [-1, 0.5]])
    region_start = start
    matrix_size = matrix.shape[0]

    # Process each loop
    loops_added = 0
    loops_skipped = 0
    for idx, loop in loops_data.iterrows():
        # Get loop anchors (use genomic coordinates, NOT bin indices!)
        anchor1_genome = (loop['start1'] + loop['end1']) / 2
        anchor2_genome = (loop['start2'] + loop['end2']) / 2

        # Check if anchors are within matrix bounds
        if not (region_start <= anchor1_genome < region_start + matrix_size * resolution and
                region_start <= anchor2_genome < region_start + matrix_size * resolution):
            loops_skipped += 1
            continue

        # Apply 45-degree rotation using GENOMIC coordinates
        coords1_reversed = np.array([anchor2_genome, anchor1_genome])
        result1 = np.dot(coords1_reversed, rotation_matrix)
        x1, y1 = result1[1], result1[0]

        coords2_reversed = np.array([anchor1_genome, anchor2_genome])
        result2 = np.dot(coords2_reversed, rotation_matrix)
        x2, y2 = result2[1], result2[0]

        # Plot the two circular markers
        ax.scatter(x1, y1, s=loop_size, c='none', marker='o',
                  edgecolors=loop_color, alpha=loop_alpha, linewidths=0.8,
                  zorder=10)
        ax.scatter(x2, y2, s=loop_size, c='none', marker='o',
                  edgecolors=loop_color, alpha=loop_alpha, linewidths=0.8,
                  zorder=10)

        loops_added += 1

    print(f"✓ Added {loops_added} loops, skipped {loops_skipped}")


# =============================================================================
# F-3, F-4, F-5: add_loops_to_heatmap + plot_heatmap_with_loops + plot_multi_heatmap_with_loops
# =============================================================================

def add_loops_to_heatmap(
    ax,
    loops_data,
    matrix_data,
    matrix_start,
    resolution,
    loop_color='blue',
    loop_alpha=0.6,
    loop_size=50
):
    """
    在热图上添加loops标注
    
    
    Args:
        ax: matplotlib轴对象
        loops_data: loops数据DataFrame
        matrix_data: 接触矩阵数据
        matrix_start: 矩阵起始位置（bp）
        resolution: 分辨率（bp）
        loop_color: loops的颜色
        loop_alpha: loops的透明度
        loop_size: loops的大小
    """
    print(f"\n开始添加loops标注...")
    print(f"总loops数量: {len(loops_data)}")
    
    # 将基因组坐标转换为矩阵bin坐标
    for idx, loop in loops_data.iterrows():
        # 计算第一个锚点的bin范围
        s_l = range(loop['start1']//resolution, int(np.ceil(loop['end1']/float(resolution))))
        # 计算第二个锚点的bin范围
        e_l = range(loop['start2']//resolution, int(np.ceil(loop['end2']/float(resolution))))
        
        # 找到信号最强的点
        si, ei = None, None
        max_signal = -1
        
        for i in s_l:
            for j in e_l:
                st = i - matrix_start//resolution
                et = j - matrix_start//resolution
                if (0 <= st < matrix_data.shape[0]) and (0 <= et < matrix_data.shape[0]):
                    if matrix_data[st,et] > max_signal:
                        si, ei = st, et
                        max_signal = matrix_data[st,et]
        
        # 如果找到了有效的点，绘制标注
        if si is not None:
            # 绘制对称的两个点
            ax.scatter(ei, si, s=loop_size, c='none', marker='o',
                      edgecolors=loop_color, alpha=loop_alpha, linewidths=0.5)
            ax.scatter(si, ei, s=loop_size, c='none', marker='o',
                      edgecolors=loop_color, alpha=loop_alpha, linewidths=0.5)
            
            # 输出loops信息
            print(f"\nLoops {idx+1}:")
            print(f"  锚点1: {loop['chrom1']}:{loop['start1']}-{loop['end1']}")
            print(f"  锚点2: {loop['chrom2']}:{loop['start2']}-{loop['end2']}")
            print(f"  矩阵坐标: ({si}, {ei})")
            print(f"  信号强度: {max_signal:.2f}")
    
    print(f"\nloops标注完成！")


def plot_heatmap_with_loops(
    mcool_path,
    loops_path,
    chrom,
    start,
    end,
    resolution,
    output_path,
    cmap='Reds',
    vmin=None,
    vmax=None,
    color_scale='linear',
    loop_color='blue',
    loop_alpha=0.6,
    loop_size=50,
    balance=False,
    dpi=1000
):
    """
    绘制带有loops标注的热图
    
    
    Args:
        mcool_path: mcool文件路径
        loops_path: loops文件路径
        chrom: 染色体
        start: 起始位置
        end: 结束位置
        resolution: 分辨率
        output_path: 输出图片路径
        cmap: 颜色映射，默认为'Reds'
        vmin: 统一的最小值，用于控制颜色范围，默认为None（自动计算）
        vmax: 统一的最大值，用于控制颜色范围，默认为None（自动计算）
        color_scale: 颜色条缩放方式，'linear'或'log'，默认为'linear'
        loop_color: loops的颜色
        loop_alpha: loops的透明度
        loop_size: loops的大小
        balance: 是否使用平衡矩阵
        dpi: 分辨率
    """
    import os
    from cfizz.viz.heatmap import read_matrix_from_cooler, plot_single_heatmap
    from cfizz.io.loops import read_loops
    from cfizz.viz.layout import setup_axes, setup_colorbar, save_figure
    from cfizz.utils.coordinates import generate_output_filename, get_matrix_range
    
    # 1. 读取接触矩阵
    matrix = read_matrix_from_cooler(
        file_path=mcool_path,
        resolution=resolution,
        chrom=chrom,
        start_pos=start,
        end_pos=end,
        balance=balance
    )
    if matrix is None:
        raise ValueError(f"无法从文件 {mcool_path} 读取矩阵")
    
    # 添加主对角线mask
    mask_width = 2  # 主对角线两侧的宽度
    n = matrix.shape[0]
    for i in range(n):
        for j in range(max(0, i-mask_width), min(n, i+mask_width+1)):
            matrix[i,j] = 0
            matrix[j,i] = 0
    
    # 2. 读取loops数据
    loops_data = read_loops(
        file_path=loops_path,
        chrom=chrom,
        start=start,
        end=end
    )
    
    # 3. 绘制热图
    fig = plot_single_heatmap(
        matrix=matrix,
        vmin=vmin,
        vmax=vmax,
        color_scale=color_scale,
        cmap=cmap,
        start_pos=start,
        end_pos=end,
        chrom=chrom,
        balance=balance
    )
    
    # 4. 获取热图轴对象
    ax = fig.axes[0]
    
    # 5. 添加loops标注
    add_loops_to_heatmap(
        ax=ax,
        loops_data=loops_data,
        matrix_data=matrix,
        matrix_start=start,
        resolution=resolution,
        loop_color=loop_color,
        loop_alpha=loop_alpha,
        loop_size=loop_size
    )
    
    # 6. 保存图片（直接用 output_path 保存 PNG + SVG）
    # 注意：output_path 是完整路径（不含扩展名），save_figure 需要支持这个
    # 先保存 PNG
    fig.savefig(f"{output_path}.png", dpi=dpi, bbox_inches='tight')
    # 再保存 SVG
    fig.savefig(f"{output_path}.svg", dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"图片已保存: {output_path}.png/.svg")


def plot_multi_heatmap_with_loops(
    mcool_paths,
    loops_paths,
    chrom,
    start,
    end,
    resolution,
    output_path,
    sample_names=None,
    cmap='Reds',
    vmin=None,
    vmax=None,
    color_scale='linear',
    loop_color='blue',
    loop_alpha=0.6,
    loop_size=50,
    balance=False,
    dpi=1000,
    plot_size=4
):
    """
    绘制多个样本的带有loops标注的热图
    
    
    Args:
        mcool_paths: mcool文件路径列表
        loops_paths: loops文件路径列表
        chrom: 染色体
        start: 起始位置
        end: 结束位置
        resolution: 分辨率
        output_path: 输出图片路径
        sample_names: 样本名称列表，默认为None（使用mcool文件名）
        cmap: 颜色映射，默认为'Reds'
        vmin: 统一的最小值，用于控制颜色范围，默认为None（自动计算）
        vmax: 统一的最大值，用于控制颜色范围，默认为None（自动计算）
        color_scale: 颜色条缩放方式，'linear'或'log'，默认为'linear'
        loop_color: loops的颜色
        loop_alpha: loops的透明度
        loop_size: loops的大小
        balance: 是否使用平衡矩阵
        dpi: 分辨率
        plot_size: 子图大小（单位:cm），默认为 4
    """
    from cfizz.viz.heatmap import read_matrix_from_cooler, plot_multi_heatmap
    from cfizz.io.loops import read_loops
    from cfizz.viz.layout import setup_axes, setup_colorbar, save_figure
    from cfizz.utils.coordinates import get_matrix_range
    
    # 1. 读取所有矩阵
    matrices = []
    for mcool_path in mcool_paths:
        matrix = read_matrix_from_cooler(
            file_path=mcool_path,
            resolution=resolution,
            chrom=chrom,
            start_pos=start,
            end_pos=end,
            balance=balance
        )
        if matrix is None:
            raise ValueError(f"无法从文件 {mcool_path} 读取矩阵")
            
        # 添加主对角线mask
        mask_width = 2  # 主对角线两侧的宽度
        n = matrix.shape[0]
        for i in range(n):
            for j in range(max(0, i-mask_width), min(n, i+mask_width+1)):
                matrix[i,j] = 0
                matrix[j,i] = 0
                
        matrices.append(matrix)
    
    # 2. 绘制多子图热图
    fig = plot_multi_heatmap(
        matrices=matrices,
        sample_names=sample_names,
        vmin=vmin,
        vmax=vmax,
        color_scale=color_scale,
        cmap=cmap,
        start_pos=start,
        end_pos=end,
        chrom=chrom,
        balance=balance,
        plot_size=plot_size
    )
    
    # 3. 为每个子图添加loops标注
    for i, (matrix, loops_path) in enumerate(zip(matrices, loops_paths)):
        # 健壮性判断：loops_path为None或空字符串时跳过loops绘制
        if loops_path is None or (isinstance(loops_path, str) and not loops_path.strip()):
            continue  # 不画loops，只画heatmap
        # 读取loops数据
        loops_data = read_loops(
            file_path=loops_path,
            chrom=chrom,
            start=start,
            end=end
        )
        # 获取对应的子图轴对象
        ax = fig.axes[i]
        # 添加loops标注
        add_loops_to_heatmap(
            ax=ax,
            loops_data=loops_data,
            matrix_data=matrix,
            matrix_start=start,
            resolution=resolution,
            loop_color=loop_color,
            loop_alpha=loop_alpha,
            loop_size=loop_size
        )
    
    # 4. 保存图片（直接保存 PNG + SVG）
    fig.savefig(f"{output_path}.png", dpi=dpi, bbox_inches='tight')
    fig.savefig(f"{output_path}.svg", dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"图片已保存: {output_path}.png/.svg")



def plot_loop_comparison(
    loops1: pd.DataFrame,
    loops2: pd.DataFrame,
    resolution: int,
    merge_distance: int = None,
    figsize: Tuple[float, float] = (12, 6),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Compare two sets of chromatin loops.
    
    This creates a Venn-like comparison showing:
    - Shared loops between two samples
    - Unique loops in each sample
    
    Parameters
    ----------
    loops1 : pd.DataFrame
        First set of loops
    loops2 : pd.DataFrame
        Second set of loops
    resolution : int
        Bin resolution (bp)
    merge_distance : int, optional
        Distance for merging nearby loops (default: resolution)
    figsize : tuple
        Figure size in inches
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    """
    if merge_distance is None:
        merge_distance = resolution
    
    # Calculate statistics
    n_loops1 = len(loops1)
    n_loops2 = len(loops2)
    
    # Find overlapping loops
    shared = 0
    for _, loop1 in loops1.iterrows():
        for _, loop2 in loops2.iterrows():
            if (abs(loop1['start1'] - loop2['start1']) < merge_distance and
                abs(loop1['start2'] - loop2['start2']) < merge_distance):
                shared += 1
                break
    
    n_shared = shared
    n_unique1 = n_loops1 - n_shared
    n_unique2 = n_loops2 - n_shared
    
    # Create figure with bar chart
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Left: Venn-like diagram using circles
    ax1 = axes[0]
    circle1 = plt.Circle((0.4, 0.5), 0.3, color='blue', alpha=0.3)
    circle2 = plt.Circle((0.6, 0.5), 0.3, color='red', alpha=0.3)
    ax1.add_patch(circle1)
    ax1.add_patch(circle2)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect('equal')
    ax1.axis('off')
    
    # Add labels
    ax1.text(0.25, 0.5, f'{n_unique1}\nunique', ha='center', va='center', fontsize=12)
    ax1.text(0.75, 0.5, f'{n_unique2}\nunique', ha='center', va='center', fontsize=12)
    ax1.text(0.5, 0.5, f'{n_shared}\nshared', ha='center', va='center', fontsize=12, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Right: Bar chart
    ax2 = axes[1]
    categories = ['Sample 1', 'Sample 2', 'Shared']
    values = [n_unique1, n_unique2, n_shared]
    colors = ['blue', 'red', 'purple']
    ax2.bar(categories, values, color=colors, alpha=0.6)
    ax2.set_ylabel('Number of loops', fontsize=10)
    
    # Add value labels on bars
    for i, v in enumerate(values):
        ax2.text(i, v + max(values) * 0.02, str(v), ha='center', fontsize=10)
    
    # Add title
    if title:
        fig.suptitle(title)
    else:
        fig.suptitle('Loop Comparison')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig
