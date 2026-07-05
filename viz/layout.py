"""
Layout calculation functions for cfizz package.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from cfizz.utils.coordinates import print_coordinate, format_coordinate
from typing import Optional
from matplotlib.colors import LinearSegmentedColormap


def setup_plot_style():
    """
    设置统一的绘图样式，包括字体和字号
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5


def setup_ratio_colormap(
    colors: list = None,
    bad_color: str = '#EEEEEE',
    name: str = 'custom_coolwarm'
) -> LinearSegmentedColormap:
    """
    设置比值矩阵（如O/E）的颜色映射
    
    参数:
        colors (list): 颜色列表，默认为 ['blue', 'white', 'red']
        bad_color (str): 无效值的颜色，默认为 '#EEEEEE'
        name (str): 颜色映射的名称，默认为 'custom_coolwarm'
    
    返回:
        LinearSegmentedColormap: 自定义的颜色映射对象
    """
    if colors is None:
        colors = ['blue', 'white', 'red']
    
    cmap = LinearSegmentedColormap.from_list(name, colors)
    cmap.set_bad(color=bad_color)
    
    return cmap


def calculate_heatmap_layout(
    n_plots: int,
    plot_size: float,
    has_bar: bool = False,
    bar_height_ratio: float = 0.3,
    margin_left: float = 0.8,
    margin_right: float = 1.6,
    margin_bottom: float = 0.8,
    n_cols: Optional[int] = None,
    n_rows: Optional[int] = None
) -> dict:
    """
    计算热图布局参数
    
    Args:
        n_plots: 子图数量
        plot_size: 单个子图的基础尺寸（单位：cm）
        has_bar: 是否包含柱状图
        bar_height_ratio: 柱状图高度与热图高度的比例
        margin_left: 左侧边距（厘米）
        margin_right: 右侧边距（厘米）
        margin_bottom: 底部边距（厘米）
        n_cols: 列数，如果为None则自动计算
        n_rows: 行数，如果为None则自动计算
    
    Returns:
        dict: 包含布局参数的字典
    """
    # 1. 计算所有绝对尺寸（单位：cm）
    heatmap_height = plot_size
    heatmap_width = plot_size
    bar_height = heatmap_height * bar_height_ratio if has_bar else 0
    gap = plot_size * 0.05
    colorbar_width = plot_size * 0.05
    
    # 2. 计算总尺寸
    if n_cols is None:
        n_cols = n_plots
    if n_rows is None:
        n_rows = (n_plots + n_cols - 1) // n_cols
    
    total_width = margin_left + n_cols * heatmap_width + (n_cols - 1) * gap + margin_right
    
    if has_bar:
        total_height = margin_bottom + n_rows * (bar_height + heatmap_height + gap) - gap
    else:
        total_height = margin_bottom + n_rows * (heatmap_height + gap) - gap
    
    # 3. 计算相对位置
    plot_width = heatmap_width / total_width
    plot_height = heatmap_height / total_height
    bar_height_relative = bar_height / total_height if has_bar else 0
    margin_left_relative = margin_left / total_width
    margin_right_relative = margin_right / total_width
    margin_bottom_relative = margin_bottom / total_height
    gap_relative = gap / total_width
    colorbar_left = (margin_left + n_cols * (heatmap_width + gap)) / total_width
    colorbar_width_relative = colorbar_width / total_width
    
    # 4. 转换为英寸
    cm_to_inch = 1/2.54
    fig_width = total_width * cm_to_inch
    fig_height = total_height * cm_to_inch
    
    return {
        'plot_size': plot_size,
        'heatmap_height': heatmap_height,
        'heatmap_width': heatmap_width,
        'bar_height': bar_height,
        'gap': gap,
        'colorbar_width': colorbar_width,
        'total_width': total_width,
        'total_height': total_height,
        'margin_left': margin_left,
        'margin_right': margin_right,
        'margin_bottom': margin_bottom,
        'plot_width': plot_width,
        'plot_height': plot_height,
        'bar_height_relative': bar_height_relative,
        'margin_left_relative': margin_left_relative,
        'margin_right_relative': margin_right_relative,
        'margin_bottom_relative': margin_bottom_relative,
        'gap_relative': gap_relative,
        'colorbar_left': colorbar_left,
        'colorbar_width_relative': colorbar_width_relative,
        'fig_width': fig_width,
        'fig_height': fig_height,
        'n_cols': n_cols,
        'n_rows': n_rows,
        'hspace': gap_relative,
        'vspace': gap_relative
    }


def calculate_heatmap_with_tracks_layout(
    n_plots: int,
    plot_size: float,
    tracks_config: list,
    track_gap: float = 0.1,
    margin_left: float = 0.8,
    margin_right: float = 1.6,
    margin_bottom: float = 0.8,
    n_cols: Optional[int] = None,
    n_rows: Optional[int] = None
) -> dict:
    """gap只用于列间，纵向无gap。track从最下方依次向上累加"""
    heatmap_height = plot_size
    heatmap_width = plot_size
    gap = plot_size * 0.05
    colorbar_width = plot_size * 0.05
    
    track_heights = []
    total_track_height = 0
    for track in tracks_config:
        if 'height_cm' in track:
            track_height = track['height_cm']
        else:
            track_ratio = track.get('height_ratio', 0.2)
            track_height = heatmap_height * track_ratio
        track_heights.append(track_height)
        total_track_height += track_height
    
    if n_cols is None:
        n_cols = min(8, n_plots)
    if n_rows is None:
        n_rows = (n_plots + n_cols - 1) // n_cols
    
    total_width = margin_left + n_cols * heatmap_width + (n_cols - 1) * gap + margin_right
    total_height = heatmap_height + margin_bottom + total_track_height
    
    plot_width = heatmap_width / total_width
    plot_height = heatmap_height / total_height
    margin_left_relative = margin_left / total_width
    margin_right_relative = margin_right / total_width
    margin_bottom_relative = margin_bottom / total_height
    colorbar_left = (margin_left + n_cols * (heatmap_width + gap)) / total_width
    colorbar_width_relative = colorbar_width / total_width
    track_heights_relative = [h / total_height for h in track_heights]
    heatmap_bottom = (total_height - heatmap_height) / total_height
    
    track_positions = []
    current_y = 0.0
    for i, track in enumerate(reversed(tracks_config)):
        track_height_rel = track_heights_relative[len(tracks_config)-1-i]
        track_positions.append({
            'left': margin_left_relative,
            'bottom': current_y,
            'width': plot_width,
            'height': track_height_rel,
            'height_cm': track_heights[len(tracks_config)-1-i],
            'type': track.get('type', 'unknown'),
            'config': track
        })
        current_y += track_height_rel
    track_positions = list(reversed(track_positions))
    
    cm_to_inch = 1/2.54
    fig_width = total_width * cm_to_inch
    fig_height = total_height * cm_to_inch
    
    return {
        'plot_size': plot_size,
        'heatmap_height': heatmap_height,
        'heatmap_width': heatmap_width,
        'gap': gap,
        'colorbar_width': colorbar_width,
        'total_width': total_width,
        'total_height': total_height,
        'margin_left': margin_left,
        'margin_right': margin_right,
        'margin_bottom': margin_bottom,
        'total_track_height': total_track_height,
        'track_gap': 0.0,
        'plot_width': plot_width,
        'plot_height': plot_height,
        'margin_left_relative': margin_left_relative,
        'margin_right_relative': margin_right_relative,
        'margin_bottom_relative': margin_bottom_relative,
        'colorbar_left': colorbar_left,
        'colorbar_width_relative': colorbar_width_relative,
        'track_heights_relative': track_heights_relative,
        'track_gap_relative': 0.0,
        'track_positions': track_positions,
        'fig_width': fig_width,
        'fig_height': fig_height,
        'n_cols': n_cols,
        'n_rows': n_rows,
        'hspace': gap / total_width,
        'vspace': gap / total_width,
        'n_tracks': len(tracks_config),
        'heatmap_bottom': heatmap_bottom
    }


def calculate_rotated_heatmap_layout(
    n_plots: int,
    plot_size: float = 4,
    margin_bottom: float = None,
    triangle_ratio: float = 1.0,
    colorbar_height_cm: float = None,
    colorbar_width_cm: float = 1.5,
    offset_ratio: float = 0.02
) -> dict:
    """计算45度旋转热图的布局参数"""
    diagonal_length = plot_size * np.sqrt(2)
    n_cols = 1
    n_rows = min(8, n_plots)
    gap = diagonal_length * 0.5 * triangle_ratio * 0.05
    
    if margin_bottom is None:
        margin_bottom = max(0.3, plot_size * 0.1)
    total_width = diagonal_length
    total_height = margin_bottom + n_rows * (diagonal_length * 0.5 * triangle_ratio) + (n_rows - 1) * gap
    
    heatmap_width = diagonal_length / total_width
    heatmap_height = (diagonal_length * 0.5 * triangle_ratio) / total_height
    margin_bottom_relative = margin_bottom / total_height
    gap_relative = gap / total_width
    
    cm_to_inch = 1/2.54
    fig_width = total_width * cm_to_inch
    fig_height = total_height * cm_to_inch
    
    subplot_positions = []
    for i in range(n_plots):
        row = i
        heatmap_left = 0
        heatmap_bottom = margin_bottom_relative + (n_rows - row - 1) * (heatmap_height + gap_relative)
        colorbar_left = 0.15
        colorbar_bottom = 0
        
        if colorbar_height_cm is None:
            colorbar_height_cm = min(0.2, (diagonal_length * 0.5 * triangle_ratio) * 0.8)
        colorbar_height_relative = colorbar_height_cm / total_height
        colorbar_width_relative = colorbar_width_cm / total_width
        
        subplot_positions.append({
            'heatmap_left': heatmap_left,
            'heatmap_bottom': heatmap_bottom,
            'heatmap_width': heatmap_width,
            'heatmap_height': heatmap_height,
            'colorbar_left': colorbar_left,
            'colorbar_bottom': colorbar_bottom,
            'colorbar_width': colorbar_width_relative,
            'colorbar_height': colorbar_height_relative,
            'triangle_ratio': triangle_ratio
        })
    
    last_heatmap_bottom = subplot_positions[-1]['heatmap_bottom']
    last_colorbar_height = colorbar_height_relative
    
    return {
        'plot_size': plot_size,
        'diagonal_length': diagonal_length,
        'total_width': total_width,
        'total_height': total_height,
        'margin_bottom': margin_bottom,
        'colorbar_height_cm': colorbar_height_cm,
        'colorbar_width_cm': colorbar_width_cm,
        'heatmap_width': heatmap_width,
        'heatmap_height': heatmap_height,
        'margin_bottom_relative': margin_bottom_relative,
        'gap_relative': gap_relative,
        'fig_width': fig_width,
        'fig_height': fig_height,
        'n_cols': n_cols,
        'n_rows': n_rows,
        'triangle_ratio': triangle_ratio,
        'subplot_positions': subplot_positions,
        'hspace': gap_relative,
        'vspace': gap_relative,
        'last_heatmap_bottom': last_heatmap_bottom,
        'colorbar_height': last_colorbar_height
    }


def log2_and_mask(matrix):
    """对矩阵取log2并mask掉无效值"""
    if matrix is None:
        raise ValueError("输入矩阵为None，无法进行处理")
    matrix = matrix.astype(float)
    safe_matrix = np.where(matrix > 0, matrix, np.nan)
    log2_m = np.log2(safe_matrix)
    log2_m[~np.isfinite(log2_m)] = np.nan
    return np.ma.masked_invalid(log2_m)


def log10_and_mask(matrix):
    """对矩阵取log10并mask掉无效值"""
    matrix = matrix.astype(float)
    log10_m = np.full_like(matrix, np.nan)
    valid_mask = (~np.isnan(matrix)) & (matrix > 0)
    log10_m[valid_mask] = np.log10(matrix[valid_mask])
    return np.ma.masked_invalid(log10_m)


def setup_axes(ax):
    """设置坐标轴和边框"""
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.tick_params(
        axis='both',
        bottom=False, top=False, left=False, right=False,
        labelbottom=False, labeltop=False, labelleft=False, labelright=False
    )
    for spine in ['right', 'top', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(0)
    return xmin, xmax, ymin, ymax


def add_coordinate_labels(ax, xmin, xmax, ymin, ymax, start_pos, end_pos, chrom, fontsize=5, offset_ratio=0.02):
    """添加坐标标签到热图"""
    offset = offset_ratio * (xmax - xmin)
    ax.text(xmin, ymin + offset, print_coordinate(start_pos), va='top', ha='left', fontsize=fontsize)
    ax.text(xmax, ymin + offset, print_coordinate(end_pos), va='top', ha='right', fontsize=fontsize)
    ax.text(-offset, ymax, print_coordinate(start_pos), 
            rotation=90, va='top', ha='right', fontsize=fontsize)
    ax.text(-offset, ymin, print_coordinate(end_pos), 
            rotation=90, va='bottom', ha='right', fontsize=fontsize)
    ax.text((xmin + xmax) / 2, ymin + 2 * offset, chrom, 
            va='top', ha='center', fontsize=fontsize)
    ax.text(-2 * offset, (ymin + ymax) / 2, chrom, 
            rotation=90, va='center', ha='right', fontsize=fontsize)


def add_rotated_coordinate_labels(ax, xmin, xmax, ymin, start_pos, end_pos, chrom, fontsize=5, offset_ratio=0.05):
    """为45度旋转的热图添加坐标标签（只添加底部标签）"""
    offset = offset_ratio * (xmax - xmin)
    ax.text(xmin, ymin - offset, print_coordinate(start_pos), 
            va='top', ha='left', fontsize=fontsize)
    ax.text(xmax, ymin - offset, print_coordinate(end_pos), 
            va='top', ha='right', fontsize=fontsize)
    ax.text((xmin + xmax) / 2, ymin - 2 * offset, chrom, va='top', ha='center', fontsize=fontsize)


def setup_colorbar(fig, sc, vmin, vmax, balance=False, 
                  color_scale='linear',
                  colorbar_left=None, colorbar_bottom=None, 
                  colorbar_width=None, colorbar_height=None,
                  label_fontsize=5, tick_fontsize=5,
                  label=None,
                  label_orientation='vertical'):
    """
    设置颜色条
    
    Parameters:
        fig: matplotlib.figure.Figure, 图形对象
        sc: matplotlib.cm.ScalarMappable, 颜色映射对象
        vmin, vmax: float, 颜色条范围
        balance: bool, 是否使用平衡矩阵
        color_scale: str, 'linear' 或 'log'，决定中间刻度的计算方式
        colorbar_left, colorbar_bottom, colorbar_width, colorbar_height: float, 颜色条位置和大小
        label_fontsize: int, 标签字体大小
        tick_fontsize: int, 刻度字体大小
        label: str, 自定义标签文本
        label_orientation: str, 'vertical' 或 'horizontal'
            - 'vertical': label 旋转 90°（适合普通正方形热图）
            - 'horizontal': label 水平放置（适合 45 度旋转热图）
    """
    # 设置全局字体为 Arial
    plt.rcParams['font.family'] = 'Arial'
    
    cax = fig.add_axes([colorbar_left, colorbar_bottom, colorbar_width, colorbar_height])
    
    if balance:
        format_str = '%.3g'
        default_label = "normalized contacts"
    else:
        format_str = '%.0f'
        default_label = "contacts"
    
    cbar = fig.colorbar(sc, cax=cax, format=format_str)
    
    # 根据 color_scale 计算中间刻度
    if color_scale == 'log':
        # Log 缩放：使用几何平均数（对数空间中点）
        mid = np.sqrt(vmin * vmax)
    else:
        # 线性缩放：使用算术平均数
        mid = (vmin + vmax) / 2
    
    # 只设置3个主刻度，关闭次刻度
    cbar.set_ticks([vmin, mid, vmax])
    cax.minorticks_off()  # 关闭次刻度
    cax.tick_params(labelsize=tick_fontsize, length=1)
    
    # 根据 label_orientation 设置 label 的位置和旋转
    label_text = label if label is not None else default_label
    if label_orientation == 'horizontal':
        # 水平 label（适合 45 度旋转热图）- 左对齐
        cax.text(0, -0.3, label_text,
                 ha='left', va='top',
                 fontsize=label_fontsize,
                 rotation=0,
                 transform=cax.transAxes)
    else:
        # 垂直 label（适合普通正方形热图）
        cax.text(0.5, -0.1, label_text,
                 ha='center', va='top',
                 fontsize=label_fontsize,
                 rotation=90,
                 transform=cax.transAxes)


def setup_horizontal_colorbar(
    fig, sc, vmin, vmax,
    colorbar_left, colorbar_bottom, colorbar_width, colorbar_height,
    fontsize=5, label="contacts", balance=False
):
    """为45度旋转的热图设置水平颜色条"""
    cax = fig.add_axes([colorbar_left, colorbar_bottom, colorbar_width, colorbar_height])
    
    if balance:
        format_str = '%.3g'
    else:
        format_str = '%.0f'

    cbar = fig.colorbar(sc, cax=cax, orientation='horizontal', format=format_str)
    cbar.outline.set_linewidth(0.1)
    
    ticks = [vmin, vmax/2, vmax]
    cbar.set_ticks(ticks)
    cbar.ax.minorticks_off()
    cbar.ax.tick_params(
        axis='both',
        bottom=True, top=False, left=True, right=True,
        labelbottom=True, labeltop=False, labelleft=True, labelright=True,
        labelsize=fontsize, length=0.2, width=0.1, pad=0.25,
        rotation=15
    )
    
    for label in cbar.ax.get_xticklabels():
        label.set_ha('right')
        label.set_va('top')
    
    cbar.set_label("")


def mask_diagonal(matrix, mask_width=0):
    """对矩阵的主对角线进行mask处理"""
    if mask_width <= 0:
        return matrix
    
    masked_matrix = matrix.copy()
    n = masked_matrix.shape[0]
    
    for i in range(n):
        for j in range(max(0, i-mask_width), min(n, i+mask_width+1)):
            masked_matrix[i, j] = 0
            masked_matrix[j, i] = 0
    
    return masked_matrix


def save_figure(fig, output_path, dpi=300):
    """保存图片"""
    import os
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=dpi)
    plt.close()


def generate_output_filename(
    output_dir: str,
    prefix: str = "plot",
    chrom: str = None,
    region: tuple = None,
    resolution: int = None,
    balance: bool = None,
    suffix: str = None,
    ext: str = "png"
) -> str:
    """
    生成有规律的文件名
    
    命名格式: {prefix}_{chrom}_{start}-{end}_{resolution}_{balance}.{ext}
    
    Examples
    --------
    >>> generate_output_filename("/output", chrom="chr1", region=(0, 50000000), resolution=1000000, balance=False)
    '/output/plot_chr1_0M-50M_1Mb_raw.png'
    
    >>> generate_output_filename("/output", chrom="chr3", region=(100000000, 150000000), resolution=1000000, balance=True)
    '/output/plot_chr3_100M-150M_1Mb_balanced.svg'
    """
    parts = [prefix]
    
    if chrom:
        # chr1 -> 1, chrX -> X
        parts.append(chrom.replace('chr', ''))
    
    if region:
        start, end = region
        # 智能格式化大小
        if end >= 1000000000:
            start_str = f"{start//1000000000}G"
            end_str = f"{end//1000000000}G"
        elif end >= 1000000:
            start_str = f"{start//1000000}M"
            end_str = f"{end//1000000}M"
        else:
            start_str = f"{start//1000}K"
            end_str = f"{end//1000}K"
        parts.append(f"{start_str}-{end_str}")
    
    if resolution:
        if resolution >= 1000000:
            parts.append(f"{resolution//1000000}Mb")
        elif resolution >= 1000:
            parts.append(f"{resolution//1000}kb")
    
    if balance is not None:
        parts.append("raw" if not balance else "balanced")
    
    if suffix:
        parts.append(suffix)
    
    filename = "_".join(parts) + f".{ext}"
    return os.path.join(output_dir, filename)


def save_figure_multi_format(
    fig,
    output_path: str,
    dpi: int = 300,
    formats: list = None
):
    """
    同时保存多种格式的图片（PNG + SVG）
    
    Parameters
    ----------
    fig : plt.Figure
        Matplotlib figure 对象
    output_path : str
        输出路径（不含扩展名）
    dpi : int
        位图格式的分辨率
    formats : list
        要保存的格式列表，默认 ["png", "svg"]
        
    Returns
    -------
    saved_files : list
        保存的文件路径列表
    """
    if formats is None:
        formats = ["png", "svg"]
    
    import os
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    
    # 分离文件名和路径
    base_path = os.path.splitext(output_path)[0]
    
    saved_files = []
    for fmt in formats:
        save_path = f"{base_path}.{fmt}"
        if fmt == "svg":
            fig.savefig(save_path, format="svg", bbox_inches='tight')
        elif fmt == "png":
            fig.savefig(save_path, format="png", dpi=dpi, bbox_inches='tight')
        elif fmt == "jpg" or fmt == "jpeg":
            fig.savefig(save_path, format="jpeg", dpi=dpi, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
        saved_files.append(save_path)
    
    plt.close(fig)
    return saved_files
