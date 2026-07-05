"""
Plot style utilities for cfizz visualization.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from typing import Optional


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
    
    Args:
        colors (list): 颜色列表，默认为 ['blue', 'white', 'red']
        bad_color (str): 无效值的颜色，默认为 '#EEEEEE'
        name (str): 颜色映射的名称，默认为 'custom_coolwarm'
    
    Returns:
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
    # 计算所有绝对尺寸（单位：cm）
    heatmap_height = plot_size
    heatmap_width = plot_size
    
    # 柱状图尺寸
    bar_height = heatmap_height * bar_height_ratio if has_bar else 0
    
    # 子图间隔
    gap = plot_size * 0.05
    
    # 颜色条尺寸
    colorbar_width = plot_size * 0.05
    
    # 计算行列数
    if n_cols is None:
        n_cols = n_plots
    if n_rows is None:
        n_rows = (n_plots + n_cols - 1) // n_cols
    
    # 计算总宽度
    total_width = margin_left + n_cols * heatmap_width + (n_cols - 1) * gap + margin_right
    
    # 计算总高度
    if has_bar:
        total_height = margin_bottom + n_rows * (bar_height + heatmap_height + gap) - gap
    else:
        total_height = margin_bottom + n_rows * (heatmap_height + gap) - gap
    
    # 计算相对位置（比例）
    plot_width = heatmap_width / total_width
    plot_height = heatmap_height / total_height
    bar_height_relative = bar_height / total_height if has_bar else 0
    margin_left_relative = margin_left / total_width
    margin_right_relative = margin_right / total_width
    margin_bottom_relative = margin_bottom / total_height
    gap_relative = gap / total_width
    colorbar_left = (margin_left + n_cols * (heatmap_width + gap)) / total_width
    colorbar_width_relative = colorbar_width / total_width
    
    # 转换为英寸
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


def setup_axes(ax):
    """
    设置坐标轴和边框
    
    参数:
        ax: matplotlib的axes对象
    """
    # 设置坐标轴
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


def print_coordinate(pos):
    """
    格式化坐标显示
    
    参数:
        pos: 位置值（bp）
    
    返回:
        格式化后的字符串，根据数值大小自动选择单位：
        - 大于等于1M时使用M为单位
        - 小于1M时使用k为单位
    """
    if pos >= 1e6:  # 大于等于1M
        if pos % 1000000 == 0:  # 如果是整M
            return f"{pos//1000000}M"
        else:  # 如果不是整M
            return f"{pos/1000000:.2f}M"
    else:  # 小于1M
        if pos % 1000 == 0:  # 如果是整k
            return f"{pos//1000}k"
        else:  # 如果不是整k
            return f"{pos/1000:.1f}k"


def format_coordinate(pos):
    """兼容字符串和数字的坐标格式化"""
    if isinstance(pos, str):
        return pos
    try:
        pos = int(pos)
        if pos % 1000000 == 0:
            return f'{pos // 1000000}M'
        return f'{pos / 1000000:.2f}M'
    except (ValueError, TypeError):
        return str(pos)


def add_coordinate_labels(ax, xmin, xmax, ymin, ymax, start, end, chrom, fontsize=5, offset_ratio=0.02):
    """添加坐标标签到热图"""
    offset = offset_ratio * (xmax - xmin)
    
    # 添加底部坐标标签
    ax.text(xmin, ymin + offset, print_coordinate(start), va='top', ha='left', fontsize=fontsize)
    ax.text(xmax, ymin + offset, print_coordinate(end), va='top', ha='right', fontsize=fontsize)
    
    # 添加左侧坐标标签
    ax.text(-offset, ymax, print_coordinate(start),
            rotation=90, va='top', ha='right', fontsize=fontsize)
    ax.text(-offset, ymin, print_coordinate(end),
            rotation=90, va='bottom', ha='right', fontsize=fontsize)
    
    # 添加染色体标签
    ax.text((xmin + xmax) / 2, ymin + 2 * offset, chrom,
            va='top', ha='center', fontsize=fontsize)
    ax.text(-2 * offset, (ymin + ymax) / 2, chrom,
            rotation=90, va='center', ha='right', fontsize=fontsize)


def setup_diff_colorbar(
    fig, sc, vmin, vmax,
    colorbar_left, colorbar_bottom, colorbar_width, colorbar_height,
    max_label, min_label, label_text, tick_fontsize=5, label_fontsize=5
):
    """设置差异矩阵的颜色条"""
    # 创建颜色条轴
    cb_ax = fig.add_axes([colorbar_left, colorbar_bottom, colorbar_width, colorbar_height])
    
    # 添加颜色条
    cbar = fig.colorbar(sc, cax=cb_ax)
    cbar.set_ticks([])
    cbar.set_ticklabels([])
    cbar.ax.tick_params(left=False, right=False, labelleft=False, labelright=False)
    cbar.set_label("")
    
    # 添加标签
    cb_ax.text(0.5, 1.1, max_label, ha='center', va='bottom', fontsize=tick_fontsize, rotation=0, transform=cb_ax.transAxes)
    cb_ax.text(0.5, -0.1, min_label, ha='center', va='top', fontsize=tick_fontsize, rotation=0, transform=cb_ax.transAxes)
    cb_ax.text(0.5, -0.5, label_text, ha='center', va='top', fontsize=label_fontsize, rotation=90, transform=cb_ax.transAxes)


def setup_colorbar(fig, sc, vmin, vmax, balance=False,
                  colorbar_left=None, colorbar_bottom=None,
                  colorbar_width=None, colorbar_height=None,
                  label_fontsize=5, tick_fontsize=5,
                  label=None):
    """设置颜色条"""
    cax = fig.add_axes([colorbar_left, colorbar_bottom, colorbar_width, colorbar_height])
    
    if balance:
        format_str = '%.3g'
        default_label = "normalized contacts"
    else:
        format_str = '%.0f'
        default_label = "contacts"
    
    cbar = fig.colorbar(sc, cax=cax, format=format_str)
    cbar.set_ticks([vmin, (vmin+vmax)/2, vmax])
    cax.tick_params(labelsize=tick_fontsize, length=1)
    
    label_text = label if label is not None else default_label
    cax.text(0.5, -0.1, label_text,
             ha='center', va='top',
             fontsize=label_fontsize,
             rotation=90,
             transform=cax.transAxes)


def save_figure(fig, output_prefix):
    """保存图片"""
    fig.savefig(f'{output_prefix}.svg', format='svg', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_prefix}.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"图片已保存: {output_prefix}.svg/.png")


def log2_and_mask(matrix):
    """对矩阵进行 log2 转换并遮罩无效值"""
    matrix = np.array(matrix, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        log_matrix = np.log2(matrix)
    log_matrix = np.ma.masked_invalid(log_matrix)
    return log_matrix


def add_bar_coordinates(ax, label="", show_labels=True, fontsize=5, xmin=None, xmax=None, ymin=None, ymax=None):
    """
    为柱状图添加 y 轴坐标和标签
    
    关键：根据数据范围计算对称的刻度值，并向上取0.5的整数倍
    """
    if xmin is None or xmax is None or ymin is None or ymax is None:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
    
    # 关闭刻度和边框
    ax.tick_params(axis='both', which='both', length=0,
                  labelbottom=False, labeltop=False,
                  labelleft=False, labelright=False)
    for spine in ['right', 'top', 'bottom']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(0.1)
    
    # 计算对称的刻度值范围
    abs_max = max(abs(ymin), abs(ymax))
    # 向上取0.5的整数倍
    abs_max = np.ceil(abs_max * 2) / 2
    
    # 只使用最大值和最小值作为刻度
    yticks = [-abs_max, abs_max]
    
    # 设置y轴范围为对称的范围
    ax.set_ylim(-abs_max, abs_max)
    
    if show_labels:
        for y in yticks:
            ax.text(0, y, f'{y:.1f}',
                    ha='right', va='center',
                    fontsize=fontsize, fontname='Arial',
                    transform=ax.get_yaxis_transform())
        if label:
            ax.text(0, 0, label,
                    ha='right', va='center',
                    fontsize=fontsize, fontname='Arial',
                    rotation=90,
                    transform=ax.get_yaxis_transform())
    
    return ymin, ymax


def add_rotated_coordinate_labels(ax, xmin, xmax, ymin, start, end, chrom, fontsize=5, offset_ratio=0.02):
    """为45度旋转的热图添加坐标标签（只添加底部标签）"""
    offset = offset_ratio * (xmax - xmin)
    
    # 添加底部坐标标签
    ax.text(xmin, ymin - offset, print_coordinate(start),
            va='top', ha='left', fontsize=fontsize)
    ax.text(xmax, ymin - offset, print_coordinate(end),
            va='top', ha='right', fontsize=fontsize)
    
    # 添加染色体标签
    ax.text((xmin + xmax) / 2, ymin - 2 * offset, chrom, va='top', ha='center', fontsize=fontsize)


def mask_diagonal(matrix, mask_width=0):
    """对矩阵的主对角线进行 mask 处理"""
    if mask_width <= 0:
        return matrix
    
    masked_matrix = matrix.copy()
    n = masked_matrix.shape[0]
    
    for i in range(n):
        for j in range(max(0, i-mask_width), min(n, i+mask_width+1)):
            masked_matrix[i, j] = np.nan
    
    return masked_matrix


def log10_and_mask(matrix):
    """对矩阵进行 log10 转换并遮罩无效值"""
    matrix = np.array(matrix, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        log_matrix = np.log10(matrix)
    log_matrix = np.ma.masked_invalid(log_matrix)
    return log_matrix


def generate_output_filename(prefix, chrom=None, start=None, end=None, resolution=None):
    """生成输出文件名"""
    parts = [prefix]
    if chrom:
        parts.append(chrom)
    if start is not None and end is not None:
        parts.append(f"{start}-{end}")
    if resolution:
        parts.append(f"{resolution}")
    return "_".join(parts)


def save_figure_multi_format(fig, output_prefix, formats=None):
    """保存多种格式的图片"""
    if formats is None:
        formats = ['svg', 'png']
    
    for fmt in formats:
        if fmt == 'png':
            fig.savefig(f'{output_prefix}.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
        elif fmt == 'svg':
            fig.savefig(f'{output_prefix}.svg', format='svg', dpi=300, bbox_inches='tight')
        elif fmt == 'pdf':
            fig.savefig(f'{output_prefix}.pdf', format='pdf', bbox_inches='tight')
    
    print(f"图片已保存: {output_prefix}.{'/'.join(formats)}")


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
    """计算带 track 的热图布局"""
    heatmap_height = plot_size
    heatmap_width = plot_size
    gap = plot_size * 0.05
    colorbar_width = plot_size * 0.05
    
    # 计算 track 高度
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
    
    # 行列数
    if n_cols is None:
        n_cols = min(8, n_plots)
    if n_rows is None:
        n_rows = (n_plots + n_cols - 1) // n_cols
    
    total_width = margin_left + n_cols * heatmap_width + (n_cols - 1) * gap + margin_right
    total_height = heatmap_height + margin_bottom + total_track_height
    
    # 相对位置
    plot_width = heatmap_width / total_width
    plot_height = heatmap_height / total_height
    margin_left_relative = margin_left / total_width
    margin_bottom_relative = margin_bottom / total_height
    
    # Track 位置
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
        'total_width': total_width,
        'total_height': total_height,
        'margin_left': margin_left,
        'margin_right': margin_right,
        'margin_bottom': margin_bottom,
        'total_track_height': total_track_height,
        'plot_width': plot_width,
        'plot_height': plot_height,
        'margin_left_relative': margin_left_relative,
        'margin_bottom_relative': margin_bottom_relative,
        'colorbar_left': (margin_left + n_cols * (heatmap_width + gap)) / total_width,
        'colorbar_width_relative': colorbar_width / total_width,
        'track_heights_relative': track_heights_relative,
        'track_positions': track_positions,
        'fig_width': fig_width,
        'fig_height': fig_height,
        'n_cols': n_cols,
        'n_rows': n_rows,
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
    """计算45度旋转热图的布局"""
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
        'last_heatmap_bottom': last_heatmap_bottom,
        'colorbar_height': last_colorbar_height
    }


def setup_horizontal_colorbar(
    fig, sc, vmin, vmax,
    ax,
    fontsize: int = 5,
    label: str = "contacts",
    colorbar_left: float = None,
    colorbar_bottom: float = None,
    colorbar_width: float = None,
    colorbar_height: float = None,
    balance: bool = False
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