"""
Saddle Plot (Compartment Cumulative Analysis) End-to-End Module.

This module provides functions for generating saddle plots from Hi-C data,
useful for visualizing chromatin compartment (A/B) interaction patterns.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker
from cytoolz import merge
import pandas as pd
import cooltools
import cooler
import pickle
from pathlib import Path
from typing import List, Dict, Optional, Any
import warnings

# 设置numpy的警告过滤器
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in divide')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in log2')

# === 从 cfizz 复制 import ===
from cfizz.viz.layout import (
    setup_plot_style,
    setup_axes,
    save_figure,
    save_figure_multi_format,
    log2_and_mask,
    setup_ratio_colormap,
)
from cfizz.viz.plot_style import setup_diff_colorbar


class MinOneMaxFormatter(ticker.LogFormatter):
    """自定义对数格式化器，用于颜色条标签"""
    def __init__(self, vmin, vmax):
        super().__init__()
        self.vmin = vmin
        self.vmax = vmax
        
    def set_locs(self, locs=None):
        self._sublabels = set([self.vmin % 10 * 10, self.vmax % 10, 1])

    def __call__(self, x, pos=None):
        if x not in [self.vmin, 1, self.vmax]:
            return ""
        else:
            return "{x:g}".format(x=x)


def calculate_saddle_strength(interaction_sum: np.ndarray, interaction_count: np.ndarray) -> np.ndarray:
    """
    计算鞍形强度

    参数:
        interaction_sum (np.ndarray): 相互作用总和矩阵
        interaction_count (np.ndarray): 相互作用计数矩阵

    返回:
        np.ndarray: 鞍形强度数组
    """
    m, n = interaction_sum.shape
    if m != n:
        raise ValueError("`saddledata` should be square.")

    ratios = np.zeros(n)
    for k in range(1, n):
        intra_sum = np.nansum(interaction_sum[0:k, 0:k]) + np.nansum(interaction_sum[n - k : n, n - k : n])
        intra_count = np.nansum(interaction_count[0:k, 0:k]) + np.nansum(interaction_count[n - k : n, n - k : n])
        intra = intra_sum / intra_count

        inter_sum = np.nansum(interaction_sum[0:k, n - k : n]) + np.nansum(interaction_sum[n - k : n, 0:k])
        inter_count = np.nansum(interaction_count[0:k, n - k : n]) + np.nansum(interaction_count[n - k : n, 0:k])
        inter = inter_sum / inter_count

        ratios[k] = intra / inter
    return ratios


def plot_saddle(
    track: pd.DataFrame,
    saddledata: np.ndarray,
    n_bins: int = 98,
    vrange: tuple = None,
    cmap: str = "coolwarm",
    scale: str = "log",
    vmin: float = 0.5,
    vmax: float = 2,
    color: str = None,
    title: str = None,
    xlabel: str = None,
    ylabel: str = None,
    clabel: str = None,
    fig: plt.Figure = None,
    fig_kws: dict = None,
    heatmap_kws: dict = None,
    margin_kws: dict = None,
    cbar_kws: dict = None,
    subplot_spec: GridSpec = None,
    output_file: str = None,
    dpi: int = 800
) -> dict:
    """
    绘制鞍形图

    参数:
        track (pd.DataFrame): 包含特征向量的数据框
        saddledata (np.ndarray): 鞍形矩阵数据
        n_bins (int): 分箱数量，默认98
        vrange (tuple, optional): 值范围
        cmap (str, optional): 颜色映射
        scale (str, optional): 缩放方式 ('log' 或 'linear')
        vmin (float, optional): 最小值
        vmax (float, optional): 最大值
        color (str, optional): 颜色
        title (str, optional): 标题
        xlabel (str, optional): x轴标签
        ylabel (str, optional): y轴标签
        clabel (str, optional): 颜色条标签
        fig (plt.Figure, optional): 图形对象
        fig_kws (dict, optional): 图形参数
        heatmap_kws (dict, optional): 热图参数
        margin_kws (dict, optional): 边缘图参数
        cbar_kws (dict, optional): 颜色条参数
        subplot_spec (GridSpec, optional): 子图规格
        output_file (str, optional): 输出文件路径
        dpi (int, optional): 图像DPI

    返回:
        dict: 包含所有子图对象的字典
    """
    # 设置绘图风格
    setup_plot_style()
    
    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    #将EV1，舍弃极大极小的2.5%数据之后，排序分组，将连续的数据转化成离散的，其实就是给全基因组排了个序。
    digitized_track, binedges = cooltools.digitize(
        track, n_bins, qrange=(0.075, 0.975)
    )
    print(digitized_track)
    print(binedges)
    # ###第一行确保我们只使用有效的分组数据;第二行计算每个分组的中心值，这些值后续会用于：# 绘制鞍形图的坐标轴 # 计算不同分组之间的相互作用频率 # 帮助解释鞍形图中的模式
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    lo, hi = 0.075, 0.975
    binedges = np.linspace(lo, hi, n_bins + 1)

    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata

    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    if subplot_spec is not None:
        gs = GridSpecFromSubplotSpec(
            nrows=2,
            ncols=3,
            width_ratios=[0.2, 1, 0.1],
            height_ratios=[0.2, 1],
            wspace=0.05,
            hspace=0.05,
            subplot_spec=subplot_spec
        )
    else:
        gs = GridSpec(
            nrows=2,
            ncols=3,
            width_ratios=[0.2, 1, 0.1],
            height_ratios=[0.2, 1],
            wspace=0.05,
            hspace=0.05,
        )

    grid = {}
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5

    if fig is None:
        fig_kws_default = dict(figsize=(5 * 1 / 2.54, 5 * 1 / 2.54))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    # 设置热图
    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

    # 设置边缘图
    margin_kws_default = dict(edgecolor=None, facecolor='gray', linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})

    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])
    plt.barh(
        binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws
    )
    plt.xlim(plt.xlim()[1], plt.xlim()[0])
    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)

    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])
    plt.bar(
        binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws
    )
    plt.xlim(lo, hi)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        cb = plt.colorbar(img, format=MinOneMaxFormatter(vmin, vmax), cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter(vmin, vmax))

    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid['ax_heatmap'].grid(False)
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    
    return grid 


def calculate_saddle_matrix(
    clr: cooler.Cooler,
    track: pd.DataFrame,
    n_bins: int = 98,
    contact_type: str = 'cis',
    view_df: pd.DataFrame = None,
    nproc: int = 8
) -> tuple:
    """
    计算鞍形矩阵

    参数:
        clr (cooler.Cooler): Hi-C数据对象
        track (pd.DataFrame): 特征向量数据框
        n_bins (int): 分箱数量，默认98
        contact_type (str): 相互作用类型，'cis'或'trans'
        view_df (pd.DataFrame): 视图数据框
        nproc (int): 并行处理的核心数，默认8

    返回:
        tuple: (interaction_sum, interaction_count) 相互作用总和矩阵和计数矩阵
    """
    # 计算期望值
    cvd = cooltools.expected_cis(
        clr=clr,
        view_df=view_df,
        nproc=nproc
    )
    
    # 计算鞍形矩阵
    interaction_sum, interaction_count = cooltools.saddle(
        clr,
        cvd,
        track,
        contact_type,
        n_bins=n_bins,
        qrange=(0.075, 0.975),  # 默认参数
        view_df=view_df
    )
    
    return interaction_sum, interaction_count


def calculate_saddle_layout(
    n_plots: int = 1,
    heatmap_size: float = 4.0,  # 热图尺寸（厘米）
    n_cols: int = None,  # 列数，如果为None则自动计算
    n_rows: int = None,  # 行数，如果为None则自动计算
) -> dict:
    """
    计算鞍形图的布局参数
    
    参数:
        n_plots: 子图数量
        heatmap_size: 热图尺寸（厘米）
        n_cols: 列数，如果为None则自动计算
        n_rows: 行数，如果为None则自动计算
    
    返回:
        dict: 包含布局参数的字典
    """
    # 1. 计算所有绝对尺寸（单位：厘米）
    # 热图尺寸
    heatmap_height = heatmap_size
    heatmap_width = heatmap_size
    
    # 柱状图尺寸（热图尺寸的0.15）
    bar_size = heatmap_size * 0.15
    
    # 子图间隔（热图尺寸的0.05）
    gap = heatmap_size * 0.05
    
    # 颜色条尺寸
    colorbar_width = heatmap_size * 0.05  # 颜色条宽度为热图边长的5%
    
    # 边距计算
    margin_left = heatmap_size * 0.2  # 左侧边距
    margin_right = heatmap_size * 0.4  # 右侧边距
    margin_top = heatmap_size * 0.2  # 顶部边距（包含柱状图高度和间隔）
    
    # 2. 计算总尺寸
    # 计算行列数
    if n_cols is None:
        n_cols = min(8, n_plots)  # 默认最多8列
    if n_rows is None:
        n_rows = (n_plots + n_cols - 1) // n_cols  # 向上取整
    
    # 计算总宽度
    # 总宽度 = 左侧边距 + (n_cols-1)个间隔 + n_cols个热图宽度 + 右侧边距
    total_width = margin_left + (n_cols - 1) * gap + n_cols * heatmap_width + margin_right
    
    # 计算总高度
    # 总高度 = 顶部边距（0.2） + 热图高度
    total_height = margin_top + heatmap_height
    
    # 3. 计算相对位置（比例）
    # 热图区域
    heatmap_width_relative = heatmap_width / total_width
    heatmap_height_relative = heatmap_height / total_height
    
    # 柱状图区域（分别计算水平和垂直方向的相对尺寸）
    bar_size_horizontal_relative = bar_size / total_width  # 左侧柱状图的相对宽度
    bar_size_vertical_relative = bar_size / total_height   # 顶部柱状图的相对高度
    
    # 边距和间隔
    margin_left_relative = margin_left / total_width
    margin_right_relative = margin_right / total_width
    margin_top_relative = margin_top / total_height
    gap_relative = gap / total_width
    
    # 颜色条位置（放在最后一个子图的右侧）
    colorbar_left = (margin_left + n_cols * heatmap_width + n_cols * gap) / total_width
    colorbar_width_relative = colorbar_width / total_width
    
    # 4. 转换为英寸（用于matplotlib）
    cm_to_inch = 1/2.54
    fig_width = total_width * cm_to_inch
    fig_height = total_height * cm_to_inch
    
    return {
        # 实际尺寸（厘米）
        'heatmap_size': heatmap_size,
        'heatmap_height': heatmap_height,
        'heatmap_width': heatmap_width,
        'bar_size': bar_size,
        'gap': gap,
        'colorbar_width': colorbar_width,
        'total_width': total_width,
        'total_height': total_height,
        'margin_left': margin_left,
        'margin_right': margin_right,
        'margin_top': margin_top,
        
        # 相对位置（比例）
        'heatmap_width_relative': heatmap_width_relative,
        'heatmap_height_relative': heatmap_height_relative,
        'bar_size_horizontal_relative': bar_size_horizontal_relative,  # 左侧柱状图的相对宽度
        'bar_size_vertical_relative': bar_size_vertical_relative,      # 顶部柱状图的相对高度
        'margin_left_relative': margin_left_relative,
        'margin_right_relative': margin_right_relative,
        'margin_top_relative': margin_top_relative,
        'gap_relative': gap_relative,
        'colorbar_left': colorbar_left,
        'colorbar_width_relative': colorbar_width_relative,
        
        # 图形尺寸（英寸）
        'fig_width': fig_width,
        'fig_height': fig_height,
        
        # 布局信息
        'n_cols': n_cols,
        'n_rows': n_rows,
        
        # 子图间距
        'hspace': gap_relative,  # 水平间距
        'vspace': gap_relative   # 垂直间距
    }


def plot_easy_saddle(
    saddledata: np.ndarray,
    title: str = None,
    xlabel: str = None,
    ylabel: str = None,
    output_file: str = None,
    heatmap_size: float = 4.0,  # 热图尺寸（厘米）
    vmin: float = -2,
    vmax: float = 2,
    cmap: str = None,
):
    """
    简化版鞍形图绘制函数
    
    参数:
        saddledata: 鞍形矩阵数据（interaction_sum/interaction_count）
        title: 图形标题
        xlabel: x轴标签
        ylabel: y轴标签
        output_file: 输出文件路径
        heatmap_size: 热图尺寸（厘米）
        vmin: 颜色映射最小值
        vmax: 颜色映射最大值
        cmap: 颜色映射
    """
    # 设置统一的绘图样式
    setup_plot_style()
    
    # 计算布局参数
    layout = calculate_saddle_layout(heatmap_size=heatmap_size)
    
    # 创建图形
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # 对saddledata做log2_and_mask处理
    masked_matrix = log2_and_mask(saddledata)
    
    # 创建主热图区域
    ax_heatmap = fig.add_axes([
        layout['margin_left_relative'],
        0,
        layout['heatmap_width_relative'],
        layout['heatmap_height_relative']
    ])
    
    cmap = setup_ratio_colormap()
    # 设置颜色映射
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    
    # 绘制热图
    im = ax_heatmap.imshow(
        masked_matrix,
        cmap=cmap,
        norm=norm,
        aspect='auto',
        interpolation='none',
    )

    _, _, _, _ = setup_axes(ax_heatmap)
    # 创建模拟的边缘图数据
    n = len(masked_matrix)
    x = np.linspace(0, 1, n)
    
    # Y轴边缘图（左侧）
    ax_margin_y = fig.add_axes([
        0,
        0,
        layout['bar_size_horizontal_relative'],
        layout['heatmap_height_relative']
    ])
    ax_margin_y.barh(x, width=x, height=1/n, align='edge', color='gray', alpha=0.6)
    ax_margin_y.set_xlim(1, 0)
    ax_margin_y.set_ylim(1, 0)
    ax_margin_y.text(0.5, 0.8, 'EV1 rank', ha='center', va='center', fontsize=5, fontname='Arial', color='black',rotation=90)
    _, _, _, _ = setup_axes(ax_margin_y)
    
    # X轴边缘图（顶部）
    ax_margin_x = fig.add_axes([
        layout['margin_left_relative'],
        layout['heatmap_height_relative'] + layout['gap_relative'],
        layout['heatmap_width_relative'],
        layout['bar_size_vertical_relative']
    ])
    ax_margin_x.bar(x, height=x, width=1/len(x), align='edge', color='gray', alpha=0.6)
    ax_margin_x.set_xlim(0, 1)
    ax_margin_x.set_ylim(0, 1)
    ax_margin_x.text(0.8, 0.5, 'EV1 rank', ha='center', va='center', fontsize=5, fontname='Arial', color='black')
    _, _, _, _ = setup_axes(ax_margin_x)
    
    # 设置颜色条
    setup_diff_colorbar(
        fig=fig,
        sc=im,
        vmin=vmin,
        vmax=vmax,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.6,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=0.15,
        max_label=str(vmax),
        min_label=str(vmin),
        label_text='log2(obs/exp)',
        tick_fontsize=5,
        label_fontsize=5
    )
    
    # 保存图形
    if output_file:
        save_figure_multi_format(fig, output_file, dpi=300, formats=['png', 'svg'])


def clear_saddle_cache(cache_dir: str = "cache"):
    """
    清除鞍形图缓存数据
    
    参数:
        cache_dir (str): 缓存目录路径
    """
    cache_dir = Path(cache_dir)
    if cache_dir.exists():
        for cache_file in cache_dir.glob("*.pkl"):
            cache_file.unlink()
        print(f"已清除缓存目录: {cache_dir}")
    else:
        print(f"缓存目录不存在: {cache_dir}")


def process_saddle_sample(args):
    """处理单个样品的saddle图计算
    
    Args:
        args: 包含(cool_file, eigenvector_file, sample_name, idx, cache_dir, n_bins, contact_type, nproc)的元组
    
    Returns:
        Dict: 处理结果字典
    """
    cool_file, eigenvector_file, sample_name, idx, cache_dir, n_bins, contact_type, nproc = args
    try:
        # 创建缓存目录
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成缓存文件名
        cache_file = cache_dir / f"{sample_name}_saddle_data.pkl"
        
        # 尝试从缓存读取数据
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    saddledata = pickle.load(f)
                print(f"从缓存加载数据成功: {cache_file}")
            except Exception as e:
                print(f"读取缓存文件失败: {e}")
                saddledata = None
        else:
            saddledata = None
        
        # 如果缓存不存在或读取失败，重新计算
        if saddledata is None:
            # 加载数据
            clr = cooler.Cooler(cool_file)
            eigenvector_track = pd.read_csv(eigenvector_file, sep='\t')[['chrom', 'start', 'end', 'E1']]
            
            # 计算鞍形矩阵
            interaction_sum, interaction_count = calculate_saddle_matrix(
                clr=clr,
                track=eigenvector_track,
                n_bins=n_bins,
                contact_type=contact_type,
                nproc=nproc
            )
            saddledata = interaction_sum / interaction_count
            
            # 保存到缓存
            try:
                with open(cache_file, 'wb') as f:
                    pickle.dump(saddledata, f)
                print(f"数据已保存到缓存: {cache_file}")
            except Exception as e:
                print(f"保存缓存文件失败: {e}")
        
        return {
            'status': 'success',
            'saddledata': saddledata,
            'sample_name': sample_name,
            'idx': idx
        }
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e),
            'sample_name': sample_name,
            'idx': idx
        }


def plot_multi_saddle(
    results: List[Dict[str, Any]],
    output_prefix: str,
    vmin: float = -2,
    vmax: float = 2,
    heatmap_size: float = 1.8,
    n_cols: Optional[int] = None,
    n_rows: Optional[int] = None
) -> None:
    """绘制多个样品的saddle图
    
    Args:
        results: 并行计算的结果列表
        output_prefix: 输出文件前缀
        vmin: 热图颜色条最小值，默认-2
        vmax: 热图颜色条最大值，默认2
        heatmap_size: 热图尺寸（厘米）
        n_cols: 列数，如果为None则自动计算
        n_rows: 行数，如果为None则自动计算
    """
    # 设置布局
    n_samples = len(results)
    layout = calculate_saddle_layout(
        n_plots=n_samples,
        heatmap_size=heatmap_size,
        n_cols=n_cols,
        n_rows=n_rows
    )
    
    # 创建图形
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # 设置统一的绘图样式
    setup_plot_style()
    
    # 创建颜色映射
    cmap = setup_ratio_colormap()
    
    # 绘制所有样品的saddle图
    last_sc = None  # 用于存储最后一个子图的sc对象
    for result in results:
        if result['status'] == 'success':
            idx = result['idx']
            sample_name = result['sample_name']
            
            # 计算子图位置
            # 水平位置 = 左侧边距 + idx * (热图宽度 + 间隔)
            left = layout['margin_left_relative'] + idx * (layout['heatmap_width_relative'] + layout['gap_relative'])
            bottom = 0
            
            # 热图区域
            ax_heatmap = fig.add_axes([
                left,
                bottom,
                layout['heatmap_width_relative'],
                layout['heatmap_height_relative']
            ])
            
            # 绘制热图
            masked_matrix = log2_and_mask(result['saddledata'])
            norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
            sc = ax_heatmap.imshow(masked_matrix, cmap=cmap, aspect='auto', interpolation='none', norm=norm)
            last_sc = sc  # 更新最后一个sc对象
            
            # 设置坐标轴
            setup_axes(ax_heatmap)
            
            # 添加样本名称
            ax_heatmap.text(0.5, -0.02, sample_name, 
                          transform=ax_heatmap.transAxes, ha='center', va='top', fontsize=5)
            
            # 只在第一个子图添加柱状图
            if idx == 0:
                # Y轴边缘图（左侧）
                ax_margin_y = fig.add_axes([
                    0,
                    0,
                    layout['bar_size_horizontal_relative'],
                    layout['heatmap_height_relative']
                ])
                x = np.linspace(0, 1, len(masked_matrix))
                ax_margin_y.barh(x, width=x, height=1/len(x), align='edge', color='gray', alpha=0.6)
                ax_margin_y.set_xlim(1, 0)
                ax_margin_y.set_ylim(1, 0)
                ax_margin_y.text(0.2, 0.8, 'EV1 rank', ha='center', va='center', fontsize=5, fontname='Arial', color='black', rotation=90)
                setup_axes(ax_margin_y)
                
                # X轴边缘图（顶部）
                ax_margin_x = fig.add_axes([
                    left,
                    layout['heatmap_height_relative'] + layout['gap_relative'],  # 添加间隔
                    layout['heatmap_width_relative'],
                    layout['bar_size_vertical_relative']
                ])
                ax_margin_x.bar(x, height=x, width=1/len(x), align='edge', color='gray', alpha=0.6)
                ax_margin_x.set_xlim(0, 1)
                ax_margin_x.set_ylim(0, 1)
                ax_margin_x.text(0.8, 0.2, 'EV1 rank', ha='center', va='center', fontsize=5, fontname='Arial', color='black')
                setup_axes(ax_margin_x)
    
    # 添加颜色条（只在最后一个子图添加）
    if last_sc is not None:
        setup_diff_colorbar(
            fig=fig,
            sc=last_sc,
            vmin=vmin,
            vmax=vmax,
            colorbar_left=layout['colorbar_left'],
            colorbar_bottom=0.6,
            colorbar_width=layout['colorbar_width_relative'],
            colorbar_height=0.15,
            max_label=str(vmax),
            min_label=str(vmin),
            label_text='log2(obs/exp)',
            tick_fontsize=5,
            label_fontsize=5
        )
    
    # 保存图片（多格式）
    save_figure_multi_format(fig, output_prefix, dpi=300, formats=['png', 'svg'])


def generate_multi_saddle(
    cool_files: List[str],
    eigenvector_files: List[str],
    output_dir: str,
    sample_names: List[str],
    cache_dir: str = "cache",
    n_bins: int = 98,
    contact_type: str = 'cis',
    heatmap_size: float = 1.8,
    vmin: float = -2,
    vmax: float = 2,
    n_cols: Optional[int] = None,
    n_rows: Optional[int] = None,
    max_workers: int = 28,
    nproc: int = 8
) -> Dict[str, Any]:
    """生成多个样品的saddle图
    
    Args:
        cool_files: Hi-C数据文件路径列表
        eigenvector_files: 特征向量文件路径列表
        output_dir: 输出目录
        sample_names: 样本名称列表
        cache_dir: 缓存目录
        n_bins: 分箱数量
        contact_type: 相互作用类型
        heatmap_size: 热图尺寸（厘米）
        vmin: 颜色映射最小值
        vmax: 颜色映射最大值
        n_cols: 列数，如果为None则自动计算
        n_rows: 行数，如果为None则自动计算
        max_workers: 最大并行进程数，默认28
        nproc: 每个进程使用的CPU核心数，默认8
    
    Returns:
        Dict[str, Any]: 处理结果字典
    """
    from concurrent.futures import ProcessPoolExecutor
    
    # 确保输出目录存在
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 生成输出文件名
    output_prefix = output_dir / "multi_saddle_plot"
    
    # 准备并行处理参数
    args = [(cool_file, eigenvector_file, sample_name, idx, cache_dir, n_bins, contact_type, nproc) 
            for idx, (cool_file, eigenvector_file, sample_name) in enumerate(zip(cool_files, eigenvector_files, sample_names))]
    
    # 并行处理所有样本
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_saddle_sample, args))
    
    # 按原始顺序排序结果
    results.sort(key=lambda x: x['idx'])
    
    # 使用新的可视化函数
    plot_multi_saddle(
        results=results,
        output_prefix=str(output_prefix),
        vmin=vmin,
        vmax=vmax,
        heatmap_size=heatmap_size,
        n_cols=n_cols,
        n_rows=n_rows
    )
    
    # 返回处理结果
    return {
        'status': 'success',
        'output_prefix': str(output_prefix),
        'results': results
    }


def generate_single_saddle(
    cool_file: str,
    eigenvector_file: str,
    output_dir: str,
    sample_name: str,
    cache_dir: str = "cache",
    n_bins: int = 98,
    contact_type: str = 'cis',
    heatmap_size: float = 4,
    vmin: float = -2,
    vmax: float = 2,
    nproc: int = 8
) -> Dict[str, Any]:
    """生成单个样品的saddle图
    
    Args:
        cool_file: Hi-C数据文件路径
        eigenvector_file: 特征向量文件路径
        output_dir: 输出目录
        sample_name: 样本名称
        cache_dir: 缓存目录
        n_bins: 分箱数量
        contact_type: 相互作用类型
        heatmap_size: 热图尺寸（厘米）
        vmin: 颜色映射最小值
        vmax: 颜色映射最大值
        nproc: 并行处理的核心数，默认8
    
    Returns:
        Dict[str, Any]: 处理结果字典
    """
    # 创建输出目录
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 生成输出文件名
    output_file = output_dir / f"{sample_name}_saddle_plot"
    
    # 处理单个样本
    result = process_saddle_sample((
        cool_file,
        eigenvector_file,
        sample_name,
        0,  # idx
        cache_dir,
        n_bins,
        contact_type,
        nproc
    ))
    
    if result['status'] == 'success':
        # 绘制鞍形图
        plot_easy_saddle(
            saddledata=result['saddledata'],
            title=sample_name,
            output_file=str(output_file),
            heatmap_size=heatmap_size,
            vmin=vmin,
            vmax=vmax
        )
        return {
            'status': 'success',
            'output_file': str(output_file),
            'sample_name': sample_name
        }
    else:
        return {
            'status': 'error',
            'error': result['error'],
            'sample_name': sample_name,
            'output_file': str(output_file)
        }
