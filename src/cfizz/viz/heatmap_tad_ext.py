"""
TAD-related heatmap functions (45-degree rotated).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.collections import QuadMesh
import itertools
from typing import Optional, List


# =============================================================================
# 45-degree rotated heatmap functions
# =============================================================================

def pcolormesh_45deg(
    ax,
    matrix_c: np.ndarray,
    start: int = 0,
    resolution: int = 1,
    norm: Optional[Normalize] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: str = 'Blues',
    alpha: float = 1.0
) -> QuadMesh:
    """绘制45度旋转的热图
    
    参数:
        ax: matplotlib轴对象
        matrix_c: 要绘制的矩阵数据
        start: 起始位置，默认0
        resolution: 分辨率，默认1
        norm: 颜色归一化对象，可选
        vmin: 颜色范围最小值
        vmax: 颜色范围最大值
        cmap: 颜色映射，默认'Reds'
        alpha: 透明度，默认1.0
    
    返回:
        QuadMesh: pcolormesh对象
    """
    start_pos_vector = [start + resolution * i for i in range(len(matrix_c) + 1)]
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(
        np.array([(i[1], i[0]) for i in itertools.product(start_pos_vector[::-1], start_pos_vector)]), 
        t
    )
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    
    # 确保vmin/vmax被提供
    if vmin is None or vmax is None:
        raise ValueError("vmin 和 vmax 必须显式传入！")

    # 只对非nan值进行clip
    processed_matrix = matrix_c.copy()
    mask = ~np.isnan(processed_matrix)
    processed_matrix[mask] = np.clip(processed_matrix[mask], vmin, vmax)
    
    im = ax.pcolormesh(x, y, np.flipud(processed_matrix), norm=norm, cmap=cmap, alpha=alpha)
    im.set_rasterized(True)
    return im


def add_tad_boundaries_to_heatmap(
    ax,
    insulation_data,
    matrix_data: np.ndarray,
    matrix_start: int,
    resolution: int,
    boundary_cmap: str = 'Greys',
    boundary_alpha: float = 0.6
) -> None:
    """在热图上添加TAD边界标注
    
    Args:
        ax: matplotlib轴对象
        insulation_data: 包含Insulation Score数据的DataFrame
        matrix_data: 接触矩阵数据
        matrix_start: 矩阵起始位置
        resolution: 分辨率
        boundary_cmap: 边界颜色映射
        boundary_alpha: 边界透明度
    """
    from cfizz.viz.heatmap import mark_boundaries_from_insulation
    
    # 生成边界标记矩阵
    boundary_matrix = mark_boundaries_from_insulation(
        insulation_data=insulation_data,
        matrix_data=matrix_data,
        matrix_start=matrix_start,
        resolution=resolution
    )
    
    # 叠加边界标记
    pcolormesh_45deg(
        ax, 
        boundary_matrix, 
        start=matrix_start, 
        resolution=resolution,
        cmap=boundary_cmap,
        alpha=boundary_alpha, 
        vmin=0, 
        vmax=1
    )


def plot_heatmap_with_tad_boundaries(
    mcool_paths: List[str],
    insulation_paths: List[str],
    chrom: str,
    start: int,
    end: int,
    resolution: int,
    output_path: str,
    window_size: int = 100000,
    cmap: str = 'Reds',
    vmin: float = None,
    vmax: float = None,
    color_scale: str = 'linear',
    plot_size: float = 4,
    triangle_ratio: float = 1,
    boundary_cmap: str = 'Greys',
    boundary_alpha: float = 0.6,
    balance: bool = False,
    dpi: int = 1000,
    sample_names: List[str] = None
) -> dict:
    """绘制带有TAD边界标注的45度旋转热图，支持多个样品
    
    Args:
        mcool_paths: mcool文件路径列表
        insulation_paths: insulation score文件路径列表
        chrom: 染色体
        start: 起始位置
        end: 结束位置
        resolution: 分辨率
        output_path: 输出图片路径
        cmap: 颜色映射，默认为'Reds'
        vmin: 最小值，默认为None（自动计算）
        vmax: 最大值，默认为None（自动计算）
        color_scale: 颜色条缩放方式，'linear'或'log'
        plot_size: 基础尺寸（厘米）
        triangle_ratio: 上三角显示比例
        window_size: 窗口大小（bp），默认100kb
        boundary_cmap: 边界颜色映射
        boundary_alpha: 边界透明度
        balance: 是否使用平衡矩阵
        dpi: 分辨率
        sample_names: 样本名称列表
    """
    from cfizz.viz.heatmap import (
        read_matrix_from_cooler,
        get_bin_index,
        mark_boundaries_from_insulation
    )
    from cfizz.utils.coordinates import get_matrix_range
    from cfizz.viz.layout import (
        calculate_rotated_heatmap_layout,
        setup_axes,
        add_rotated_coordinate_labels,
        setup_horizontal_colorbar,
        save_figure_multi_format
    )
    import os
    
    # 检查输入列表长度是否一致
    if len(mcool_paths) != len(insulation_paths):
        raise ValueError("mcool_paths和insulation_paths的长度必须相同")
    
    n_plots = len(mcool_paths)
    
    # 如果没有提供样本名称，使用默认名称
    if sample_names is None:
        sample_names = [f"Sample {i+1}" for i in range(n_plots)]
    
    # 1. 加载所有接触矩阵数据
    data_list = []
    visible_matrices = []
    for mcool_path in mcool_paths:
        data = read_matrix_from_cooler(
            file_path=mcool_path,
            resolution=resolution,
            chrom=chrom,
            start_pos=start,
            end_pos=end,
            balance=balance
        )
        if data is None:
            raise ValueError(f"无法从文件 {mcool_path} 读取矩阵")
        data_list.append(data)
        
        # 计算实际可视化的矩阵范围
        matrix_size = data.shape[0]
        diagonal_extension = int(matrix_size * triangle_ratio)
        visible_matrix = np.zeros_like(data)
        for j in range(matrix_size):
            for k in range(matrix_size):
                if abs(j - k) <= diagonal_extension:
                    visible_matrix[j, k] = data[j, k]
        visible_matrices.append(visible_matrix)
    
    # 2. 加载所有Insulation Score数据
    from cfizz.io.insulation import read_insulation_scores
    insulation_data_list = []
    for insulation_path in insulation_paths:
        insulation_data = read_insulation_scores(
            file_path=insulation_path,
            windows=window_size,
            chrom=chrom,
            start=start,
            end=end
        )
        insulation_data_list.append(insulation_data)
    
    # 3. 计算布局
    layout = calculate_rotated_heatmap_layout(
        n_plots=n_plots,
        plot_size=plot_size,
        triangle_ratio=triangle_ratio
    )
    
    # 4. 创建图形
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # 5. 计算所有矩阵的最大最小值
    all_values = np.concatenate([m[~np.isnan(m)] for m in visible_matrices])
    vmin, vmax = get_matrix_range(all_values, vmin, vmax)
    if color_scale == 'log':
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    # 6. 遍历每个样本，创建子图
    for i in range(n_plots):
        pos = layout['subplot_positions'][i]
        
        # 创建热图轴
        ax = fig.add_axes([
            pos['heatmap_left'],
            pos['heatmap_bottom'],
            pos['heatmap_width'],
            pos['heatmap_height']
        ])
        
        # 绘制热图
        im = pcolormesh_45deg(
            ax, data_list[i], start=start, resolution=resolution, 
            norm=norm, cmap=cmap, vmin=vmin, vmax=vmax
        )
        
        # 添加TAD边界标注（read_insulation_scores 返回 DataFrame）
        add_tad_boundaries_to_heatmap(
            ax=ax,
            insulation_data=insulation_data_list[i],
            matrix_data=data_list[i],
            matrix_start=start,
            resolution=resolution,
            boundary_cmap=boundary_cmap,
            boundary_alpha=boundary_alpha
        )
        
        # 设置坐标轴和标签
        ax.set_aspect(0.5)
        ax.set_ylim(0, (end - start) * triangle_ratio)
        xmin, xmax, ymin, ymax = setup_axes(ax)
        
        # 只在最后一个子图添加坐标标签
        if i == n_plots - 1:
            add_rotated_coordinate_labels(ax, xmin, xmax, ymin, start, end, chrom)
        
        # 添加样本名称到左侧
        ax.text(0, 0.5, sample_names[i],
                ha='right', va='center',
                fontsize=5, fontname='Arial',
                rotation=90,
                transform=ax.get_yaxis_transform())
        
        # 设置颜色条（只在最后一个子图添加）
        if i == n_plots - 1:
            colorbar_bottom = layout['last_heatmap_bottom'] - layout['colorbar_height'] - 0.01
            setup_horizontal_colorbar(
                fig=fig,
                sc=im,
                vmin=vmin,
                vmax=vmax,
                fontsize=5,
                label="contacts",
                colorbar_left=pos['colorbar_left'],
                colorbar_bottom=colorbar_bottom,
                colorbar_width=pos['colorbar_width'],
                colorbar_height=pos['colorbar_height'],
                balance=balance
            )
    
    # 7. 保存图片
    save_figure_multi_format(fig, output_path, dpi=dpi)
    
    return {
        'status': 'success',
        'output_file': output_path,
        'n_samples': n_plots
    }
