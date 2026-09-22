"""
APA (Aggregate Peak Analysis) end-to-end visualization module.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, Dict, Optional
from scipy.special import ndtr

from cfizz.io.loops import read_loops
from cfizz.viz.heatmap import (
    plot_single_heatmap, 
    read_matrix_from_cooler, 
    plot_multi_heatmap,
    calculate_heatmap_layout,
    get_matrix_range,
    setup_axes,
    setup_colorbar
)
from cfizz.viz.layout import (
    setup_plot_style,
    save_figure_multi_format
)
# 注意: setup_axes 在 heatmap 和 layout 都有, 用 heatmap 的(跟 L-2 一致)


def extract_apa_submatrix(matrix: np.ndarray, positions: List[Tuple[int, int]], window: int = 5) -> List[np.ndarray]:
    """提取APA子矩阵

    Args:
        matrix: 接触矩阵
        positions: 位置列表，每个元素为(start, end)坐标对
        window: 窗口大小

    Returns:
        子矩阵列表
    """
    submatrices = []
    for si, ei in positions:
        if si > ei:
            si, ei = ei, si
        submatrix = matrix[si-window:si+window+1, ei-window:ei+window+1]
        if submatrix.shape == (2*window+1, 2*window+1):
            submatrices.append(submatrix)
    return submatrices



def _default_apa_label(balance: bool) -> str:
    """
    根据 balance 状态返回 APA colorbar 的默认 label.

    与 5_3 TAD pileup 用的 cbar_label 文案保持完全一致, 让两个累计图视觉对齐.

    Parameters
    ----------
    balance : bool
        True → 使用平衡矩阵 (normalized contacts), 文案 'mean normalized contacts'.
        False → 使用 raw counts, 文案 'mean contacts'.

    Returns
    -------
    label : str
    """
    return 'mean normalized contacts' if balance else 'mean contacts'


def _apply_apa_color_scale(matrix, color_scale):
    """
    对 APA 累计矩阵应用 color_scale 变换.

    与 compute_apa_oe_pileup_multi 内部的 color_scale 处理语义一致 —
    linear 直接显示, log10/log2 对数压缩.

    Notes
    -----
    与 viz/pileup.py 同样的 filled(1e-10) 兜底策略 — 0 值在 log 变换后变 -10.
    raw + sum + log 路径 caller 需自行传 vmin=0 规避.
    """
    if color_scale == 'linear':
        return matrix
    if color_scale == 'log2':
        return np.log2(np.ma.masked_less_equal(matrix, 0).filled(1e-10))
    if color_scale in ('log10', 'log10_linear'):
        return np.log10(np.ma.masked_less_equal(matrix, 0).filled(1e-10))
    raise ValueError(
        f"_apply_apa_color_scale: color_scale 必须是 'linear' / 'log2' / 'log10' / 'log10_linear', "
        f"实得 {color_scale!r}"
    )


def analyze_apa(
    submatrices: List[np.ndarray],
    window: int = 5,
    corner_size: int = 3,
    method: str = 'mean'
) -> Tuple[np.ndarray, float, float, float, float]:
    """
    分析APA结果（主流定义，兼容旧代码）

    Args:
        submatrices: 子矩阵列表
        window: 中心点坐标（如5）
        corner_size: 角落区域大小（如3）
        method: 子矩阵累积方式, 'mean' / 'median' / 'sum', 默认 'mean'.

    Returns:
        avg: 累计信号矩阵（按 method 聚合）
        score: APA分数（中心点/左下角均值）
        z: z分数
        p: p值
        maxi: 最大值（右上角均值*5）

    Notes
    -----
    1/99% 异常值截断逻辑保留: per-loop 子矩阵均值先做 p1/p99 截断再聚合。
    score/z/p/maxi 仍然基于聚合后的 avg 矩阵计算（语义不变）。
    """
    if not submatrices:
        return None, 0, 0, 1, 0

    method = method.lower()
    if method not in ('mean', 'median', 'sum'):
        raise ValueError(
            f"analyze_apa: method 必须是 'mean' / 'median' / 'sum', 实得 {method!r}"
        )

    # 去除异常值
    mean_arr = np.array([np.mean(arr) for arr in submatrices])
    p99 = np.percentile(mean_arr, 99)
    p1 = np.percentile(mean_arr, 1)
    mask = (mean_arr < p99) & (mean_arr > p1)
    stacked = np.array(submatrices)[mask]

    if method == 'mean':
        avg = np.mean(stacked, axis=0)
    elif method == 'median':
        avg = np.median(stacked, axis=0)
    else:  # method == 'sum'
        avg = np.sum(stacked, axis=0)

    # 背景（左下角）
    lowerpart = avg[-corner_size:, :corner_size]
    # 右上角
    upperpart = avg[:corner_size, -corner_size:]
    # 中心点
    center_value = avg[window, window]

    # APA score
    score = center_value / lowerpart.mean() if lowerpart.mean() != 0 else 0
    # z-score
    z = (center_value - lowerpart.mean()) / lowerpart.std() if lowerpart.std() != 0 else 0
    # p-value
    p = 1 - ndtr(z)
    # maxi
    maxi = upperpart.mean() * 5

    return avg, score, z, p, maxi

def extract_submatrix_for_loop(args):
    """为单个loop提取子矩阵的辅助函数
    
    Args:
        args: 包含所有必要参数的元组
            (mcool_path, loop, resolution, window, min_distance, balance)
    
    Returns:
        子矩阵或None（如果提取失败）
    """
    mcool_path, loop, resolution, window, min_distance, balance = args
    
    # 计算bin坐标
    x, y = loop['start1']//resolution, loop['start2']//resolution
    if abs(y-x) < min_distance:
        return None
        
    # 读取该loop周围的矩阵区域
    matrix = read_matrix_from_cooler(
        file_path=mcool_path,
        resolution=resolution,
        chrom=loop['chrom1'],
        start_pos=loop['start1'] - window*resolution,
        end_pos=loop['end2'] + window*resolution,
        balance=balance
    )
    if matrix is None:
        return None
        
    # 提取子矩阵
    si = x - (loop['start1'] - window*resolution)//resolution
    ei = y - (loop['start1'] - window*resolution)//resolution
    if si > ei:
        si, ei = ei, si
        
    try:
        submatrix = matrix[si-window:si+window+1, ei-window:ei+window+1]
        if submatrix.shape == (2*window+1, 2*window+1):
            return submatrix
    except IndexError:
        return None
    
    return None

def plot_apa_heatmap_visualization(
    matrix: np.ndarray,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    dpi: int = 200,
    balance: bool = False,
    cbar_label: Optional[str] = None,
    color_scale: str = 'linear',
) -> plt.Figure:
    """绘制APA热图

    Args:
        matrix: APA分析得到的平均信号矩阵
        vmin: 颜色条最小值
        vmax: 颜色条最大值
        cmap: 颜色映射名称
        dpi: 图片分辨率
        balance: 是否使用平衡矩阵
        cbar_label: 颜色条标签. None → 按 balance 自动选
                    (True → 'mean normalized contacts', False → 'mean contacts')
        color_scale: 颜色变换, 'linear' / 'log2' / 'log10' / 'log10_linear', 默认 'linear'.
                    在 imshow 前对 matrix 应用, 影响显示数值范围.

    Returns:
        matplotlib.figure.Figure: 包含热图的图形对象
    """
    # 处理cmap参数
    if isinstance(cmap, str):
        cmap = plt.colormaps[cmap]

    # 计算布局参数
    layout = calculate_heatmap_layout(n_plots=1, plot_size=4.0)

    # 创建图形
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))

    # 添加热图区域
    ax = fig.add_axes([layout['margin_left']/layout['total_width'],
                      layout['margin_bottom']/layout['total_height'],
                      layout['plot_width'],
                      layout['plot_height']])

    # 应用 color_scale 变换 (linear / log2 / log10) — 先变换再 derive vmin/vmax
    matrix_disp = _apply_apa_color_scale(matrix, color_scale)

    # 获取矩阵范围 — 必须从 matrix_disp 派生, 否则 vmin/vmax 与显示数据尺度不一致
    vmin, vmax = get_matrix_range(matrix_disp, vmin, vmax)

    # 绘制热图
    sc = ax.imshow(matrix_disp, cmap=cmap, aspect='auto', interpolation='none',
                  vmax=vmax, vmin=vmin)

    # 设置坐标轴和边框
    setup_axes(ax)

    # 添加颜色条
    setup_colorbar(
        fig=fig,
        sc=sc,
        vmin=vmin,
        vmax=vmax,
        balance=balance,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.75,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=0.15,
        label_fontsize=6,
        tick_fontsize=5,
        label=cbar_label if cbar_label is not None else _default_apa_label(balance)
    )

    return fig

def plot_apa_heatmap(
    mcool_path: str,
    loops_path: str,
    output_path: str,
    window: int = 5,
    corner_size: int = 3,
    min_distance: int = 10,
    resolution: int = 10000,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    dpi: int = 200,
    balance: bool = True,
    n_processes: int = None,
    method: str = 'mean',
    cbar_label: Optional[str] = None,
    color_scale: str = 'linear',
):
    """绘制APA热图

    Args:
        mcool_path: mcool文件路径
        loops_path: loops文件路径
        output_path: 输出图片路径
        window: 窗口大小
        corner_size: 角落区域大小
        min_distance: 最小距离（bin数）
        resolution: 分辨率
        vmin: 颜色条最小值
        vmax: 颜色条最大值
        cmap: 颜色映射名称
        dpi: 图片分辨率
        balance: 是否使用平衡矩阵
        n_processes: 并行处理的进程数，默认为None（使用CPU核心数）
        method: 子矩阵累积方式 'mean'/'median'/'sum', 默认 'mean'.
        cbar_label: 颜色条标签. None → 按 balance 自动选.
        color_scale: 颜色变换 'linear'/'log2'/'log10'/'log10_linear', 默认 'linear'.
    """
    # 1. 读取loops数据
    loops_data = read_loops(loops_path)

    # 2. 并行提取子矩阵
    if n_processes is None:
        n_processes = cpu_count()

    print(f"\n开始并行提取子矩阵...")
    print(f"使用进程数: {n_processes}")
    print(f"总loops数量: {len(loops_data)}")

    # 准备参数
    args_list = [(mcool_path, loop, resolution, window, min_distance, balance)
                 for _, loop in loops_data.iterrows()]

    # 使用进程池并行处理
    with Pool(n_processes) as pool:
        results = pool.map(extract_submatrix_for_loop, args_list)

    # 过滤掉None结果
    submatrices = [matrix for matrix in results if matrix is not None]

    print(f"成功提取子矩阵数量: {len(submatrices)}")

    # 3. 分析APA
    avg, score, z, p, maxi = analyze_apa(submatrices, window, corner_size, method=method)

    # 4. 绘制热图
    fig = plot_apa_heatmap_visualization(
        matrix=avg,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        dpi=dpi,
        balance=balance,
        cbar_label=cbar_label,
        color_scale=color_scale,
    )

    # 5. 保存图片
    save_figure_multi_format(fig, output_path, dpi=dpi)

def calculate_multi_apa(
    mcool_paths: List[str],
    loops_paths: List[str],
    window: int = 7,
    corner_size: int = 3,
    min_distance: int = 10,
    resolution: int = 10000,
    balance: bool = True,
    n_processes: int = None,
    sample_names: List[str] = None,
    method: str = 'mean',
) -> Tuple[List[np.ndarray], List[float], List[str], List[int]]:
    """计算多个样本的APA分析结果

    Args:
        mcool_paths: mcool文件路径列表
        loops_paths: loops文件路径列表
        window: 窗口大小
        corner_size: 角落区域大小
        min_distance: 最小距离（bin数）
        resolution: 分辨率
        balance: 是否使用平衡矩阵
        n_processes: 并行处理的进程数，默认为None（使用CPU核心数）
        sample_names: 样本名称列表（可选）
        method: 子矩阵累积方式 'mean'/'median'/'sum', 默认 'mean'.

    Returns:
        all_avg_matrices: 所有样本的累计信号矩阵列表
        all_scores: 所有样本的APA分数列表
        sample_names: 样本名称列表
        all_loop_counts: 所有样本成功提取的loops数量列表
    """
    # 生成样本名称 - 优先使用外部传入的
    if sample_names is None:
        sample_names = []
        for path in mcool_paths:
            filename = os.path.basename(path)
            size = filename.split('.')[-2]  # 获取倒数第二个部分
            sample_names.append(size)

    n_samples = len(mcool_paths)

    # 存储所有样本的APA结果
    all_avg_matrices = []
    all_scores = []
    all_loop_counts = []

    # 处理每个样本
    for i, (mcool_path, loops_path, sample_name) in enumerate(zip(mcool_paths, loops_paths, sample_names)):
        print(f"\n处理样本 {i+1}/{n_samples}: {sample_name}")

        # 1. 读取loops数据
        loops_data = read_loops(loops_path)

        # 2. 并行提取子矩阵
        if n_processes is None:
            n_processes = cpu_count()

        print(f"开始并行提取子矩阵...")
        print(f"使用进程数: {n_processes}")
        print(f"总loops数量: {len(loops_data)}")

        # 准备参数
        args_list = [(mcool_path, loop, resolution, window, min_distance, balance)
                     for _, loop in loops_data.iterrows()]

        # 调试输出：显示传递给 extract_submatrix_for_loop 的前3个参数示例
        if i == 0:  # 只为第一个样本打印详细调试信息
            print(f"\n调试信息 (样本 {sample_name}):")
            print(f"  传递给 extract_submatrix_for_loop 的参数结构:")
            print(f"  参数1: mcool_path (字符串)")
            print(f"  参数2: loop (字典，包含chrom1, start1, start2, end2)")
            print(f"  参数3: resolution (整数) = {resolution}")
            print(f"  参数4: window (整数) = {window}")
            print(f"  参数5: min_distance (整数) = {min_distance}")
            print(f"  参数6: balance (布尔值) = {balance}")
            if len(args_list) > 0:
                print(f"\n  第一个loop的参数详情:")
                arg_example = args_list[0]
                print(f"    mcool_path: {arg_example[0]}")
                print(f"    loop: {arg_example[1]}")
                print(f"    resolution: {arg_example[2]}")
                print(f"    window: {arg_example[3]}")
                print(f"    min_distance: {arg_example[4]}")
                print(f"    balance: {arg_example[5]}")

        # 使用进程池并行处理
        with Pool(n_processes) as pool:
            results = pool.map(extract_submatrix_for_loop, args_list)

        # 过滤掉None结果
        submatrices = [matrix for matrix in results if matrix is not None]

        print(f"成功提取子矩阵数量: {len(submatrices)}")

        # 3. 分析APA
        avg, score, z, p, maxi = analyze_apa(submatrices, window, corner_size, method=method)
        all_avg_matrices.append(avg)
        all_scores.append(score)
        all_loop_counts.append(len(submatrices))

    return all_avg_matrices, all_scores, sample_names, all_loop_counts

def visualize_multi_apa(
    all_avg_matrices: List[np.ndarray],
    all_scores: List[float],
    sample_names: List[str],
    output_path: str,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    dpi: int = 200,
    balance: bool = True,
    plot_size: float = 2,
    loop_counts: List[int] = None,
    cbar_label: Optional[str] = None,
    color_scale: str = 'linear',
) -> None:
    """可视化多个样本的APA分析结果

    Args:
        all_avg_matrices: 所有样本的累计信号矩阵列表
        all_scores: 所有样本的APA分数列表
        sample_names: 样本名称列表
        output_path: 输出图片路径
        vmin: 颜色条最小值
        vmax: 颜色条最大值
        cmap: 颜色映射名称
        dpi: 图片分辨率
        balance: 是否使用平衡矩阵
        plot_size: 每个子图的大小
        loop_counts: 每个样本成功提取的loops数量列表（可选）
        cbar_label: 颜色条标签. None → 按 balance 自动选.
        color_scale: 颜色变换 'linear'/'log2'/'log10'/'log10_linear', 默认 'linear'.
    """
    # 处理cmap参数
    if isinstance(cmap, str):
        cmap = plt.colormaps[cmap]

    # 计算布局参数
    setup_plot_style()
    n_samples = len(sample_names)
    n_cols = n_samples
    # n_cols = min(8, n_samples)  # 每行最多8个子图
    layout = calculate_heatmap_layout(n_plots=n_cols, plot_size=plot_size)

    # 创建图形
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))

    # 计算所有矩阵的统一颜色范围 — 先应用 color_scale 变换再 derive, 避免 vmin/vmax 与显示数据尺度不一致
    all_disps = []
    for matrix in all_avg_matrices:
        if matrix is not None and not (isinstance(matrix, np.ndarray) and np.all(np.isnan(matrix))):
            all_disps.append(_apply_apa_color_scale(matrix, color_scale))
    all_values = np.concatenate([m.flatten() for m in all_disps]) if all_disps else np.array([])
    vmin, vmax = get_matrix_range(all_values, vmin, vmax)
    print(f"统一颜色范围: vmin={vmin:.4f}, vmax={vmax:.4f}")

    # 绘制所有子图
    sc = None  # 初始化 sc 以防没有有效矩阵
    for i, (avg, score, sample_name) in enumerate(zip(all_avg_matrices, all_scores, sample_names)):
        # 获取loops数量（如果有的话）
        loop_count = loop_counts[i] if loop_counts and i < len(loop_counts) else None

        # 计算子图位置（相对位置）
        col = i % n_cols
        row = i // n_cols

        # 计算子图位置，使用与plot_multi_heatmap相同的逻辑
        left = layout['margin_left']/layout['total_width'] + col * (layout['plot_width'] + layout['hspace'])
        bottom = layout['margin_bottom']/layout['total_height'] + (layout['n_rows'] - 1 - row) * (layout['plot_height'] + layout['vspace'])

        # 添加子图
        ax = fig.add_axes([left, bottom, layout['plot_width'], layout['plot_height']])

        # 判断是否为无loops样本
        if avg is None or (isinstance(avg, np.ndarray) and np.all(np.isnan(avg))):
            # 空子图，不画热图
            ax.axis('off')
            title_text = f"{sample_name}\nno loops"
            if loop_count is not None:
                title_text += f"\n({loop_count} loops)"
            ax.set_title(title_text, fontsize=5, pad=2)
            if i == 0:
                ax.text(0.0, 1.01, "APA\nscore:", transform=ax.transAxes, fontsize=5, ha='left', va='bottom')
            continue

        # 正常绘制热图
        # 应用 color_scale 变换 (linear / log2 / log10)
        avg_disp = _apply_apa_color_scale(avg, color_scale)
        sc = ax.imshow(avg_disp, cmap=cmap, aspect='auto', interpolation='none',
                      vmax=vmax, vmin=vmin)

        # 设置坐标轴和边框
        setup_axes(ax)

        # 添加样本名称和分数
        title_text = f"{sample_name}\n{score:.4g}"
        if loop_count is not None:
            title_text += f"\n({loop_count} loops)"
        ax.set_title(title_text, fontsize=5, pad=2)  # 减小标题与图片的间距
        if i == 0:  # 只在第一个子图添加APA score标签
            ax.text(0.0, 1.01, "APA\nscore:", transform=ax.transAxes, fontsize=5, ha='left', va='bottom')

    # 添加颜色条
    if sc is not None:
        setup_colorbar(
            fig=fig,
            sc=sc,
            vmin=vmin,
            vmax=vmax,
            balance=balance,
            colorbar_left=layout['colorbar_left'],
            colorbar_bottom=0.75,
            colorbar_width=layout['colorbar_width_relative'],
            colorbar_height=0.15,
            label_fontsize=6,
            tick_fontsize=5,
            label=cbar_label if cbar_label is not None else _default_apa_label(balance)
        )

    # 保存图片
    save_figure_multi_format(fig, output_path, dpi=dpi)

def plot_multi_apa_heatmap(
    mcool_paths: List[str],
    loops_paths: List[str],
    output_path: str,
    sample_names: List[str] = None,
    window: int = 5,
    corner_size: int = 3,
    min_distance: int = 10,
    resolution: int = 10000,
    vmin: float = None,
    vmax: float = None,
    cmap: str = 'Reds',
    dpi: int = 200,
    balance: bool = True,
    n_processes: int = None,
    plot_size: float = 2,
    method: str = 'mean',
    cbar_label: Optional[str] = None,
    color_scale: str = 'linear',
) -> Dict[str, float]:
    """绘制多个样本的APA热图

    Args:
        mcool_paths: mcool文件路径列表
        loops_paths: loops文件路径列表
        output_path: 输出图片路径
        sample_names: 样本名称列表，默认为None（使用mcool文件名）
        window: 窗口大小
        corner_size: 角落区域大小
        min_distance: 最小距离（bin数）
        resolution: 分辨率
        vmin: 颜色条最小值
        vmax: 颜色条最大值
        cmap: 颜色映射名称
        dpi: 图片分辨率
        balance: 是否使用平衡矩阵
        n_processes: 并行处理的进程数，默认为None（使用CPU核心数）
        plot_size: 每个子图的大小
        method: 子矩阵累积方式 'mean'/'median'/'sum', 默认 'mean'.
        cbar_label: 颜色条标签. None → 按 balance 自动选.
        color_scale: 颜色变换 'linear'/'log2'/'log10'/'log10_linear', 默认 'linear'.

    Returns:
        Dict[str, float]: 样本名称到APA分数的映射
    """
    # 1. 计算APA结果
    all_avg_matrices, all_scores, used_sample_names, all_loop_counts = calculate_multi_apa(
        mcool_paths=mcool_paths,
        loops_paths=loops_paths,
        window=window,
        corner_size=corner_size,
        min_distance=min_distance,
        resolution=resolution,
        balance=balance,
        n_processes=n_processes,
        sample_names=sample_names,
        method=method,
    )

    # 2. 可视化结果
    # 调试输出：显示传递给 extract_submatrix_for_loop 的参数示例
    print("\n=== 调试信息：传递给 extract_submatrix_for_loop 的参数示例 ===")
    if len(mcool_paths) > 0:
        sample_idx = 0  # 显示第一个样本的参数作为示例
        loops_data_sample = read_loops(loops_paths[sample_idx]) if loops_paths else pd.DataFrame()
        if len(loops_data_sample) > 0:
            # 显示第一个loop的参数
            first_loop = loops_data_sample.iloc[0]
            debug_args = (
                mcool_paths[sample_idx],
                {
                    'chrom1': first_loop['chrom1'],
                    'start1': first_loop['start1'],
                    'start2': first_loop['start2'],
                    'end2': first_loop['end2']
                },
                resolution,
                window,
                min_distance,
                balance
            )
            print(f"样本 {used_sample_names[sample_idx]} 的第一个loop参数:")
            print(f"  - mcool路径: {debug_args[0]}")
            print(f"  - loop信息: {debug_args[1]}")
            print(f"  - 分辨率: {debug_args[2]}")
            print(f"  - 窗口大小: {debug_args[3]}")
            print(f"  - 最小距离: {debug_args[4]}")
            print(f"  - 是否平衡: {debug_args[5]}")
            print(f"总共将处理 {len(loops_data_sample)} 个loops")
            print("=" * 60)

    visualize_multi_apa(
        all_avg_matrices=all_avg_matrices,
        all_scores=all_scores,
        sample_names=used_sample_names,
        output_path=output_path,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        dpi=dpi,
        balance=balance,
        plot_size=plot_size,
        loop_counts=all_loop_counts,
        cbar_label=cbar_label,
        color_scale=color_scale,
    )

    # 3. 返回APA分数
    return dict(zip(used_sample_names, all_scores))


# =============================================================================
# O/E APA pileup (阶段 2 新增)
# =============================================================================
# 设计思路:
#   1. 子矩阵抽取走 cooltools.pileup (原生支持 BEDPE features_df + expected_df)
#      — 与 viz/pileup.py:cooltools_tad_pileup 同一套基础设施
#   2. 累积方式自写: np.nanmean / np.nanmedian / np.nansum (跟 TAD pileup.method 对齐)
#   3. 颜色变换 log2(obs/exp) (与 TAD O/E pileup 对齐)
#   4. 复用 cfizz.analyze.compartment.get_view_df (跟 6_3 / 5_3 O/E pileup 一致)
#   5. 不复用现有 extract_submatrix_for_loop / analyze_apa (O/E 模式语义不同,
#      1/99% 异常值截断和 APA score 在 O/E 模式下不适用)


def _read_loops_for_cooltools_pileup(loops_path: str) -> pd.DataFrame:
    """
    读 loops 并转成 cooltools.pileup 期望的 BEDPE 格式 features_df.

    cooltools.pileup 接受列:
        - bed:    chrom, start, end
        - bedpe:  chrom1, start1, end1, chrom2, start2, end2
    Loops 天然是 bedpe (两 anchor), 所以走 bedpe 路径.

    Returns
    -------
    features_df : pd.DataFrame
        列: chrom1, start1, end1, chrom2, start2, end2
    """
    loops_data = read_loops(loops_path)
    # read_loops 已经返回 6 列 BEDPE (chrom1/start1/end1/chrom2/start2/end2)
    # 但保险起见, 显式取这几列
    required = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    missing = [c for c in required if c not in loops_data.columns]
    if missing:
        raise ValueError(
            f"loops 文件 {loops_path!r} 缺少 cooltools 期望的 BEDPE 列: {missing}. "
            f"实际列: {list(loops_data.columns)}"
        )
    return loops_data[required].copy()


def _aggregate_pileup_stack(
    stack: np.ndarray,
    method: str = 'mean'
) -> np.ndarray:
    """
    把 cooltools.pileup 返回的 (n, D, D) stack 按 method 聚合为 (D, D) 矩阵.

    与 viz/pileup.py:cooltools_tad_pileup (L579-586) 完全相同的逻辑 —
    TAD pileup 和 APA O/E 共享同一套聚合语义.
    """
    method = method.lower()
    if method == 'mean':
        return np.nanmean(stack, axis=0)
    if method == 'median':
        return np.nanmedian(stack, axis=0)
    if method == 'sum':
        return np.nansum(stack, axis=0)
    raise ValueError(
        f"_aggregate_pileup_stack: method 必须是 'mean' / 'median' / 'sum', 实得 {method!r}"
    )


def compute_apa_oe_pileup_single(
    mcool_path: str,
    loops_path: str,
    expected_df: pd.DataFrame,
    output_path: str,
    sample_name: Optional[str] = None,
    flank: int = 70_000,
    resolution: int = 10_000,
    balance: bool = True,
    method: str = 'mean',
    vmin: float = -1.0,
    vmax: float = 1.0,
    cmap: str = 'coolwarm',
    color_scale: str = 'log2',
    cbar_label: Optional[str] = None,
    dpi: int = 1000,
    plot_size: float = 4.0,
) -> str:
    """
    单 sample O/E APA pileup (走 cooltools.pileup 原生 obs/exp 路径).

    与 compute_tad_oe_pileup (5_3 用) 设计对称 — 不同点:
      - features_df 来自 loops 文件 (BEDPE 格式)
      - flank 是 window 中心到边界的 bp 数, 对应 APA 的中心圈尺寸 (window=flank/resolution)
      - 累积方式 (method) 与 TAD O/E pileup 一致

    Parameters
    ----------
    mcool_path : str
        .mcool 路径 (本函数会拼接 ::resolutions/{resolution}).
    loops_path : str
        loops 文件路径 (BEDPE 6 列格式: chrom1/start1/end1/chrom2/start2/end2).
    expected_df : pd.DataFrame
        cooltools.expected_cis 返回的 P(s) DataFrame, 全染色体 per-sample.
        建议由 cfizz.analyze.oe.compute_expected_cis_dataframe() 生成.
    output_path : str
        输出图片前缀 (无扩展名), 会落 .png + .svg.
    sample_name : str, optional
        标题用的样本名. None → 从 mcool_path 派生.
    flank : int
        cooltools.pileup 的 flank 参数, 默认 70_000 (= 7 bins @ 10kb resolution,
        对应 APA window=7).
    resolution : int
        Hi-C 分辨率 (bp), 默认 10_000.
    balance : bool
        是否使用 balance 权重 (clr_weight_name='weight' if balance else None).
        注意 O/E 模式天然依赖 balance 矩阵, 这里 balance 只是控制 weight 列名.
    method : str
        子矩阵累积方式 'mean' / 'median' / 'sum', 默认 'mean'.
    vmin, vmax : float
        颜色条范围, 默认 ±1 (log2(obs/exp) 合理范围).
    cmap : str
        颜色映射, 默认 'coolwarm' (发散型).
    color_scale : str
        颜色变换: 'log2' / 'log10' / 'linear', 默认 'log2'.
    cbar_label : str, optional
        颜色条标签. None → 按 color_scale 派生 ('log2(obs/exp)' / 'log10(obs/exp)' / 'obs/exp').
    dpi : int
        图片分辨率, 默认 1000.
    plot_size : float
        热图尺寸 (cm), 默认 4.0.

    Returns
    -------
    output_path : str
        输出图片路径 (无扩展名).
    """
    try:
        import cooler as cooler_lib
        import cooltools
    except ImportError:
        raise ImportError(
            "compute_apa_oe_pileup_single 需要 cooler + cooltools"
        )

    from cfizz.analyze.compartment import get_view_df
    from cfizz.viz.layout import save_figure_multi_format, calculate_heatmap_layout

    if sample_name is None:
        sample_name = os.path.basename(mcool_path).split('.')[0]

    # 1) 读 loops → BEDPE features_df
    loops_bedpe = _read_loops_for_cooltools_pileup(loops_path)
    print(f"  → {sample_name}: {len(loops_bedpe)} loops loaded for O/E pileup")

    # 2) cooltools.pileup
    clr = cooler_lib.Cooler(f"{mcool_path}::resolutions/{resolution}")
    view_df = get_view_df(clr)

    stack = cooltools.pileup(
        clr=clr,
        features_df=loops_bedpe,
        view_df=view_df,
        expected_df=expected_df,
        flank=flank,
        min_diag='auto',
        clr_weight_name='weight' if balance else None,
    )
    print(f"  → {sample_name}: stack shape {stack.shape}")

    # 3) 按 method 聚合
    mtx = _aggregate_pileup_stack(stack, method=method)

    # 4) 颜色变换
    if color_scale == 'log2':
        mtx_disp = np.log2(np.ma.masked_less_equal(mtx, 0).filled(1e-10))
        auto_label = 'log2(obs/exp)'
    elif color_scale == 'log10':
        mtx_disp = np.log10(np.ma.masked_less_equal(mtx, 0).filled(1e-10))
        auto_label = 'log10(obs/exp)'
    else:  # 'linear'
        mtx_disp = mtx
        auto_label = 'obs/exp'

    final_cbar_label = cbar_label if cbar_label is not None else auto_label

    # 5) 画图
    layout = calculate_heatmap_layout(n_plots=1, plot_size=plot_size)
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    ax = fig.add_axes([
        layout['margin_left']/layout['total_width'],
        layout['margin_bottom']/layout['total_height'],
        layout['plot_width'],
        layout['plot_height']
    ])
    sc = ax.imshow(mtx_disp, cmap=cmap, aspect='auto', interpolation='none',
                   vmin=vmin, vmax=vmax)

    from cfizz.viz.heatmap import setup_axes
    setup_axes(ax)
    # 计算 APA score (基于未变换的 obs/exp 矩阵 mtx)
    # score = center_value / corner_mean, 反映 loop 中心相对角落的富集倍数
    # corner_size=3 匹配 analyze_apa 默认值
    D = mtx.shape[0]
    window_px = D // 2
    apa_corner_size = 3
    center_val = mtx[window_px, window_px]
    corner_mean = mtx[-apa_corner_size:, :apa_corner_size].mean()
    apa_score = center_val / corner_mean if corner_mean > 0 else 0.0

    # Title 风格 (模仿 visualize_multi_apa L557-563):
    # sample_name\nscore\n(N loops)
    title_text = f"{sample_name}\n{apa_score:.4g}"
    if len(loops_bedpe) > 0:
        title_text += f"\n({len(loops_bedpe)} loops)"
    ax.set_title(title_text, fontsize=5, pad=2)
    # "APA\nscore:" 标签 (single 永远是 i=0, 也加)
    ax.text(0.0, 1.01, "APA\nscore:", transform=ax.transAxes, fontsize=5, ha='left', va='bottom')

    # 6) colorbar
    from cfizz.viz.heatmap import setup_colorbar
    setup_colorbar(
        fig=fig, sc=sc, vmin=vmin, vmax=vmax, balance=balance,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.75,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=0.15,
        label_fontsize=6,
        tick_fontsize=5,
        label=final_cbar_label,
    )

    # 7) 保存
    save_figure_multi_format(fig, output_path, dpi=dpi)
    return output_path


def compute_apa_oe_pileup_multi(
    mcool_paths: List[str],
    loops_paths: List[str],
    expected_dfs: List[pd.DataFrame],
    output_path: str,
    sample_names: Optional[List[str]] = None,
    flank: int = 70_000,
    resolution: int = 10_000,
    balance: bool = True,
    method: str = 'mean',
    vmin: float = -1.0,
    vmax: float = 1.0,
    cmap: str = 'coolwarm',
    color_scale: str = 'log2',
    cbar_label: Optional[str] = None,
    dpi: int = 1000,
    plot_size: float = 4.0,
) -> Dict[str, str]:
    """
    多 sample O/E APA pileup 并排, 走 cooltools.pileup 原生 obs/exp 路径.

    与 compute_tad_oe_pileup (5_3 用) 设计对称 — per-sample 各自跑
    cooltools.pileup, 各自按 method 聚合, 然后共享一个 vmin/vmax + 共享 colorbar.

    Parameters
    ----------
    mcool_paths : list of str
        per-sample mcool 路径列表.
    loops_paths : list of str
        per-sample loops 文件路径列表 (顺序与 mcool_paths 对齐).
    expected_dfs : list of pd.DataFrame
        per-sample 全染色体 P(s) DataFrame (与 mcool_paths 对齐).
        建议由 cfizz.analyze.oe.compute_expected_cis_per_sample() 批量生成.
    output_path : str
        输出图片前缀 (无扩展名).
    sample_names : list of str, optional
        标题用的样本名列表. None → 从 mcool_paths 各自派生.
    其他参数: 同 compute_apa_oe_pileup_single.

    Returns
    -------
    output_paths : dict[str, str]
        {sample_name: output_path} 映射.
    """
    try:
        import cooler as cooler_lib
        import cooltools
    except ImportError:
        raise ImportError(
            "compute_apa_oe_pileup_multi 需要 cooler + cooltools"
        )

    from cfizz.analyze.compartment import get_view_df
    from cfizz.viz.layout import save_figure_multi_format, calculate_heatmap_layout

    assert len(mcool_paths) == len(loops_paths) == len(expected_dfs), \
        f"mcool_paths ({len(mcool_paths)}) / loops_paths ({len(loops_paths)}) / expected_dfs ({len(expected_dfs)}) 长度必须一致"

    n_samples = len(mcool_paths)
    if sample_names is None:
        sample_names = [os.path.basename(p).split('.')[0] for p in mcool_paths]

    # 1) per-sample 跑 cooltools.pileup + 聚合
    all_matrices = []
    all_loop_counts = []  # 跟踪每个 sample 的 loop 数 (title 显示用)
    for i in range(n_samples):
        mcool_path = mcool_paths[i]
        loops_path = loops_paths[i]
        expected_df = expected_dfs[i]
        sample_name = sample_names[i]

        loops_bedpe = _read_loops_for_cooltools_pileup(loops_path)
        print(f"  → {sample_name}: {len(loops_bedpe)} loops loaded for O/E pileup")
        all_loop_counts.append(len(loops_bedpe))

        clr = cooler_lib.Cooler(f"{mcool_path}::resolutions/{resolution}")
        view_df = get_view_df(clr)

        stack = cooltools.pileup(
            clr=clr,
            features_df=loops_bedpe,
            view_df=view_df,
            expected_df=expected_df,
            flank=flank,
            min_diag='auto',
            clr_weight_name='weight' if balance else None,
        )
        mtx = _aggregate_pileup_stack(stack, method=method)
        all_matrices.append(mtx)

    # 2) 颜色变换 (per-matrix)
    if color_scale == 'log2':
        all_matrices_disp = [
            np.log2(np.ma.masked_less_equal(m, 0).filled(1e-10))
            for m in all_matrices
        ]
        auto_label = 'log2(obs/exp)'
    elif color_scale == 'log10':
        all_matrices_disp = [
            np.log10(np.ma.masked_less_equal(m, 0).filled(1e-10))
            for m in all_matrices
        ]
        auto_label = 'log10(obs/exp)'
    else:  # 'linear'
        all_matrices_disp = all_matrices
        auto_label = 'obs/exp'

    final_cbar_label = cbar_label if cbar_label is not None else auto_label

    # 3) 多 panel 画图
    n_cols = min(8, n_samples)
    layout = calculate_heatmap_layout(n_plots=n_cols, plot_size=plot_size)
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5

    from cfizz.viz.heatmap import setup_axes, setup_colorbar
    sc = None
    # 用 zip 同时遍历 mtx_disp (变换后, 用于显示), mtx (原始 obs/exp, 用于算 score), 和 n_loops
    for i, (mtx_disp, mtx, sample_name, n_loops) in enumerate(
        zip(all_matrices_disp, all_matrices, sample_names, all_loop_counts)
    ):
        col = i % n_cols
        row = i // n_cols
        left = layout['margin_left']/layout['total_width'] + col * (layout['plot_width'] + layout['hspace'])
        bottom = layout['margin_bottom']/layout['total_height'] + (layout['n_rows'] - 1 - row) * (layout['plot_height'] + layout['vspace'])
        ax = fig.add_axes([left, bottom, layout['plot_width'], layout['plot_height']])

        sc = ax.imshow(mtx_disp, cmap=cmap, aspect='auto', interpolation='none',
                       vmin=vmin, vmax=vmax)
        setup_axes(ax)
        # 计算 APA score (基于未变换的 obs/exp 矩阵 mtx)
        # score = center_value / corner_mean, 反映 loop 中心相对角落的富集倍数
        # corner_size=3 匹配 analyze_apa 默认值
        D = mtx.shape[0]
        window_px = D // 2
        apa_corner_size = 3
        center_val = mtx[window_px, window_px]
        corner_mean = mtx[-apa_corner_size:, :apa_corner_size].mean()
        apa_score = center_val / corner_mean if corner_mean > 0 else 0.0

        # Title 风格 (模仿 visualize_multi_apa L557-563):
        # sample_name\nscore\n(N loops)
        title_text = f"{sample_name}\n{apa_score:.4g}"
        if n_loops > 0:
            title_text += f"\n({n_loops} loops)"
        ax.set_title(title_text, fontsize=5, pad=2)
        # "APA\nscore:" 标签只在第一个子图加
        if i == 0:
            ax.text(0.0, 1.01, "APA\nscore:", transform=ax.transAxes, fontsize=5, ha='left', va='bottom')

    # 4) 共享 colorbar
    if sc is not None:
        setup_colorbar(
            fig=fig, sc=sc, vmin=vmin, vmax=vmax, balance=balance,
            colorbar_left=layout['colorbar_left'],
            colorbar_bottom=0.75,
            colorbar_width=layout['colorbar_width_relative'],
            colorbar_height=0.15,
            label_fontsize=6,
            tick_fontsize=5,
            label=final_cbar_label,
        )

    # 5) 保存
    save_figure_multi_format(fig, output_path, dpi=dpi)
    return {sample_names[i]: output_path for i in range(n_samples)}
