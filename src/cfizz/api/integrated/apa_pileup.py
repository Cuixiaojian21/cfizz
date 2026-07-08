"""
APA (Aggregate Peak Analysis) end-to-end visualization module.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, Dict
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



def analyze_apa(submatrices: List[np.ndarray], window: int = 5, corner_size: int = 3) -> Tuple[np.ndarray, float, float, float, float]:
    """
    分析APA结果（主流定义，兼容旧代码）

    Args:
        submatrices: 子矩阵列表
        window: 中心点坐标（如5）
        corner_size: 角落区域大小（如3）

    Returns:
        avg: 平均信号矩阵
        score: APA分数（中心点/左下角均值）
        z: z分数
        p: p值
        maxi: 最大值（右上角均值*5）
    """
    if not submatrices:
        return None, 0, 0, 1, 0

    # 去除异常值
    mean_arr = np.array([np.mean(arr) for arr in submatrices])
    p99 = np.percentile(mean_arr, 99)
    p1 = np.percentile(mean_arr, 1)
    mask = (mean_arr < p99) & (mean_arr > p1)
    avg = np.mean(np.array(submatrices)[mask], axis=0)

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
    balance: bool = False
) -> plt.Figure:
    """绘制APA热图
    
    Args:
        matrix: APA分析得到的平均信号矩阵
        vmin: 颜色条最小值
        vmax: 颜色条最大值
        cmap: 颜色映射名称
        dpi: 图片分辨率
        balance: 是否使用平衡矩阵
        
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
    
    # 获取矩阵范围
    vmin, vmax = get_matrix_range(matrix, vmin, vmax)
    
    # 绘制热图
    sc = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='none', 
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
        label=''  # 不显示标签
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
    n_processes: int = None
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
    avg, score, z, p, maxi = analyze_apa(submatrices, window, corner_size)
    
    # 4. 绘制热图
    fig = plot_apa_heatmap_visualization(
        matrix=avg,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        dpi=dpi,
        balance=balance
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
    sample_names: List[str] = None
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
        
    Returns:
        all_avg_matrices: 所有样本的平均信号矩阵列表
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
        avg, score, z, p, maxi = analyze_apa(submatrices, window, corner_size)
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
    loop_counts: List[int] = None
) -> None:
    """可视化多个样本的APA分析结果
    
    Args:
        all_avg_matrices: 所有样本的平均信号矩阵列表
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
    
    # 计算所有矩阵的统一颜色范围
    all_values = np.concatenate([
        matrix.flatten() for matrix in all_avg_matrices
        if matrix is not None and not (isinstance(matrix, np.ndarray) and np.all(np.isnan(matrix)))
    ])
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
        sc = ax.imshow(avg, cmap=cmap, aspect='auto', interpolation='none', 
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
            label=''  # 不显示标签
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
    plot_size: float = 2
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
        sample_names=sample_names
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
        loop_counts=all_loop_counts
    )
    
    # 3. 返回APA分数
    return dict(zip(used_sample_names, all_scores))
