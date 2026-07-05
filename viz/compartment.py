"""
Compartment visualization module for cfizz package.
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize

from .plot_style import (
    setup_plot_style,
    calculate_heatmap_layout,
    setup_axes,
    add_coordinate_labels,
    add_bar_coordinates,
    setup_diff_colorbar,
    log2_and_mask,
    save_figure,
)

from cfizz.analyze.compartment import process_compartment


def calculate_chromosome_data(
    mcool_path: str,
    chrom: str,
    resolution: int,
    fasta_path: str,
    output_dir: str,
    sample_name: str = "sample",
    n_eigs: int = 3,
    force_recompute: bool = False,
    load_only: bool = False,
) -> Tuple[pd.DataFrame, np.ndarray]:
    res_str = "1.0M" if resolution == 1_000_000 else f"{resolution // 1000}k"
    eig_file = os.path.join(output_dir, f'1_1.{sample_name}.{res_str}.E1.tsv')
    prefix = f"{sample_name}_{res_str}_{chrom}_whole"
    oe_file = os.path.join(output_dir, f"1_3.{prefix}.observed_expected_normalized.txt")
    
    # load_only 模式：只读取，不计算
    if load_only:
        if os.path.exists(eig_file) and os.path.exists(oe_file):
            print(f"\n[load_only] 加载已有计算结果: {output_dir}")
            try:
                eig_df = pd.read_csv(eig_file, sep='\t')
                oe_matrix = np.loadtxt(oe_file)
                return eig_df, oe_matrix
            except Exception as e:
                raise FileNotFoundError(f"load_only 模式：无法加载 {eig_file} 或 {oe_file}: {e}")
        else:
            raise FileNotFoundError(f"load_only 模式：文件不存在 - {eig_file} 或 {oe_file}")
    
    # 检查缓存
    if os.path.exists(eig_file) and os.path.exists(oe_file) and not force_recompute:
        print(f"\n发现已存在的计算结果，直接加载...")
        try:
            eig_df = pd.read_csv(eig_file, sep='\t')
            oe_matrix = np.loadtxt(oe_file)
            return eig_df, oe_matrix
        except Exception as e:
            print(f"加载已有结果时出错: {str(e)}")
            print("将重新计算...")
    
    # 计算 compartment
    print(f"\n开始计算 compartment...")
    eig_df = process_compartment(
        mcool_path=mcool_path,
        resolution=resolution,
        fasta_path=fasta_path,
        output_dir=output_dir,
        n_eigs=n_eigs
    )
    
    # 计算 O/E 矩阵（带 npy 缓存）
    print(f"\n开始计算 O/E 矩阵...")
    from cfizz.analyze.oe import load_or_compute_oe_matrix
    oe_matrix, oe_metadata = load_or_compute_oe_matrix(
        mcool_path=mcool_path,
        chrom=chrom,
        resolution=resolution,
        output_dir=output_dir,
        sample_name=sample_name,
        balance=False,  # 默认使用原始矩阵
        force_recompute=force_recompute
    )
    
    print("\n=== O/E 标准化矩阵统计信息 ===")
    for key, value in oe_metadata.get('stats', {}).items():
        if isinstance(value, (int, float)):
            print(f"{key}: {value:.6f}")
        else:
            print(f"{key}: {value}")
    
    return eig_df, oe_matrix


def plot_heatmap_with_e1(
    oe_matrix: np.ndarray,
    e1_df: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    output_prefix: str,
    vmin: float = -2,
    vmax: float = 2,
    cmap: Optional[LinearSegmentedColormap] = None,
    plot_size: float = 4.0,
    bar_height_ratio: float = 0.3
):
    """
    绘制 O/E 热图 + E1 柱状图
    
    Args:
        oe_matrix: O/E 矩阵
        e1_df: E1 主成分 DataFrame
        chrom: 染色体名称
        start: 起始位置
        end: 结束位置
        output_prefix: 输出文件前缀
        vmin: 颜色范围最小值
        vmax: 颜色范围最大值
        cmap: 颜色映射
        plot_size: 单个子图的基础尺寸（单位：cm）
        bar_height_ratio: 柱状图高度与热图高度的比例
    """
    # 数据检查
    print("\n=== 数据检查信息 ===")
    print(f"O/E 矩阵形状: {oe_matrix.shape}")
    print(f"E1 数据框形状: {e1_df.shape}")
    
    # 检查 E1 数据中的空值
    na_count = e1_df['E1'].isna().sum()
    if na_count > 0:
        print(f"警告：E1 数据中存在 {na_count} 个空值")
        e1_df['E1'] = e1_df['E1'].fillna(0)
        print("已将空值填充为 0")
    
    # 确保数据对齐
    if len(e1_df) != oe_matrix.shape[0]:
        print(f"警告：O/E 矩阵行数({oe_matrix.shape[0]})与 E1 数据长度({len(e1_df)})不匹配")
        min_len = min(len(e1_df), oe_matrix.shape[0])
        e1_df = e1_df.iloc[:min_len]
        oe_matrix = oe_matrix[:min_len, :min_len]
        print(f"已截取前{min_len}个 bin 的数据")
    
    # 设置统一的绘图样式
    setup_plot_style()
    
    if cmap is None:
        cmap = LinearSegmentedColormap.from_list('custom_coolwarm', ['blue', 'white', 'red'])
        cmap.set_bad(color='#EEEEEE')
    
    # 设置布局
    layout = calculate_heatmap_layout(
        n_plots=1,
        plot_size=plot_size,
        has_bar=True,
        bar_height_ratio=bar_height_ratio,
        margin_bottom=0.4
    )
    
    # 创建图形
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # 热图区域
    ax_heat = fig.add_axes([
        layout['margin_left_relative'],
        layout['margin_bottom_relative'] + layout['bar_height_relative'] + layout['gap_relative'],
        layout['plot_width'],
        layout['plot_height']
    ])
    
    # 绘制热图
    masked_matrix = log2_and_mask(oe_matrix)
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    sc = ax_heat.imshow(masked_matrix, cmap=cmap, aspect='auto', interpolation='none', norm=norm)
    
    # 设置坐标轴
    xmin, xmax, ymin, ymax = setup_axes(ax_heat)
    
    # 添加坐标标签
    if all(v is not None for v in [start, end, chrom]):
        add_coordinate_labels(ax_heat, xmin, xmax, ymin, ymax, start, end, chrom, fontsize=5)
    
    # 添加颜色条
    setup_diff_colorbar(
        fig=fig,
        sc=sc,
        vmin=vmin,
        vmax=vmax,
        colorbar_left=layout['colorbar_left'],
        colorbar_bottom=0.8,
        colorbar_width=layout['colorbar_width_relative'],
        colorbar_height=0.15,
        max_label=str(vmax),
        min_label=str(vmin),
        label_text='log2(obs/exp)',
        tick_fontsize=5,
        label_fontsize=5
    )
    
    # E1 柱状图
    ax_bar = fig.add_axes([
        layout['margin_left_relative'],
        0,
        layout['plot_width'],
        layout['bar_height_relative']
    ])
    
    # 绘制 E1 柱状图
    x = np.arange(len(e1_df))
    colors = ['red' if v > 0 else 'blue' for v in e1_df['E1']]
    ax_bar.bar(x, e1_df['E1'], width=1.0, color=colors, edgecolor='none')
    
    # 调整 x 轴范围
    ax_bar.set_xlim(-0.5, len(e1_df) - 0.5)
    
    # 设置 E1 柱状图的坐标轴
    setup_axes(ax_bar)
    
    # 添加柱状图坐标标签
    add_bar_coordinates(ax_bar, show_labels=True, fontsize=5, label="EV1")
    
    # 保存图片
    save_figure(fig, output_prefix)


def plot_compartment(
    eig_tsv_path: str,
    oe_npy_path: str,
    output_dir: str,
    chrom: str,
    start_pos: int = 0,
    end_pos: Optional[int] = None,
    resolution: int = 100000,
    vmin: float = -2,
    vmax: float = 2,
    sample_name: str = "sample",
    plot_size: float = 4.0,
    bar_height_ratio: float = 0.3
) -> Dict[str, Any]:
    """
    生成 compartment 分析的可视化(纯可视化函数,读已计算产物)。
    
    宪法对齐: §3 P3.1 算/可视分离。**此函数 = 多 sample 一条龙 generate_multi_compartment 的特殊情况**(单 sample)。
    **不调任何计算函数**(calculate_chromosome_data / process_compartment / load_or_compute_oe_matrix)。
    只读预计算的 eigenvector.tsv + oe.npy,做区域过滤 + 拼图 + save。
    
    Parameters
    ----------
    eig_tsv_path : str
        预计算的 eigenvector.{res}.tsv 路径(从 2_1_compartment_calculation 产物)
    oe_npy_path : str
        预计算的 {sample}_{chrom}_{res}.oe.npy 路径
    output_dir : str
        输出目录(自动 mkdir)
    chrom : str
        染色体号
    start_pos : int
        起始位置(bp),默认 0
    end_pos : int
        结束位置(bp),默认 None 表示整条 chr
    resolution : int
        分辨率(bp),用于 bin 计算
    vmin, vmax : float
        颜色条范围
    sample_name : str
        用于文件名
    plot_size : float
        图形尺寸(cm),默认 4.0
    bar_height_ratio : float
        E1 柱状图高度与热图高度的比例
    
    Returns
    -------
    Dict[str, str] : {'status', 'output_file', ...}
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if end_pos is None:
        import numpy as np
        first_oe = np.load(oe_npy_path)
        end_pos = first_oe.shape[0] * resolution
    
    # 直接调 slice_compartment_region(纯计算,单 sample)
    results = [slice_compartment_region((
        eig_tsv_path, oe_npy_path, sample_name, 0,  # idx=0
        chrom, resolution, start_pos, end_pos, output_dir, 3
    ))]
    
    # 拼图 + save(只调可视化函数)
    from .plot_style import calculate_heatmap_layout
    layout = calculate_heatmap_layout(
        n_plots=1, plot_size=plot_size, has_bar=True,
        bar_height_ratio=bar_height_ratio, margin_bottom=0.4
    )
    # 实际生产文件名
    res_str = "1.0M" if resolution == 1_000_000 else f"{resolution//1000}k"
    start_str = f"{start_pos//1_000_000}M" if start_pos >= 1_000_000 else f"{start_pos//1000}k"
    end_str = f"{end_pos//1_000_000}M" if end_pos >= 1_000_000 else f"{end_pos//1000}k"
    output_prefix = os.path.join(output_dir, f"{sample_name}_{chrom}_{res_str}_{start_str}-{end_str}")
    
    result = results[0]
    if result['status'] == 'success':
        plot_heatmap_with_e1(
            oe_matrix=result['oe_matrix'],
            e1_df=result['eig_df'],
            chrom=chrom,
            start=start_pos,
            end=end_pos,
            output_prefix=output_prefix,
            vmin=vmin, vmax=vmax,
            plot_size=plot_size,
            bar_height_ratio=bar_height_ratio
        )
    
    return {
        'status': result['status'],
        'output_file': f"{output_prefix}.png",
        'sample_name': sample_name,
        'chrom': chrom,
        'region': f"{start_pos}-{end_pos}",
        'resolution': resolution
    }


# =============================================================================
# 以下是 cfizz 特有的接口，保持原有风格（返回 Figure）
# =============================================================================

def plot_eigenvector(
    eigenvector_data: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    resolution: int,
    eig_column: str = "E1",
    color_up: str = "#E41A1C",
    color_down: str = "#377EB8",
    figsize: Tuple[float, float] = (10, 2),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    绘制特征向量轨道图
    
    Parameters
    ----------
    eigenvector_data : pd.DataFrame
        DataFrame with columns [chrom, start, end, E1, ...]
    chrom : str
        染色体名称
    start : int
        起始位置（bp）
    end : int
        结束位置（bp）
    resolution : int
        分辨率（bp）
    eig_column : str
        特征向量列名
    color_up : str
        正值颜色（A compartment）
    color_down : str
        负值颜色（B compartment）
    figsize : tuple
        图形大小
    title : str, optional
        标题
    save_path : str, optional
        保存路径
    dpi : int
        分辨率
        
    Returns
    -------
    fig : plt.Figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # 提取区域数据
    region_data = eigenvector_data[
        (eigenvector_data['chrom'] == chrom) &
        (eigenvector_data['start'] >= start) &
        (eigenvector_data['end'] <= end)
    ].copy()
    
    if len(region_data) == 0:
        raise ValueError(f"No data found for {chrom}:{start}-{end}")
    
    region_data[eig_column] = region_data[eig_column].fillna(0)
    
    # 创建条形图
    x = np.arange(len(region_data))
    values = region_data[eig_column].values
    colors = np.where(values > 0, color_up, color_down)
    
    ax.bar(x, values, width=1.0, color=colors, edgecolor='none')
    ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
    
    ax.set_xlim(-0.5, len(region_data) - 0.5)
    ax.set_xlabel(f'{chrom} position (bins)', fontsize=10)
    ax.set_ylabel('E1', fontsize=10)
    
    if title:
        ax.set_title(title)
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


# =============================================================================
# T-3.2/T-3.4 函数: slice_compartment_region / plot_multi_compartment / generate_multi_compartment
# =============================================================================

def slice_compartment_region(args):
    """
    对单 sample 的计算结果做区域过滤 + O/E 切片(多 sample 并行友好的纯计算函数)。
    
    设计哲学(算/可视分离宪法 P3.1):
      - 纯计算,无 plt/imshow/savefig
      - **不调任何计算函数**(calculate_chromosome_data)
      - **不调任何 hicviz 函数**
      - 输入: 
          - eig_tsv_path: 预计算的 eigenvector.{res}.tsv 路径(从 calculate_chromosome_data 产物)
          - oe_npy_path: 预计算的 {sample}_{chrom}_{res}.oe.npy 路径(从 calculate_chromosome_data 产物)
      - 输出: 切片后的 region_eig_df + region_oe_matrix
      - 附加: idx 用于多 sample 并行结果排序
      - 附加: try/except 包装,单 sample 失败不挂整个并行批次
    """
    # 兼容旧调用方式(10 个参数)
    if len(args) == 10:
        eig_tsv_path, oe_npy_path, sample_name, idx, chrom, resolution, start_pos, end_pos, output_dir, n_eigs = args
    else:
        # 新 11 个参数(包含 force_recompute 但已废弃)
        eig_tsv_path, oe_npy_path, sample_name, idx, chrom, resolution, start_pos, end_pos, output_dir, n_eigs, _force_recompute_unused = args
    
    try:
        # 读 eigenvector.tsv(全基因组)
        eig_df = pd.read_csv(eig_tsv_path, sep='\t')
        
        # 读 O/E npy(单 chr 全矩阵)
        oe_matrix = np.load(oe_npy_path)
        
        # 提取对应区域的 eig_df
        region_eig_df = eig_df[
            (eig_df['chrom'] == chrom) & 
            (eig_df['start'] >= start_pos) & 
            (eig_df['end'] <= end_pos)
        ].copy()
        
        # 提取对应区域的 O/E 矩阵
        start_bin = start_pos // resolution
        end_bin = end_pos // resolution
        region_oe_matrix = oe_matrix[start_bin:end_bin, start_bin:end_bin]
        
        # 确保数据对齐
        if len(region_eig_df) != region_oe_matrix.shape[0]:
            min_len = min(len(region_eig_df), region_oe_matrix.shape[0])
            region_eig_df = region_eig_df.iloc[:min_len]
            region_oe_matrix = region_oe_matrix[:min_len, :min_len]
        
        return {
            'status': 'success',
            'eig_df': region_eig_df,
            'oe_matrix': region_oe_matrix,
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


def plot_multi_compartment(
    results: List[Dict[str, Any]],
    output_prefix: str,
    vmin: float = -2,
    vmax: float = 2,
    plot_size: float = 3.0,
    bar_height_ratio: float = 0.3,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    chrom: Optional[str] = None,
    group_name: Optional[str] = None
) -> None:
    """
    绘制多个样本的 compartment 可视化(O/E 热图 + E1 柱状图,水平拼接)。
    
    设计哲学: 单/多样本统一。接收 `results: List[Dict]`,每个 Dict 来自 slice_compartment_region。
    **不调任何 hicviz 函数**,内部用本地 calculate_heatmap_layout / setup_diff_colorbar / save_figure。
    """
    # 设置布局(用本地 plot_style.calculate_heatmap_layout)
    n_samples = len(results)
    layout = calculate_heatmap_layout(
        n_plots=n_samples,
        plot_size=plot_size,
        has_bar=True,
        bar_height_ratio=bar_height_ratio,
        margin_bottom=0.3,
        n_cols=n_samples,
        n_rows=1
    )
    
    # 创建图形
    fig = plt.figure(figsize=(layout['fig_width'], layout['fig_height']))
    
    # 设置统一的绘图样式(本地 plot_style)
    setup_plot_style()
    
    # 创建颜色映射
    cmap = LinearSegmentedColormap.from_list('custom_coolwarm', ['blue', 'white', 'red'])
    cmap.set_bad(color='#EEEEEE')
    
    # 计算全局 E1 范围
    e1_max = None
    for result in results:
        if result['status'] == 'success':
            current_max = max(abs(result['eig_df']['E1'].min()), abs(result['eig_df']['E1'].max()))
            if e1_max is None or current_max > e1_max:
                e1_max = current_max
    
    # 向上取 0.5 的整数倍
    if e1_max is not None:
        e1_max = np.ceil(e1_max * 2) / 2
    
    # 绘制所有样本的热图和 E1 柱状图
    sc = None
    for result in results:
        if result['status'] == 'success':
            idx = result['idx']
            sample_name = result['sample_name']
            
            col = idx
            
            # 热图区域
            ax_heat = fig.add_axes([
                layout['margin_left_relative'] + col * (layout['plot_width'] + layout['gap_relative']),
                layout['margin_bottom_relative'] + layout['bar_height_relative'] + layout['gap_relative'],
                layout['plot_width'],
                layout['plot_height']
            ])
            
            # 绘制热图
            masked_matrix = log2_and_mask(result['oe_matrix'])
            norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
            sc = ax_heat.imshow(masked_matrix, cmap=cmap, aspect='auto', interpolation='none', norm=norm)
            
            # 设置坐标轴
            xmin, xmax, ymin, ymax = setup_axes(ax_heat)
            
            # 只在第一个子图添加坐标标签
            if idx == 0 and all(v is not None for v in [start_pos, end_pos, chrom]):
                add_coordinate_labels(ax_heat, xmin, xmax, ymin, ymax, start_pos, end_pos, chrom, fontsize=5)
            
            # 添加样本名称
            ax_heat.text(0.5, 1.02, sample_name, 
                        transform=ax_heat.transAxes, ha='center', va='bottom', fontsize=5)
            
            # E1 柱状图
            ax_bar = fig.add_axes([
                layout['margin_left_relative'] + col * (layout['plot_width'] + layout['gap_relative']),
                0,
                layout['plot_width'],
                layout['bar_height_relative']
            ])
            
            # 绘制 E1 柱状图
            x = np.arange(len(result['eig_df']))
            colors = ['red' if v > 0 else 'blue' for v in result['eig_df']['E1']]
            ax_bar.bar(x, result['eig_df']['E1'], width=1.0, color=colors, edgecolor='none')
            
            ax_bar.set_xlim(-0.5, len(result['eig_df']) - 0.5)
            setup_axes(ax_bar)
            
            if e1_max is not None:
                ax_bar.set_ylim(-e1_max, e1_max)
            
            add_bar_coordinates(ax_bar, show_labels=(idx==0), fontsize=5, label="EV1")
    
    # 添加颜色条
    if sc is not None:
        setup_diff_colorbar(
            fig=fig,
            sc=sc,
            vmin=vmin,
            vmax=vmax,
            colorbar_left=layout['colorbar_left'],
            colorbar_bottom=0.8,
            colorbar_width=layout['colorbar_width_relative'],
            colorbar_height=0.15,
            max_label=str(vmax),
            min_label=str(vmin),
            label_text='log2(obs/exp)',
            tick_fontsize=5,
            label_fontsize=5
        )
    
    # 保存图片
    if group_name and group_name not in output_prefix:
        save_figure(fig, f"{output_prefix}_{group_name}")
    else:
        save_figure(fig, output_prefix)


def generate_multi_compartment(
    eig_tsv_paths: List[str],
    oe_npy_paths: List[str],
    output_dir: str,
    sample_names: List[str],
    chrom: str = "chr8",
    resolution: int = 100000,
    start_pos: Optional[int] = None,
    end_pos: Optional[int] = None,
    n_eigs: int = 3,
    vmin: float = -2,
    vmax: float = 2,
    plot_size: float = 3.0,
    bar_height_ratio: float = 0.3,
    group_name: Optional[str] = None,
    max_workers: int = 1
) -> Dict[str, Any]:
    """
    生成多个样本的 compartment 可视化(纯可视化函数,读已计算产物)。
    
    宪法对齐: §3 P3.1 算/可视分离 + 单/多样本统一。
    **不调任何计算函数**(calculate_chromosome_data / process_compartment)。
    只读预计算的 eigenvector.tsv + oe.npy,做并行切片 + 拼图 + save。
    
    Parameters
    ----------
    eig_tsv_paths : List[str]
        预计算的 eigenvector.{res}.tsv 路径列表(每个 sample 一个)
    oe_npy_paths : List[str]
        预计算的 npy 路径列表(每个 sample 一个)
    output_dir : str
        输出目录
    sample_names : List[str]
        样本名称列表(长度 = len(eig_tsv_paths))
    chrom : str
        染色体号
    resolution : int
        分辨率(bp)
    start_pos, end_pos : int
        区域范围(默认 None = 整条 chr)
    n_eigs : int
        保留兼容,实际不再使用
    vmin, vmax : float
        颜色条范围
    plot_size : float
        多 sample 折中尺寸(cm),默认 3.0
    bar_height_ratio : float
        E1 柱状图高度比例
    group_name : str
        输出文件名前缀(用于多 sample 区分)
    max_workers : int
        保留兼容,切片是纯 IO,几乎不需要并行
    
    Returns
    -------
    Dict[str, Any] : {'status', 'output_prefix', 'results'}
    """
    from concurrent.futures import ThreadPoolExecutor
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 默认 region
    if start_pos is None:
        start_pos = 0
    if end_pos is None:
        first_oe = np.load(oe_npy_paths[0])
        end_pos = first_oe.shape[0] * resolution
    
    # 准备并行切片参数(只 IO + 切片,无需 ProcessPool)
    args = [
        (eig_tsv_paths[i], oe_npy_paths[i], sample_names[i], i, chrom, resolution,
         start_pos, end_pos, output_dir, n_eigs)
        for i in range(len(eig_tsv_paths))
    ]
    
    # ThreadPoolExecutor(纯 IO + 切片,Thread 即可)
    if max_workers == 1:
        results = [slice_compartment_region(a) for a in args]
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(slice_compartment_region, args))
    
    # 按 idx 排序
    results.sort(key=lambda x: x['idx'])
    
    # 生成输出 prefix
    res_str = "1.0M" if resolution == 1_000_000 else f"{resolution//1000}k"
    start_str = f"{start_pos//1_000_000}M" if start_pos >= 1_000_000 else f"{start_pos//1000}k"
    end_str = f"{end_pos//1_000_000}M" if end_pos >= 1_000_000 else f"{end_pos//1000}k"
    group_prefix = f"{group_name}_" if group_name else ""
    output_prefix = os.path.join(output_dir, f'{group_prefix}multi_compartment_{res_str}_{chrom}_{start_str}-{end_str}')
    
    # 调 plot_multi_compartment(纯拼图 + save)
    plot_multi_compartment(
        results=results,
        output_prefix=output_prefix,
        vmin=vmin, vmax=vmax,
        plot_size=plot_size,
        bar_height_ratio=bar_height_ratio,
        start_pos=start_pos, end_pos=end_pos,
        chrom=chrom, group_name=group_name
    )
    
    return {
        'status': 'success',
        'output_prefix': output_prefix,
        'results': results
    }


