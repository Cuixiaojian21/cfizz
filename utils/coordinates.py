"""
Coordinate and utility functions for cfizz package.

"""

import os
import numpy as np


def ensure_dir(path):
    """确保目录存在"""
    os.makedirs(path, exist_ok=True)


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


def get_res_str(value, unit='bp'):
    """
    将数值转换为带单位的字符串
    
    Parameters:
        value: int or float, 要转换的数值
        unit: str, 单位类型，可选值：
            - 'bp': 碱基对，转换为kb/Mb
            - 'pos': 位置，转换为kb/Mb
            - 'other': 其他单位，保持原样
    
    Returns:
        str: 格式化的字符串
    """
    if unit == 'bp':
        if value >= 1e6:
            return f"{value/1e6:.0f}Mb"
        elif value >= 1e3:
            return f"{value/1e3:.0f}kb"
        else:
            return f"{value}bp"
    elif unit == 'pos':
        if value >= 1e6:
            return f"{value/1e6:.1f}Mb"
        elif value >= 1e3:
            return f"{value/1e3:.1f}kb"
        else:
            return f"{value}bp"
    else:
        return str(value)


def get_matrix_range(matrix, vmin=None, vmax=None):
    """
    获取矩阵的值范围，确保使用非零值
    
    参数:
        matrix: np.ndarray, 矩阵
        vmin, vmax: 如果提供，则直接使用；否则自动计算
    
    返回:
        vmin, vmax: 矩阵的值范围
    """
    if vmin is None or vmax is None:
        # 获取所有非零值
        nonzero = matrix[np.nonzero(matrix)]
        
        if vmin is None:
            # 使用非零值的最小值
            vmin = np.min(nonzero) if nonzero.size > 0 else 0
            
        if vmax is None:
            if nonzero.size > 0:
                # 使用93百分位数
                perc = np.percentile(nonzero, 93)
                # 找到最接近perc的实际值
                vmax = nonzero[np.abs(nonzero - perc).argmin()]
            else:
                vmax = 1
                
    return vmin, vmax


def generate_output_filename(sample_name, resolution, chrom, start_pos, end_pos, prefix='heatmap', balance=False):
    """
    生成统一的输出文件名
    
    Parameters:
        sample_name: str, 样品名称
        resolution: int, 分辨率（bp）
        chrom: str, 染色体号
        start_pos: int, 起始位置（bp）
        end_pos: int, 结束位置（bp）
        prefix: str, 文件名前缀，默认为'heatmap'
        balance: bool, 是否使用balance矩阵，默认为False
    
    Returns:
        str: 格式化的文件名（不含扩展名）
        格式：{prefix}_{sample_name}_{resolution}_{chrom}_{start}-{end}
        例如：single_Nalm6.400M_100kb_chr2_100.0Mb-110.0Mb
    """
    res_str = get_res_str(resolution, unit='bp')
    start_str = get_res_str(start_pos, unit='pos')
    end_str = get_res_str(end_pos, unit='pos')
    balance_str = 'balanced' if balance else 'raw'
    return f"{prefix}_{sample_name}_{res_str}_{chrom}_{start_str}-{end_str}_{balance_str}"