#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
动态计算可视化范围的工具函数

从 script/a_8 / a_9 / a_10 移植,全部是纯算法,不依赖 cfizz 其他模块。
"""

from typing import Tuple


def calc_symmetric_range(
    center: int,
    half_span_bins: int,
    resolution: int
) -> Tuple[int, int]:
    """
    对称扩展可视化范围(以 center 为中心,上下游各 N × resolution)
    
    锚点: a_8 L228-245 calculate_visualization_range(算法 1:1,单位改成 bins)
    
    参数:
        center: 中心位置 (bp)
        half_span_bins: 上下游各扩展的 bin 数(N)
        resolution: 分辨率 (bp)
    
    返回:
        (start, end) — start 不小于 0
    
    行为特征:
        1. 输入: center (bp) + half_span_bins (整数 N) + resolution (bp)
        2. 输出: (start, end) 整数元组
        3. 跟 6_2 关系: 给 compartment 差异区域 + a_10 anchor 距离对称用
        4. 跟 a_8 算法关系: 1:1 算法,但单位从 bp 改 bins(更通用,换 1M 分辨率也合理)
        5. 跟其他函数关系: 跟 calc_anchor_range 配合(后者内部调它)
        6. 错误处理: start < 0 → max(0, ...);center < 0 抛 ValueError
    
    示例:
        >>> calc_symmetric_range(50_000_000, half_span_bins=25, resolution=100_000)
        (47_500_000, 52_500_000)  # 上下游各 2.5Mb @ 100kb
        >>> calc_symmetric_range(50_000_000, half_span_bins=25, resolution=1_000_000)
        (25_000_000, 75_000_000)  # 上下游各 25Mb @ 1Mb(换分辨率自动适配)
    """
    if center < 0:
        raise ValueError(f"center must be non-negative, got {center}")
    if half_span_bins < 0:
        raise ValueError(f"half_span_bins must be non-negative, got {half_span_bins}")
    if resolution <= 0:
        raise ValueError(f"resolution must be positive, got {resolution}")
    half_span_bp = half_span_bins * resolution
    start = max(0, center - half_span_bp)
    end = center + half_span_bp
    return int(start), int(end)


def calc_smart_union_range(
    s1_start: int, s1_end: int,
    s2_start: int, s2_end: int,
    s1_dist_up: float, s1_dist_down: float,
    s2_dist_up: float, s2_dist_down: float,
    default_flank: int,
    resolution: int
) -> Tuple[int, int]:
    """
    智能取并集(a_9 TAD 边界差异专用,通用)
    
    锚点: script/a_9 L166-223 calculate_smart_range
    
    逻辑:
        - sample1: 上游 = s1_start - s1_dist_up (inf → s1_start - default_flank)
                   下游 = s1_end + s1_dist_down (inf → s1_end + default_flank)
        - sample2: 同上
        - plot_start = min(s1_upstream, s2_upstream) 后 max(0, ...) 兜底
        - plot_end = max(s1_downstream, s2_downstream)
        - ± 1 bin (resolution) 安全边距
    
    参数: 8 个值 + default_flank + resolution
    返回: (plot_start, plot_end)
    
    行为特征:
        1. 输入: sample1/sample2 各自的 (start, end) + (dist_up, dist_down) + default_flank + resolution
        2. 输出: (plot_start, plot_end) 整数元组
        3. 跟 6_2 关系: 给 a_9 TAD 边界差异区域用(6_1 没生成 visualization_regions 文件,6_2 自包含)
        4. 跟 a_9 算法关系: 1:1 复刻 L177-223
        5. 跟其他函数关系: 独立
        6. 错误处理: np.isinf(...) 兜底;plot_start < 0 → max(0, ...)
    """
    import numpy as np
    
    s1_upstream = s1_start - s1_dist_up if not np.isinf(s1_dist_up) else s1_start - default_flank
    s1_downstream = s1_end + s1_dist_down if not np.isinf(s1_dist_down) else s1_end + default_flank
    s2_upstream = s2_start - s2_dist_up if not np.isinf(s2_dist_up) else s2_start - default_flank
    s2_downstream = s2_end + s2_dist_down if not np.isinf(s2_dist_down) else s2_end + default_flank
    
    plot_start = min(s1_upstream, s2_upstream)
    plot_end = max(s1_downstream, s2_downstream)
    
    plot_start = max(0, plot_start)
    plot_start = plot_start - resolution
    plot_end = plot_end + resolution
    
    return int(plot_start), int(plot_end)


def calc_anchor_range(start1: int, end1: int, start2: int, end2: int) -> Tuple[int, int]:
    """
    Anchor 距离对称扩展(loop / APA / pileup 通用)
    
    锚点: script/a_10 L120-130 calculate_visualization_range
    
    算法:
        center1 = (start1 + end1) / 2
        center2 = (start2 + end2) / 2
        anchor_distance = center2 - center1
        start = max(0, center1 - anchor_distance)
        end = center2 + anchor_distance
    
    参数: 4 个 BEDPE 边界值
    返回: (start, end) 整数元组
    
    行为特征:
        1. 输入: BEDPE 两端的 (start, end)
        2. 输出: (start, end) 整数元组
        3. 跟 6_2 关系: 给 a_10 loop 差异区域用
        4. 跟 a_10 算法关系: 1:1 复刻 L120-130
        5. 跟其他函数关系: 内部调 calc_symmetric_range(逻辑等价)
        6. 错误处理: start < 0 → max(0, ...);anchor_distance < 0 抛 ValueError
    """
    center1 = (start1 + end1) / 2
    center2 = (start2 + end2) / 2
    anchor_distance = center2 - center1
    if anchor_distance < 0:
        raise ValueError(f"start2/end2 must be after start1/end1, got anchor_distance={anchor_distance}")
    
    start = max(0, int(center1 - anchor_distance))
    end = int(center2 + anchor_distance)
    return start, end
