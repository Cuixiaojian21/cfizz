#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BEDPE 格式工具函数

从 script/a_10 移植,用于生成单 loop 标注文件。
"""

from typing import Any


def make_loop_region_id(row: Any) -> str:
    """
    生成 loop 的 region id: {chrom}_{start1}_{end1}_{start2}_{end2}
    
    锚点: script/a_10 L115-117 make_loop_region_id
    
    参数:
        row: 含有 chrom1/start1/end1/start2/end2 属性的对象(DataFrame row 或 dict)
    
    返回:
        region id 字符串
    
    行为特征:
        1. 输入: 1 个 row(含 chrom1/start1/end1/start2/end2)
        2. 输出: region id 字符串
        3. 跟 6_2 关系: 给 plot_multi_heatmap_with_loops 生成文件名用
        4. 跟 a_10 算法关系: 1:1 复刻
        5. 跟其他函数关系: 跟 write_single_loop_bedpe 配合(后者用这个做文件名)
        6. 错误处理: 缺列时抛 KeyError/AttributeError
    """
    return f"{row['chrom1']}_{int(row['start1'])}_{int(row['end1'])}_{int(row['start2'])}_{int(row['end2'])}"


def write_single_loop_bedpe(row: Any, output_path: str) -> str:
    """
    写单 loop 的 BEDPE 16 列标注文件
    
    锚点: script/a_10 L133-142 generate_single_loop_file
    
    16 列格式(HiCCUPS 标准):
        chrom1 / start1 / end1 / chrom2 / start2 / end2 / . / . / . / .
        / donut_fe / pvalue1 / donut_q / ll_fe / pvalue2 / ll_q
    
    参数:
        row: 含有 chrom1/start1/end1/chrom2/start2/end2/donut_fe/ll_fe 的对象
        output_path: 输出文件路径
    
    返回:
        输出文件路径
    
    行为特征:
        1. 输入: 1 行 loop(含 chrom1/start1/end1/chrom2/start2/end2/donut_fe/ll_fe)
        2. 输出: 文件路径
        3. 跟 6_2 关系: 给 plot_multi_heatmap_with_loops 喂单 loop 标注用
        4. 跟 a_10 算法关系: 1:1 复刻 L138-141
        5. 跟其他函数关系: 跟 make_loop_region_id 配合(后者生成文件名)
        6. 错误处理: 缺列时用 .get(..., 0) 兜底
    """
    output_path = str(output_path)
    with open(output_path, 'w') as f:
        f.write(
            f"{row['chrom1']}\t{row['start1']}\t{row['end1']}\t"
            f"{row['chrom2']}\t{row['start2']}\t{row['end2']}\t"
            f".\t.\t.\t.\t"
            f"{row.get('donut_fe', 0)}\t0\t0\t"
            f"{row.get('ll_fe', 0)}\t0\t0\n"
        )
    return output_path
