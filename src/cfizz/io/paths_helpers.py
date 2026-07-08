"""
便捷 Helper 函数(for 8_1 等 Example 快速使用)

每个 helper = 组合 StagePath + FeaturePath + make_filename,
暴露业务友好的属性 API。

用法:
    >>> from cfizz.io.paths_helpers import compartment, tad, loop
    >>> p = compartment("50_1", "/tmp/8_1/1_computation/2_compartment", "chr1", 1_000_000)
    >>> p.eig_tsv
    '/tmp/8_1/1_computation/2_compartment/1_2_eigenvector.50_1.chr1.1.0M.tsv'
"""

import os as _os
from .paths import StagePath, FeaturePath, ComputeFeature


# =============================================================================
# 通用 res_str 计算
# =============================================================================

def _res_str(resolution: int) -> str:
    """计算分辨率字符串(1Mb → '1.0M', 100kb → '100k')"""
    if resolution == 1_000_000:
        return "1.0M"
    elif resolution >= 1_000_000:
        return f"{resolution // 1_000_000}.0M"
    else:
        return f"{resolution // 1000}k"


def _window_str(window_bp: int) -> str:
    """计算 window 字符串(100000 → '100kb')"""
    return f"{window_bp // 1000}kb"


# =============================================================================
# Compartment Helper
# =============================================================================

def compartment(
    sample_name: str,
    output_dir: str,
    chrom: str,
    resolution: int,
) -> dict:
    """
    Stage 1 compartment 产物路径 helper
    
    产物命名(process_compartment + load_or_compute_oe_matrix 实际输出):
    - eigenvector: {output_dir}/eigenvector.{res_str}.tsv
                    例: eigenvector.1.0M.tsv
    - gc_cov:       {output_dir}/gc_cov.{res_str}.tsv
                    例: gc_cov.1.0M.tsv
    - oe_npy:       {output_dir}/{sample}_{chrom}_{res_str_m}.oe.npy
                    例: 50_1_chr1_1M.oe.npy  (注意:oe 用 "1M" 不是 "1.0M")
    - oe_meta:      {output_dir}/{sample}_{chrom}_{res_str_m}.oe_meta.json
    其中 res_str = "1.0M" / "100k", res_str_m = "1M" / "100k"
    
    Returns:
        dict with keys: eig_tsv, gc_cov_tsv, oe_npy, oe_meta_json
    """
    res_str = _res_str(resolution)          # "1.0M" / "100k"
    res_str_m = f"{resolution // 1_000_000}M" if resolution >= 1_000_000 else f"{resolution // 1000}k"
    
    eig_tsv = _os.path.join(output_dir, f"eigenvector.{res_str}.tsv")
    gc_cov_tsv = _os.path.join(output_dir, f"gc_cov.{res_str}.tsv")
    # oe 用 sample_chrom_1M 格式
    oe_npy = _os.path.join(output_dir, f"{sample_name}_{chrom}_{res_str_m}.oe.npy")
    oe_meta_json = _os.path.join(output_dir, f"{sample_name}_{chrom}_{res_str_m}.oe_meta.json")
    
    return {
        "eig_tsv": eig_tsv,
        "gc_cov_tsv": gc_cov_tsv,
        "oe_npy": oe_npy,
        "oe_meta_json": oe_meta_json,
    }


# =============================================================================
# TAD Helper
# =============================================================================

def tad(
    sample_name: str,
    output_dir: str,
    resolution: int,
    windows: list,
    mcool_path: str = None,
) -> dict:
    """
    Stage 1 TAD 产物路径 helper
    
    产物命名(process_tads 实际输出):
    - insulation: {output_dir}/{basename}/1_0.{basename}.{resolution}.insulation.tsv
                   例: 1_0.50_1_1000.10000.insulation.tsv
                   (basename = mcool_path.split('.')[0])
    - boundaries: {output_dir}/{basename}/2_0.{basename}.{resolution}.{window}kb.boundaries.tsv
                  例: 2_0.50_1_1000.10000.100kb.boundaries.tsv
    - tads:       {output_dir}/{basename}/2_1.{basename}.{resolution}.{window}kb.tads.tsv
    
    mcool_path 可选:若提供,basename = mcool_path.split('.')[0]
                    若不提供,basename = sample_name
    """
    basename = sample_name if mcool_path is None else _os.path.basename(mcool_path).split('.')[0]
    res_str = str(resolution)      # "10000"
    
    # insulation → 直接在 output_dir 下(无子目录)
    ins_tsv = _os.path.join(output_dir, f"1_0.{basename}.{res_str}.insulation.tsv")
    
    # boundaries + tads → 同目录
    def boundaries_tsv(window_bp: int) -> str:
        ws = _window_str(window_bp)
        return _os.path.join(output_dir, f"2_0.{basename}.{res_str}.{ws}.boundaries.tsv")
    
    def tads_tsv(window_bp: int) -> str:
        ws = _window_str(window_bp)
        return _os.path.join(output_dir, f"2_1.{basename}.{res_str}.{ws}.tads.tsv")
    
    return {
        "basename": basename,
        "insulation_tsv": ins_tsv,
        "boundaries_tsv": boundaries_tsv,   # callable: boundaries_tsv(window_bp)
        "tads_tsv": tads_tsv,               # callable: tads_tsv(window_bp)
    }


# =============================================================================
# Loop Helper
# =============================================================================

def loop(
    basename: str,
    output_dir: str,
    resolution: int,
) -> dict:
    """
    Stage 1 Loop 产物路径 helper
    
    产物命名(process_loops 实际输出):
    - loops: {output_dir}/{basename}.{res_k}k.loops.txt
             例: 50_1_1000.10k.loops.txt
             (basename = mcool_path.split('.')[0],无子目录)
    """
    res_k = resolution // 1000    # 10_000 → 10
    loops_txt = _os.path.join(output_dir, f"{basename}.{res_k}k.loops.txt")
    
    return {
        "basename": basename,
        "loops_txt": loops_txt,
    }
