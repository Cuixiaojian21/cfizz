"""
Path management for cfizz outputs.

集中管理产物路径生成,避免 example / 上层代码手写路径导致 bug。

设计哲学(见 DESIGN_PHILOSOPHY §4):
- P4.1: x_y 编号系统,stage(1-4) + 子步骤 y(从 1 起)
- P4.2: 文件名格式 = `stage_y_特征.限定.扩展`
- P4.3: 目录格式 = `stage_x_特征名/`
- P4.4: Stage 3 累计数据 + 图同目录
- P4.5: Stage 2 viz 文件前缀用 `2_y` 即使在 3_aggregation/
- P4.6: 三层抽象 StagePath → FeaturePath → make_filename

使用示例:
    >>> from cfizz.io.paths import StagePath, FeaturePath, ComputeFeature, make_filename
    >>> sp = StagePath(output_dir="/tmp/8_1")
    >>> fp = FeaturePath(stage_dir=sp.computation_dir, feature=ComputeFeature.COMPARTMENT)
    >>> fp.file_for(file_feature="eigenvector", qualifiers=["50_1", "1Mb"], ext="tsv")
    '/tmp/8_1/1_computation/2_compartment/1_2_eigenvector.50_1.1Mb.tsv'
"""

from enum import IntEnum
from dataclasses import dataclass
from typing import Union as _Union, Protocol
import os as _os
from pathlib import Path
import inspect


# =============================================================================
# Stage + Feature 枚举
# =============================================================================

class Stage(IntEnum):
    """4 个 stage 大类"""
    COMPUTE = 1       # 数据生产
    VISUALIZE = 2     # 单特征 viz
    AGGREGATE = 3     # 累计分析
    DIFFERENTIAL = 4  # 差异分析(未来)


class ComputeFeature(IntEnum):
    """Stage 1 数据生产:4 个特征"""
    DISTANCE = 1
    COMPARTMENT = 2
    TAD = 3
    LOOP = 4


class VizFeature(IntEnum):
    """Stage 2 viz:4 个 viz 类型"""
    HEATMAP = 1
    COMPARTMENT = 2
    TAD = 3
    LOOP = 4


class AggregateFeature(IntEnum):
    """Stage 3 累计:3 个累计类型"""
    SADDLE = 1
    TAD_PILEUP = 2
    APA = 3


# =============================================================================
# 类型别名(模块内使用)
# =============================================================================

# Union of all feature enums for type annotation
Feature = ComputeFeature | VizFeature | AggregateFeature  # type: ignore[operator, misc]


# =============================================================================
# x 编号查表(每 stage 内部,特征名 string → x int)
# 用 string key 而非 enum key,避免 VizFeature(ComputeFeature.COMPARTMENT) 互转报错
# key 必须跟 self.feature.name.lower() 完全对齐(DISTANCE → "distance")
# =============================================================================

_COMPUTE_X = {
    "distance": 1, "compartment": 2,
    "tad": 3, "loop": 4,
}
_VIZ_X = {
    "heatmap": 1, "compartment": 2,
    "tad": 3, "loop": 4,
}
_AGGREGATE_X = {
    "saddle": 1, "tad_pileup": 2, "apa": 3,
}

# Stage → 查表字典 的映射(单源真相)
_STAGE_TO_X_MAP = {
    Stage.COMPUTE: _COMPUTE_X,
    Stage.VISUALIZE: _VIZ_X,
    Stage.AGGREGATE: _AGGREGATE_X,
}

# Stage → 该 stage 支持的特征名集合(运行时 assert 用)
_EXPECTED_FEATURE_NAMES = {
    Stage.COMPUTE: {"distance", "compartment", "tad", "loop"},
    Stage.VISUALIZE: {"heatmap", "compartment", "tad", "loop"},
    Stage.AGGREGATE: {"saddle", "tad_pileup", "apa"},
}


# =============================================================================
# 静态查表(feature → sub 编号)
# (Stage, feature_name, file_feature) → sub
# =============================================================================

_FILE_SUB_MAP = {
    # Stage 1: compartment 子产物
    (Stage.COMPUTE, "compartment", "gc_cov"): 1,
    (Stage.COMPUTE, "compartment", "eigenvector"): 2,
    (Stage.COMPUTE, "compartment", "oe"): 3,
    (Stage.COMPUTE, "compartment", "oe_meta"): 4,
    # Stage 1: tad 子产物
    (Stage.COMPUTE, "tad", "insulation"): 1,
    (Stage.COMPUTE, "tad", "boundaries"): 2,
    (Stage.COMPUTE, "tad", "tads"): 3,
    # Stage 1: loop 子产物
    (Stage.COMPUTE, "loop", "loops"): 1,
    # Stage 1: distance 子产物
    (Stage.COMPUTE, "distance", "distance"): 1,
    # Stage 2: heatmap viz
    (Stage.VISUALIZE, "heatmap", "single_heatmap"): 1,
    (Stage.VISUALIZE, "heatmap", "multi_heatmap"): 2,
    (Stage.VISUALIZE, "heatmap", "45deg_heatmap"): 3,
    # Stage 2: compartment viz
    (Stage.VISUALIZE, "compartment", "compartment_heatmap"): 1,
    (Stage.VISUALIZE, "compartment", "compartment_heatmap_45deg"): 2,
    # Stage 2: tad viz
    (Stage.VISUALIZE, "tad", "tad_heatmap"): 1,
    (Stage.VISUALIZE, "tad", "tad_heatmap_45deg"): 2,
    # Stage 2: loop viz
    (Stage.VISUALIZE, "loop", "loop_heatmap"): 1,
    # Stage 3: saddle
    (Stage.AGGREGATE, "saddle", "saddle_matrix"): 1,
    (Stage.AGGREGATE, "saddle", "saddle_heatmap"): 1,
    # Stage 3: tad_pileup
    (Stage.AGGREGATE, "tad_pileup", "tad_pileup"): 1,
    (Stage.AGGREGATE, "tad_pileup", "boundary_pileup_heatmap"): 1,
    # Stage 3: apa
    (Stage.AGGREGATE, "apa", "apa"): 1,
    (Stage.AGGREGATE, "apa", "apa_heatmap"): 1,
}


# =============================================================================
# Stage 大目录类
# =============================================================================

@dataclass(frozen=True)
class StagePath:
    """
    4 个 stage 大目录的根路径
    
    Attributes:
        output_dir: 根目录(用户给的)
    
    Properties:
        computation_dir: 1_computation/ 绝对路径
        visualization_dir: 2_visualization/
        aggregation_dir: 3_aggregation/
        differential_dir: 4_differential/
    
    Example:
        >>> sp = StagePath(output_dir="/tmp/8_1")
        >>> sp.computation_dir
        '/tmp/8_1/1_computation'
    """
    output_dir: str
    
    @property
    def computation_dir(self) -> str:
        return f"{self.output_dir}/1_computation"
    
    @property
    def visualization_dir(self) -> str:
        return f"{self.output_dir}/2_visualization"
    
    @property
    def aggregation_dir(self) -> str:
        return f"{self.output_dir}/3_aggregation"
    
    @property
    def differential_dir(self) -> str:
        return f"{self.output_dir}/4_differential"
    
    def dir_of(self, stage: Stage) -> str:
        """给定 stage,返回对应大目录路径"""
        return {
            Stage.COMPUTE: self.computation_dir,
            Stage.VISUALIZE: self.visualization_dir,
            Stage.AGGREGATE: self.aggregation_dir,
            Stage.DIFFERENTIAL: self.differential_dir,
        }[stage]


# =============================================================================
# 文件名生成器
# =============================================================================

def make_filename(stage: int, sub: int, feature: str,
                  qualifiers: list, ext: str) -> str:
    """
    生成文件名 stage_y_特征.限定.扩展
    
    Args:
        stage: 1 / 2 / 3 / 4
        sub: 该 stage 内的子步骤编号(从 1 起)
        feature: 特征名(eigenvector / oe / saddle_matrix / ...)
        qualifiers: 限定列表 [sample, res, chrom, window, ...](按顺序拼)
        ext: 扩展(tsv / npy / png / pkl / json)
    
    Returns:
        文件名字符串(stage_y_特征.限定.扩展)
    
    Example:
        >>> make_filename(stage=1, sub=2, feature="eigenvector",
        ...              qualifiers=["50_1", "1Mb"], ext="tsv")
        '1_2_eigenvector.50_1.1Mb.tsv'
        >>> make_filename(stage=2, sub=1, feature="compartment_heatmap",
        ...              qualifiers=["50_1", "chr1", "1Mb"], ext="png")
        '2_1_compartment_heatmap.50_1.chr1.1Mb.png'
    """
    parts = [f"{stage}_{sub}", feature] + [str(q) for q in qualifiers if q is not None]
    return ".".join(parts) + f".{ext}"



# =============================================================================
# 特征二级目录类(核心)
# =============================================================================

@dataclass(frozen=True)
class FeaturePath:
    """
    特征二级目录(给定 stage 大目录 + 特征 enum → 二级目录路径)
    
    Attributes:
        stage_dir: 上层的 stage 大目录路径(如 /tmp/8_1/1_computation)
        feature: 特征 enum(ComputeFeature / VizFeature / AggregateFeature 之一,IntEnum 类型)
        stage: Stage 枚举(默认 COMPUTE)
    
    Properties:
        x: 特征 x 编号(查表)
        subdir: 完整二级目录路径(x_特征名,**不带 stage 前缀**,stage 由父目录承担)
        feature_name: 特征名小写字符串(从 enum.name.lower() 推导)
    
    Methods:
        file_for(file_feature, qualifiers, ext): 自动查表得 sub,生成完整文件路径
        file(sub, file_feature, qualifiers, ext): 显式指定 sub(高级用法)
    
    Example:
        >>> sp = StagePath(output_dir="/tmp/8_1")
        >>> fp = FeaturePath(stage_dir=sp.computation_dir, feature=ComputeFeature.COMPARTMENT)
        >>> fp.subdir
        '/tmp/8_1/1_computation/2_compartment'              ← x=2 (compartment),不带 stage 前缀
        >>> fp.file_for(file_feature="eigenvector", qualifiers=["50_1", "1Mb"], ext="tsv")
        '/tmp/8_1/1_computation/2_compartment/1_2_eigenvector.50_1.1Mb.tsv'
    """
    stage_dir: str
    feature: Feature
    stage: Stage = Stage.COMPUTE
    
    @property
    def feature_name(self) -> str:
        """特征名小写字符串(从 enum.name.lower() 推导,如 COMPARTMENT → 'compartment')"""
        return self.feature.name.lower()
    
    @property
    def x(self) -> int:
        """特征 x 编号(根据 stage + feature.name 查表)"""
        fname = self.feature.name.lower()  # "compartment" / "tad" / ...
        # 运行时 assert:检查 feature 是否跟 stage 匹配
        expected = _EXPECTED_FEATURE_NAMES.get(self.stage)
        if expected is None:
            raise ValueError(f"Stage {self.stage} (x 查表) 待 T-9.3 后续扩展")
        if fname not in expected:
            raise ValueError(
                f"feature '{fname}' 与 stage '{self.stage.name}' 不匹配。"
                f"Stage {self.stage.name} 支持: {sorted(expected)}"
            )
        # 查表(全部用 string key,无 enum 互转问题)
        mapping = _STAGE_TO_X_MAP[self.stage]
        return mapping[fname]
    
    @property
    def subdir(self) -> str:
        """完整二级目录路径 = stage_dir / x_特征名(不带 stage 前缀)"""
        return f"{self.stage_dir}/{self.x}_{self.feature_name}"
    
    def file_for(self, file_feature: str, qualifiers: list, ext: str) -> str:
        """
        静态查表:给定 file_feature → 自动得 sub → 生成完整文件路径
        
        静态映射(本任务定):
        compartment 子产物: gc_cov=1, eigenvector=2, oe=3, oe_meta=4
        tad 子产物: insulation=1, boundaries=2, tads=3
        loop 子产物: loops=1
        heatmap viz: single_heatmap=1, multi_heatmap=2, 45deg_heatmap=3
        compartment viz: compartment_heatmap=1, compartment_heatmap_45deg=2
        tad viz: tad_heatmap=1, tad_heatmap_45deg=2
        loop viz: loop_heatmap=1
        saddle aggregate: saddle_matrix=1, saddle_heatmap=1
        tad_pileup aggregate: tad_pileup=1, boundary_pileup_heatmap=1
        apa aggregate: apa=1, apa_heatmap=1
        """
        sub = _FILE_SUB_MAP.get((self.stage, self.feature_name, file_feature))
        if sub is None:
            raise KeyError(
                f"未定义 ({self.stage}, '{self.feature_name}', '{file_feature}') 的 sub 编号。"
                f"请在 _FILE_SUB_MAP 加映射"
            )
        return f"{self.subdir}/{make_filename(int(self.stage), sub, file_feature, qualifiers, ext)}"
    
    def file(self, sub: int, file_feature: str, qualifiers: list, ext: str) -> str:
        """
        显式指定 sub(高级用法,绕过 file_for 的查表)
        """
        return f"{self.subdir}/{make_filename(int(self.stage), sub, file_feature, qualifiers, ext)}"


# =============================================================================
# 动态 sub 注册器(高级用法)
# =============================================================================

class SubstepRegistry:
    """
    动态追踪某目录内"下一个可用 sub 编号"
    
    扫描目录下现有的 `stage_y_xxx.{ext}` 文件,返回 max(sub) + 1
    
    高级用法(适用于"sub 不固定,按调用顺序"场景)
    默认推荐用 FeaturePath.file_for()(静态查表,更稳)
    
    Example:
        >>> reg = SubstepRegistry("/tmp/8_1/1_computation/1_2_compartment")
        >>> reg.next_sub  # 若目录已有 1_0_gc_cov.tsv, 1_2_eigenvector.tsv, 返回 4
    """
    def __init__(self, subdir: str):
        self.subdir = subdir
    
    @property
    def next_sub(self) -> int:
        """扫描目录下所有 stage_y_xxx.ext 文件,返回 max(sub) + 1;空目录返回 1"""
        max_y = 0
        if not _os.path.isdir(self.subdir):
            return 1
        for fname in _os.listdir(self.subdir):
            # 解析 "stage_y_xxx.ext" 中的 y
            parts = fname.split("_", 2)
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                max_y = max(max_y, int(parts[1]))
        return max_y + 1 if max_y > 0 else 1


# =============================================================================
# 函数签名反射(高级用法)
# =============================================================================

def extract_qualifiers(func, locals_dict: dict) -> dict:
    """
    从函数签名 + 调用参数自动提取 qualifier 字典
    
    用于"qualifier 自动从函数参数提取"的场景
    
    Args:
        func: 函数对象(用 func.__name__ 或直接传函数)
        locals_dict: 函数的 locals() 字典
    
    Returns:
        dict: {参数名: str(value)} 排除 self / args / kwargs
    
    Example:
        >>> def compute_eigenvector(mcool_path, resolution, sample_name):
        ...     qualifiers = extract_qualifiers(compute_eigenvector, locals())
        ...     return qualifiers
        >>> compute_eigenvector("50_1/50_1.mcool", 1_000_000, "50_1")
        # qualifiers = {"mcool_path": "50_1/50_1.mcool", "resolution": "1000000",
        #                "sample_name": "50_1"}
    """
    sig = inspect.signature(func)
    params = sig.parameters
    result = {}
    for name, value in locals_dict.items():
        if name in ("self", "cls", "args", "kwargs"):
            continue
        if name in params:
            result[name] = str(value)
    return result


# =============================================================================
# PathBuilder Protocol(扩展性)
# =============================================================================

# =============================================================================
# 便捷 Helper 函数(for 8_1 等 Example 快速使用)
# 每个 helper = 组合 StagePath + FeaturePath + make_filename,暴露业务友好的 dict API
# 用法: from cfizz.io.paths import compartment, tad, loop
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
    - oe_npy:       {output_dir}/{sample}_{chrom}_{res_str_m}.oe.npy
                    例: 50_1_chr1_1M.oe.npy  (注意:oe 用 "1M" 不是 "1.0M")
    - oe_meta:      {output_dir}/{sample}_{chrom}_{res_str_m}.oe_meta.json
    
    Returns:
        dict: {"eig_tsv", "gc_cov_tsv", "oe_npy", "oe_meta_json"}
    """
    res_str = _res_str(resolution)          # "1.0M" / "100k"
    res_str_m = f"{resolution // 1_000_000}M" if resolution >= 1_000_000 else f"{resolution // 1000}k"
    
    return {
        "eig_tsv": _os.path.join(output_dir, f"eigenvector.{res_str}.tsv"),
        "gc_cov_tsv": _os.path.join(output_dir, f"gc_cov.{res_str}.tsv"),
        "oe_npy": _os.path.join(output_dir, f"{sample_name}_{chrom}_{res_str_m}.oe.npy"),
        "oe_meta_json": _os.path.join(output_dir, f"{sample_name}_{chrom}_{res_str_m}.oe_meta.json"),
    }


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
    - boundaries: {output_dir}/{basename}/2_0.{basename}.{resolution}.{window}kb.boundaries.tsv
    
    mcool_path 可选:若提供,basename = mcool_path.split('.')[0];
                    若不提供,basename = sample_name
    
    Returns:
        dict: {"basename", "insulation_tsv", "boundaries_tsv", "tads_tsv"}
        boundaries_tsv 和 tads_tsv 是 callable(window_bp) → str
    """
    basename = sample_name if mcool_path is None else _os.path.basename(mcool_path).split('.')[0]
    res_str = str(resolution)
    
    def boundaries_tsv(window_bp: int) -> str:
        ws = _window_str(window_bp)
        return _os.path.join(output_dir, f"2_0.{basename}.{res_str}.{ws}.boundaries.tsv")
    
    def tads_tsv(window_bp: int) -> str:
        ws = _window_str(window_bp)
        return _os.path.join(output_dir, f"2_1.{basename}.{res_str}.{ws}.tads.tsv")
    
    return {
        "basename": basename,
        "insulation_tsv": _os.path.join(output_dir, f"1_0.{basename}.{res_str}.insulation.tsv"),
        "boundaries_tsv": boundaries_tsv,
        "tads_tsv": tads_tsv,
    }


def loop(
    basename: str,
    output_dir: str,
    resolution: int,
) -> dict:
    """
    Stage 1 Loop 产物路径 helper
    
    产物命名(process_loops 实际输出):
    - loops: {output_dir}/{basename}.{res_k}k.loops.txt
    
    basename = mcool_path.split('.')[0](如 '50_1_1000'),不是 sample_name
    
    Returns:
        dict: {"basename", "loops_txt"}
    """
    res_k = resolution // 1000
    return {
        "basename": basename,
        "loops_txt": _os.path.join(output_dir, f"{basename}.{res_k}k.loops.txt"),
    }


class PathBuilder(Protocol):
    """
    所有 path 类的契约(扩展性预留)
    
    未来加新 path 类(如 DifferentialPaths),只要实现 paths() 返回 dict,
    调用方就可以通过 dict 模式访问路径,无需修改现有代码
    
    Example:
        >>> class DifferentialPaths:
        ...     def paths(self) -> dict[str, str]:
        ...         return {"diff_eigenvector": self.diff_eigenvector_tsv, ...}
        >>> dp = DifferentialPaths(...)
        >>> dp.paths()["diff_eigenvector"]
        '/tmp/8_1/4_differential/4_1_diff_eigenvector/diff_eigenvector.50_1.tsv'
    """
    def paths(self) -> dict[str, str]: ...
