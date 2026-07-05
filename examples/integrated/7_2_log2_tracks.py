#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
7_2: 多组学整合可视化 + log2 fold change track(基于 7_1)

在 7_1 基础上,把最后 2 条 RNA track 换成 log2 fold change:
- ATAC-Seq log2(hiPSC_var/hiPSC_nor) bw(idx 6)
- RNA log2(hiPSC_var/hiPSC_nor) bw(idx 7)

log2 文件特性:
- 数据范围有正值(实验组 > 对照组)和负值(实验组 < 对照组)
- 不是单 sample,是 "实验组 vs 对照组 差异"
- 本次不指定 min/max,看 log2 数据自然范围(自动算 + 20% headroom)

复刻 l_2_multi_gene_visualization.py (T-9.14 已迁到 demo 数据集)
的所有业务逻辑,但完全基于 cfizz()。

数据:iPSC 多组学(hiPSC_nor + hiPSC_var 三组学) — T-9.14 切到 demo 数据集
- Hi-C: cfizz/demo/data/hiPSC_{nor,var}_chr17.mcool
- GTF: cfizz/demo/data/{gene}.gtf
- BigWig(全基因组): cfizz/demo/data/hiPSC_{nor,var}_{ATAC-Seq,CUTTag-CTCF,CUTTag-H3K27ac,chr17}_mean.bw
- Log2 bw(区域限定 75.4-76.34M): cfizz/demo/data/{ATAC-Seq,RNA}_hiPSC_var_div_hiPSC_nor_log2.75.4-76.34M.bw
- 默认基因 = FOXJ1(5kb,在 940kb 范围内)
- 默认区域 = chr17:75400000-76340000(940kb,统一区域)

使用:
    # 直接跑(无 CLI,参数硬编码)
    python 7_2_log2_tracks.py
    # 产出:cfizz/demo/output/7_2_log2_tracks/

设计:
- 基于 7_1 模板,8 条 BW 顺序:ATAC×2 → H3K27ac×2 → CTCF×2 → log2(ATAC+RNA)×2
- USE_TRACK_YLIM = False:不指定 min/max,看 log2 自然范围
- import 改为 cfizz.api.integrated.{quick_plot_integrated, GenomeRange}
- 默认 width=8 / left_margin=1.0 / right_margin=2.0(跟 l_2)
- 产物目录:examples_output/7_2/
"""

import os
import sys
import tempfile
from pathlib import Path

# T-9.14: 让 python 7_2_log2_tracks.py 直接跑(不需要 PYTHONPATH)
# 相对路径:7_× 在 cfizz/examples/integrated/<name>.py → 上 4 层到项目根目录
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from cfizz.api.integrated import quick_plot_integrated
from cfizz.api.integrated.heatmap_tracks import GenomeRange



# ============================================================
# T-7.8 verify: 测试负值 min_value 接口能力
# ============================================================
# - 前 6 条原 BW min=0(不变,ATAC/H3K27AC/CTCF 不该有负值)
# - 后 2 条 log2 min=-0.3(ATAC)/-1.5(RNA)— 测试负 min_value 能不能生效
# ============================================================
USE_TRACK_YLIM = True
TRACK_YLIM = {
    "max_value": [0.5, 0.5, 0.200, 0.200, 0.5, 0.5, 0.1, 0.6],  # 8 条全设,后 2 条给 log2 留够空间
    "min_value": [0, 0, 0, 0, 0, 0, -0.1, -0.1],               # 后 2 条设负值
}

# ============================================================================
# T-9.14: 配置区域(相对路径 demo 数据集)
# ============================================================================
DEMO_DATA = Path(__file__).resolve().parent.parent.parent / "demo" / "data"

# Hi-C (.mcool) 文件
HIC_FILES = {
    "hiPSC_nor": DEMO_DATA / "hiPSC_nor_chr17.mcool",  # 原 C7
    "hiPSC_var": DEMO_DATA / "hiPSC_var_chr17.mcool",  # 原 D11
}

# GTF 目录(已经按基因切好)
GENE_GTF_DIR = DEMO_DATA  # FOXJ1.gtf / 19_genes_subset.gtf 等直接放在这里

# Enhancers BED 文件(T-9.14: 用区域子集)
ENHANCER_BED_FILE = DEMO_DATA / "concordant_enhancer.chr17_75.4-76.34M.bed"

# CUTTag/ATAC BigWig 目录(同 l_2:chr17 数据)
MERGED_BW_DIR = DEMO_DATA

# RNA-Seq BigWig 目录(同 l_2:chr17 数据)
RNA_BW_DIR = DEMO_DATA

# 输出目录(T-9.14: 相对 demo/output)
OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "demo" / "output" / "7_2_log2_tracks"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 区域扩展配置(bp)
SCALE_CONFIG = {
    "small": 100_000,   # 100kb
    "large": 1_000_000,  # 1Mb
}

# 染色体最大长度映射
CHR_MAX_LENGTH = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569,
}

# 默认 flank 配置
DEFAULT_FLANK = 500_000  # 500kb

# 样本颜色配置(T-9.14: 跟 sample name 同步)
SAMPLE_COLORS = {
    "hiPSC_nor": "#666666",  # 原 C7 灰色(对照组)
    "hiPSC_var": "#C0392B",  # 原 D11 深红(实验组)
}

# 比较组配置(T-9.14: 跟 sample name + comparison key 同步)
# 顺序:hiPSC_nor 在前(对照组),hiPSC_var 在后(实验组)
COMPARISONS = {
    "hiPSC_var--hiPSC_nor": {
        "hic": ["hiPSC_nor", "hiPSC_var"],
        "atac": ["hiPSC_nor", "hiPSC_var"],
        "h3k27ac": ["hiPSC_nor", "hiPSC_var"],
        "ctcf": ["hiPSC_nor", "hiPSC_var"],
        "rna": ["hiPSC_nor", "hiPSC_var"],
    },
}


# ============================================================================
# 辅助函数(同 l_2 风格)
# ============================================================================
def get_hic_path(sample: str) -> str:
    """获取 Hi-C mcool 路径"""
    hic_file = HIC_FILES[sample]
    return str(hic_file) if hic_file.exists() else None


def get_bw_path(sample: str, assay_type: str, chrom: str) -> str:
    """获取 BigWig 文件路径(同 l_2 文件命名约定)"""
    if assay_type == "rna":
        bw_file = RNA_BW_DIR / f"{sample}_{chrom}_mean.bw"
        if bw_file.exists():
            return str(bw_file)
        return None

    if assay_type == "atac":
        bw_file = MERGED_BW_DIR / f"{sample}_ATAC-Seq_{chrom}_mean.bw"
    elif assay_type == "h3k27ac":
        bw_file = MERGED_BW_DIR / f"{sample}_CUTTag-H3K27ac_{chrom}_mean.bw"
    elif assay_type == "ctcf":
        bw_file = MERGED_BW_DIR / f"{sample}_CUTTag-CTCF_{chrom}_mean.bw"
    else:
        return None

    if bw_file.exists():
        return str(bw_file)
    return None


# T-7.8: 新加 log2 bw 路径函数(独立于 get_bw_path,因为 log2 文件不需要 sample+chrom)
def get_log2_bw_path(assay_type: str) -> str:
    """T-7.8: 获取 log2 fold change bw 文件路径(hiPSC_var vs hiPSC_nor)

    log2 bw 文件特性:
    - 全基因组一份(不按 chrom 分片)
    - 不分 sample(是 hiPSC_var vs hiPSC_nor 差异)
    - 值范围有正有负(hiPSC_var > hiPSC_nor 为正,hiPSC_var < hiPSC_nor 为负)
    """
    if assay_type == "atac_log2":
        bw_file = DEMO_DATA / "ATAC-Seq_hiPSC_var_div_hiPSC_nor_log2.75.4-76.34M.bw"  # T-9.14
    elif assay_type == "rna_log2":
        bw_file = DEMO_DATA / "RNA_hiPSC_var_div_hiPSC_nor_log2.75.4-76.34M.bw"  # T-9.14
    else:
        return None
    if bw_file.exists():
        return str(bw_file)
    return None


def get_gene_gtf_path(gene_name: str) -> str:
    """获取基因 GTF 路径(已切好)"""
    gene_gtf = GENE_GTF_DIR / f"{gene_name}.gtf"
    return str(gene_gtf) if gene_gtf.exists() else None


def get_enhancers_bed_path() -> str:
    """获取 Enhancers BED 文件路径"""
    if ENHANCER_BED_FILE.exists():
        return str(ENHANCER_BED_FILE)
    return None


def extract_gene_region(gene_name: str) -> dict:
    """从 GTF 提取基因坐标"""
    gtf_path = GENE_GTF_DIR / f"{gene_name}.gtf"
    if not gtf_path.exists():
        return {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue
            attrs = {}
            for kv in parts[8].strip(";").split("; "):
                if " " in kv:
                    k, v = kv.split(" ", 1)
                    attrs[k.replace('"', "")] = v.replace('"', "")
            if attrs.get("gene_name") == gene_name:
                return {
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "gene_name": gene_name,
                }
    return {}


def parse_region(region_str: str) -> tuple:
    """解析 chr:start-end,限制在染色体边界"""
    parts = region_str.split(":")
    if len(parts) != 2:
        raise ValueError(f"区域格式错误: {region_str}")
    chrom, start_end = parts
    start, end = map(int, start_end.split("-"))

    chr_max = CHR_MAX_LENGTH.get(chrom, 83257441)
    start = max(1, start)
    end = min(chr_max, end)
    if start >= end:
        start = max(1, end - 100000)

    return chrom, start, end


def calculate_region(gene_info: dict, flank: int) -> str:
    """根据基因信息计算区域"""
    chrom = gene_info["chrom"]
    start = max(1, gene_info["start"] - flank)
    end = gene_info["end"] + flank
    return f"{chrom}:{start}-{end}"


def merge_gtf_files(gene_names: list, output_path: str) -> bool:
    """合并多个基因的 GTF 文件(同 l_2)"""
    all_lines = []
    for gene_name in gene_names:
        gtf_path = get_gene_gtf_path(gene_name)
        if gtf_path and Path(gtf_path).exists():
            with open(gtf_path, "r") as f:
                all_lines.extend(f.readlines())
            print(f"  [OK] 添加GTF: {gene_name}")
    if all_lines:
        with open(output_path, "w") as f:
            f.writelines(all_lines)
        print(f"  [合并] 已合并 {len(gene_names)} 个基因的GTF到: {output_path}")
        return True
    return False




# ============================================================================
# 主可视化函数
# ============================================================================
def run_multi_gene_visualization(
    gene_names: list, comparison: str, region: str, **kwargs
):
    """多基因多组学可视化(同 l_2)"""
    chrom, start, end = parse_region(region)
    genes_str = (
        "_".join(gene_names[:3])
        if len(gene_names) <= 3
        else f"{gene_names[0]}_etc{len(gene_names)}"
    )

    print(f"\n{'=' * 60}")
    print(f"多基因可视化: {', '.join(gene_names)} ({region})")
    print(f"比较组: {comparison}")
    print(f"{'=' * 60}")

    samples = COMPARISONS[comparison]

    # ============ 准备 Hi-C configs ============
    hics = []
    for i, sample in enumerate(samples["hic"]):
        hic_path = get_hic_path(sample)
        if hic_path:
            hics.append({
                "file": hic_path,
                "triangle_ratio": 1,
                "cmap": "Reds",
                "color_scale": "linear",
                "balance": True,
                "name": sample,
                "flip_vertical": (i == 1),  # hiPSC_var 翻转(跟 l_2 一致)
            })
            print(f"  [OK] Hi-C {sample}")
        else:
            print(f"  [警告] Hi-C {sample}: 文件不存在")

    # ============ 准备 Tracks(同 l_2 顺序) ============
    tracks = []
    temp_gtf = None

    # 1. GTF track(贴近 Hi-C,第一个)
    if len(gene_names) == 1:
        # 单基因,直接用现成 gtf
        gtf_path = get_gene_gtf_path(gene_names[0])
        if gtf_path:
            tracks.append({
                "file": gtf_path,
                "color": "#666666",
                "name": None,
                "height_cm": kwargs.get("gtf_height", 0.5),
                "gtf_style": "flybase",
                "labels": True,
                "fontsize": kwargs.get("fontsize", 5),
                "color_utr": "blue",
                "border_color": "black",
                "color_backbone": "black",
                "line_width": 1.0,
            })
            print(f"  [OK] GTF: {gene_names[0]}")
        else:
            print(f"  [警告] GTF {gene_names[0]} 不存在")
    else:
        # 多基因,合并 GTF
        temp_gtf = tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False)
        temp_gtf.close()
        if merge_gtf_files(gene_names, temp_gtf.name):
            tracks.append({
                "file": temp_gtf.name,
                "color": "#666666",
                "name": None,
                "height_cm": kwargs.get("gtf_height", 0.5) * min(len(gene_names), 7),
                "gtf_style": "flybase",
                "labels": True,
                "fontsize": kwargs.get("fontsize", 5),
                "color_utr": "blue",
                "border_color": "black",
                "color_backbone": "black",
                "line_width": 1.0,
            })

    # 2. Enhancers BED track(同 l_2)
    enhancer_bed_path = get_enhancers_bed_path()
    if enhancer_bed_path:
        tracks.append({
            "file": enhancer_bed_path,
            "color": "purple",
            "name": "Enhancers",
            "height_cm": 0.8,
            "labels": False,
            "fontsize": 5,
            "line_width": 0.5,
        })
        print(f"  [OK] Enhancers BED")
    else:
        print(f"  [警告] Enhancers BED: 文件不存在")

    # 3. 前 6 条原 BW(ATAC + H3K27ac + CTCF × hiPSC_nor/hiPSC_var)+ 后 2 条 log2(ATAC + RNA × hiPSC_var/hiPSC_nor)
    # T-7.8: 不再走 rna mean,后 2 条改 log2 fold change(可负值)
    bw_assay_order = ["atac", "h3k27ac", "ctcf"]

    bw_track_idx = 0
    for assay in bw_assay_order:
        for sample in samples[assay]:
            bw_path = get_bw_path(sample, assay, chrom)
            if not bw_path:
                print(f"  [警告] {assay.upper()} {sample}: 文件不存在")
                continue

            track_dict = {
                "file": bw_path,
                "color": SAMPLE_COLORS[sample],
                "name": f"{sample} {assay.upper()}",
                "height_cm": kwargs.get("track_height", 0.8),
                "min_value": 0,
            }

            # T-7.5b: 注入 max_value
            # T-7.6: 用脚本开头 TRACK_YLIM dict 注入每条 BW 的 ymin/ymax(可短不可长)
            injected = False
            if TRACK_YLIM:
                if "max_value" in TRACK_YLIM and bw_track_idx < len(TRACK_YLIM["max_value"]):
                    track_dict["max_value"] = TRACK_YLIM["max_value"][bw_track_idx]
                    injected = True
                if "min_value" in TRACK_YLIM and bw_track_idx < len(TRACK_YLIM["min_value"]):
                    track_dict["min_value"] = TRACK_YLIM["min_value"][bw_track_idx]
                    injected = True
                if injected:
                    print(f"  [T-7.6] 注入 BW track {bw_track_idx}: max={track_dict.get('max_value')}, min={track_dict.get('min_value')}")
                else:
                    print(f"  [T-7.6] BW track {bw_track_idx}: TRACK_YLIM 已设但无该 idx,保持默认 (min={track_dict.get('min_value')}, max=自动)")
            else:
                print(f"  [T-7.6] BW track {bw_track_idx}: TRACK_YLIM=None,保持默认 (min={track_dict.get('min_value')}, max=自动)")

            tracks.append(track_dict)
            bw_track_idx += 1
            print(f"  [OK] {assay.upper()} {sample}")

    # T-7.8: 加 2 条 log2 fold change track(idx 6: ATAC log2, idx 7: RNA log2)
    # 关键:不设 min_value=0(让 matplotlib 自动算 min,可能为负)
    # T-9.14: log2 display name 跟 sample name 同步(D11/C7 → hiPSC_var/hiPSC_nor)
    log2_assays = [
        ("atac_log2", "ATAC log2(hiPSC_var/hiPSC_nor)", "black"),
        ("rna_log2", "RNA log2(hiPSC_var/hiPSC_nor)", "darkred"),
    ]
    for log2_assay, log2_name, log2_color in log2_assays:
        bw_path = get_log2_bw_path(log2_assay)
        if not bw_path:
            print(f"  [警告] {log2_assay}: 文件不存在")
            continue

        track_dict = {
            "file": bw_path,
            "color": log2_color,
            "name": log2_name,
            "height_cm": kwargs.get("track_height", 0.8),
            # 不设 min_value=0 —— log2 数据需要看负值范围
        }

        # T-7.6 注入逻辑保留(跟原 6 条一致)
        injected = False
        if TRACK_YLIM:
            if "max_value" in TRACK_YLIM and bw_track_idx < len(TRACK_YLIM["max_value"]):
                track_dict["max_value"] = TRACK_YLIM["max_value"][bw_track_idx]
                injected = True
            if "min_value" in TRACK_YLIM and bw_track_idx < len(TRACK_YLIM["min_value"]):
                track_dict["min_value"] = TRACK_YLIM["min_value"][bw_track_idx]
                injected = True
            if injected:
                print(f"  [T-7.6] 注入 BW track {bw_track_idx}: max={track_dict.get('max_value')}, min={track_dict.get('min_value')}")
            else:
                print(f"  [T-7.6] BW track {bw_track_idx}: TRACK_YLIM 已设但无该 idx,保持默认 (min={track_dict.get('min_value')}, max=自动)")
        else:
            print(f"  [T-7.6] BW track {bw_track_idx}: TRACK_YLIM=None,保持默认 (min={track_dict.get('min_value')}, max=自动)")

        tracks.append(track_dict)
        bw_track_idx += 1
        print(f"  [OK] {log2_assay}: {log2_name}")

    if not hics:
        print("错误: 没有有效的Hi-C文件")
        if temp_gtf and Path(temp_gtf.name).exists():
            os.unlink(temp_gtf.name)
        return

    # 输出文件名
    output_dir = Path(kwargs.get("output_dir", OUTPUT_DIR))
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = str(
        output_dir / f"MultiGene_{chrom}_{start}_{end}_{genes_str}_{comparison}"
    )
    print(f"  输出: {output_prefix}")

    # 执行可视化(完全用 cfizz 的 quick_plot_integrated)
    quick_plot_integrated(
        hics=hics,
        tracks=tracks,
        region=GenomeRange(chrom, start, end),
        output=output_prefix,
        width_cm=kwargs.get("width", 8),         # l_2 默认 8
        gap_cm=kwargs.get("gap", 0.1),
        left_margin_cm=kwargs.get("left_margin", 1.0),
        right_margin_cm=kwargs.get("right_margin", 2.0),
        dpi=kwargs.get("dpi", 300),
    )

    print(f"  完成: {output_prefix}.png")

    # 清理临时文件
    if temp_gtf and Path(temp_gtf.name).exists():
        os.unlink(temp_gtf.name)


# ============================================================================
# T-9.14: 直接入口(硬编码参数,无 argparse)
# ============================================================================
def main():
    # 7_× 统一范围(940kb)
    GENES = ["FOXJ1"]  # 7_2 单基因
    REGION = "chr17:75400000-76340000"  # 940kb / 94 bin 整数倍
    COMPARISON = "hiPSC_var--hiPSC_nor"

    print("=" * 60)
    print("7_2: 多组学整合可视化 + log2 fold change track(基于 7_1)")
    print(f"  genes: {GENES}")
    print(f"  region: {REGION}")
    print(f"  comparison: {COMPARISON}")
    print(f"  output_dir: {OUTPUT_DIR}")
    print("=" * 60)

    run_multi_gene_visualization(
        gene_names=GENES,
        comparison=COMPARISON,
        region=REGION,
        output_dir=str(OUTPUT_DIR),
        width=8, gap=0.1, track_height=0.8, gtf_height=0.5,
        fontsize=5, left_margin=1.0, right_margin=2.0, dpi=300,
    )

    print("\n" + "=" * 60)
    print("完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
