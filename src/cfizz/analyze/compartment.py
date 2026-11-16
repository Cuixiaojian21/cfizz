"""
A/B compartment analysis for Hi-C data.

The module computes compartment eigenvectors from balanced cis contact
matrices with :func:`cooltools.eigs_cis`, assigns A/B labels, and provides
quality-control and differential-analysis helpers.

Design notes
------------
1. **E1 selection.** ``cooltools`` sorts eigenvectors by ``abs(eigenvalue)``.
   On a near-degenerate ``+a / -a`` pair this may select the negative-λ
   vector, which describes an "anti-compartment" pattern rather than A/B.
   ``eig_select="positive"`` promotes the largest positive-λ vector to E1.

2. **Region definition.** A whole-chromosome decomposition can put the
   p-arm/q-arm split on PC1 instead of the intra-arm checkerboard. With
   ``view="arm"`` each chromosome is split at its centromere (UCSC
   cytoBandIdeo ``acen`` midpoint) and p/q arms are decomposed separately.
   Arms with fewer than ``min_arm_bins`` valid bins stay as one region.

3. **Quality control.** :func:`compartment_qc` reports per-region
   checkerboard structure, GC correlation, block counts and arm metrics,
   so that silent failures are visible.

4. **Differential analysis.** :func:`calculate_compartment_diff` can ignore
   near-zero E1 sign flips with ``min_abs_e1`` and classify them as
   ``Uncertain`` instead of A/B switches.

The centromere annotation is read from ``demo/data/hg38.cytoBand.ideo.txt.gz``
in the source tree. Pass ``cytoband_file=`` or set ``CFIZZ_CYTOBAND`` when the
file lives elsewhere (e.g. after ``pip install``).
"""

from __future__ import annotations

import gzip
import os
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.linalg import eigh

warnings.filterwarnings('ignore', category=RuntimeWarning)

# =============================================================================
# Cytoband / centromere
# =============================================================================

CYTOBAND_FILENAME = "hg38.cytoBand.ideo.txt.gz"


def resolve_cytoband_path(
    cytoband_file: Optional[os.PathLike | str] = None,
) -> Optional[Path]:
    """Resolve the cytoBandIdeo path.

    Resolution order:

    1. explicit ``cytoband_file``;
    2. ``CFIZZ_CYTOBAND`` environment variable;
    3. ``<repo>/demo/data/hg38.cytoBand.ideo.txt.gz`` (source checkout);
    4. ``<cwd>/demo/data/hg38.cytoBand.ideo.txt.gz``.

    Returns ``None`` when no file is found.
    """
    candidates: List[Path] = []
    if cytoband_file is not None:
        path = Path(cytoband_file)
        if not path.is_file():
            raise FileNotFoundError(f"cytoBandIdeo file not found: {path}")
        return path
    env_path = os.environ.get("CFIZZ_CYTOBAND")
    if env_path:
        path = Path(env_path)
        if not path.is_file():
            raise FileNotFoundError(
                f"CFIZZ_CYTOBAND points to a missing file: {path}"
            )
        return path

    # src/cfizz/analyze/compartment.py -> repo root is parents[3]
    repo_root = Path(__file__).resolve().parents[3]
    candidates.append(repo_root / "demo" / "data" / CYTOBAND_FILENAME)
    candidates.append(Path.cwd() / "demo" / "data" / CYTOBAND_FILENAME)

    for path in candidates:
        if path.is_file():
            return path
    return None


def _format_resolution(resolution: int) -> str:
    """File-name token for a resolution: ``100k`` / ``1.0M``."""
    if resolution < 1_000_000:
        return f"{resolution // 1000}k"
    return f"{resolution / 1_000_000:.1f}M"



# =============================================================================
# =============================================================================

def calc_gc_cov(
    bins: pd.DataFrame,
    fasta_path: str,
    out_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Calculate GC content for genomic bins using bioframe.

    This function uses bioframe.frac_gc() to compute the GC content
    for each genomic bin in the provided bins DataFrame.

    Parameters
    ----------
    bins : pd.DataFrame
        DataFrame with columns: chrom, start, end (bin coordinates)
    fasta_path : str
        Path to reference genome FASTA file
    out_path : str, optional
        Path to save GC content results

    Returns
    -------
    gc_cov : pd.DataFrame
        DataFrame with columns: chrom, start, end, gc (GC content)

    Examples
    --------
    >>> clr = cooler.Cooler("sample.mcool::resolutions/100000")
    >>> bins = clr.bins()[:]
    >>> gc_cov = calc_gc_cov(bins, "hg38.fa")
    >>> print(f"GC range: {gc_cov['gc'].min():.3f} to {gc_cov['gc'].max():.3f}")
    """
    try:
        import bioframe
    except ImportError:
        raise ImportError("bioframe required: pip install bioframe")

    # Load genome FASTA
    hg_genome = bioframe.load_fasta(fasta_path)

    # Calculate GC content
    gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg_genome)

    # Save if path provided
    if out_path:
        gc_cov.to_csv(out_path, index=False, sep='\t')

    return gc_cov


def get_view_df(clr) -> pd.DataFrame:
    """
    Get chromosome view DataFrame from a cooler object.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object

    Returns
    -------
    view_df : pd.DataFrame
        DataFrame with columns: chrom, start, end, name
    """
    return pd.DataFrame({
        'chrom': clr.chromnames,
        'start': 0,
        'end': clr.chromsizes.values,
        'name': clr.chromnames
    })


# =============================================================================
# Centromere coordinates and arm-level view
# =============================================================================

def _chrom_key(chrom: str) -> str:
    """Normalize ``chr1`` / ``1`` for cross-reference matching."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _load_cytoband(
    path: Path,
) -> Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, int]]:
    """Read acen bands and chromosome sizes from a UCSC cytoBandIdeo file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    bands: Dict[str, List[Tuple[int, int]]] = {}
    sizes: Dict[str, int] = {}
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, start, end, _band, giemsa = (
                parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            )
            sizes[chrom] = max(sizes.get(chrom, 0), end)
            if "acen" in giemsa:
                bands.setdefault(chrom, []).append((start, end))
    return bands, sizes


def _cytoband_matches_clr(
    band_sizes: Dict[str, int],
    clr,
    rel_tol: float = 0.002,
) -> bool:
    """Check that cytoband and cooler describe the same assembly."""
    matched = 0
    checked = 0
    by_key = {_chrom_key(chrom): chrom for chrom in band_sizes}
    for chrom, size in clr.chromsizes.items():
        key = by_key.get(_chrom_key(chrom))
        if key is None:
            continue
        checked += 1
        cb_size = band_sizes[key]
        tol = max(100_000, rel_tol * int(size))
        if abs(cb_size - int(size)) <= tol:
            matched += 1
    return checked > 0 and matched >= max(1, checked // 2)


def load_centromere_positions(
    cytoband_file: Optional[os.PathLike | str] = None,
) -> Dict[str, int]:
    """Read centromere positions from a UCSC cytoBandIdeo file.

    The centromere is the midpoint of the chromosome's ``acen`` band
    interval(s).

    Parameters
    ----------
    cytoband_file : path-like, optional
        cytoBandIdeo file (gzip or plain text). When omitted, the resolver
        checks ``CFIZZ_CYTOBAND`` and then ``demo/data/hg38.cytoBand.ideo.txt.gz``.

    Returns
    -------
    dict
        ``{chrom: centromere_position_bp}``
    """
    path = resolve_cytoband_path(cytoband_file)
    if path is None:
        raise FileNotFoundError(
            "cytoBandIdeo file not found. Pass cytoband_file=..., set "
            "CFIZZ_CYTOBAND, or download it with:\n"
            "  curl -o demo/data/hg38.cytoBand.ideo.txt.gz \\\n"
            "    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz"
        )

    bands, _sizes = _load_cytoband(path)
    return {
        chrom: (min(s for s, _ in v) + max(e for _, e in v)) // 2
        for chrom, v in bands.items()
    }


def build_arm_view_df(
    clr,
    cytoband_file: Optional[os.PathLike | str] = None,
    min_arm_bins: Optional[int] = None,
) -> pd.DataFrame:
    """Build a view DataFrame split at each chromosome's centromere.

    Every chromosome becomes two regions (``<chrom>_p`` / ``<chrom>_q``).
    A chromosome stays as one region when its centromere is unknown, the
    assembly does not match the cytoband, or either arm has fewer than
    ``min_arm_bins`` valid bins (e.g. acrocentric short arms).

    Parameters
    ----------
    clr : cooler.Cooler
        Balanced cooler (needs a ``weight`` column).
    cytoband_file : path-like, optional
        UCSC cytoBandIdeo file. See :func:`resolve_cytoband_path`.
    min_arm_bins : int, optional
        Minimum number of valid bins per arm. When omitted, it is derived
        from the resolution (10 Mb worth of bins, at least 10).

    Returns
    -------
    pd.DataFrame
        Columns: ``chrom, start, end, name``.
    """
    path = resolve_cytoband_path(cytoband_file)
    if path is None:
        warnings.warn(
            "cytoBandIdeo not found; falling back to whole-chromosome view. "
            "Pass cytoband_file=... or set CFIZZ_CYTOBAND to enable arm view.",
            UserWarning,
        )
        return get_view_df(clr)

    bands, band_sizes = _load_cytoband(path)
    if not _cytoband_matches_clr(band_sizes, clr):
        warnings.warn(
            "cytoBandIdeo chromosome sizes do not match the cooler; "
            "falling back to whole-chromosome view.",
            UserWarning,
        )
        return get_view_df(clr)

    centromeres = {
        chrom: (min(s for s, _ in v) + max(e for _, e in v)) // 2
        for chrom, v in bands.items()
    }

    bins = clr.bins()[:]
    has_weight = "weight" in bins.columns
    valid = bins["weight"].notna() if has_weight else pd.Series(True, index=bins.index)

    if min_arm_bins is None:
        if len(bins) > 1:
            resolution = int(np.median(np.diff(bins["start"].to_numpy())))
            resolution = resolution if resolution > 0 else 1
        else:
            resolution = 1
        min_arm_bins = max(10, int(np.ceil(10_000_000 / resolution)))

    rows: List[Dict[str, Any]] = []
    for chrom in clr.chromnames:
        size = int(clr.chromsizes[chrom])
        cen = centromeres.get(chrom)
        if cen is None:
            cen = centromeres.get(_chrom_key(chrom))
        if cen is None:
            cen = centromeres.get(f"chr{_chrom_key(chrom)}")

        if cen is None or not (0 < cen < size):
            rows.append({"chrom": chrom, "start": 0, "end": size, "name": chrom})
            continue

        sub = bins[bins["chrom"] == chrom]
        n_p = int(valid[sub.index][sub["start"] < cen].sum())
        n_q = int(valid[sub.index][sub["start"] >= cen].sum())
        if n_p < min_arm_bins or n_q < min_arm_bins:
            rows.append({"chrom": chrom, "start": 0, "end": size, "name": chrom})
            continue

        rows.append({"chrom": chrom, "start": 0, "end": cen, "name": f"{chrom}_p"})
        rows.append({"chrom": chrom, "start": cen, "end": size, "name": f"{chrom}_q"})

    return pd.DataFrame(rows)


# =============================================================================
# Eigenvector selection
# =============================================================================

def select_compartment_eigenvector(
    eigenvalues: pd.DataFrame,
    eigenvector_df: pd.DataFrame,
    eig_select: str = "positive",
    n_eigs: int = 3,
) -> pd.DataFrame:
    """选择哪个特征向量作为 E1 (A/B 区室模式)。

    A/B 区室模式在 (O/E − 1) 矩阵上必须对应**正特征值**: 同区室 bin 之间的
    接触高于期望 ⇒ 二次型为正 ⇒ 特征值为正。cooltools 按 |λ| 排序是符号盲的,
    在 ± 近简并对上会选到负 λ 的那个 —— 那是"反 A/B 模式", 棋盘强度为负。

    Parameters
    ----------
    eigenvalues : pd.DataFrame
        cooltools.eigs_cis 返回的 eigvals_table (每个 region 一行, 含 eigval1..N)
    eigenvector_df : pd.DataFrame
        cooltools.eigs_cis 返回的 eigvec_table (每个 bin 一行, 含 E1..EN)
    eig_select : str
        "positive" —— 在 λ > 0 的候选中取 λ 最大者提为 E1 (推荐)
        "abs"      —— 保持 cooltools 原行为 (|λ| 最大者), 供 A/B 对照
    n_eigs : int
        特征向量个数

    Returns
    -------
    pd.DataFrame
        重排后的 eigenvector_df (E1..EN 仍是同一组正交基的置换)

    Notes
    -----
    cooltools 的 eigvals_table 与 eigvec_table 严格一一对应 (eigval_k ↔ E_k):
    `_phase_eigs` 在 sort_metric=None 时只翻符号不重排, sort_metric 非 None 时
    eigvals 与 eigvecs 一起重排。本函数已用独立稠密 eigh 对照验证 (r = ±1.000000)。
    """
    if eig_select not in ("positive", "abs"):
        raise ValueError(f"未知 eig_select: {eig_select!r} (可选: 'positive' / 'abs')")
    if eig_select == "abs" or n_eigs <= 1:
        return eigenvector_df

    val_cols = [f"eigval{i + 1}" for i in range(n_eigs)]
    vec_cols = [f"E{i + 1}" for i in range(n_eigs)]

    out = eigenvector_df.copy()
    for _, row in eigenvalues.iterrows():
        lams = np.array([row[c] for c in val_cols], dtype=float)
        pos = np.where(lams > 0)[0]
        if len(pos) == 0:
            # 没有正 λ → 保持原顺序 (退化情况由 QC 表报警)
            continue
        k = int(pos[np.argmax(lams[pos])])
        if k == 0:
            continue

        # 取出该 region 的 bin (按坐标区间)
        sel = (
            (out["chrom"] == row["chrom"])
            & (out["start"] >= int(row["start"]))
            & (out["start"] < int(row["end"]))
        )
        if not sel.any():
            continue
        # E1 <- E_{k+1}, 其余依次前移 (保持正交基的置换)
        order = [k] + [i for i in range(n_eigs) if i != k]
        block = out.loc[sel, vec_cols].to_numpy(copy=True)
        out.loc[sel, vec_cols] = block[:, order]

    return out


# =============================================================================
# 质量体检 (QC)
# =============================================================================

def _count_blocks(e: np.ndarray, starts: np.ndarray,
                  min_abs: float = 0.05) -> Tuple[int, int, float]:
    """统计连续 A 块数、连续 B 块数、最长单块占该 region 长度的比例。

    |E1| <= min_abs 的 bin 视为"未定", 用最近的确定 bin 填充, 避免近零抖动
    产生虚假的小块。
    """
    if len(e) == 0:
        return 0, 0, 0.0
    s = np.sign(e)
    solid = np.abs(e) > min_abs
    if not solid.any():
        return 0, 0, 0.0
    idx = np.where(solid)[0]
    filled = np.full(len(s), np.nan)
    for j, i in enumerate(idx):
        lo = idx[j - 1] if j > 0 else 0
        filled[lo:i + 1] = s[i]
    f = filled[np.isfinite(filled)]
    if len(f) == 0:
        return 0, 0, 0.0
    changes = np.diff(f) != 0
    bounds = np.where(np.concatenate([[True], changes]))[0]
    ends = np.concatenate([bounds[1:], [len(f)]])
    n_a = n_b = 0
    block_lens: List[float] = []
    for b, en in zip(bounds, ends):
        if en <= b:
            continue
        if f[b] > 0:
            n_a += 1
        else:
            n_b += 1
        block_lens.append(float(starts[min(en, len(starts)) - 1] - starts[b]))
    span = float(starts[-1] - starts[0]) if len(starts) > 1 else 0.0
    max_frac = (max(block_lens) / span) if (span > 0 and block_lens) else 0.0
    return n_a, n_b, max_frac


def compartment_qc(
    eigenvector_df: pd.DataFrame,
    gc_cov: Optional[pd.DataFrame] = None,
    centromere_positions: Optional[Dict[str, int]] = None,
    view_df: Optional[pd.DataFrame] = None,
    n_eigs: int = 3,
) -> pd.DataFrame:
    """Run per-region quality checks on a computed eigenvector table.

    This is the "is every chromosome sane?" diagnostic. It reports:

    - **gc_corr > 0**: the A side must be the GC-rich side (phasing check).
    - **n_A_blocks / n_B_blocks / n_runs**: number of contiguous A/B runs.
      Few runs with a large ``max_block_frac`` indicate a whole-arm or
      whole-region split rather than an intra-arm checkerboard.
    - **pctA / absE1_median**: compartment balance and signal amplitude.
    - **arm_purity / arm_contrast**: only meaningful for a region spanning a
      whole chromosome; a purity near 1 plus high contrast indicates the two
      arms were called as different compartments (the 2009 Science artefact).

    Parameters
    ----------
    eigenvector_df : pd.DataFrame
        Eigenvector table (chrom, start, end, E1, ...).
    gc_cov : pd.DataFrame, optional
        GC table (chrom, start, end, gc/GC).
    centromere_positions : dict, optional
        ``{chrom: position}``. Defaults to the resolved cytoBand file.
    view_df : pd.DataFrame, optional
        Region table used for the decomposition. Required to separate
        adjacent arm regions (their coordinates are contiguous).
    n_eigs : int
        Number of eigenvectors (kept for API compatibility).

    Returns
    -------
    pd.DataFrame
        One row per (chrom, region) with QC metrics.
    """
    import scipy.stats as st

    if centromere_positions is None:
        try:
            centromere_positions = load_centromere_positions()
        except FileNotFoundError:
            centromere_positions = {}

    gc_map = None
    if gc_cov is not None:
        # bioframe.frac_gc 的返回列名在不同版本是 'gc' 或 'GC'
        gc_col = "gc" if "gc" in gc_cov.columns else (
            "GC" if "GC" in gc_cov.columns else None)
        if gc_col is not None:
            gc_map = gc_cov.set_index(["chrom", "start", "end"])[gc_col]

    rows: List[Dict[str, Any]] = []

    if view_df is None:
        raise ValueError(
            "compartment_qc 需要 view_df 才能区分 region —— "
            "拆臂后两臂坐标连续, 无法从特征向量表本身还原边界。")

    for _, rv in view_df.iterrows():
        chrom, r_start, r_end = rv["chrom"], int(rv["start"]), int(rv["end"])
        region_name = rv["name"] if "name" in view_df.columns else chrom
        sub = eigenvector_df[
            (eigenvector_df["chrom"] == chrom)
            & (eigenvector_df["start"] >= r_start)
            & (eigenvector_df["start"] < r_end)
        ].sort_values("start")
        if len(sub) == 0:
            continue

        e = sub["E1"].to_numpy(dtype=float)
        valid = np.isfinite(e)
        e = e[valid]
        starts = sub["start"].to_numpy()[valid]
        if len(e) == 0:
            continue

        n_a, n_b, max_frac = _count_blocks(e, starts)
        row: Dict[str, Any] = {
            "chrom": chrom,
            "region": region_name,
            "start": int(starts[0]),
            "end": int(starts[-1]),
            "n_bins_valid": int(len(e)),
            "n_A_blocks": n_a,
            "n_B_blocks": n_b,
            "n_runs": n_a + n_b,
            "max_block_frac": max_frac,
            "pctA": float(100 * (e > 0).mean()),
            "absE1_median": float(np.median(np.abs(e))),
        }

        # GC 相关
        if gc_map is not None:
            gcv = gc_map.reindex(
                pd.MultiIndex.from_arrays(
                    [sub["chrom"].to_numpy()[valid],
                     sub["start"].to_numpy()[valid],
                     sub["end"].to_numpy()[valid]]
                )
            ).to_numpy()
            m = np.isfinite(gcv) & np.isfinite(e)
            row["gc_corr"] = float(st.spearmanr(gcv[m], e[m])[0]) if m.sum() > 10 else np.nan
        else:
            row["gc_corr"] = np.nan

        # 臂内均一性 / 臂间对比度
        # 仅当该 region 覆盖整条染色体时才有意义 —— 拆臂后每个 region 只有一条臂,
        # 此时"臂内均一性"应当按 region 自身算 (即该 region 内的 A/B 混杂程度)。
        cen = centromere_positions.get(chrom)
        chrom_end = int(eigenvector_df.loc[eigenvector_df["chrom"] == chrom, "end"].max())
        spans_whole = (r_start == 0) and (r_end >= chrom_end)
        if cen is not None and spans_whole:
            p_mask = starts < cen
            q_mask = ~p_mask
            if p_mask.sum() > 10 and q_mask.sum() > 10:
                pctA_p = float(100 * (e[p_mask] > 0).mean())
                pctA_q = float(100 * (e[q_mask] > 0).mean())
                row["pctA_p"] = pctA_p
                row["pctA_q"] = pctA_q
                row["arm_contrast"] = abs(pctA_p - pctA_q) / 2
                row["arm_purity"] = max(
                    max(pctA_p, 100 - pctA_p), max(pctA_q, 100 - pctA_q)
                ) / 100
            else:
                row["pctA_p"] = np.nan
                row["pctA_q"] = np.nan
                row["arm_contrast"] = np.nan
                row["arm_purity"] = np.nan
        elif cen is not None:
            # 单臂 region: 用该臂自身的 A/B 混杂度代替 arm_purity
            pctA_r = float(100 * (e > 0).mean())
            row["pctA_p"] = np.nan
            row["pctA_q"] = np.nan
            row["arm_contrast"] = np.nan
            row["arm_purity"] = max(pctA_r, 100 - pctA_r) / 100
        else:
            row["pctA_p"] = np.nan
            row["pctA_q"] = np.nan
            row["arm_contrast"] = np.nan
            row["arm_purity"] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# =============================================================================

def compute_eigenvector_cooltools(
    clr,
    gc_cov: pd.DataFrame,
    view_df: pd.DataFrame,
    n_eigs: int = 3,
    eig_select: str = "positive",
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Compute eigenvector using cooltools.eigs_cis.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object
    gc_cov : pd.DataFrame
        GC content DataFrame
    view_df : pd.DataFrame
        Chromosome view DataFrame
    n_eigs : int
        Number of eigenvectors to compute
    eig_select : str
        E1 选取规则 (见 select_compartment_eigenvector):
        "positive" —— 取 λ>0 中最大者 (推荐, 修 chr13 反号)
        "abs"      —— cooltools 原行为

    Returns
    -------
    eigenvalues : np.ndarray
        Eigenvalues for each chromosome
    eigenvector_df : pd.DataFrame
        DataFrame with columns: chrom, start, end, weight, E1, E2, ...
    """
    import cooltools

    # Compute eigenvectors using cooltools
    eigenvalues, eigenvector_df = cooltools.eigs_cis(
        clr,
        gc_cov,
        view_df=view_df,
        n_eigs=n_eigs
    )

    # 选出真正的 A/B 区室模式作为 E1
    eigenvector_df = select_compartment_eigenvector(
        eigenvalues, eigenvector_df, eig_select=eig_select, n_eigs=n_eigs
    )

    return eigenvalues, eigenvector_df


# =============================================================================
# Complete Process Functions
# =============================================================================

def process_compartment(
    mcool_path: str,
    resolution: int,
    fasta_path: Optional[str] = None,
    output_dir: str = ".",
    n_eigs: int = 3,
    eig_select: str = "positive",
    view: str = "arm",
    cytoband_file: Optional[str] = None,
    min_arm_bins: Optional[int] = None,
    write_qc: bool = True,
    sample_name: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compute compartment eigenvectors for one .mcool.

    Output files
    ------------
    - ``eigenvector.{sample}.{res}.tsv``
    - ``gc_cov.{res}.tsv`` (written when ``fasta_path`` is provided)
    - ``compartment_qc.{sample}.{res}.tsv`` (when ``write_qc=True``)

    Parameters
    ----------
    mcool_path : str
        Path to .mcool file.
    resolution : int
        Resolution in bp (e.g. 100000).
    fasta_path : str, optional
        Reference genome FASTA. Without it, GC is set to 0.5.
    output_dir : str
        Output directory.
    n_eigs : int
        Number of eigenvectors to compute.
    eig_select : str
        ``"positive"`` promotes the largest positive-λ vector to E1.
        ``"abs"`` keeps cooltools' absolute-λ ordering.
    view : str
        ``"arm"`` decomposes p/q arms separately at the centromere.
        ``"chromosome"`` uses one region per chromosome.
    cytoband_file : path-like, optional
        UCSC cytoBandIdeo file. Defaults to ``CFIZZ_CYTOBAND`` or
        ``demo/data/hg38.cytoBand.ideo.txt.gz``.
    min_arm_bins : int, optional
        Minimum valid bins per arm. Defaults to 10 Mb worth of bins.
    write_qc : bool
        Write ``compartment_qc.{sample}.{res}.tsv``.
    sample_name : str, optional
        Sample label used in output filenames. Defaults to the first
        dot-separated token of the .mcool filename.

    Returns
    -------
    pd.DataFrame
        Eigenvector table with columns ``chrom, start, end, weight, E1...``.
    """
    import cooler

    os.makedirs(output_dir, exist_ok=True)
    if sample_name is None:
        sample_name = Path(mcool_path).name.split(".")[0]

    print(f"\n[process_compartment] Computing compartment...")
    print(f"  mcool: {mcool_path}")
    print(f"  sample: {sample_name}   resolution: {resolution}")
    print(f"  eig_select: {eig_select}   view: {view}")

    clr = cooler.Cooler(f'{mcool_path}::resolutions/{resolution}')
    bins = clr.bins()[:]

    res_str = _format_resolution(resolution)
    gc_cov_path = os.path.join(output_dir, f'gc_cov.{res_str}.tsv')

    if fasta_path is None:
        print("  No fasta provided, using default GC=0.5")
        gc_cov = pd.DataFrame({
            'chrom': bins['chrom'],
            'start': bins['start'],
            'end': bins['end'],
            'gc': 0.5,
        })
    elif os.path.exists(gc_cov_path):
        print(f"  Loading existing GC content: {gc_cov_path}")
        gc_cov = pd.read_csv(gc_cov_path, sep='\t')
    else:
        print("  Calculating GC content with bioframe.frac_gc...")
        gc_cov = calc_gc_cov(bins, fasta_path, gc_cov_path)
        print(f"  GC content saved to: {gc_cov_path}")

    if view == "arm":
        view_df = build_arm_view_df(
            clr, cytoband_file=cytoband_file, min_arm_bins=min_arm_bins
        )
        n_split = int((view_df["name"] != view_df["chrom"]).sum())
        print(f"  View: arm-level ({n_split} arm regions, "
              f"{len(view_df)} regions total)")
    elif view == "chromosome":
        view_df = get_view_df(clr)
        print(f"  View: whole chromosome ({len(view_df)} regions)")
    else:
        raise ValueError(f"未知 view: {view!r} (可选: 'arm' / 'chromosome')")

    print(f"  Computing eigenvectors (n_eigs={n_eigs})...")
    eigenvalues, eigenvector_df = compute_eigenvector_cooltools(
        clr, gc_cov, view_df, n_eigs=n_eigs, eig_select=eig_select
    )

    eig_path = os.path.join(
        output_dir, f'eigenvector.{sample_name}.{res_str}.tsv'
    )
    eigenvector_df.to_csv(eig_path, sep='\t', index=False)
    print(f"  Eigenvector saved to: {eig_path}")

    if write_qc:
        try:
            centromeres = load_centromere_positions(cytoband_file)
        except FileNotFoundError:
            centromeres = None
        qc = compartment_qc(
            eigenvector_df,
            gc_cov=gc_cov,
            view_df=view_df,
            centromere_positions=centromeres,
            n_eigs=n_eigs,
        )
        qc_path = os.path.join(
            output_dir, f'compartment_qc.{sample_name}.{res_str}.tsv'
        )
        qc.to_csv(qc_path, sep='\t', index=False)
        print(f"  QC table saved to: {qc_path}")

    return eigenvector_df


# =============================================================================
# Original functions (kept for reference/backward compatibility)
# =============================================================================

def compute_compartment(
    oe_matrix: np.ndarray,
    n_components: int = 1
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute A/B compartment using eigenvector decomposition (custom method).

    The first eigenvector (E1) reflects the chessboard pattern of Hi-C data,
    where positive values indicate A compartments (open chromatin) and
    negative values indicate B compartments (closed chromatin).

    Note: This is a custom implementation. For cooler data, use
    process_compartment() which uses cooltools.eigs_cis for better results.
    """
    # Replace NaN values with 1 (neutral O/E value)
    oe_matrix = np.nan_to_num(oe_matrix, nan=1.0)

    # Ensure symmetric matrix
    if not np.allclose(oe_matrix, oe_matrix.T):
        oe_matrix = (oe_matrix + oe_matrix.T) / 2

    # Calculate correlation matrix
    n = oe_matrix.shape[0]
    corr_matrix = np.ones((n, n))

    for i in range(n):
        for j in range(n):
            if i != j:
                corr_matrix[i, j] = oe_matrix[i, j] * oe_matrix[j, i]

    # Replace diagonal with average of off-diagonal
    mask = ~np.eye(n, dtype=bool)
    valid_values = oe_matrix[mask]
    diag_value = np.nanmean(valid_values) if len(valid_values) > 0 else 1.0
    np.fill_diagonal(corr_matrix, diag_value)

    # Handle any remaining NaN/Inf values
    corr_matrix = np.nan_to_num(corr_matrix, nan=1.0, posinf=1.0, neginf=-1.0)

    # Eigenvalue decomposition
    eigenvalues, eigenvectors = eigh(corr_matrix)

    # Sort by eigenvalue (descending)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return eigenvectors[:, :n_components], eigenvalues[:n_components]


def compute_eigenvector(
    matrix: np.ndarray,
    n_components: int = 1,
    method: str = "classic"
) -> pd.DataFrame:
    """
    Compute eigenvectors for compartment analysis.
    """
    if method == "classic":
        eigenvectors, eigenvalues = compute_compartment(matrix, n_components)

    n = matrix.shape[0]
    data = {
        'start': np.arange(n),
        'end': np.arange(n) + 1
    }

    for i in range(n_components):
        data[f'E{i+1}'] = eigenvectors[:, i]


def assign_compartments(
    eigenvector: np.ndarray,
    gc_content: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, str]:
    """
    Assign A/B compartment labels based on eigenvector.
    """
    if gc_content is not None:
        correlation = np.corrcoef(eigenvector, gc_content)[0, 1]
        if correlation < 0:
            eigenvector = -eigenvector
            orientation = "inverted"
        else:
            orientation = "direct"
    else:
        orientation = "unknown"

    compartment_labels = np.where(eigenvector > 0, 1, -1)
    return compartment_labels, orientation


def calculate_compartment_strength(
    matrix: np.ndarray,
    compartment_labels: np.ndarray,
    n_bins: int = 10
) -> Dict[str, float]:
    """
    Calculate compartment strength metrics.
    """
    n = matrix.shape[0]

    within_a = within_b = between = 0
    count_a = count_b = count_between = 0

    for i in range(n):
        for j in range(i + 1, n):
            if compartment_labels[i] == 1 and compartment_labels[j] == 1:
                within_a += matrix[i, j]
                count_a += 1
            elif compartment_labels[i] == -1 and compartment_labels[j] == -1:
                within_b += matrix[i, j]
                count_b += 1
            else:
                between += matrix[i, j]
                count_between += 1

    avg_within_a = within_a / count_a if count_a > 0 else 0
    avg_within_b = within_b / count_b if count_b > 0 else 0
    avg_between = between / count_between if count_between > 0 else 0

    return {
        'avg_within_A': avg_within_a,
        'avg_within_B': avg_within_b,
        'avg_between': avg_between,
        'strength_A': avg_within_a / avg_between if avg_between > 0 else 0,
        'strength_B': avg_within_b / avg_between if avg_between > 0 else 0
    }


# =============================================================================
# =============================================================================

def calculate_saddle_matrix(
    matrix: np.ndarray,
    track: pd.DataFrame,
    n_bins: int = 98,
    qrange: Tuple[float, float] = (0.05, 0.95)
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate saddle score matrix using digitization and aggregation.

    Parameters
    ----------
    matrix : np.ndarray
        O/E normalized contact matrix
    track : pd.DataFrame
        DataFrame with compartment scores. Must have 'chrom', 'start', 'end', and E1 column.
    n_bins : int
        Number of bins for saddle
    qrange : tuple
        Quantile range for digitization

    Returns
    -------
    saddledata : np.ndarray
        Saddle matrix (n_bins x n_bins)
    track_values : np.ndarray
        Original track values (filtered)
    binedges : np.ndarray
        Bin edges for the track values
    """
    track_values = track['E1'].values
    n = min(len(track_values), matrix.shape[0])

    track_values = track_values[:n]
    matrix_cropped = matrix[:n, :n]

    binedges = np.percentile(track_values, np.linspace(0, 100, n_bins + 1))
    digitized = np.digitize(track_values, binedges[1:-1])

    saddledata = np.zeros((n_bins, n_bins))
    counts = np.zeros((n_bins, n_bins))

    for i in range(n):
        for j in range(n):
            di = digitized[i]
            dj = digitized[j]
            if 0 <= di < n_bins and 0 <= dj < n_bins:
                saddledata[di, dj] += matrix_cropped[i, j]
                counts[di, dj] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        saddledata = np.where(counts > 0, saddledata / counts, np.nan)

    return saddledata, track_values, binedges


def saddle_score(
    matrix: np.ndarray,
    track: pd.DataFrame,
    n_bins: int = 98,
    qrange: Tuple[float, float] = (0.075, 0.975)
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate saddle score matrix (legacy function).
    """
    return calculate_saddle_matrix(matrix, track, n_bins, qrange)


# =============================================================================
# Complete Analysis Pipeline
# =============================================================================

def calculate_compartment_analysis(
    mcool_path: str,
    chrom: str,
    resolution: int,
    n_bins_saddle: int = 10,
    fasta_path: Optional[str] = None,
    output_dir: str = ".",
    view: str = "chromosome",
    eig_select: str = "positive",
    sample_name: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Perform a complete single-chromosome compartment analysis.

    Workflow: compute E1 with :func:`process_compartment`, build the saddle
    matrix from the chromosome-wide O/E matrix, and summarize A/B statistics.

    ``view`` defaults to ``"chromosome"`` here because the saddle matrix is
    computed on one chromosome-wide O/E matrix. Set ``view="arm"`` when
    arm-level E1 values are desired.

    Returns
    -------
    dict
        Keys: ``eigenvector_df``, ``saddle_matrix``, ``saddle_stats``,
        ``compartment_stats``.
    """
    from cfizz.analyze.oe import calculate_oe_matrix_cooltools

    # Step 1: Calculate O/E matrix
    print(f"[calculate_compartment_analysis] Step 1: O/E calculation...")
    OE, _ = calculate_oe_matrix_cooltools(
        mcool_path=mcool_path,
        chrom=chrom,
        resolution=resolution,
        balance=True,
        smooth=True
    )

    # Step 2: Compute eigenvectors
    print(f"[calculate_compartment_analysis] Step 2: Eigenvector computation...")
    eig_df = process_compartment(
        mcool_path=mcool_path,
        resolution=resolution,
        fasta_path=fasta_path,
        output_dir=output_dir,
        n_eigs=2,
        eig_select=eig_select,
        view=view,
        sample_name=sample_name,
    )

    # Filter for our chromosome
    eig_df = eig_df[eig_df['chrom'] == chrom].copy()

    # Step 3: Calculate saddle matrix
    print(f"[calculate_compartment_analysis] Step 3: Saddle matrix calculation...")
    saddle_matrix, _, _ = calculate_saddle_matrix(
        matrix=OE,
        track=eig_df,
        n_bins=n_bins_saddle
    )

    # Step 4: Calculate statistics
    E1 = eig_df['E1'].values[:OE.shape[0]]
    n_a = np.sum(E1 > 0)
    n_b = np.sum(E1 <= 0)

    E1_pos = E1 > 0
    E1_neg = E1 <= 0

    strength_a = np.nanmean(OE[np.ix_(E1_pos, E1_pos)]) if np.sum(E1_pos) > 0 else np.nan
    strength_b = np.nanmean(OE[np.ix_(E1_neg, E1_neg)]) if np.sum(E1_neg) > 0 else np.nan

    results = {
        'eigenvector_df': eig_df,
        'saddle_matrix': saddle_matrix,
        'saddle_stats': {
            'n_bins': n_bins_saddle,
            'mean_aa': np.nanmean(saddle_matrix[:n_bins_saddle//2, :n_bins_saddle//2]),
            'mean_bb': np.nanmean(saddle_matrix[n_bins_saddle//2:, n_bins_saddle//2:])
        },
        'compartment_stats': {
            'n_A_bins': n_a,
            'n_B_bins': n_b,
            'pct_A': 100 * n_a / len(E1) if len(E1) > 0 else 0,
            'pct_B': 100 * n_b / len(E1) if len(E1) > 0 else 0,
            'strength_A': strength_a,
            'strength_B': strength_b
        }
    }

    return results


# =============================================================================
# Differential Analysis Functions (移植自 script/a_4_compartment_diff.py)
# Source: a_4 L48, L143, L188, L240, L258, L271, L368
# =============================================================================

import logging

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math

# 配置字体和绘图参数 (a_4 L61-65)
FONT_CONFIG = {
    'family': 'Arial',
    'weight': 'regular',
    'size': 5,
}

CM = 1 / 2.54  # cm转inch (a_4 L48)

# 配色方案 (a_4 L69)
COLOR_PALETTE = ['#264653', '#27736F', '#299D92', '#8AB17C', '#E8C56B', '#F2A361', '#E66F51']

# Compartment状态颜色映射 (a_4 L72-77)
COMPARTMENT_COLORS = {
    'Stable_A': COLOR_PALETTE[0],   # 深青色 - A compartment稳定
    'Stable_B': COLOR_PALETTE[1],   # 青色 - B compartment稳定
    'A_to_B': COLOR_PALETTE[5],     # 橙色 - A转B
    'B_to_A': COLOR_PALETTE[6],     # 红橙色 - B转A
}

# 兼容性别名 (旧代码用 COLORS_DICT)
COLORS_DICT = COMPARTMENT_COLORS

def max_to_top(max_value: float) -> float:
    """计算绘图时的最大值上限 (a_4 L113-121)"""
    all_weishu = len(str(int(max_value)))
    top_value = math.ceil(max_value / (10 ** (all_weishu - 1))) * (10 ** (all_weishu - 1))
    if top_value < max_value * 1.2:
        return top_value
    else:
        top_value = math.ceil(max_value / (10 ** (all_weishu - 2))) * (10 ** (all_weishu - 2))
        return top_value


def parse_comparison_func(comparison_str: str) -> Tuple[str, str]:
    """解析比较组字符串，返回 (treatment, control)"""
    if '--' in comparison_str:
        parts = comparison_str.split('--')
        return parts[0], parts[1]
    return comparison_str, comparison_str


def load_compartment_data(file_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    加载compartment数据文件。
    Source: a_4 L143

    Parameters
    ----------
    file_path : str
        E1 TSV文件路径，列: chrom, start, end, E1
    logger : logging.Logger
        日志记录器

    Returns
    -------
    pd.DataFrame
        清理后的DataFrame，只包含 chrom, start, end, E1 列

    行为特征 6 问:
    1. 输入: TSV文件，列名可能有大小写差异
    2. 输出: DataFrame [chrom, start, end, E1]
    3. 跟5_1产物关系: 读eigenvector.{sample}.100k.tsv，不重算
    4. 跟a_4算法关系: 完全一致，L143-185
    5. 调谁: pd.read_csv, dropna
    6. 错误处理: 列名映射，dropna处理
    """
    logger.info(f"正在读取文件: {file_path}")

    try:
        df = pd.read_csv(file_path, sep='\t')

        expected_cols = {'chrom', 'start', 'end', 'E1'}
        actual_cols = set(df.columns)

        if not expected_cols.issubset(actual_cols):
            logger.warning(f"列名不匹配。期望: {expected_cols}, 实际: {actual_cols}")
            col_mapping = {}
            for col in actual_cols:
                col_lower = col.lower()
                if 'chr' in col_lower or 'chrom' in col_lower:
                    col_mapping['chrom'] = col
                elif 'start' in col_lower:
                    col_mapping['start'] = col
                elif 'end' in col_lower:
                    col_mapping['end'] = col
                elif 'e1' in col_lower:
                    col_mapping['E1'] = col

            if col_mapping:
                df = df.rename(columns=col_mapping)

        df_clean = df.dropna(subset=['E1']).copy()

        for col in ['start', 'end']:
            if col in df_clean.columns:
                df_clean.loc[:, col] = pd.to_numeric(df_clean[col], errors='coerce')

        df_clean = df_clean.dropna(subset=['start', 'end']).copy()

        logger.info(f"  -> 有效数据行数: {len(df_clean)}")
        logger.info(f"  -> 染色体数量: {df_clean['chrom'].nunique()}")
        logger.info(f"  -> E1值范围: [{df_clean['E1'].min():.3f}, {df_clean['E1'].max():.3f}]")

        return df_clean[['chrom', 'start', 'end', 'E1']].copy()

    except Exception as e:
        logger.error(f"读取文件失败 {file_path}: {str(e)}")
        raise


def calculate_compartment_diff(
    df_control: pd.DataFrame,
    df_treatment: pd.DataFrame,
    treatment_name: str,
    logger: logging.Logger,
    min_abs_e1: float = 0.0,
) -> pd.DataFrame:
    """
    计算两个样本间的E1差异，逐bin AB状态对比。
    Source: a_4 L188

    核心函数，内含 classify_compartment_change 内部函数。
    分类4类: Stable_A / Stable_B / A_to_B / B_to_A

    min_abs_e1 幅度门槛。
        |E1| 很小的 bin (两个样本中任一方) 的符号翻转通常是噪声而非真实区室变化。
        min_abs_e1 > 0 时, 这类 bin 归入 'Uncertain' 类, 不计入 A/B 翻转。
        min_abs_e1 = 0.0 (默认) 时关闭该门槛。

    Parameters
    ----------
    df_control : pd.DataFrame
        对照组E1数据
    df_treatment : pd.DataFrame
        处理组E1数据
    treatment_name : str
        处理组名称
    logger : logging.Logger
        日志记录器
    min_abs_e1 : float
        幅度门槛; 0.0 = 关闭 (原版行为)

    Returns
    -------
    pd.DataFrame
        差异数据，包含 E1_control, E1_{treatment_name}, compartment_change_{treatment_name}

    行为特征 6 问:
    1. 输入: 两个E1 DataFrame，通过chrom/start/end合并
    2. 输出: 合并后的DataFrame，新增dE1和compartment_change列
    3. 跟5_1产物关系: 读两个eigenvector.tsv，不重算
    4. 跟a_4算法关系: 完全一致，L188-237 (min_abs_e1=0 时)
    5. 调谁: pd.merge, numpy操作
    6. 错误处理: 空DataFrame时返回空DataFrame
    """
    logger.info(f"\n计算 {treatment_name} 与 control 的差异...")

    merged = pd.merge(df_control, df_treatment, on=['chrom', 'start', 'end'],
                      suffixes=('_control', f'_{treatment_name}'),
                      how='inner')

    if len(merged) == 0:
        logger.warning(f"  -> 警告: 没有找到匹配的坐标区域")
        return pd.DataFrame()

    logger.info(f"  -> 匹配区域数: {len(merged)}")

    e1_control = merged['E1_control']
    e1_sample = merged[f'E1_{treatment_name}']

    dE1_col = f'dE1_{treatment_name}_vs_control'
    merged[dE1_col] = e1_sample - e1_control

    abs_dE1_col = f'abs_dE1_{treatment_name}_vs_control'
    rel_dE1_col = f'rel_dE1_{treatment_name}_vs_control'

    merged[abs_dE1_col] = np.abs(merged[dE1_col])

    merged[rel_dE1_col] = np.where(
        np.abs(e1_control) > 1e-10,
        np.abs(merged[dE1_col]) / np.abs(e1_control),
        0
    )

    def classify_compartment_change(e1_ctrl, e1_sample):
        """内部函数: 分类E1变化"""
        if e1_ctrl > 0 and e1_sample > 0:
            return 'Stable_A'
        elif e1_ctrl < 0 and e1_sample < 0:
            return 'Stable_B'
        elif e1_ctrl > 0 and e1_sample < 0:
            return 'A_to_B'
        elif e1_ctrl < 0 and e1_sample > 0:
            return 'B_to_A'
        else:
            return 'Boundary'

    merged[f'compartment_change_{treatment_name}'] = merged.apply(
        lambda row: classify_compartment_change(row['E1_control'], row[f'E1_{treatment_name}']),
        axis=1
    )

    # 幅度门槛 —— 近零 E1 的符号抖动归入 Uncertain
    if min_abs_e1 > 0:
        col = f'compartment_change_{treatment_name}'
        uncertain = (merged['E1_control'].abs() < min_abs_e1) | \
                    (merged[f'E1_{treatment_name}'].abs() < min_abs_e1)
        n_unc = int(uncertain.sum())
        merged.loc[uncertain, col] = 'Uncertain'
        logger.info(f"  -> |E1| < {min_abs_e1} 归入 Uncertain: {n_unc} bin "
                    f"({100 * n_unc / len(merged):.1f}%)")

    return merged


def save_plot_data(
    df_diff: pd.DataFrame,
    treatment_name: str,
    output_path: Path,
    logger: logging.Logger
) -> Path:
    """
    保存绘图专用数据文件。
    Source: a_4 L240

    Parameters
    ----------
    df_diff : pd.DataFrame
        差异数据
    treatment_name : str
        处理组名称
    output_path : Path
        输出目录
    logger : logging.Logger
        日志记录器

    Returns
    -------
    Path
        保存的绘图数据文件路径
    """
    plot_columns = [
        'chrom', 'start', 'end',
        'E1_control',
        f'E1_{treatment_name}',
        f'compartment_change_{treatment_name}'
    ]

    df_plot = df_diff[plot_columns].copy()

    plot_data_file = output_path / f"compartment_plot_data.tsv"
    df_plot.to_csv(plot_data_file, sep='\t', index=False)
    logger.info(f"  -> 绘图数据: {plot_data_file}")

    return plot_data_file


def load_plot_data(plot_data_file: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    加载绘图专用数据文件。
    Source: a_4 L258

    Parameters
    ----------
    plot_data_file : Path
        绘图数据文件路径
    logger : logging.Logger
        日志记录器

    Returns
    -------
    pd.DataFrame
        绘图数据
    """
    logger.info(f"读取绘图数据: {plot_data_file}")

    if not plot_data_file.exists():
        raise FileNotFoundError(f"绘图数据文件不存在: {plot_data_file}\n请先运行 'compute' 或 'all' 模式生成数据")

    df = pd.read_csv(plot_data_file, sep='\t')
    logger.info(f"  -> 数据行数: {len(df)}")

    return df


def plot_compartment_scatter(df: pd.DataFrame, treatment_name: str, comparison: str,
                             output_dir: Path, logger: logging.Logger) -> None:
    """绘制compartment差异散点图（带边缘分布） (a_4 L271-359)"""
    logger.info(f"\n生成 {comparison} 散点图...")

    figure = plt.figure(figsize=(8 * CM, 6.4 * CM), dpi=600)
    figure.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.95, wspace=0.1, hspace=0.1)

    grid = plt.GridSpec(6, 6, wspace=0.05, hspace=0.05, figure=figure)
    ax1 = figure.add_subplot(grid[0, 0:5])
    ax2 = figure.add_subplot(grid[1:6, 5])
    ax3 = figure.add_subplot(grid[1:6, 0:5])

    e1_control = df['E1_control']
    e1_sample = df[f'E1_{treatment_name}']
    max_val = max(e1_control.max(), e1_sample.max())
    top_value = max_to_top(max_val)
    bottom_value = top_value * 0.1 * -1
    bottom_value = bottom_value if bottom_value < -1.5 else -1.5

    sns.kdeplot(data=df, x='E1_control', color=COLOR_PALETTE[2], fill=True,
                common_norm=False, legend=False, alpha=0.5, linewidth=0.5, cut=0, ax=ax1)
    sns.kdeplot(data=df, y=f'E1_{treatment_name}', color=COLOR_PALETTE[5], fill=True,
                common_norm=False, legend=False, alpha=0.5, linewidth=0.5, cut=0, ax=ax2)

    x_line = np.linspace(bottom_value, top_value, 1000)
    plt.plot(x_line, x_line, color="#bbbbbb", linewidth=0.1)

    change_col = f'compartment_change_{treatment_name}'
    for status, color in COMPARTMENT_COLORS.items():
        subset = df[df[change_col] == status]
        if len(subset) > 0:
            ax3.scatter(x=subset['E1_control'], y=subset[f'E1_{treatment_name}'],
                       s=1, alpha=0.5, color=color, edgecolors=color, linewidths=0.5)

    total = len(df)

    legend_elements = []
    for status, color in COMPARTMENT_COLORS.items():
        count = len(df[df[change_col] == status])
        if count > 0:
            pct = count / total * 100
            legend_elements.append(mlines.Line2D([0], [0], marker='o', color='w',
                                              markerfacecolor=color, markersize=4,
                                              label=f'{status}\nn={count:,}\n({pct:.1f}%)\n'))

    if legend_elements:
        ax2.legend(handles=legend_elements, prop=FONT_CONFIG, labelspacing=0.4,
                  handleheight=1.5, handletextpad=0.2, loc=(1.01,0.3),
                  frameon=False, title='')

    ax1.spines[:].set_linewidth(0)
    ax1.tick_params(width=0.6, length=2.5, labelsize=6)
    ax1.set_xticks([])
    ax1.set_xlabel("")
    ax1.set_yticks([])
    ax1.set_ylabel("")

    ax2.spines[:].set_linewidth(0)
    ax2.tick_params(width=0.6, length=2.5, labelsize=6)
    ax2.set_xticks([])
    ax2.set_xlabel("")
    ax2.set_yticks([])
    ax2.set_ylabel("")

    ax3.spines[:].set_linewidth(0.4)
    ax3.tick_params(width=0.6, length=2.5, labelsize=7)

    _, control_name = parse_comparison_func(comparison)
    ax3.set_xlabel(f'Compartment of {control_name}', fontsize=7, x=0.55)
    ax3.set_ylabel(f'Compartment of {treatment_name}', fontsize=7, y=0.55)

    ax1.set_xlim(bottom_value, top_value)
    ax2.set_ylim(bottom_value, top_value)
    ax3.set_xlim(bottom_value, top_value)
    ax3.set_ylim(bottom_value, top_value)

    ax3.axhline(y=0, color='gray', linestyle='--', linewidth=0.3, alpha=0.5)
    ax3.axvline(x=0, color='gray', linestyle='--', linewidth=0.3, alpha=0.5)

    output_svg = output_dir / f"compartment_{comparison}_scatter.svg"
    output_png = output_dir / f"compartment_{comparison}_scatter.png"

    figure.savefig(output_svg, format='svg')
    figure.savefig(output_png, format='png', dpi=600)
    plt.close()

    logger.info(f"  -> 散点图已保存: {output_svg}")
    logger.info(f"  -> 散点图已保存: {output_png}")


def merge_diff_regions(
    diff_df: pd.DataFrame,
    treatment_name: str,
    change_type: str = 'A_to_B',
    min_region_length: int = 3,
    resolution: int = 100000
) -> pd.DataFrame:
    """
    合并连续差异区域。
    Source: a_8 L182-225 (注释说明)

    用于 A_to_B / B_to_A 找连续bin合并成region块。

    Parameters
    ----------
    diff_df : pd.DataFrame
        差异数据DataFrame
    treatment_name : str
        处理组名称
    change_type : str
        变化类型 ('A_to_B' 或 'B_to_A')
    min_region_length : int
        最小区域长度(bin数)
    resolution : int
        分辨率(bp)，默认100kb

    Returns
    -------
    pd.DataFrame
        合并后的差异区域
    """
    change_col = f'compartment_change_{treatment_name}'

    filtered = diff_df[diff_df[change_col] == change_type].copy()

    if len(filtered) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'n_bins', 'change_type'])

    filtered = filtered.sort_values(['chrom', 'start']).reset_index(drop=True)

    regions = []
    current_chrom = None
    current_start = None
    current_end = None
    current_count = 0

    for _, row in filtered.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']

        if current_chrom is None:
            current_chrom = chrom
            current_start = start
            current_end = end
            current_count = 1
        elif chrom != current_chrom or start > current_end + resolution:
            if current_count >= min_region_length:
                regions.append({
                    'chrom': current_chrom,
                    'start': current_start,
                    'end': current_end,
                    'n_bins': current_count,
                    'change_type': change_type
                })
            current_chrom = chrom
            current_start = start
            current_end = end
            current_count = 1
        else:
            current_end = end
            current_count += 1

    if current_count >= min_region_length:
        regions.append({
            'chrom': current_chrom,
            'start': current_start,
            'end': current_end,
            'n_bins': current_count,
            'change_type': change_type
        })

    return pd.DataFrame(regions)


def analyze_single_comparison(
    comparison: str,
    output_root: str,
    run_mode: str,
    control_e1_path: str = None,
    treatment_e1_path: str = None,
    min_abs_e1: float = 0.0,
):
    """
    分析单个比较组的Compartment差异。
    Source: a_4 L368

    主流程函数，整合上面所有函数，支持 run_mode='all'/'compute'/'plot'。

    Parameters
    ----------
    comparison : str
        比较组字符串，格式 "处理组--对照组"
    output_root : str
        输出根目录
    run_mode : str
        运行模式: 'all', 'compute', 'plot'
    control_e1_path : str, optional
        对照组E1文件路径，如果不提供则使用默认推导
    treatment_e1_path : str, optional
        处理组E1文件路径，如果不提供则使用默认推导
    min_abs_e1 : float, optional
        |E1| below this threshold (in either sample) is classified as
        ``Uncertain`` instead of an A/B switch. Default 0 disables the filter.

    行为特征 6 问:
    1. 输入: comparison字符串，E1文件路径
    2. 输出: diff.tsv, plot_data.tsv, scatter图
    3. 跟5_1产物关系: 读eigenvector.tsv，不重算
    4. 跟a_4算法关系: 完全一致，L368-420
    5. 调谁: load_compartment_data, calculate_compartment_diff等
    6. 错误处理: 空数据时打印错误并返回
    """
    def setup_logger(output_path, comparison):
        """设置日志记录器"""
        log_dir = output_path / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        log_file = log_dir / f"{comparison}.log"
        logger = logging.getLogger(f"compartment_diff_{comparison}")
        logger.setLevel(logging.INFO)

        if logger.handlers:
            logger.handlers.clear()

        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)

        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        logger.addHandler(fh)
        logger.addHandler(ch)

        return logger

    treatment, control = parse_comparison_func(comparison)

    output_path = Path(output_root) / "compartment"
    output_path.mkdir(parents=True, exist_ok=True)

    logger = setup_logger(output_path, comparison)

    logger.info("=" * 70)
    logger.info(f"Compartment差异分析: {comparison}")
    logger.info(f"处理组: {treatment}")
    logger.info(f"对照组: {control}")
    logger.info(f"运行模式: {run_mode}")
    logger.info("=" * 70)

    try:
        if run_mode in ['all', 'compute']:
            logger.info("\n1. 加载compartment数据")
            logger.info("-" * 70)

            if control_e1_path is None or treatment_e1_path is None:
                logger.warning("请提供 control_e1_path 和 treatment_e1_path 参数")
                return

            logger.info(f"  对照组({control})文件: {control_e1_path}")
            logger.info(f"  处理组({treatment})文件: {treatment_e1_path}")

            df_control = load_compartment_data(control_e1_path, logger)
            df_treatment = load_compartment_data(treatment_e1_path, logger)

            logger.info("\n2. 计算compartment差异")
            logger.info("-" * 70)

            df_diff = calculate_compartment_diff(
                df_control, df_treatment, treatment, logger,
                min_abs_e1=min_abs_e1,
            )

            if len(df_diff) == 0:
                logger.error("没有有效的差异数据，跳后续分析")
                return

            logger.info("\n3. 保存差异数据")
            logger.info("-" * 70)

            diff_output = output_path / f"compartment_{comparison}_diff.tsv"
            df_diff.to_csv(diff_output, sep='\t', index=False)
            logger.info(f"  -> 完整差异数据: {diff_output}")

            save_plot_data(df_diff, treatment, output_path, logger)

        if run_mode in ['all', 'plot']:
            logger.info("\n4. 生成统计图")
            logger.info("-" * 70)

            plot_data_file = output_path / "compartment_plot_data.tsv"

            if not plot_data_file.exists():
                logger.error(f"绘图数据文件不存在: {plot_data_file}")
                return

            df_plot = load_plot_data(plot_data_file, logger)
            plot_compartment_scatter(df_plot, treatment, comparison, output_path, logger)

        logger.info("\n" + "=" * 70)
        logger.info(f"分析完成！")
        logger.info(f"输出目录: {output_path}")
        logger.info("=" * 70)

    except Exception as e:
        logger.error(f"分析失败: {str(e)}")
        raise
