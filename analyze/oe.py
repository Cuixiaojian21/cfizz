"""
Observed/Expected (O/E) normalization module.

Functions for calculating O/E normalized matrices.
"""

import numpy as np
import os
from typing import Tuple, Dict, Optional, List, Union
from pathlib import Path
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor

from ..io.cooler import CoolerReader, extract_contact_matrix


def compute_expected(M: np.ndarray) -> np.ndarray:
    """
    Compute expected matrix based on distance-dependent decay.
    
    The expected matrix represents the expected contact frequency
    as a function of genomic distance, calculated as the average
    contact frequency along each diagonal (constant distance).
    
    Parameters
    ----------
    M : np.ndarray
        Contact matrix (symmetric, square)
    
    Returns
    -------
    E : np.ndarray
        Expected matrix of the same shape as M
        
    Examples
    --------
    >>> expected = compute_expected(observed_matrix)
    """
    n, m = M.shape
    E = np.zeros_like(M, dtype=float)
    
    # Calculate the average for each diagonal (constant distance)
    diagonal_averages = {}
    for offset in range(-n + 1, n):
        diagonal = np.diag(M, offset)
        if diagonal.size > 0:
            avg = diagonal.mean()
            diagonal_averages[offset] = avg
    
    # Fill the expected matrix
    for offset, avg in diagonal_averages.items():
        if offset >= 0:
            for i in range(n - offset):
                E[i, i + offset] = avg
        else:
            for i in range(m + offset):
                E[i - offset, i] = avg
    
    return E


def calculate_oe_matrix(
    O: np.ndarray,
    E: Optional[np.ndarray] = None,
    mask_diagonal: bool = True,
    n_diag: int = 1
) -> Tuple[np.ndarray, Dict]:
    """
    Calculate O/E (Observed/Expected) normalized matrix.

    The O/E normalization removes the distance-dependent decay
    from the contact matrix, making it easier to identify
    structural features like compartments, TADs, and loops.

    Parameters
    ----------
    O : np.ndarray
        Observed contact matrix
    E : np.ndarray, optional
        Expected matrix. If None, computes from O.
    mask_diagonal : bool
        Whether to mask the main diagonal (and the `n_diag` diagonals
        on each side) as NaN. Default True — useful for visualization
        and downstream analysis to avoid spurious values along the
        main diagonal.
    n_diag : int
        Number of diagonals on each side of the main diagonal to mask
        when `mask_diagonal` is True. Default 1 (just the main diagonal).

    Returns
    -------
    OE : np.ndarray
        O/E normalized matrix
    stats : dict
        Statistics about the normalization

    Examples
    --------
    >>> OE, stats = calculate_oe_matrix(observed)
    >>> OE, stats = calculate_oe_matrix(observed, expected)
    >>> OE, stats = calculate_oe_matrix(observed, mask_diagonal=True, n_diag=2)
    """
    n, m = O.shape

    # Compute expected if not provided
    if E is None:
        E = compute_expected(O)

    OE = np.full_like(O, np.nan, dtype=float)

    # Calculate O/E ratio for each element
    for offset in range(-n + 1, n):
        mask = np.diag(E, offset) != 0
        diagonal_O = np.diag(O, offset)
        diagonal_E = np.diag(E, offset)

        if diagonal_E.size > 0:
            for idx, valid in enumerate(mask):
                if valid:
                    i = idx if offset >= 0 else idx - offset
                    j = idx + offset if offset >= 0 else idx
                    if E[i, j] != 0:
                        if O[i, j] == 0:
                            OE[i, j] = np.nan
                        else:
                            OE[i, j] = O[i, j] / E[i, j]

    # Fill diagonal
    np.fill_diagonal(OE, OE.diagonal())

    # Fill lower triangle (mirror upper triangle)
    lower_tri_mask = np.tril(np.ones_like(OE, dtype=bool), -1)
    OE[lower_tri_mask] = OE.T[lower_tri_mask]

    # Mask the main diagonal and n_diag diagonals on each side as NaN
    # (the diagonal itself is already in `lower_tri_mask` for off-diag,
    # but we apply a wider band to suppress noisy values near the diagonal)
    if mask_diagonal and n_diag > 0:
        for k in range(-n_diag, n_diag + 1):
            np.fill_diagonal(OE[:, k:] if k >= 0 else OE[:k, -k:], np.nan)
        # The above np.fill_diagonal with slicing is fragile; use a clean approach:
        # build a band mask via |i-j| <= n_diag
        i_idx, j_idx = np.ogrid[:n, :m]
        band_mask = np.abs(i_idx - j_idx) <= n_diag
        OE[band_mask] = np.nan

    # Calculate statistics
    valid_oe = OE[~np.isnan(OE)]
    stats = {
        "n_elements": len(valid_oe),
        "min": np.min(valid_oe) if len(valid_oe) > 0 else np.nan,
        "max": np.max(valid_oe) if len(valid_oe) > 0 else np.nan,
        "mean": np.mean(valid_oe) if len(valid_oe) > 0 else np.nan,
        "median": np.median(valid_oe) if len(valid_oe) > 0 else np.nan,
        "std": np.std(valid_oe) if len(valid_oe) > 0 else np.nan,
        "p25": np.percentile(valid_oe, 25) if len(valid_oe) > 0 else np.nan,
        "p75": np.percentile(valid_oe, 75) if len(valid_oe) > 0 else np.nan,
        "p90": np.percentile(valid_oe, 90) if len(valid_oe) > 0 else np.nan,
        "p95": np.percentile(valid_oe, 95) if len(valid_oe) > 0 else np.nan,
        "n_nan": np.sum(np.isnan(OE)),
        "method": "naive"  # differentiate from cooltools method
    }

    return OE, stats


def extract_contact_matrix_from_file(
    mcool_path: str,
    chrom: str,
    resolution: int = 10000,
    balance: bool = True
) -> Tuple[np.ndarray, int]:
    """
    Extract contact matrix from a cooler file.
    
    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    chrom : str
        Chromosome name
    resolution : int
        Resolution in bp
    balance : bool
        Whether to return balanced matrix
        
    Returns
    -------
    matrix : np.ndarray
        Contact matrix
    chrom_length : int
        Chromosome length
    """
    reader = CoolerReader(mcool_path, resolution=resolution)
    return extract_contact_matrix(reader, chrom, balance=balance)


def process_oe_normalization(
    mcool_path: str,
    chrom: str,
    resolution: int = 10000,
    output_dir: Optional[str] = None,
    sample_name: Optional[str] = None,
    force_recompute: bool = False,
    balance: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
    """
    Process full O/E normalization workflow.
    
    Loads contact matrix, computes expected matrix,
    calculates O/E normalized matrix, and optionally saves results.
    
    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    chrom : str
        Chromosome name
    resolution : int
        Resolution in bp
    output_dir : str, optional
        Directory to save results
    sample_name : str, optional
        Sample name for file naming
    force_recompute : bool
        Whether to recompute even if cached files exist
    balance : bool
        Whether to use balanced matrix
        
    Returns
    -------
    observed : np.ndarray
        Observed contact matrix
    expected : np.ndarray
        Expected matrix
    oe : np.ndarray
        O/E normalized matrix
    stats : dict
        Normalization statistics
        
    Examples
    --------
    >>> M, E, OE, stats = process_oe_normalization(
    ...     "sample.mcool", "chr1", resolution=10000
    ... )
    """
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        basename = os.path.basename(mcool_path).split('.')[0]
        res_str = f"{resolution // 1000}k" if resolution < 1e6 else f"{resolution // int(1e6)}M"
        sample_prefix = f"{sample_name}." if sample_name else f"{basename}."
        # 文件命名规则: {序}.{sample}.{resolution}.{chrom}.{type}.txt
        # 数字前缀保留用户的层递检查习惯,统一分隔符便于排序
        prefix = f"{sample_prefix}{res_str}.{chrom}"

        observed_file = os.path.join(output_dir, f"1_1.{prefix}.observed.txt")
        expected_file = os.path.join(output_dir, f"1_2.{prefix}.expected.txt")
        oe_file = os.path.join(output_dir, f"1_3.{prefix}.oe.txt")

        if (not force_recompute and
            os.path.exists(observed_file) and
            os.path.exists(expected_file) and
            os.path.exists(oe_file)):
            # Load cached results
            M = np.loadtxt(observed_file)
            E = np.loadtxt(expected_file)
            OE = np.loadtxt(oe_file)
            # 不读 stats.txt（已删除），重新计算 stats 保证 key 一致
            _, stats = calculate_oe_matrix(M, E)
            return M, E, OE, stats

    # Load contact matrix
    reader = CoolerReader(mcool_path, resolution=resolution)
    M, chrom_length = extract_contact_matrix(reader, chrom, balance=balance)

    # Compute expected matrix
    E = compute_expected(M)

    # Calculate O/E normalized matrix
    OE, stats = calculate_oe_matrix(M, E)

    # Save results if output directory specified
    if output_dir:
        np.savetxt(observed_file, M, fmt='%.6f', delimiter='\t')
        np.savetxt(expected_file, E, fmt='%.6f', delimiter='\t')
        np.savetxt(oe_file, OE, fmt='%.6f', delimiter='\t')

    return M, E, OE, stats


def batch_compute_oe(
    mcool_paths: List[str],
    chroms: List[str],
    resolution: int,
    output_dir: str,
    n_processes: int = 8
) -> Dict[Tuple[str, str], str]:
    """
    Batch compute O/E matrices for multiple samples and chromosomes.
    
    Parameters
    ----------
    mcool_paths : list
        List of .mcool file paths
    chroms : list
        List of chromosome names
    resolution : int
        Resolution in bp
    output_dir : str
        Cache output directory
    n_processes : int
        Number of parallel processes
        
    Returns
    -------
    cache_paths : dict
        Dictionary mapping (sample, chrom) to cache file path
        
    Examples
    --------
    >>> cache = batch_compute_oe(
    ...     ["sample1.mcool", "sample2.mcool"],
    ...     ["chr1", "chr2"],
    ...     resolution=10000,
    ...     output_dir="./cache"
    ... )
    """
    os.makedirs(output_dir, exist_ok=True)
    tasks = []
    cache_paths = {}
    
    for mcool_path in mcool_paths:
        sample = os.path.basename(mcool_path).split('.')[0]
        for chrom in chroms:
            cache_path = os.path.join(
                output_dir, 
                f"{sample}.{chrom}.{resolution}bp.oe.npy"
            )
            cache_paths[(sample, chrom)] = cache_path
            tasks.append((mcool_path, chrom, resolution, cache_path))
    
    def compute_and_cache(args):
        mcool_path, chrom, resolution, cache_path = args
        if os.path.exists(cache_path):
            return cache_path
        try:
            _, _, oe, _ = process_oe_normalization(
                mcool_path=mcool_path,
                chrom=chrom,
                resolution=resolution,
                output_dir=None
            )
            np.save(cache_path, oe)
            return cache_path
        except Exception as e:
            print(f"O/E computation failed: {mcool_path} {chrom}: {e}")
            return None
    
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        list(executor.map(compute_and_cache, tasks))
    
    return cache_paths


def batch_load_oe(
    cache_paths: Dict[tuple, str],
    samples: List[str],
    chroms: List[str]
) -> Dict[tuple, np.ndarray]:
    """
    Load cached O/E matrices.
    
    Parameters
    ----------
    cache_paths : dict
        Dictionary mapping (sample, chrom) to file path
    samples : list
        List of sample names
    chroms : list
        List of chromosome names
        
    Returns
    -------
    matrices : dict
        Dictionary mapping (sample, chrom) to O/E matrix
    """
    matrices = {}
    for sample in samples:
        for chrom in chroms:
            path = cache_paths.get((sample, chrom))
            if path and os.path.exists(path):
                matrices[(sample, chrom)] = np.load(path)
            else:
                print(f"O/E cache not found: {sample} {chrom}")
    return matrices


# Convenience function for quick O/E normalization
def quick_oe(
    mcool_path: str,
    chrom: str,
    resolution: int = 10000,
    balance: bool = True
) -> Tuple[np.ndarray, Dict]:
    """
    Quick O/E normalization for a single region.
    
    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    chrom : str
        Chromosome name
    resolution : int
        Resolution in bp
    balance : bool
        Whether to use balanced matrix
        
    Returns
    -------
    oe : np.ndarray
        O/E normalized matrix
    stats : dict
        Normalization statistics
        
    Examples
    --------
    >>> OE, stats = quick_oe("sample.mcool", "chr1")
    """
    _, _, OE, stats = process_oe_normalization(
        mcool_path=mcool_path,
        chrom=chrom,
        resolution=resolution,
        output_dir=None,
        balance=balance
    )
    return OE, stats


def calculate_decay(
    matrix: np.ndarray,
    resolution: int,
    max_distance: int = 20000000
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate distance decay curve.
    
    This is a simple wrapper around the distance module function
    for convenience.
    
    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    resolution : int
        Bin resolution (bp)
    max_distance : int
        Maximum distance to include (bp)
        
    Returns
    -------
    distances : np.ndarray
        Distance values (bp)
    contacts : np.ndarray
        Average contact frequency at each distance
    """
    from .distance import calculate_distance_decay
    return calculate_distance_decay(matrix, resolution, max_distance)


def log2_and_mask(
    matrix: np.ndarray,
    mask_value: float = 0
) -> np.ndarray:
    """
    Apply log2 transformation and mask invalid values.
    
    Parameters
    ----------
    matrix : np.ndarray
        Input matrix
    mask_value : float
        Value to use for masked elements
        
    Returns
    -------
    result : np.ndarray
        Log2 transformed matrix
    """
    result = np.log2(matrix + 1e-10)
    result = np.where(np.isfinite(result), result, mask_value)
    return result


# =============================================================================
# Cooltools-based O/E calculation
# =============================================================================

def compute_expected_cooltools(
    mcool_path: str,
    chrom: str,
    resolution: int = 1000000,
    balance: bool = True,
    smooth: bool = True,
    ignore_diags: int = 2
) -> Tuple[np.ndarray, int]:
    """
    Compute expected matrix using cooltools.expected_cis.
    
    This function uses the cooltools library to calculate the expected
    interaction frequencies as a function of genomic distance.
    
    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    chrom : str
        Chromosome name
    resolution : int
        Bin resolution in bp
    balance : bool
        Whether to use balanced matrix
    smooth : bool
        Apply smoothing to expected curve
    ignore_diags : int
        Number of diagonals to ignore at beginning
        
    Returns
    -------
    expected_curve : np.ndarray
        1D array with expected value for each distance (bin)
    n_bins : int
        Number of bins in the chromosome
        
    Examples
    --------
    >>> expected, n = compute_expected_cooltools(
    ...     "sample.mcool", "chr1", resolution=1000000
    ... )
    """
    try:
        import cooltools
        import cooler
    except ImportError:
        raise ImportError("cooltools and cooler are required: pip install cooltools cooler")
    
    # Open cooler file
    clr = cooler.Cooler(f"{mcool_path}::/resolutions/{resolution}")
    
    # Calculate expected using cooltools
    expected_df = cooltools.expected_cis(
        clr,
        view_df=None,
        intra_only=True,
        smooth=smooth,
        aggregate_smoothed=True,
        smooth_sigma=0.1,
        clr_weight_name='weight' if balance else None,
        ignore_diags=ignore_diags
    )
    
    # Filter for the specific chromosome (column name is 'region1' in cooltools)
    expected_df = expected_df[expected_df['region1'] == chrom]
    
    if len(expected_df) == 0:
        raise ValueError(f"No expected data found for chromosome {chrom}")
    
    # Get chromosome length
    chrom_length = clr.chromsizes[chrom]
    n_bins = chrom_length // resolution
    
    # Get expected curve (use balanced.avg if available, otherwise count.avg)
    if 'balanced.avg' in expected_df.columns:
        expected_curve = expected_df['balanced.avg'].values
    else:
        expected_curve = expected_df['count.avg'].values
    
    # Get distances in bins
    distances = expected_df['dist'].values
    
    # Pad to full length if needed
    if len(expected_curve) < n_bins:
        # Extend with the last expected value
        expected_curve = np.pad(expected_curve, (0, n_bins - len(expected_curve)), 
                                mode='edge')
    
    return expected_curve[:n_bins], n_bins


def calculate_oe_matrix_cooltools(
    mcool_path: str,
    chrom: str,
    resolution: int = 1000000,
    balance: bool = False,
    view_df=None,
    smooth: bool = True,
    ignore_diags: int = 2,
    nproc: int = 1
) -> Tuple[np.ndarray, Dict]:
    """
    Calculate O/E matrix using cooltools.
    
    This is the preferred method for O/E calculation as it uses
    the well-tested cooltools library.
    
    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    chrom : str
        Chromosome name
    resolution : int
        Bin resolution in bp
    balance : bool
        Whether to use balanced matrix (ICE normalized).
        Default is False (use raw counts for both matrix and expected curve).
    view_df : DataFrame, optional
        ViewFrame for specifying genomic regions
    smooth : bool
        Apply smoothing to expected curve
    ignore_diags : int
        Number of diagonals to ignore
    nproc : int
        Number of processes for calculation
        
    Returns
    -------
    oe_matrix : np.ndarray
        O/E normalized matrix
    stats : dict
        Statistics about the normalization
        
    Examples
    --------
    >>> OE, stats = calculate_oe_matrix_cooltools(
    ...     "sample.mcool", "chr1", resolution=1000000
    ... )
    """
    try:
        import cooltools
        import cooler
    except ImportError:
        raise ImportError("cooltools and cooler are required: pip install cooltools cooler")
    
    from ..io.cooler import CoolerReader, extract_contact_matrix
    
    # Open cooler file
    clr = cooler.Cooler(f"{mcool_path}::/resolutions/{resolution}")
    
    # Get chromosome length and matrix
    chrom_length = clr.chromsizes[chrom]
    n_bins = chrom_length // resolution
    
    # Fetch the contact matrix
    reader = CoolerReader(mcool_path, resolution=resolution)
    O, _ = extract_contact_matrix(reader, chrom, balance=balance)
    
    # Ensure matrix is square
    O = O[:n_bins, :n_bins]
    
    # Calculate expected using cooltools
    expected_df = cooltools.expected_cis(
        clr,
        view_df=view_df,
        intra_only=True,
        smooth=smooth,
        aggregate_smoothed=True,
        smooth_sigma=0.1,
        clr_weight_name='weight' if balance else None,
        ignore_diags=ignore_diags,
        nproc=nproc
    )
    
    # Filter for the specific chromosome (column name is 'region1' in cooltools)
    expected_df = expected_df[expected_df['region1'] == chrom]
    
    # Build expected matrix from curve
    n = O.shape[0]
    E = np.zeros((n, n), dtype=float)
    
    # Get expected values by distance
    if 'balanced.avg' in expected_df.columns:
        col_name = 'balanced.avg'
    else:
        col_name = 'count.avg'
    
    for _, row in expected_df.iterrows():
        dist = int(row['dist'])
        if 0 <= dist < n:
            expected_val = row[col_name]
            if np.isfinite(expected_val):
                # Fill this diagonal (upper triangle for positive distance)
                for i in range(n - dist):
                    E[i, i + dist] = expected_val
    
    # Mirror to lower triangle
    E = E + E.T - np.diag(np.diag(E))
    
    # Calculate O/E matrix
    with np.errstate(divide='ignore', invalid='ignore'):
        OE = np.where(E > 0, O / E, np.nan)
    
    # Calculate statistics
    valid_oe = OE[~np.isnan(OE)]
    stats = {
        'n_elements': len(valid_oe),
        'min': np.min(valid_oe) if len(valid_oe) > 0 else np.nan,
        'max': np.max(valid_oe) if len(valid_oe) > 0 else np.nan,
        'mean': np.mean(valid_oe) if len(valid_oe) > 0 else np.nan,
        'median': np.median(valid_oe) if len(valid_oe) > 0 else np.nan,
        'std': np.std(valid_oe) if len(valid_oe) > 0 else np.nan,
        'method': 'cooltools.expected_cis'
    }
    
    return OE, stats


def quick_oe_cooltools(
    mcool_path: str,
    chrom: str,
    resolution: int = 1000000,
    balance: bool = False
) -> Tuple[np.ndarray, Dict]:
    """
    Quick O/E calculation using cooltools (convenience function).
    
    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
    chrom : str
        Chromosome name
    resolution : int
        Bin resolution in bp
    balance : bool
        Whether to use balanced matrix
        
    Returns
    -------
    oe_matrix : np.ndarray
        O/E normalized matrix
    stats : dict
        Statistics about the normalization
        
    Examples
    --------
    >>> OE, stats = quick_oe_cooltools("sample.mcool", "chr1", resolution=1000000)
    """
    return calculate_oe_matrix_cooltools(
        mcool_path=mcool_path,
        chrom=chrom,
        resolution=resolution,
        balance=balance
    )


def save_oe_matrix_npy(
    oe_matrix: np.ndarray,
    output_path: str,
    metadata: dict = None
) -> str:
    """
    保存 O/E 矩阵为 npy 格式
    
    Parameters
    ----------
    oe_matrix : np.ndarray
        O/E 矩阵
    output_path : str
        输出文件路径
    metadata : dict, optional
        元数据信息 (chrom, resolution, sample_name 等)
        
    Returns
    -------
    str : 保存的文件路径
        
    Examples
    --------
    >>> save_oe_matrix_npy(oe_matrix, "chr2_1M.oe.npy", {"chrom": "chr2", "resolution": 1000000})
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # 使用 float32 节省空间
    np.save(output_path, oe_matrix.astype(np.float32))
    
    # 如果提供了元数据，同时保存
    if metadata is not None:
        import json
        meta_path = output_path.replace('.npy', '_meta.json')
        with open(meta_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"元数据已保存: {meta_path}")
    
    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"O/E 矩阵已保存: {output_path} ({size_mb:.2f} MB)")
    
    return output_path


def load_oe_matrix_npy(
    npy_path: str,
    metadata_path: str = None
) -> Tuple[np.ndarray, dict]:
    """
    从 npy 文件加载 O/E 矩阵
    
    Parameters
    ----------
    npy_path : str
        npy 文件路径
    metadata_path : str, optional
        元数据文件路径
        
    Returns
    -------
    Tuple[np.ndarray, dict] : (O/E 矩阵, 元数据)
        
    Examples
    --------
    >>> oe_matrix, metadata = load_oe_matrix_npy("chr2_1M.oe.npy")
    """
    if not os.path.exists(npy_path):
        raise FileNotFoundError(f"文件不存在: {npy_path}")
    
    oe_matrix = np.load(npy_path)
    print(f"O/E 矩阵已加载: {npy_path} (shape: {oe_matrix.shape})")
    
    metadata = {}
    if metadata_path is None:
        metadata_path = npy_path.replace('.npy', '_meta.json')
    
    if os.path.exists(metadata_path):
        import json
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        print(f"元数据已加载: {metadata_path}")
    
    return oe_matrix, metadata


def get_npy_cache_path(
    output_dir: str,
    sample_name: str,
    chrom: str,
    resolution: int
) -> str:
    """
    生成 O/E 矩阵 npy 缓存文件的路径
    
    Parameters
    ----------
    output_dir : str
        输出目录
    sample_name : str
        样本名称
    chrom : str
        染色体名称
    resolution : int
        分辨率 (bp)
        
    Returns
    -------
    str : 缓存文件路径
    """
    res_str = f"{resolution // 1000}k" if resolution < 1e6 else f"{resolution // int(1e6)}M"
    filename = f"{sample_name}_{chrom}_{res_str}.oe.npy"
    return os.path.join(output_dir, filename)


def load_or_compute_oe_matrix(
    mcool_path: str,
    chrom: str,
    resolution: int,
    output_dir: str,
    sample_name: str = None,
    balance: bool = False,
    force_recompute: bool = False
) -> Tuple[np.ndarray, dict]:
    """
    加载或计算 O/E 矩阵（带缓存）
    
    Parameters
    ----------
    mcool_path : str
        mcool 文件路径
    chrom : str
        染色体名称
    resolution : int
        分辨率 (bp)
    output_dir : str
        输出目录
    sample_name : str, optional
        样本名称
    balance : bool
        是否使用平衡矩阵
    force_recompute : bool
        是否强制重新计算
        
    Returns
    -------
    Tuple[np.ndarray, dict] : (O/E 矩阵, 元数据)
    """
    if sample_name is None:
        sample_name = os.path.basename(mcool_path).split('.')[0]
    
    npy_path = get_npy_cache_path(output_dir, sample_name, chrom, resolution)
    
    # 检查缓存
    if not force_recompute and os.path.exists(npy_path):
        print(f"📦 [缓存命中] 从 npy 加载 O/E 矩阵: {npy_path}")
        return load_oe_matrix_npy(npy_path)
    
    # 计算 O/E 矩阵
    print(f"🧮 计算 O/E 矩阵: {chrom} @ {resolution} bp")
    oe_matrix, stats = calculate_oe_matrix_cooltools(
        mcool_path=mcool_path,
        chrom=chrom,
        resolution=resolution,
        balance=balance
    )
    
    # 获取染色体长度（用于计算 start_pos）
    import cooler
    clr = cooler.Cooler(f'{mcool_path}::/resolutions/{resolution}')
    chrom_length = clr.chromsizes[chrom]
    
    # 构建元数据（确保值为 Python 原生类型）
    metadata = {
        'mcool_path': mcool_path,
        'chrom': chrom,
        'resolution': int(resolution),
        'start_pos': 0,
        'end_pos': int(chrom_length),
        'sample_name': sample_name,
        'balance': balance,
        'shape': list(oe_matrix.shape),
        'stats': stats
    }
    
    # 保存为 npy
    save_oe_matrix_npy(oe_matrix, npy_path, metadata)
    
    return oe_matrix, metadata
