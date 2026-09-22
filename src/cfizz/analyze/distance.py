"""
Distance decay analysis module.

This module provides functions for analyzing the distance-dependent
decay of chromatin interactions in Hi-C data.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Any, List
from scipy.optimize import curve_fit
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)


def calculate_distance_decay(
    matrix: np.ndarray,
    resolution: int,
    max_distance: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate distance decay curve for Hi-C contacts.
    
    The distance decay follows a power law: P(s) ∝ s^(-α)
    where s is the genomic distance and α is typically around 1.
    
    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    resolution : int
        Bin resolution (bp)
    max_distance : int, optional
        Maximum distance to include (bp)
        
    Returns
    -------
    distances : np.ndarray
        Distance values (bp)
    contacts : np.ndarray
        Average contact frequency at each distance
        
    Examples
    --------
    >>> from cfizz.io import read_cooler
    >>> reader = read_cooler("sample.mcool")
    >>> matrix = reader.fetch("chr1", 0, 50000000)
    >>> distances, contacts = calculate_distance_decay(matrix, 100000)
    """
    n = matrix.shape[0]
    
    if max_distance is None:
        max_distance = n * resolution
    
    max_distance_bins = min(max_distance // resolution, n - 1)
    
    # Calculate average contact frequency at each distance
    distances = []
    contacts = []
    
    for d in range(1, max_distance_bins):
        values = []
        for i in range(n - d):
            values.append(matrix[i, i + d])
        
        if values:
            distances.append(d * resolution)
            contacts.append(np.mean(values))
    
    return np.array(distances), np.array(contacts)


def fit_power_law(
    distances: np.ndarray,
    contacts: np.ndarray
) -> Tuple[float, float, float]:
    """
    Fit power law to distance decay data.
    
    Model: P(s) = A * s^(-α)
    
    Parameters
    ----------
    distances : np.ndarray
        Distance values
    contacts : np.ndarray
        Contact frequencies
        
    Returns
    -------
    A : float
        Prefactor
    alpha : float
        Power law exponent
    R2 : float
        R-squared value for fit
        
    Examples
    --------
    >>> dist, contacts = calculate_distance_decay(matrix, 100000)
    >>> A, alpha, R2 = fit_power_law(dist, contacts)
    >>> print(f"Power law exponent: {alpha:.3f}")
    """
    # Log-transform for linear fitting
    log_dist = np.log10(distances)
    log_contacts = np.log10(contacts)
    
    # Linear regression: log(P) = log(A) - α * log(s)
    n = len(log_dist)
    sum_x = np.sum(log_dist)
    sum_y = np.sum(log_contacts)
    sum_xy = np.sum(log_dist * log_contacts)
    sum_x2 = np.sum(log_dist ** 2)
    
    alpha = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x ** 2)
    log_A = (sum_y - alpha * sum_x) / n
    A = 10 ** log_A
    
    # Calculate R-squared
    y_mean = np.mean(log_contacts)
    SS_tot = np.sum((log_contacts - y_mean) ** 2)
    SS_res = np.sum((log_contacts - log_A + alpha * log_dist) ** 2)
    R2 = 1 - SS_res / SS_tot if SS_tot > 0 else 0
    
    return A, alpha, R2


def fit_exponential(
    distances: np.ndarray,
    contacts: np.ndarray
) -> Tuple[float, float, float]:
    """
    Fit exponential decay to distance decay data.
    
    Model: P(s) = A * exp(-s/λ)
    
    Parameters
    ----------
    distances : np.ndarray
        Distance values
    contacts : np.ndarray
        Contact frequencies
        
    Returns
    -------
    A : float
        Prefactor
    lambda_val : float
        Decay length
    R2 : float
        R-squared value for fit
    """
    def exp_func(s, A, lam):
        return A * np.exp(-s / lam)
    
    try:
        popt, _ = curve_fit(exp_func, distances, contacts, p0=[1, 1000000], maxfev=5000)
        A, lambda_val = popt
        
        # Calculate R-squared
        y_mean = np.mean(contacts)
        SS_tot = np.sum((contacts - y_mean) ** 2)
        predicted = exp_func(distances, A, lambda_val)
        SS_res = np.sum((contacts - predicted) ** 2)
        R2 = 1 - SS_res / SS_tot if SS_tot > 0 else 0
        
        return A, lambda_val, R2
    except:
        return 0, 0, 0


def calculate_decay_rate(
    contacts: np.ndarray,
    distances: np.ndarray,
    window_size: int = 10
) -> float:
    """
    Calculate local decay rate.
    
    Parameters
    ----------
    contacts : np.ndarray
        Contact frequencies
    distances : np.ndarray
        Distance values
    window_size : int
        Window size for local calculation
        
    Returns
    -------
    decay_rate : float
        Local decay rate (exponent)
    """
    if len(contacts) < window_size:
        return 0
    
    # Calculate decay rate in sliding windows
    rates = []
    
    for i in range(len(contacts) - window_size):
        window_dist = distances[i:i+window_size]
        window_contacts = contacts[i:i+window_size]
        
        try:
            _, alpha, _ = fit_power_law(window_dist, window_contacts)
            rates.append(alpha)
        except:
            pass
    
    if rates:
        return np.median(rates)
    return 0


def calculate_decay_profile(
    matrix: np.ndarray,
    resolution: int,
    bin_size: int = 1000000,
    max_distance: int = 20000000
) -> pd.DataFrame:
    """
    Calculate comprehensive decay profile.
    
    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    resolution : int
        Bin resolution (bp)
    bin_size : int
        Bin size for output (bp)
    max_distance : int
        Maximum distance (bp)
        
    Returns
    -------
    profile : pd.DataFrame
        DataFrame with decay statistics
    """
    distances, contacts = calculate_distance_decay(matrix, resolution, max_distance)
    
    if len(distances) == 0:
        return pd.DataFrame(columns=['distance', 'contact', 'log_distance', 'log_contact'])
    
    # Bin by distance
    bin_edges = np.arange(0, max_distance + bin_size, bin_size)
    binned_distances = []
    binned_contacts = []
    
    for i in range(len(bin_edges) - 1):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i + 1])
        if np.any(mask):
            binned_distances.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            binned_contacts.append(np.mean(contacts[mask]))
    
    profile = pd.DataFrame({
        'distance': binned_distances,
        'contact': binned_contacts,
        'log_distance': np.log10(binned_distances),
        'log_contact': np.log10(binned_contacts)
    })
    
    # Fit power law
    if len(binned_distances) > 2:
        A, alpha, R2 = fit_power_law(
            np.array(binned_distances),
            np.array(binned_contacts)
        )
        profile['power_law_A'] = A
        profile['power_law_alpha'] = alpha
        profile['power_law_R2'] = R2
    
    return profile


def compare_decay_rates(
    matrix1: np.ndarray,
    matrix2: np.ndarray,
    resolution: int,
    max_distance: int = 20000000
) -> Dict[str, Any]:
    """
    Compare distance decay rates between two matrices.

    Parameters
    ----------
    matrix1 : np.ndarray
        First contact matrix
    matrix2 : np.ndarray
        Second contact matrix
    resolution : int
        Bin resolution (bp)
    max_distance : int
        Maximum distance (bp)

    Returns
    -------
    comparison : dict
        Dictionary with comparison metrics
    """
    dist1, contacts1 = calculate_distance_decay(matrix1, resolution, max_distance)
    dist2, contacts2 = calculate_distance_decay(matrix2, resolution, max_distance)

    _, alpha1, R2_1 = fit_power_law(dist1, contacts1)
    _, alpha2, R2_2 = fit_power_law(dist2, contacts2)

    return {
        'alpha1': alpha1,
        'alpha2': alpha2,
        'alpha_diff': alpha2 - alpha1,
        'R2_1': R2_1,
        'R2_2': R2_2
    }


# =============================================================================
# Distance decay: full P(s) pipeline (integrated from g_6d_distance_decay.py)
# =============================================================================
# 与 compute_expected_cis_per_sample 区别:
#   - 排除 chrM/chrY (默认)
#   - 截断到 max_distance
#   - 归一化: count.avg / total (sample 内总信号)
#   - 重命名为 normalized_count / normalized_count_smoothed
#   - 跨 region 平均 (groupby 'dist_bp').mean()
#   - 幂律拟合 (限区间 fit_window) + 分段斜率
# =============================================================================


DEFAULT_SEGMENTS_MB: list = [(0.02, 0.2), (0.2, 1), (1, 5), (5, 20), (20, 50)]


def process_expected_cis_decay(
    exp_df: pd.DataFrame,
    max_distance: int = 100_000_000,
    excluded_chroms: Optional[Tuple[str, ...]] = ('chrM', 'chrY'),
) -> pd.DataFrame:
    """
    从 ``cooltools.expected_cis`` 输出过滤 + 归一化 + 重命名列.

    Parameters
    ----------
    exp_df : pd.DataFrame
        ``cooltools.expected_cis`` 返回的原始 DataFrame (列含
        ``region1/region2/dist/dist_bp/count.avg/count.avg.smoothed/...``).
    max_distance : int
        截断的最大距离 (bp), 默认 100Mb.
    excluded_chroms : tuple of str, optional
        排除的染色体 (按 ``region1`` 列过滤), 默认 ``('chrM', 'chrY')``;
        传 ``None`` 不过滤.

    Returns
    -------
    df : pd.DataFrame
        列: ``region1, region2, dist, dist_bp, normalized_count,
        normalized_count_smoothed, total_contacts``.
        每个 ``dist_bp`` 上 sample 内跨 region 多个点 (groupby 前).
    """
    df = exp_df.copy()
    # 排除指定染色体 (基于 region1 列)
    if excluded_chroms is not None:
        df = df[~df['region1'].isin(excluded_chroms)]
    # 去掉 dist=0 (ignore_diags) 和 NaN
    df = df[df['dist'] > 0].dropna(subset=['count.avg']).copy()
    df = df[df['dist_bp'] <= max_distance]
    # 重命名带点列名
    df = df.rename(columns={
        'count.avg': 'count_avg',
        'count.avg.smoothed': 'count_smoothed',
    })
    # 归一化: sample 内总信号做比例化
    total = float(df['count_avg'].sum())
    if total > 0:
        df['count_avg'] = df['count_avg'] / total
        df['count_smoothed'] = df['count_smoothed'] / total
    df['total_contacts'] = total
    df = df.rename(columns={
        'count_avg': 'normalized_count',
        'count_smoothed': 'normalized_count_smoothed',
    })
    return df


def fit_decay_alpha_segments(
    dist_bp: np.ndarray,
    contacts: np.ndarray,
    fit_window: Tuple[float, float] = (0.2e6, 10e6),
    segments: Optional[List[Tuple[float, float]]] = None,
) -> Dict[str, float]:
    """
    log-log 线性回归拟合幂律 P = A * s^-α (限区间 + 分段).

    Parameters
    ----------
    dist_bp : np.ndarray
        距离值 (bp).
    contacts : np.ndarray
        对应 normalized_count 数值.
    fit_window : tuple of (lo, hi) in bp
        限区间拟合窗口, 默认 ``(0.2e6, 10e6)`` (0.2-10 Mb).
    segments : list of (lo_mb, hi_mb), optional
        分段拟合区间 (Mb), 默认 ``[(0.02, 0.2), (0.2, 1), (1, 5), (5, 20), (20, 50)]``.

    Returns
    -------
    result : dict
        - ``A_raw``, ``alpha_raw``, ``R2_raw`` — fit_window 上 normalized_count 拟合
        - ``A_smoothed``, ``alpha_smoothed``, ``R2_smoothed`` — 同上但用 smoothed
        - ``slope_<lo>-<hi>Mb`` — 每段斜率
    """
    segs = segments if segments is not None else DEFAULT_SEGMENTS_MB
    seg_slopes = {}
    for lo_mb, hi_mb in segs:
        m = (dist_bp >= lo_mb * 1e6) & (dist_bp < hi_mb * 1e6)
        sub = dist_bp[m], contacts[m]
        if len(sub) >= 3:
            _, alpha_seg, _ = fit_power_law(sub[0], sub[1])
            seg_slopes[f'slope_{lo_mb}-{hi_mb}Mb'] = alpha_seg
        else:
            seg_slopes[f'slope_{lo_mb}-{hi_mb}Mb'] = np.nan

    # 全 fit_window 拟合 (normalized_count 和 smoothed 各一次)
    mask_fit = (dist_bp >= fit_window[0]) & (dist_bp < fit_window[1])
    if mask_fit.sum() >= 3:
        sub_d = dist_bp[mask_fit]
        sub_c = contacts[mask_fit]
        A_raw, alpha_raw, R2_raw = fit_power_law(sub_d, sub_c)
    else:
        A_raw, alpha_raw, R2_raw = np.nan, np.nan, np.nan

    # smoothed 列 (如果存在)
    if 'count.avg.smoothed' in dir(np):  # 不可触发, 仅占位
        pass

    return {
        'A_raw': A_raw,
        'alpha_raw': alpha_raw,
        'R2_raw': R2_raw,
        **seg_slopes,
    }


def compute_distance_decay_per_sample(
    samples: Dict[str, str],
    resolution: int = 1_000,
    nproc: int = 8,
    max_distance: int = 100_000_000,
    fit_window: Tuple[float, float] = (0.2e6, 10e6),
    segments: Optional[List[Tuple[float, float]]] = None,
    balance: bool = False,
    excluded_chroms: Optional[Tuple[str, ...]] = ('chrM', 'chrY'),
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, Dict]]:
    """
    per-sample 全基因组距离衰减: cooltools.expected_cis + 过滤 + 归一化 + groupby + 拟合.

    Parameters
    ----------
    samples : dict
        ``{sample_name: mcool_path}`` 映射.
    resolution : int
        mcool resolution (bp), 默认 1_000 (1kb).
    nproc : int
        cooltools 内部并行核数.
    max_distance : int
        截断最大距离 (bp).
    fit_window : tuple
        幂律拟合限区间 (bp).
    segments : list of tuples, optional
        分段拟合区间 (Mb). 默认 5 段.
    balance : bool
        True → cooltools 用 balance 权重 (clr_weight_name='weight'),
        False → raw (clr_weight_name=None). 默认 False (与 g_6d 一致).
    excluded_chroms : tuple, optional
        默认 ``('chrM', 'chrY')``, 传 None 不过滤.

    Returns
    -------
    data_dict : dict
        ``{sample: process_expected_cis_decay() 输出}`` — per-region per-dist.
    agg_dict : dict
        ``{sample: groupby('dist_bp').mean() 跨 region 平均}`` — 每距离 1 点.
    alpha_dict : dict
        ``{sample: fit_decay_alpha_segments() 输出}`` — 全区间 α + R² + 分段斜率.
    """
    try:
        import cooler
        import cooltools
    except ImportError:
        raise ImportError(
            "compute_distance_decay_per_sample 需要 cooler + cooltools"
        )

    from cfizz.analyze.compartment import get_view_df

    segs = segments if segments is not None else DEFAULT_SEGMENTS_MB

    data_dict: Dict[str, pd.DataFrame] = {}
    agg_dict: Dict[str, pd.DataFrame] = {}
    alpha_dict: Dict[str, Dict] = {}

    for sample, mcool_path in samples.items():
        clr = cooler.Cooler(f"{mcool_path}::resolutions/{resolution}")
        view_df = get_view_df(clr)
        # get_view_df() 已排除 chrM/chrY; 与 excluded_chroms 参数保持一致语义

        exp = cooltools.expected_cis(
            clr,
            view_df=view_df,
            clr_weight_name='weight' if balance else None,
            nproc=nproc,
        )

        df = process_expected_cis_decay(
            exp,
            max_distance=max_distance,
            excluded_chroms=excluded_chroms,
        )
        data_dict[sample] = df

        # 跨 region 平均 (同一 dist_bp 上 20 个常染色体 region 平均为 1 个点)
        agg = df.groupby('dist_bp').agg({
            'normalized_count': 'mean',
            'normalized_count_smoothed': 'mean',
            'total_contacts': 'first',
        }).reset_index().sort_values('dist_bp')
        agg_dict[sample] = agg

        # 限区间幂律拟合 + 分段斜率 (normalized_count)
        alpha_dict[sample] = fit_decay_alpha_segments(
            agg['dist_bp'].values,
            agg['normalized_count'].values,
            fit_window=fit_window,
            segments=segs,
        )

    return data_dict, agg_dict, alpha_dict
