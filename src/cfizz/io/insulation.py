"""
Insulation score computation for TAD identification.

This module provides functions to calculate insulation scores
and identify TAD boundaries from Hi-C data using cooltools.

Key functions:
  - compute_insulation_from_cooler(): Compute insulation scores from mcool file
  - calculate_insulation_score(): Compute from matrix (for loaded data)

"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
import warnings
import os
from pathlib import Path

warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# Cooler-based functions (using cooltools)
# =============================================================================

def compute_insulation_from_cooler(
    cooler_path: str,
    resolution: int,
    windows: List[int],
    nproc: int = 1,
    chunksize: int = 50_000_000,
    ignore_diags: int = 2,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Compute insulation scores from a cooler file using cooltools.
    
    This is the main function for computing insulation scores from mcool files.
    
    Parameters
    ----------
    cooler_path : str
        Path to mcool file (e.g., "sample.mcool")
    resolution : int
        Resolution in bp (e.g., 10000 for 10kb, 100000 for 100kb)
    windows : list
        List of window sizes in bp (e.g., [50000, 100000, 500000])
    nproc : int
        Number of parallel processes
    chunksize : int
        Data chunk size for parallel processing
    ignore_diags : int
        Number of diagonals to ignore
    verbose : bool
        Print progress information
        
    Returns
    -------
    pd.DataFrame
        Insulation scores DataFrame with columns:
        - chrom, start, end, binsize
        - insulation_{window} for each window size
        - boundary_strength_{window} for each window size
        - is_boundary_{window} for each window size
        
    Examples
    --------
    >>> df = compute_insulation_from_cooler(
    ...     "sample.mcool",
    ...     resolution=10000,
    ...     windows=[50000, 100000, 500000],
    ...     nproc=1
    ... )
    >>> df.head()
    """
    from cooltools.api.insulation import insulation
    
    # Load cooler
    clr = cooler.Cooler(f"{cooler_path}::/resolutions/{resolution}")
    
    # Compute insulation scores
    insulation_table = insulation(
        clr,
        window_bp=sorted(windows),
        nproc=nproc,
        chunksize=chunksize,
        ignore_diags=ignore_diags,
        verbose=verbose
    )
    
    return insulation_table


def compute_insulation_and_save(
    cooler_path: str,
    resolution: int,
    windows: List[int],
    output_file: str,
    nproc: int = 1,
    chunksize: int = 50_000_000,
    ignore_diags: int = 2,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Compute insulation scores and save to file.
    
    Parameters
    ----------
    cooler_path : str
        Path to mcool file
    resolution : int
        Resolution in bp
    windows : list
        List of window sizes in bp
    output_file : str
        Output TSV file path
    nproc : int
        Number of parallel processes
    chunksize : int
        Data chunk size
    ignore_diags : int
        Number of diagonals to ignore
    verbose : bool
        Print progress information
        
    Returns
    -------
    pd.DataFrame
        Insulation scores DataFrame
    """
    df = compute_insulation_from_cooler(
        cooler_path=cooler_path,
        resolution=resolution,
        windows=windows,
        nproc=nproc,
        chunksize=chunksize,
        ignore_diags=ignore_diags,
        verbose=verbose
    )
    
    # Save to file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, sep='\t')
    
    return df


# =============================================================================
# Matrix-based functions (for loaded data)
# =============================================================================

def calculate_insulation_score(
    matrix: np.ndarray,
    window_sizes: List[int] = [5, 10, 20],
    ignore_diagonals: int = 2
) -> pd.DataFrame:
    """
    Calculate insulation score for TAD boundary identification from matrix.
    
    This function works on loaded matrix data (not directly from cooler).
    For cooler files, use compute_insulation_from_cooler() instead.

    Insulation score measures how isolated a genomic region is from
    its neighbors. Low insulation scores indicate potential TAD boundaries.

    The formula is:
    IS(b) = log2( sum(contacts in window) / n^2 )

    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    window_sizes : list
        List of window sizes (in bins) to use
    ignore_diagonals : int
        Number of diagonals to ignore

    Returns
    -------
    insulation_df : pd.DataFrame
        DataFrame with columns [start, end, insulation_score]

    Examples
    --------
    >>> from cfizz.io import read_cooler
    >>> reader = read_cooler("data.mcool::/resolutions/100000")
    >>> matrix = reader.fetch("chr1", 0, 10000000)
    >>> ins_df = calculate_insulation_score(matrix, window_sizes=[10])
    >>> boundaries = ins_df[ins_df['is_boundary']].head()
    """
    n_bins = matrix.shape[0]

    # Create bin coordinates (assuming 1-indexed bins)
    starts = np.arange(n_bins)
    ends = starts + 1

    # Calculate insulation scores for each window size
    all_scores = {}

    for window_size in window_sizes:
        scores = np.zeros(n_bins)

        for i in range(n_bins):
            # Define the square window around bin i
            window_start = max(0, i - window_size)
            window_end = min(n_bins, i + window_size + 1)

            # Sum contacts in the window (excluding the diagonal)
            window_sum = 0
            count = 0

            for r in range(window_start, window_end):
                for c in range(window_start, window_end):
                    if abs(r - c) >= ignore_diagonals:
                        window_sum += matrix[r, c]
                        count += 1

            # Calculate insulation score
            if count > 0:
                scores[i] = window_sum / count
            else:
                scores[i] = np.nan

        # Log2 transform
        scores = np.log2(scores + 1)
        all_scores[f'insulation_{window_size}'] = scores

    # Create DataFrame
    insulation_df = pd.DataFrame({
        'start': starts,
        'end': ends,
        **all_scores
    })

    # Add default insulation score (use first window size)
    if window_sizes:
        insulation_df['insulation_score'] = all_scores[f'insulation_{window_sizes[0]}']

    # Identify boundaries using local minima
    insulation_df['is_boundary'] = _find_local_minima(
        insulation_df['insulation_score'].values,
        threshold=0.1
    )

    return insulation_df


# =============================================================================
# Helper functions
# =============================================================================

def _find_local_minima(
    values: np.ndarray,
    threshold: float = 0.1,
    min_distance: int = 5
) -> np.ndarray:
    """
    Find local minima in an array.

    Parameters
    ----------
    values : np.ndarray
        Array of values
    threshold : float
        Minimum relative drop to qualify as boundary
    min_distance : int
        Minimum distance between boundaries

    Returns
    -------
    is_boundary : np.ndarray
        Boolean array indicating boundary positions
    """
    n = len(values)
    is_boundary = np.zeros(n, dtype=bool)

    if n < 3:
        return is_boundary

    # Calculate the global median
    median_val = np.nanmedian(values)

    for i in range(1, n - 1):
        # Check if it's a local minimum
        if (values[i] <= values[i - 1] and
            values[i] <= values[i + 1] and
            np.isfinite(values[i])):

            # Check relative drop
            avg_neighbor = (values[i - 1] + values[i + 1]) / 2
            if avg_neighbor - values[i] > threshold * abs(avg_neighbor):
                is_boundary[i] = True

    # Remove boundaries that are too close
    boundary_indices = np.where(is_boundary)[0]
    if len(boundary_indices) > 1:
        keep = [boundary_indices[0]]
        for idx in boundary_indices[1:]:
            if idx - keep[-1] >= min_distance:
                keep.append(idx)

        is_boundary[:] = False
        is_boundary[keep] = True

    return is_boundary


def calculate_boundary_strength(
    insulation_df: pd.DataFrame,
    matrix: np.ndarray,
    resolution: int,
    window_size: int = 10
) -> pd.DataFrame:
    """
    Calculate boundary strength for each TAD boundary.

    Boundary strength is measured as the ratio of insulation
    at the boundary to the average insulation in the surrounding region.

    Parameters
    ----------
    insulation_df : pd.DataFrame
        DataFrame with insulation scores
    matrix : np.ndarray
        Hi-C contact matrix
    resolution : int
        Bin resolution (bp)
    window_size : int
        Window size for strength calculation

    Returns
    -------
    boundary_df : pd.DataFrame
        DataFrame with boundary positions and strengths
    """
    boundaries = insulation_df[insulation_df['is_boundary']].copy()

    strengths = []
    for _, row in boundaries.iterrows():
        pos = row['start']
        pos_bin = pos // resolution

        # Calculate average insulation in surrounding region
        window_start = max(0, pos_bin - window_size)
        window_end = min(len(insulation_df), pos_bin + window_size + 1)

        local_insulation = insulation_df.iloc[window_start:window_end]['insulation_score'].mean()
        boundary_insulation = row['insulation_score']

        # Boundary strength = local_avg / boundary_value
        if boundary_insulation > 0:
            strength = local_insulation / boundary_insulation
        else:
            strength = 1.0

        strengths.append(strength)

    boundaries = boundaries.copy()
    boundaries['boundary_strength'] = strengths

    return boundaries


# =============================================================================
# Utility functions
# =============================================================================

def read_insulation_scores(
    file_path: str,
    windows: int,
    chrom: str = None,
    start: int = None,
    end: int = None
) -> pd.DataFrame:
    """
    Read Insulation Score data for a specific window size.
    
    This function reads data from cooltools output insulation files
    and extracts data for a specific window size.
    
    Parameters
    ----------
    file_path : str
        Path to insulation score file (e.g., "1_0.sample.10kb.insulation.tsv")
    windows : int
        Window size in bp (e.g., 50000, 100000, 500000)
    chrom : str, optional
        Chromosome name to filter
    start : int, optional
        Start position to filter
    end : int, optional
        End position to filter
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, insulation_score, is_boundary
    
    Examples
    --------
    >>> ins_data = read_insulation_scores(
    ...     file_path="1_0.sample.10kb.insulation.tsv",
    ...     windows=100000,
    ...     chrom="chr2",
    ...     start=100_000_000,
    ...     end=110_000_000
    ... )
    """
    # Read the file
    df = pd.read_csv(file_path, sep='\t')
    
    # Build column names based on window size
    score_col = f'log2_insulation_score_{windows}'
    boundary_col = f'is_boundary_{windows}'
    
    # Check if columns exist
    if score_col not in df.columns:
        # Try alternative naming (without log2_)
        alt_score_cols = [c for c in df.columns if str(windows) in c and 'insulation' in c.lower()]
        alt_boundary_cols = [c for c in df.columns if str(windows) in c and 'boundary' in c.lower()]
        if alt_score_cols:
            score_col = alt_score_cols[0]
        if alt_boundary_cols:
            boundary_col = alt_boundary_cols[0]
    
    # Extract relevant columns
    result = pd.DataFrame({
        'chrom': df['chrom'],
        'start': df['start'],
        'end': df['end'],
        'insulation_score': df[score_col] if score_col in df.columns else df.get(score_col.replace('log2_', ''), df.iloc[:, 3]),
        'is_boundary': df.get(boundary_col, False) if boundary_col in df.columns else df.get(boundary_col.replace('is_', ''), False)
    })
    
    # Filter by chromosome
    if chrom is not None:
        result = result[result['chrom'] == chrom]
    
    # Filter by position
    if start is not None:
        result = result[result['start'] >= start]
    
    if end is not None:
        result = result[result['end'] <= end]
    
    return result


def read_insulation_from_file(file_path: str) -> pd.DataFrame:
    """
    Read insulation scores from a TSV file.
    
    Parameters
    ----------
    file_path : str
        Path to insulation TSV file
        
    Returns
    -------
    pd.DataFrame
        Insulation scores DataFrame
    """
    return pd.read_csv(file_path, sep='\t')


def get_available_chroms(cooler_path: str, resolution: int) -> list:
    """
    Get available chromosomes from a cooler file.
    
    Parameters
    ----------
    cooler_path : str
        Path to mcool file
    resolution : int
        Resolution in bp
        
    Returns
    -------
    list
        List of chromosome names
    """
    clr = cooler.Cooler(f"{cooler_path}::/resolutions/{resolution}")
    return list(clr.chromnames)


def format_resolution(resolution: int) -> str:
    """
    Format resolution as human-readable string.
    
    Parameters
    ----------
    resolution : int
        Resolution in bp
        
    Returns
    -------
    str
        Formatted string (e.g., "100kb" or "1Mb")
    """
    if resolution >= 1_000_000:
        return f"{resolution // 1_000_000}Mb"
    else:
        return f"{resolution // 1000}kb"


# Import cooler for type hints
import cooler