"""
TAD (Topologically Associating Domain) I/O and utilities.

This module provides functions to read/write TAD coordinates
and perform basic TAD operations.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any


def read_tads_bed(bed_path: str) -> pd.DataFrame:
    """
    Read TAD coordinates from BED file.
    
    BED format: chrom start end [score]
    
    Parameters
    ----------
    bed_path : str
        Path to BED file
        
    Returns
    -------
    tad_df : pd.DataFrame
        DataFrame with columns [chrom, start, end, score]
        
    Examples
    --------
    >>> tads = read_tads_bed("tads.bed")
    >>> print(f"Found {len(tads)} TADs")
    """
    try:
        df = pd.read_csv(
            bed_path, 
            sep='\t', 
            header=None, 
            names=['chrom', 'start', 'end', 'score'],
            usecols=range(4)
        )
    except ValueError:
        df = pd.read_csv(bed_path, sep='\t', header=None, 
                        names=['chrom', 'start', 'end'])
        df['score'] = 1.0
    
    return df


def write_tads_bed(tad_df: pd.DataFrame, output_path: str) -> None:
    """
    Write TAD coordinates to BED file.
    
    Parameters
    ----------
    tad_df : pd.DataFrame
        DataFrame with TAD coordinates
    output_path : str
        Output file path
        
    Examples
    --------
    >>> write_tads_bed(tad_df, "output_tads.bed")
    """
    required_cols = ['chrom', 'start', 'end']
    
    for col in required_cols:
        if col not in tad_df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    output_df = tad_df[required_cols].copy()
    
    if 'score' in tad_df.columns:
        output_df['score'] = tad_df['score']
    else:
        output_df['score'] = 1.0
    
    output_df.to_csv(output_path, sep='\t', header=False, index=False)


class TADSet:
    """
    Container for TAD calls.
    
    Provides utilities for filtering, merging, and analyzing TAD sets.
    
    Examples
    --------
    >>> tad_set = TADSet.from_bed("tads.bed")
    >>> filtered = tad_set.filter_by_size(min_size=5000000)
    >>> filtered.save("filtered_tads.bed")
    """
    
    def __init__(self, tad_df: pd.DataFrame):
        """
        Initialize TADSet.
        
        Parameters
        ----------
        tad_df : pd.DataFrame
            DataFrame with TAD coordinates
        """
        self.tads = tad_df.copy()
    
    @classmethod
    def from_bed(cls, path: str) -> 'TADSet':
        """Create TADSet from BED file."""
        return cls(read_tads_bed(path))
    
    def save(self, path: str) -> None:
        """Save to BED file."""
        write_tads_bed(self.tads, path)
    
    def filter_by_size(
        self, 
        min_size: int = None, 
        max_size: int = None
    ) -> 'TADSet':
        """
        Filter TADs by size.
        
        Parameters
        ----------
        min_size : int, optional
            Minimum TAD size (bp)
        max_size : int, optional
            Maximum TAD size (bp)
            
        Returns
        -------
        filtered : TADSet
        """
        filtered = self.tads.copy()
        
        sizes = filtered['end'] - filtered['start']
        
        if min_size is not None:
            filtered = filtered[sizes >= min_size]
            sizes = sizes[sizes >= min_size]
        
        if max_size is not None:
            filtered = filtered[sizes <= max_size]
        
        return TADSet(filtered)
    
    def filter_by_score(self, threshold: float) -> 'TADSet':
        """
        Filter TADs by score.
        
        Parameters
        ----------
        threshold : float
            Minimum score
            
        Returns
        -------
        filtered : TADSet
        """
        if 'score' not in self.tads.columns:
            return self
        
        filtered = self.tads[self.tads['score'] >= threshold].copy()
        return TADSet(filtered)
    
    def get_sizes(self) -> np.ndarray:
        """Get TAD sizes."""
        return (self.tads['end'] - self.tads['start']).values
    
    def get_statistics(self) -> Dict[str, float]:
        """
        Get TAD statistics.
        
        Returns
        -------
        stats : dict
            Dictionary with statistics
        """
        sizes = self.get_sizes()
        return {
            'n_tads': len(self.tads),
            'mean_size': np.mean(sizes),
            'median_size': np.median(sizes),
            'min_size': np.min(sizes),
            'max_size': np.max(sizes),
            'std_size': np.std(sizes)
        }
    
    def __len__(self) -> int:
        return len(self.tads)
    
    def __repr__(self) -> str:
        return f"TADSet(n={len(self.tads)})"


def tads_to_boundaries(tad_df: pd.DataFrame) -> List[int]:
    """
    Convert TADs to list of boundary positions.
    
    Parameters
    ----------
    tad_df : pd.DataFrame
        DataFrame with TAD coordinates
        
    Returns
    -------
    boundaries : list
        List of boundary positions
    """
    boundaries = []
    for _, tad in tad_df.iterrows():
        boundaries.append(tad['start'])
        boundaries.append(tad['end'])
    
    return sorted(set(boundaries))


def boundaries_to_tads(
    boundaries: List[int],
    min_size: int = 0
) -> pd.DataFrame:
    """
    Convert boundary positions to TAD coordinates.
    
    Parameters
    ----------
    boundaries : list
        List of boundary positions (sorted)
    min_size : int
        Minimum TAD size (bp)
        
    Returns
    -------
    tad_df : pd.DataFrame
        DataFrame with TAD coordinates
    """
    if len(boundaries) < 2:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'score'])
    
    tads = []
    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1]
        size = end - start
        
        if size >= min_size:
            tads.append({
                'chrom': 'chr1',  # Will need chrom info
                'start': start,
                'end': end,
                'size': size
            })
    
    tad_df = pd.DataFrame(tads)
    if len(tad_df) > 0:
        tad_df['score'] = 1.0
    
    return tad_df


def merge_tads(
    tad_df: pd.DataFrame,
    min_size: int = 0
) -> pd.DataFrame:
    """
    Merge overlapping or adjacent TADs.
    
    Parameters
    ----------
    tad_df : pd.DataFrame
        DataFrame with TAD coordinates
    min_size : int
        Minimum TAD size (bp)
        
    Returns
    -------
    merged : pd.DataFrame
        Merged TAD coordinates
    """
    if len(tad_df) == 0:
        return tad_df
    
    # Sort by start position
    sorted_tads = tad_df.sort_values('start').reset_index(drop=True)
    
    merged = []
    current = sorted_tads.iloc[0].to_dict()
    
    for i in range(1, len(sorted_tads)):
        row = sorted_tads.iloc[i]
        
        # Check if overlapping or adjacent
        if row['start'] <= current['end']:
            # Merge by extending current
            current['end'] = max(current['end'], row['end'])
            
            if 'score' in current and 'score' in row:
                current['score'] = max(current['score'], row['score'])
        else:
            # Add current and move to next
            if current['end'] - current['start'] >= min_size:
                merged.append(current)
            current = row.to_dict()
    
    # Add the last TAD
    if current['end'] - current['start'] >= min_size:
        merged.append(current)
    
    return pd.DataFrame(merged)
