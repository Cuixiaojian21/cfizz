"""
Chromatin loops I/O and utilities.

This module provides functions to read/write chromatin loop calls
and perform basic loop operations.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any


def read_loops_bedpe(bedpe_path: str) -> pd.DataFrame:
    """
    Read chromatin loops from BEDPE file.
    
    BEDPE format: chrom1 start1 end1 chrom2 start2 end2 [score]
    
    Parameters
    ----------
    bedpe_path : str
        Path to BEDPE file
        
    Returns
    -------
    loops_df : pd.DataFrame
        DataFrame with columns [chrom1, start1, end1, chrom2, start2, end2, score]
        
    Examples
    --------
    >>> loops = read_loops_bedpe("loops.bedpe")
    >>> print(f"Found {len(loops)} loops")
    """
    columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    
    try:
        df = pd.read_csv(
            bedpe_path, 
            sep='\t', 
            header=None, 
            names=columns + ['score'],
            usecols=range(7)
        )
    except ValueError:
        df = pd.read_csv(bedpe_path, sep='\t', header=None, names=columns)
        df['score'] = 1.0
    
    return df


def write_loops_bedpe(loops_df: pd.DataFrame, output_path: str) -> None:
    """
    Write chromatin loops to BEDPE file.
    
    Parameters
    ----------
    loops_df : pd.DataFrame
        DataFrame with loop coordinates
    output_path : str
        Output file path
        
    Examples
    --------
    >>> write_loops_bedpe(loops_df, "output_loops.bedpe")
    """
    required_cols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    
    for col in required_cols:
        if col not in loops_df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    output_df = loops_df[required_cols].copy()
    
    if 'score' in loops_df.columns:
        output_df['score'] = loops_df['score']
    else:
        output_df['score'] = 1.0
    
    output_df.to_csv(output_path, sep='\t', header=False, index=False)


class LoopSet:
    """
    Container for chromatin loop calls.
    
    Provides utilities for filtering, merging, and analyzing loop sets.
    
    Examples
    --------
    >>> loop_set = LoopSet.from_bedpe("loops.bedpe")
    >>> filtered = loop_set.filter_by_score(threshold=0.5)
    >>> filtered.save("filtered_loops.bedpe")
    """
    
    def __init__(self, loops_df: pd.DataFrame):
        """
        Initialize LoopSet.
        
        Parameters
        ----------
        loops_df : pd.DataFrame
            DataFrame with loop coordinates
        """
        self.loops = loops_df.copy()
    
    @classmethod
    def from_bedpe(cls, path: str) -> 'LoopSet':
        """Create LoopSet from BEDPE file."""
        return cls(read_loops_bedpe(path))
    
    def save(self, path: str) -> None:
        """Save to BEDPE file."""
        write_loops_bedpe(self.loops, path)
    
    def filter_by_score(self, threshold: float) -> 'LoopSet':
        """
        Filter loops by score.
        
        Parameters
        ----------
        threshold : float
            Minimum score
            
        Returns
        -------
        filtered : LoopSet
        """
        if 'score' not in self.loops.columns:
            return self
        
        filtered = self.loops[self.loops['score'] >= threshold].copy()
        return LoopSet(filtered)
    
    def filter_by_distance(
        self, 
        min_distance: int = None, 
        max_distance: int = None
    ) -> 'LoopSet':
        """
        Filter loops by genomic distance.
        
        Parameters
        ----------
        min_distance : int, optional
            Minimum loop distance (bp)
        max_distance : int, optional
            Maximum loop distance (bp)
            
        Returns
        -------
        filtered : LoopSet
        """
        distances = self.loops['start2'] - self.loops['end1']
        
        filtered = self.loops.copy()
        
        if min_distance is not None:
            filtered = filtered[distances >= min_distance]
            distances = distances[distances >= min_distance]
        
        if max_distance is not None:
            filtered = filtered[distances <= max_distance]
        
        return LoopSet(filtered)
    
    def merge_loops(self, merge_distance: int = 100000) -> 'LoopSet':
        """
        Merge nearby loops.
        
        Parameters
        ----------
        merge_distance : int
            Maximum distance to merge (bp)
            
        Returns
        -------
        merged : LoopSet
        """
        if len(self.loops) == 0:
            return self
        
        # Sort by position
        sorted_loops = self.loops.sort_values(['chrom1', 'start1']).reset_index(drop=True)
        
        merged = []
        current = sorted_loops.iloc[0].to_dict()
        
        for i in range(1, len(sorted_loops)):
            row = sorted_loops.iloc[i]
            
            # Check if within merge distance
            if (row['chrom1'] == current['chrom1'] and
                abs(row['start1'] - current['start1']) <= merge_distance and
                abs(row['start2'] - current['start2']) <= merge_distance):
                
                # Merge by averaging
                current['start1'] = (current['start1'] + row['start1']) // 2
                current['end1'] = (current['end1'] + row['end1']) // 2
                current['start2'] = (current['start2'] + row['start2']) // 2
                current['end2'] = (current['end2'] + row['end2']) // 2
                
                if 'score' in current and 'score' in row:
                    current['score'] = max(current['score'], row['score'])
            else:
                merged.append(current)
                current = row.to_dict()
        
        merged.append(current)
        
        return LoopSet(pd.DataFrame(merged))
    
    def get_anchor_positions(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get anchor positions as arrays.
        
        Returns
        -------
        anchor1 : np.ndarray
            First anchor positions (midpoints)
        anchor2 : np.ndarray
            Second anchor positions (midpoints)
        """
        anchor1 = (self.loops['start1'] + self.loops['end1']) / 2
        anchor2 = (self.loops['start2'] + self.loops['end2']) / 2
        return anchor1.values, anchor2.values
    
    def get_distances(self) -> np.ndarray:
        """Get loop distances."""
        return (self.loops['start2'] - self.loops['end1']).values
    
    def __len__(self) -> int:
        return len(self.loops)
    
    def __repr__(self) -> str:
        return f"LoopSet(n={len(self.loops)})"


# =============================================================================
# =============================================================================

def read_loops(
    file_path: str,
    chrom: str = None,
    start: int = None,
    end: int = None
) -> pd.DataFrame:
    """
    读取loops文件
    
    
    Args:
        file_path: loops文件路径
        chrom: 染色体号，如果为None则读取所有染色体
        start: 起始位置（bp），如果为None则从染色体起始位置开始
        end: 结束位置（bp），如果为None则到染色体结束位置
        
    Returns:
        pd.DataFrame: 包含loops信息的DataFrame，列包括：
            - chrom1: 第一个锚点的染色体
            - start1: 第一个锚点的起始位置
            - end1: 第一个锚点的结束位置
            - chrom2: 第二个锚点的染色体
            - start2: 第二个锚点的起始位置
            - end2: 第二个锚点的结束位置
    """
    # 健壮性处理：file_path为None或空字符串时直接返回空DataFrame
    if file_path is None or (isinstance(file_path, str) and not file_path.strip()):
        return pd.DataFrame(columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
    # 读取loops文件，只读取前6列
    df = pd.read_csv(file_path, sep='\t', header=None, usecols=range(6), dtype={0: str, 3: str})
    
    # 设置列名
    df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    
    # 如果指定了染色体和位置范围，进行过滤
    if chrom is not None:
        df = df[
            ((df['chrom1'] == chrom) & (df['chrom2'] == chrom)) &
            ((df['start1'] >= start) & (df['start1'] <= end) |
             (df['start2'] >= start) & (df['start2'] <= end))
        ]
    
    return df


def loops_to_anchor_pairs(loops_df: pd.DataFrame) -> List[Tuple[int, int]]:
    """
    Convert loops to list of anchor pairs.
    
    Parameters
    ----------
    loops_df : pd.DataFrame
        DataFrame with loops
        
    Returns
    -------
    pairs : list
        List of (anchor1_pos, anchor2_pos) tuples
    """
    pairs = []
    for _, row in loops_df.iterrows():
        anchor1 = (row['start1'] + row['end1']) // 2
        anchor2 = (row['start2'] + row['end2']) // 2
        pairs.append((anchor1, anchor2))
    return pairs
