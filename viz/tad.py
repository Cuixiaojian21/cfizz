"""
TAD (Topologically Associating Domain) visualization module.

Functions for plotting TAD analysis results:
- Insulation score tracks
- TAD boundary visualization
- TAD pileup

Real implementations:
- plot_heatmap_with_tad_boundaries() -> cfizz/viz/heatmap_tad_ext.py
- plot_multi_tad_boundary_pileup() -> cfizz/viz/pileup.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
import warnings
import os

# Suppress warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# Visualization functions
# =============================================================================

def plot_insulation_track(
    insulation_data: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    window_size: int = 100000,
    color: str = "#1f77b4",
    color_boundary: str = "#d62728",
    figsize: Tuple[float, float] = (10, 3),
    show_boundaries: bool = True,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Plot insulation score track with TAD boundaries.
    
    Parameters
    ----------
    insulation_data : pd.DataFrame
        DataFrame with insulation scores
    chrom : str
        Chromosome name
    start, end : int
        Genomic region
    window_size : int
        Window size (bp)
    color : str
        Line color
    color_boundary : str
        Boundary marker color
    figsize : tuple
        Figure size
    show_boundaries : bool
        Show boundary markers
    title : str
        Plot title
    save_path : str
        Output file path
    dpi : int
        Output DPI
        
    Returns
    -------
    fig : plt.Figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Filter region data
    region_data = insulation_data[
        (insulation_data['chrom'] == chrom) &
        (insulation_data['start'] >= start) &
        (insulation_data['end'] <= end)
    ].copy()
    
    if len(region_data) == 0:
        raise ValueError(f"No data found for {chrom}:{start}-{end}")
    
    # Get positions and values
    positions = (region_data['start'] + region_data['end']) / 2
    values = region_data['insulation_score'].values
    
    # Fill NaN
    values = np.nan_to_num(values, nan=0)
    
    # Plot
    ax.plot(positions / 1e6, values, color=color, linewidth=1)
    ax.fill_between(positions / 1e6, values, alpha=0.3, color=color)
    
    # Mark boundaries
    if show_boundaries:
        strength_col = f'boundary_strength_{window_size}'
        is_boundary_col = f'is_boundary_{window_size}'
        
        if is_boundary_col in region_data.columns:
            boundaries = region_data[region_data[is_boundary_col]]
            ax.scatter(
                (boundaries['start'] + boundaries['end']) / 2 / 1e6,
                boundaries['insulation_score'],
                color=color_boundary,
                s=50,
                zorder=5,
                marker='|',
                linewidth=2
            )
    
    ax.set_xlabel(f'{chrom} position (Mb)')
    ax.set_ylabel('Insulation Score')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title(f'Insulation Score (window={window_size//1000}kb)')
    
    ax.grid(True, alpha=0.3)
    ax.set_xlim(start / 1e6, end / 1e6)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


def plot_tad_heatmap_with_boundaries(
    matrix: np.ndarray,
    boundaries: pd.DataFrame,
    resolution: int,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    vmin: float = 0,
    vmax: float = 1,
    cmap: str = "Reds",
    color_boundary: str = "#d62728",
    boundary_width: int = 2,
    figsize: Tuple[float, float] = (8, 8),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Plot heatmap with TAD boundaries marked.
    
    Parameters
    ----------
    matrix : np.ndarray
        Contact matrix
    boundaries : pd.DataFrame
        TAD boundaries with chrom, start, end columns
    resolution : int
        Bin resolution (bp)
    chrom, start, end : str, int
        Region info for filtering boundaries
    vmin, vmax : float
        Color range
    cmap : str
        Colormap
    color_boundary : str
        Boundary line color
    boundary_width : int
        Line width
    figsize : tuple
        Figure size
    title : str
        Plot title
    save_path : str
        Output file
    dpi : int
        Output DPI
        
    Returns
    -------
    fig : plt.Figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', 
                   interpolation='none', vmin=vmin, vmax=vmax)
    
    # Mark boundaries
    if chrom is not None and start is not None:
        region_bounds = boundaries[
            (boundaries['chrom'] == chrom) &
            (boundaries['start'] >= start) &
            (boundaries['end'] <= end)
        ].copy()
        
        n = matrix.shape[0]
        for _, row in region_bounds.iterrows():
            idx = (row['start'] - start) // resolution
            if 0 <= idx < n:
                for i in range(boundary_width):
                    if idx + i < n:
                        ax.axhline(y=idx + i, color=color_boundary, 
                                   linewidth=0.5, alpha=0.8)
                        ax.axvline(x=idx + i, color=color_boundary, 
                                   linewidth=0.5, alpha=0.8)
    
    ax.set_xlabel('Bin')
    ax.set_ylabel('Bin')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title('TAD Boundaries')
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Contact frequency')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


# =============================================================================
# Utilities
# =============================================================================

def get_tad_boundary_positions(
    boundaries: pd.DataFrame,
    chrom: str,
    start: int,
    end: int
) -> List[int]:
    """
    Get TAD boundary positions for a region.
    
    Parameters
    ----------
    boundaries : pd.DataFrame
        TAD boundaries DataFrame
    chrom : str
        Chromosome
    start, end : int
        Region boundaries
        
    Returns
    -------
    list
        List of boundary positions (bin indices)
    """
    region_bounds = boundaries[
        (boundaries['chrom'] == chrom) &
        (boundaries['start'] >= start) &
        (boundaries['end'] <= end)
    ]
    
    positions = []
    for _, row in region_bounds.iterrows():
        pos = (row['start'] + row['end']) // 2
        positions.append(pos)
    
    return sorted(positions)


def create_tad_summary(
    tads: pd.DataFrame,
    statistics: Dict[str, Any]
) -> str:
    """
    Create a summary string for TAD analysis.
    
    Parameters
    ----------
    tads : pd.DataFrame
        TAD DataFrame
    statistics : dict
        Statistics dictionary
        
    Returns
    -------
    str
        Summary string
    """
    summary = []
    summary.append("=" * 60)
    summary.append("TAD Analysis Summary")
    summary.append("=" * 60)
    summary.append(f"Total TADs: {len(tads)}")
    
    if len(tads) > 0:
        summary.append(f"\nSize Statistics:")
        summary.append(f"  Mean: {statistics.get('mean_size', 0)/1000:.1f} kb")
        summary.append(f"  Median: {statistics.get('median_size', 0)/1000:.1f} kb")
        summary.append(f"  Min: {statistics.get('min_size', 0)/1000:.1f} kb")
        summary.append(f"  Max: {statistics.get('max_size', 0)/1000:.1f} kb")
        summary.append(f"  Std: {statistics.get('std_size', 0)/1000:.1f} kb")
    
    summary.append("=" * 60)
    return "\n".join(summary)