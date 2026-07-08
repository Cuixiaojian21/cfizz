"""
Genome track visualization module.

Functions for plotting genomic tracks alongside Hi-C heatmaps:
- Gene tracks
- Expression tracks
- Custom bed tracks
- Multi-track integration

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
from typing import Optional, Tuple, List, Dict, Any
import warnings

# Suppress warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


def plot_tracks(
    chrom: str,
    start: int,
    end: int,
    tracks: List[Dict[str, Any]],
    track_height: float = 0.5,
    figsize: Tuple[float, float] = (12, 8),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    """
    Plot multiple genomic tracks.
    
    This function creates a multi-track visualization showing
    various genomic data types (genes, peaks, etc.) aligned
    along the genomic coordinate.
    
    Parameters
    ----------
    chrom : str
        Chromosome name
    start : int
        Start position (bp)
    end : int
        End position (bp)
    tracks : list
        List of track dictionaries with keys:
        - 'type': 'gene', 'bed', 'bedgraph' or 'interval'
        - 'data': Data for the track
        - 'name': Track name (optional)
        - 'color': Track color (optional)
    track_height : float
        Height of each track in inches
    figsize : tuple
        Figure size in inches
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    dpi : int
        Resolution for saved figure
        
    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
        
    Examples
    --------
    >>> tracks = [
    ...     {'type': 'gene', 'data': gene_df, 'name': 'Genes'},
    ...     {'type': 'bedgraph', 'data': signal_df, 'name': 'Signal'}
    ... ]
    >>> fig = plot_tracks("chr1", 0, 10000000, tracks)
    """
    n_tracks = len(tracks)
    fig_height = n_tracks * track_height + 1
    
    fig, axes = plt.subplots(n_tracks, 1, figsize=(figsize[0], fig_height), sharex=True)
    
    if n_tracks == 1:
        axes = [axes]
    
    # Add title
    if title:
        fig.suptitle(title, fontsize=12, y=0.98)
    
    for i, track in enumerate(tracks):
        ax = axes[i]
        track_type = track.get('type', 'unknown')
        track_name = track.get('name', f'Track {i+1}')
        track_color = track.get('color', '#333333')
        
        if track_type == 'gene':
            _plot_gene_track(ax, track['data'], start, end, track_color)
        elif track_type == 'bed':
            _plot_bed_track(ax, track['data'], start, end, track_color)
        elif track_type == 'bedgraph':
            _plot_bedgraph_track(ax, track['data'], start, end, track_color)
        elif track_type == 'interval':
            _plot_interval_track(ax, track['data'], start, end, track_color)
        
        ax.set_ylabel(track_name, fontsize=8, rotation=0, ha='right', va='center')
        ax.set_xlim(start, end)
        
        if i < n_tracks - 1:
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])
    
    axes[-1].set_xlabel(f'{chrom} position (bp)', fontsize=10)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig


def _plot_gene_track(
    ax: plt.Axes,
    gene_data: pd.DataFrame,
    start: int,
    end: int,
    color: str = '#333333'
):
    """Plot gene track from BED12 or similar format."""
    for _, gene in gene_data.iterrows():
        # Get gene coordinates
        if 'thickStart' in gene:
            cds_start = gene['thickStart']
            cds_end = gene['thickEnd']
        else:
            cds_start = gene['start']
            cds_end = gene['end']
        
        gene_start = gene['start']
        gene_end = gene['end']
        
        # Draw gene body
        ax.plot([gene_start, gene_end], [0.5, 0.5], color=color, linewidth=2)
        
        # Draw CDS regions (thick regions)
        ax.plot([cds_start, cds_end], [0.5, 0.5], color=color, linewidth=5)
        
        # Draw arrows for direction
        if 'strand' in gene:
            strand = gene['strand']
            mid = (gene_start + gene_end) / 2
            
            if strand == '+':
                ax.annotate('', xy=(mid + 0.01 * (end - start), 0.5), 
                           xytext=(mid, 0.5),
                           arrowprops=dict(arrowstyle='->', color=color, lw=1))
            else:
                ax.annotate('', xy=(mid - 0.01 * (end - start), 0.5), 
                           xytext=(mid, 0.5),
                           arrowprops=dict(arrowstyle='->', color=color, lw=1))
    
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)


def _plot_bed_track(
    ax: plt.Axes,
    bed_data: pd.DataFrame,
    start: int,
    end: int,
    color: str = '#333333'
):
    """Plot BED feature track."""
    y_pos = 0.5
    for _, feature in bed_data.iterrows():
        feature_start = feature['start']
        feature_end = feature['end']
        ax.plot([feature_start, feature_end], [y_pos, y_pos], 
                color=color, linewidth=3)
    
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)


def _plot_bedgraph_track(
    ax: plt.Axes,
    bedgraph_data: pd.DataFrame,
    start: int,
    end: int,
    color: str = '#333333'
):
    """Plot bedGraph track (continuous signal)."""
    positions = []
    values = []
    
    for _, row in bedgraph_data.iterrows():
        positions.extend([row['start'], row['end']])
        values.extend([row['value'], row['value']])
    
    if positions:
        ax.fill_between(positions, values, alpha=0.5, color=color)
        ax.plot(positions, values, color=color, linewidth=0.5)
    
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)


def _plot_interval_track(
    ax: plt.Axes,
    interval_data: pd.DataFrame,
    start: int,
    end: int,
    color: str = '#333333'
):
    """Plot interval track."""
    y_pos = 0.5
    for _, interval in interval_data.iterrows():
        interval_start = interval['start']
        interval_end = interval['end']
        
        # Draw rectangle
        rect = Rectangle(
            (interval_start, y_pos - 0.2),
            interval_end - interval_start,
            0.4,
            facecolor=color,
            edgecolor='none',
            alpha=0.7
        )
        ax.add_patch(rect)
    
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)


def add_tracks_to_heatmap(
    fig: plt.Figure,
    heatmap_ax: plt.Axes,
    tracks: List[Dict[str, Any]],
    position: str = 'top',
    track_height: float = 0.5
) -> List[plt.Axes]:
    """
    Add genomic tracks to an existing heatmap figure.
    
    Parameters
    ----------
    fig : plt.Figure
        Existing figure
    heatmap_ax : plt.Axes
        Existing heatmap axes
    tracks : list
        List of track dictionaries
    position : str
        'top' or 'bottom' - where to place tracks
    track_height : float
        Height of each track (in figure units)
        
    Returns
    -------
    track_axes : list
        List of created track axes
    """
    # Get heatmap position
    heatmap_pos = heatmap_ax.get_position()
    
    if position == 'top':
        # Create tracks above heatmap
        track_bottom = heatmap_pos.y1 + 0.02
        track_height_norm = track_height / fig.get_figheight()
    else:
        # Create tracks below heatmap
        track_height_norm = track_height / fig.get_figheight()
        track_bottom = heatmap_pos.y0 - track_height_norm - 0.02
    
    track_axes = []
    current_bottom = track_bottom
    
    for track in tracks:
        ax = fig.add_axes([
            heatmap_pos.x0,
            current_bottom,
            heatmap_pos.width,
            track_height_norm
        ])
        track_axes.append(ax)
        current_bottom += track_height_norm + 0.01
    
    return track_axes
