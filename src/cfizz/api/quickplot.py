"""
Quick plotting API.

High-level functions for common visualization tasks.
"""

from typing import Optional, List, Dict


def quick_plot_integrated(
    hic_file: str,
    region: str,
    output: str,
    tracks: Optional[List[str]] = None,
    resolution: int = 10000,
    cmap: str = "Reds",
    dpi: int = 300
) -> None:
    """
    Quick plot function for integrated Hi-C heatmap and tracks.
    
    This is a convenience wrapper that handles layout calculation and visualization
    in a single function call.
    
    Parameters
    ----------
    hic_file : str
        Path to Hi-C file (.mcool)
    region : str
        Genomic region (e.g., "chr1:1000000-2000000")
    output : str
        Output file path (without extension)
    tracks : list, optional
        List of track file paths
    resolution : int
        Resolution in bp
    cmap : str
        Colormap name
    dpi : int
        Output resolution
        
    Examples
    --------
    >>> from cfizz.api import quick_plot_integrated
    >>> quick_plot_integrated(
    ...     hic_file="sample.mcool",
    ...     region="chr1:1000000-2000000",
    ...     output="output/my_region",
    ...     tracks=["CTCF.bw", "RNAseq.bw"]
    ... )
    """
    # Placeholder - to be implemented
    raise NotImplementedError("Use cfizz.api.integrated.quick_plot_integrated instead")
