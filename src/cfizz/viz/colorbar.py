"""
Colorbar utilities for integrated visualizations.

Provides unified colorbar setup with horizontal/vertical orientation support.
"""

import matplotlib.pyplot as plt
from typing import Optional


def setup_colorbar(
    fig,
    sc,
    vmin,
    vmax,
    balance: bool = False,
    color_scale: str = 'linear',
    colorbar_left: Optional[float] = None,
    colorbar_bottom: Optional[float] = None,
    colorbar_width: Optional[float] = None,
    colorbar_height: Optional[float] = None,
    label_fontsize: int = 5,
    tick_fontsize: int = 5,
    label: Optional[str] = None,
    label_orientation: str = 'vertical',
    orientation: str = 'vertical',
):
    """
    Setup colorbar for heatmap visualization.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object
    sc : matplotlib.cm.ScalarMappable
        Scalar mappable object (from imshow)
    vmin, vmax : float
        Colorbar range
    balance : bool, default False
        Whether to use balanced matrix formatting
    colorbar_left : float, optional
        Left position (normalized coordinates 0-1)
    colorbar_bottom : float, optional
        Bottom position (normalized coordinates 0-1)
    colorbar_width : float, optional
        Width (normalized coordinates 0-1)
    colorbar_height : float, optional
        Height (normalized coordinates 0-1)
    label_fontsize : int, default 5
        Label font size
    tick_fontsize : int, default 5
        Tick label font size
    label : str, optional
        Custom label text, None for default
    orientation : str, default 'vertical'
        Colorbar orientation: 'vertical' or 'horizontal'
    """
    # Set Arial font globally for this function
    plt.rcParams['font.family'] = 'Arial'

    # Create colorbar axes
    cax = fig.add_axes([
        colorbar_left,
        colorbar_bottom,
        colorbar_width,
        colorbar_height
    ])

    # Set format string based on balance parameter
    if balance:
        format_str = '%.3g'  # 3 significant digits for balanced matrix
        default_label = "normalized contacts"
    else:
        format_str = '%.3g'  # 3 significant digits for raw matrix
        default_label = "contacts"

    # Add colorbar
    # Note: sc should already have the correct norm from the heatmap
    # Always use vertical orientation for colorbar
    cbar = fig.colorbar(sc, cax=cax, format=format_str)

    # Set ticks at meaningful positions
    import numpy as np
    ticks = np.linspace(vmin, vmax, 3)
    cbar.set_ticks(ticks)

    # Set tick parameters
    cax.tick_params(
        labelsize=tick_fontsize,
        length=1,
        width=0.5
    )

    # Set label - adjust position/rotation based on label_orientation
    label_text = label if label is not None else default_label
    if label_orientation == 'horizontal':
        # Below the axis, horizontal
        cax.text(
            0, -0.2, label_text,
            ha='left', va='top', fontsize=label_fontsize,
            fontname='Arial', rotation=0, transform=cax.transAxes
        )
    else:
        # Vertical: right of axis, rotated 90°
        cax.text(
            0.5, -0.1, label_text,
            ha='center', va='top', fontsize=label_fontsize,
            fontname='Arial', rotation=90, transform=cax.transAxes
        )

    return cbar
