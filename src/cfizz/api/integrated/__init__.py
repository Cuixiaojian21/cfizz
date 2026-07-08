"""
Integrated visualization API.

Provides high-level interface for creating integrated Hi-C + tracks figures.
"""

from cfizz.api.integrated.layout import (
    IntegratedLayout,
    calculate_integrated_layout,
    validate_layout,
)

from cfizz.api.integrated.heatmap_tracks import (
    HeatmapTracks,
    GenomeRange,
    _plot_rotated_heatmap,
    _add_coordinate_labels,
)

from cfizz.api.integrated.quick_plot import quick_plot_integrated

__all__ = [
    # Classes
    'IntegratedLayout',
    'HeatmapTracks',
    'GenomeRange',
    # Functions
    'calculate_integrated_layout',
    'validate_layout',
    'quick_plot_integrated',
    # Internal (for advanced usage)
    '_plot_rotated_heatmap',
    '_add_coordinate_labels',
]