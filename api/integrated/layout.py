"""
Layout calculation for integrated visualizations.


This module provides precise layout calculation for combining Hi-C heatmaps
with genomic tracks in a single figure.
"""

from typing import List, Tuple
from dataclasses import dataclass


@dataclass
class IntegratedLayout:
    """
    Precise layout configuration for integrated heatmap + tracks visualization.

    Attributes
    ----------
    total_width_cm : float
        Total figure width in centimeters (including margins)
    core_width_cm : float
        Core width of heatmap and tracks (genomic coordinate alignment region)
    total_height_cm : float
        Total figure height in centimeters
    heatmap_axes : List[float]
        Matplotlib axes position for single heatmap [left, bottom, width, height]
        (for backward compatibility)
    hic_axes_list : List[List[float]]
        List of matplotlib axes positions for each Hi-C heatmap (when using multiple)
    tracks_axes : List[List[float]]
        List of matplotlib axes positions for each track
    coordinate_label_area : List[float]
        Area reserved for coordinate labels
    gap_cm : float
        Gap between heatmap and tracks in centimeters
    triangle_ratio : float
        Ratio of triangle display (0.0 to 1.0)
    left_margin_cm : float
        Left margin in centimeters
    right_margin_cm : float
        Right margin in centimeters
    n_hics : int
        Number of Hi-C heatmaps (default 1)
    hic_gap_cm : float
        Gap between multiple Hi-C heatmaps in centimeters
    """

    total_width_cm: float
    core_width_cm: float
    total_height_cm: float
    heatmap_axes: List[float]
    hic_axes_list: List[List[float]]
    tracks_axes: List[List[float]]
    coordinate_label_area: List[float]
    gap_cm: float
    triangle_ratio: float
    left_margin_cm: float
    right_margin_cm: float
    n_hics: int
    hic_gap_cm: float
    # Additional attributes (not in dataclass for flexibility)
    n_tracks: int = 0


def calculate_integrated_layout(
    n_tracks: int,
    track_heights_cm: List[float],
    width_cm: float,
    gap_cm: float = 0.3,
    coordinate_label_height_cm: float = 0.3,
    triangle_ratio: float = 1.0,
    left_margin_cm: float = 1.0,
    right_margin_cm: float = 2.0,
    n_hics: int = 1,
    hic_gap_cm: float = 0.1
) -> IntegratedLayout:
    """
    Calculate precise layout for integrated heatmap + tracks visualization.

    This function computes the exact positions and sizes of all elements
    in an integrated figure, ensuring pixel-perfect alignment and
    publication-quality dimensions.

    Parameters
    ----------
    n_tracks : int
        Number of tracks to display
    track_heights_cm : List[float]
        Heights of each track in centimeters (length should equal n_tracks)
    width_cm : float
        Core width of the visualization components. For tracks, this is the actual
        drawing width. For the 45-degree rotated heatmap, this equals the diagonal
        length of the heatmap square. The heatmap's axes width will be width_cm/√2.
        Heatmap height is automatically calculated as width_cm/√2 to maintain
        a square aspect ratio.
        This ensures perfect genomic coordinate alignment between heatmap and tracks.
    gap_cm : float, default 0.3
        Gap between heatmap and tracks in centimeters
    coordinate_label_height_cm : float, default 0.3
        Height reserved for coordinate labels at the bottom
    triangle_ratio : float, default 1.0
        Ratio of triangle display (0.0 to 1.0). Controls how much of the
        upper triangle to show in the heatmap.
    left_margin_cm : float, default 1.0
        Left margin outside the width_cm core region (for y-axis labels/titles)
    right_margin_cm : float, default 2.0
        Right margin outside the width_cm core region
    n_hics : int, default 1
        Number of Hi-C heatmaps to display (supports multiple for comparison)
    hic_gap_cm : float, default 0.1
        Gap between multiple Hi-C heatmaps in centimeters

    Returns
    -------
    IntegratedLayout
        Complete layout configuration with precise positions
    """

    # Validate inputs
    if n_tracks < 0:
        raise ValueError("Number of tracks must be non-negative")
    if len(track_heights_cm) != n_tracks:
        raise ValueError(
            f"track_heights_cm length ({len(track_heights_cm)}) "
            f"must equal n_tracks ({n_tracks})"
        )
    if width_cm <= 0:
        raise ValueError("width_cm must be positive")
    if triangle_ratio < 0.0 or triangle_ratio > 1.0:
        raise ValueError("triangle_ratio must be between 0.0 and 1.0")
    if left_margin_cm < 0:
        raise ValueError("left_margin_cm must be non-negative")
    if right_margin_cm < 0:
        raise ValueError("right_margin_cm must be non-negative")
    if n_hics < 1:
        raise ValueError("n_hics must be at least 1")
    if hic_gap_cm < 0:
        raise ValueError("hic_gap_cm must be non-negative")

    # Calculate total dimensions
    # width_cm is the core width for heatmap and tracks
    core_width_cm = width_cm
    total_width_cm = left_margin_cm + core_width_cm + right_margin_cm

    # Calculate heatmap dimensions
    # Single heatmap height (each)
    # **KEY FIX**: For 45-degree rotation, the diagonal = width_cm, so height = width_cm / sqrt(2) ≈ width_cm / 2
    individual_heatmap_height_cm = core_width_cm / 2 * triangle_ratio

    # Total heatmap area height (all heatmaps + gaps between them)
    total_heatmap_height_cm = n_hics * individual_heatmap_height_cm + (n_hics - 1) * hic_gap_cm

    tracks_height_cm = sum(track_heights_cm)
    total_height_cm = total_heatmap_height_cm + gap_cm + tracks_height_cm + coordinate_label_height_cm

    # Calculate axes positions
    heatmap_left = left_margin_cm / total_width_cm
    heatmap_width = core_width_cm / total_width_cm

    # Calculate position of the bottom-most heatmap
    heatmap_bottom = (tracks_height_cm + coordinate_label_height_cm + gap_cm) / total_height_cm

    # Calculate axes for each Hi-C heatmap (stacked vertically)
    hic_axes_list = []
    for i in range(n_hics):
        # Each heatmap is positioned above the previous one
        individual_height_ratio = individual_heatmap_height_cm / total_height_cm
        bottom_position = heatmap_bottom + i * (individual_height_ratio + hic_gap_cm / total_height_cm)

        hic_axes = [
            heatmap_left,                       # left edge (with left margin)
            bottom_position,                    # bottom position
            heatmap_width,                      # width (core region)
            individual_height_ratio             # height
        ]
        hic_axes_list.append(hic_axes)

    # For backward compatibility: heatmap_axes is the first (or only) heatmap
    heatmap_axes = hic_axes_list[0] if n_hics > 0 else [0, 0, 0, 0]

    # Calculate tracks axes positions (bottom portion, stacked)
    # Each track spans the core width region (NOT divided by sqrt(2))
    tracks_width_normalized = core_width_cm / total_width_cm
    tracks_axes = []
    current_bottom = coordinate_label_height_cm / total_height_cm

    for track_height_cm in track_heights_cm:
        track_axes = [
            heatmap_left,                       # left edge (with left margin)
            current_bottom,                     # bottom position
            tracks_width_normalized,            # width (full core region, not /sqrt(2))
            track_height_cm / total_height_cm   # height
        ]
        tracks_axes.append(track_axes)
        current_bottom += track_height_cm / total_height_cm

    # Coordinate label area (bottom strip)
    coordinate_label_area = [
        heatmap_left,                           # left edge (with left margin)
        0.0,                                    # bottom edge
        tracks_width_normalized,                 # width (full core region, not /sqrt(2))
        coordinate_label_height_cm / total_height_cm  # height
    ]

    layout = IntegratedLayout(
        total_width_cm=total_width_cm,
        core_width_cm=core_width_cm,
        total_height_cm=total_height_cm,
        heatmap_axes=heatmap_axes,
        hic_axes_list=hic_axes_list,
        tracks_axes=tracks_axes,
        coordinate_label_area=coordinate_label_area,
        gap_cm=gap_cm,
        triangle_ratio=triangle_ratio,
        left_margin_cm=left_margin_cm,
        right_margin_cm=right_margin_cm,
        n_hics=n_hics,
        hic_gap_cm=hic_gap_cm
    )
    # Additional attribute (not in dataclass)
    layout.n_tracks = n_tracks
    
    return layout


def validate_layout(layout: IntegratedLayout) -> Tuple[bool, List[str]]:
    """
    Validate a layout configuration.

    Parameters
    ----------
    layout : IntegratedLayout
        Layout to validate

    Returns
    -------
    Tuple[bool, List[str]]
        (is_valid, list_of_errors)
    """
    errors = []

    # Check total dimensions
    if layout.total_width_cm <= 0:
        errors.append("total_width_cm must be positive")

    if layout.total_height_cm <= 0:
        errors.append("total_height_cm must be positive")

    # Check Hi-C axes list
    if len(layout.hic_axes_list) != layout.n_hics:
        errors.append(f"hic_axes_list length ({len(layout.hic_axes_list)}) must equal n_hics ({layout.n_hics})")

    for i, axes in enumerate(layout.hic_axes_list):
        if len(axes) != 4:
            errors.append(f"hic_axes_list[{i}] must have 4 elements [left, bottom, width, height]")

        left, bottom, width, height = axes
        if left < 0 or left > 1:
            errors.append(f"hic_axes_list[{i}]: left must be between 0 and 1, got {left}")
        if bottom < 0 or bottom > 1:
            errors.append(f"hic_axes_list[{i}]: bottom must be between 0 and 1, got {bottom}")
        if width <= 0 or width > 1:
            errors.append(f"hic_axes_list[{i}]: width must be between 0 and 1, got {width}")
        if height <= 0 or height > 1:
            errors.append(f"hic_axes_list[{i}]: height must be between 0 and 1, got {height}")

    # Check tracks axes positions
    for i, axes in enumerate(layout.tracks_axes):
        if len(axes) != 4:
            errors.append(f"tracks_axes[{i}] must have 4 elements [left, bottom, width, height]")

        left, bottom, width, height = axes
        if left < 0 or left > 1:
            errors.append(f"tracks_axes[{i}]: left must be between 0 and 1, got {left}")
        if bottom < 0 or bottom > 1:
            errors.append(f"tracks_axes[{i}]: bottom must be between 0 and 1, got {bottom}")
        if width <= 0 or width > 1:
            errors.append(f"tracks_axes[{i}]: width must be between 0 and 1, got {width}")
        if height <= 0 or height > 1:
            errors.append(f"tracks_axes[{i}]: height must be between 0 and 1, got {height}")

    # Check heatmap axes (for backward compatibility)
    if len(layout.heatmap_axes) != 4:
        errors.append("heatmap_axes must have 4 elements [left, bottom, width, height]")

    return len(errors) == 0, errors
