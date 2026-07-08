#!/usr/bin/env python3
"""
Quick Plot Convenience Functions for Integrated Visualization.


This module provides convenient one-stop functions for creating integrated
Hi-C heatmap and genomic tracks visualizations with minimal configuration.
"""

from typing import List, Optional
from cfizz.api.integrated.heatmap_tracks import HeatmapTracks, GenomeRange
from cfizz.api.integrated.layout import calculate_integrated_layout


def quick_plot_integrated(
    hics: Optional[List[dict]] = None,
    hic: Optional[dict] = None,
    hic_file: Optional[str] = None,
    tracks: Optional[List[dict]] = None,
    track_files: Optional[List[str]] = None,
    region: Optional[GenomeRange] = None,
    output: Optional[str] = None,
    n_tracks: int = 2,
    track_heights_cm: Optional[List[float]] = None,
    width_cm: float = 5.0,
    gap_cm: float = 0.1,
    left_margin_cm: float = 1.0,
    right_margin_cm: float = 2.0,
    dpi: int = 300,
    triangle_ratio: float = 1.0,
    hic_cmap: str = 'Reds',
    hic_color_scale: str = 'linear',
    balance: bool = True,
    resolution: int = 10000,
    track_colors: Optional[List[str]] = None,
    track_names: Optional[List[str]] = None,
    min_value: Optional[List[float]] = None,
    max_value: Optional[List[float]] = None,   # T-7.4 新加(跟 min_value 对称)
    **kwargs
) -> None:
    """
    Quick plot function for integrated Hi-C heatmap and tracks visualization.

    Supports three input formats:
    1. Multiple Hi-C: hics list + tracks list
    2. Single Hi-C: hic dict + tracks list
    3. Traditional: separate parameters

    Parameters
    ----------
    hics : List[dict], optional
        Multiple Hi-C configurations as dict list.
        Example: [{'file': 'sample1.mcool', 'cmap': 'Reds', 'name': 'Normal'}, ...]
        Supports: file, triangle_ratio, cmap, color_scale, balance, resolution, name,
        insulation_path, window_size, boundary_cmap, boundary_alpha
    hic : dict, optional
        Single Hi-C configuration dict (backward compatible).
    hic_file : str, optional
        Traditional format: single Hi-C file path.
    tracks : List[dict], optional
        Tracks configuration dict list.
        Example: [{'file': 'track1.bw', 'color': 'blue', 'name': 'CTCF', 'height_cm': 1.0}, ...]
    track_files : List[str], optional
        Traditional format: track file paths list.
    region : GenomeRange
        Genomic region to visualize.
    output : str
        Output file path (without extension).
    n_tracks : int, default 2
        Number of tracks to display.
    track_heights_cm : List[float], optional
        Height of each track in centimeters.
    width_cm : float, default 5.0
        Core width in centimeters.
    gap_cm : float, default 0.1
        Gap between heatmap and tracks.
    left_margin_cm : float, default 1.0
        Left margin in centimeters.
    right_margin_cm : float, default 2.0
        Right margin in centimeters.
    dpi : int, default 300
        Output resolution.
    triangle_ratio : float, default 1.0
        Triangle display ratio (0.0-1.0).
    hic_cmap : str, default 'Reds'
        Hi-C colormap.
    hic_color_scale : str, default 'linear'
        Hi-C color scale ('linear' or 'log').
    balance : bool, default True
        Whether to balance Hi-C data.
    resolution : int, default 10000
        Hi-C resolution in base pairs.
    track_colors : List[str], optional
        Colors for each track.
    track_names : List[str], optional
        Names for each track.
    min_value : List[float], optional
        Minimum y-axis values.
    """
    # 1. Validate required parameters
    if region is None:
        raise ValueError("region is required")
    if output is None:
        raise ValueError("output is required")

    # 2. Process Hi-C parameters
    hics_files = None
    extracted_hic_cmaps = []
    extracted_triangle_ratios = []
    extracted_hic_color_scales = []
    extracted_balances = []
    extracted_resolutions = []
    extracted_flip_vertical = []
    extracted_hic_names = []
    extracted_insulation_paths = []
    extracted_window_sizes = []
    extracted_boundary_cmaps = []
    extracted_boundary_alphas = []
    extracted_loops_paths = []
    extracted_loop_colors = []
    extracted_loop_alphas = []
    extracted_loop_sizes = []

    if hics is not None:
        # New format: hics dict list
        n_hics = len(hics)
        hics_files = []
        for hic_dict in hics:
            if 'file' not in hic_dict:
                raise ValueError("Each Hi-C dict must have 'file' parameter")
            hics_files.append(hic_dict['file'])
            extracted_hic_cmaps.append(hic_dict.get('cmap', 'Reds'))
            extracted_triangle_ratios.append(hic_dict.get('triangle_ratio', 1.0))
            extracted_hic_color_scales.append(hic_dict.get('color_scale', 'linear'))
            extracted_balances.append(hic_dict.get('balance', True))
            extracted_resolutions.append(hic_dict.get('resolution', 10000))
            extracted_flip_vertical.append(hic_dict.get('flip_vertical', False))
            extracted_hic_names.append(hic_dict.get('name', None))
            extracted_insulation_paths.append(hic_dict.get('insulation_path', None))
            extracted_window_sizes.append(hic_dict.get('window_size', 100000))
            extracted_boundary_cmaps.append(hic_dict.get('boundary_cmap', 'Blues_r'))
            extracted_boundary_alphas.append(hic_dict.get('boundary_alpha', 0.6))
            extracted_loops_paths.append(hic_dict.get('loops_path', None))
            extracted_loop_colors.append(hic_dict.get('loop_color', 'blue'))
            extracted_loop_alphas.append(hic_dict.get('loop_alpha', 0.6))
            extracted_loop_sizes.append(hic_dict.get('loop_size', 1.0))

    elif hic is not None:
        # Single Hi-C dict format (backward compatible)
        if 'file' not in hic:
            raise ValueError("hic dict must have 'file' parameter")
        hics_files = [hic['file']]
        extracted_hic_cmaps = [hic.get('cmap', 'Reds')]
        extracted_triangle_ratios = [hic.get('triangle_ratio', 1.0)]
        extracted_hic_color_scales = [hic.get('color_scale', 'linear')]
        extracted_balances = [hic.get('balance', True)]
        extracted_resolutions = [hic.get('resolution', 10000)]
        extracted_flip_vertical = [hic.get('flip_vertical', False)]
        extracted_hic_names = [hic.get('name', None)]
        extracted_insulation_paths = [hic.get('insulation_path', None)]
        extracted_window_sizes = [hic.get('window_size', 100000)]
        extracted_boundary_cmaps = [hic.get('boundary_cmap', 'Blues_r')]
        extracted_boundary_alphas = [hic.get('boundary_alpha', 0.6)]
        extracted_loops_paths = [hic.get('loops_path', None)]
        extracted_loop_colors = [hic.get('loop_color', 'blue')]
        extracted_loop_alphas = [hic.get('loop_alpha', 0.6)]
        extracted_loop_sizes = [hic.get('loop_size', 1.0)]

    elif hic_file is not None:
        # Traditional format: single Hi-C file path
        hics_files = [hic_file]
        extracted_hic_cmaps = [hic_cmap]
        extracted_triangle_ratios = [triangle_ratio]
        extracted_hic_color_scales = [hic_color_scale]
        extracted_balances = [balance]
        extracted_resolutions = [resolution]
        extracted_flip_vertical = [False]
        extracted_hic_names = [None]
        extracted_insulation_paths = [None]
        extracted_window_sizes = [100000]
        extracted_boundary_cmaps = ['Blues_r']
        extracted_boundary_alphas = [0.6]
        extracted_loops_paths = [None]
        extracted_loop_colors = ['blue']
        extracted_loop_alphas = [0.6]
        extracted_loop_sizes = [1.0]
    else:
        raise ValueError("hics, hic, or hic_file parameter is required")

    # Validate Hi-C file paths
    if not hics_files or hics_files[0] is None:
        raise ValueError("Hi-C file path cannot be empty")

    # 3. Process tracks dict format vs traditional format
    all_kwargs = dict(kwargs)
    tracks_kwargs = None

    if tracks is not None:
        # New format: tracks dict list
        track_files_out = [t['file'] for t in tracks]
        n_tracks = len(tracks)

        # Extract parameters
        extracted_colors = []
        extracted_names = []
        extracted_min_values = []
        extracted_max_values = []    # T-7.4 新加
        extracted_heights = []
        tracks_kwargs = []

        for t in tracks:
            if 'color' in t:
                extracted_colors.append(t['color'])
            if 'name' in t:
                extracted_names.append(t['name'])
            # T-7.8 fix: 总是 append(用 None 占位),保持 idx 跟 track_files 一致
            # - GTF/BED track 没 min_value/max_value 字段,改前会被 skip 导致 idx 错位
            # - 改后用 None 占位,后续 _plot_tracks L367-368 走 tk.get 兜底逻辑
            extracted_min_values.append(t.get('min_value', None))
            extracted_max_values.append(t.get('max_value', None))
            if 'height_cm' in t:
                extracted_heights.append(t['height_cm'])

            track_kwargs = {}
            for key, value in t.items():
                if key not in ['file', 'color', 'name', 'min_value', 'height_cm']:
                    track_kwargs[key] = value
            tracks_kwargs.append(track_kwargs)

        # Use extracted parameters
        if track_colors is None and extracted_colors:
            track_colors = extracted_colors
        if track_names is None and extracted_names:
            track_names = extracted_names
        if min_value is None and extracted_min_values:
            min_value = extracted_min_values
        if track_heights_cm is None and extracted_heights:
            track_heights_cm = extracted_heights

        all_kwargs['tracks_kwargs'] = tracks_kwargs

    elif track_files is not None:
        # Traditional format
        track_files_out = track_files
        n_tracks = len(track_files)
    else:
        # Support no tracks (only Hi-C heatmap)
        track_files_out = []
        n_tracks = 0

    # 4. Validate track count
    if len(track_files_out) != n_tracks:
        raise ValueError(
            f"Number of track files ({len(track_files_out)}) must match n_tracks ({n_tracks})"
        )

    # 5. Set default track heights
    if track_heights_cm is None:
        track_heights_cm = [1.0] * n_tracks
    elif len(track_heights_cm) != n_tracks:
        raise ValueError(
            f"Length of track_heights_cm ({len(track_heights_cm)}) must match n_tracks ({n_tracks})"
        )

    # 6. Calculate layout
    n_hics = len(hics_files)
    triangle_ratio_for_layout = extracted_triangle_ratios[0] if extracted_triangle_ratios else triangle_ratio

    layout = calculate_integrated_layout(
        n_tracks=n_tracks,
        track_heights_cm=track_heights_cm[::-1] if n_tracks > 0 else [],
        width_cm=width_cm,
        gap_cm=gap_cm,
        triangle_ratio=triangle_ratio_for_layout,
        left_margin_cm=left_margin_cm,
        right_margin_cm=right_margin_cm,
        n_hics=n_hics,
        hic_gap_cm=0.1
    )

    # 7. Create visualization
    HeatmapTracks.plot(
        hic_file=hics_files[0],
        track_files=track_files_out,
        region=region,
        layout=layout,
        output=output,
        dpi=dpi,
        hic_cmap=extracted_hic_cmaps[0],
        hic_color_scale=extracted_hic_color_scales[0],
        balance=extracted_balances[0],
        resolution=extracted_resolutions[0],
        track_colors=track_colors,
        track_names=track_names,
        min_value=min_value,
        max_value=max_value,    # T-7.4 新加
        hic_files=hics_files,
        hic_cmaps=extracted_hic_cmaps,
        hic_labels=extracted_hic_names,
        flip_vertical_list=extracted_flip_vertical,
        triangle_ratio_list=extracted_triangle_ratios,
        insulation_paths=extracted_insulation_paths,
        window_size_list=extracted_window_sizes,
        boundary_cmap_list=extracted_boundary_cmaps,
        boundary_alpha_list=extracted_boundary_alphas,
        loops_paths=extracted_loops_paths,
        loop_color_list=extracted_loop_colors,
        loop_alpha_list=extracted_loop_alphas,
        loop_size_list=extracted_loop_sizes,
        **all_kwargs
    )