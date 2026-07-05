"""
Integrated Heatmap + Tracks Visualization.


This module provides a high-level API for creating publication-quality
integrated figures combining Hi-C heatmaps with genomic tracks.
"""

from typing import List, Optional, Tuple, Dict, Any
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

from cfizz.io import read_cooler
from cfizz.io.insulation import read_insulation_scores
from cfizz.viz.layout import (
    setup_plot_style,
    setup_axes,
    calculate_heatmap_layout,
    setup_colorbar,
)
from cfizz.utils.coordinates import get_matrix_range, print_coordinate
from cfizz.api.integrated.layout import calculate_integrated_layout


class HeatmapTracks:
    """
    Create integrated visualizations combining Hi-C heatmaps with genomic tracks.
    
    This class provides a streamlined interface for generating publication-quality
    figures that combine 45-degree rotated Hi-C heatmaps with genomic feature tracks.
    """
    
    @staticmethod
    def plot(
        hic_file: str,
        track_files: List[str],
        region: 'GenomeRange',
        layout: 'IntegratedLayout',
        output: str,
        dpi: int = 300,
        # Hi-C specific parameters
        hic_cmap: str = 'Reds',
        hic_color_scale: str = 'linear',
        balance: bool = True,
        resolution: int = 10000,
        # Multiple Hi-C support
        hic_files: Optional[List[str]] = None,
        hic_cmaps: Optional[List[str]] = None,
        hic_labels: Optional[List[str]] = None,
        flip_vertical_list: Optional[List[bool]] = None,
        triangle_ratio_list: Optional[List[float]] = None,
        # TAD boundary support
        insulation_paths: Optional[List[str]] = None,
        window_size_list: Optional[List[int]] = None,
        boundary_cmap_list: Optional[List[str]] = None,
        boundary_alpha_list: Optional[List[float]] = None,
        # Loops annotation support
        loops_paths: Optional[List[str]] = None,
        loop_color_list: Optional[List[str]] = None,
        loop_alpha_list: Optional[List[float]] = None,
        loop_size_list: Optional[List[float]] = None,
        # Track parameters
        track_colors: Optional[List[str]] = None,
        track_names: Optional[List[str]] = None,
        min_value: Optional[List[float]] = None,
        max_value: Optional[List[float]] = None,   # T-7.4 新加(跟 min_value 对称)
        **kwargs
    ) -> None:
        """
        Plot integrated Hi-C heatmap with genomic tracks.
        
        Parameters
        ----------
        hic_file : str
            Primary Hi-C file path
        track_files : list
            List of track file paths
        region : GenomeRange
            Genomic region to visualize
        layout : IntegratedLayout
            Layout configuration
        output : str
            Output file path
        dpi : int
            Resolution for PNG
        hic_cmap : str
            Colormap for Hi-C heatmap
        hic_color_scale : str
            Color scale ('linear' or 'log')
        balance : bool
            Whether to balance the matrix
        resolution : int
            Resolution in bp
        hic_files : list, optional
            Multiple Hi-C file paths
        hic_cmaps : list, optional
            Colormaps for each Hi-C
        hic_labels : list, optional
            Labels for each Hi-C
        flip_vertical_list : list, optional
            Flip each matrix vertically
        triangle_ratio_list : list, optional
            Triangle ratio for each heatmap
        insulation_paths : list, optional
            TAD boundary detection files
        window_size_list : list, optional
            TAD window sizes
        boundary_cmap_list : list, optional
            TAD boundary colormaps
        boundary_alpha_list : list, optional
            TAD boundary alphas
        loops_paths : list, optional
            Loops annotation files
        loop_color_list : list, optional
            Loops colors
        loop_alpha_list : list, optional
            Loops alphas
        loop_size_list : list, optional
            Loops sizes
        track_colors : list, optional
            Track colors
        track_names : list, optional
            Track names
        min_value : list, optional
            Minimum y-axis values
        """
        # Determine if using multiple Hi-C
        use_multiple_hics = hic_files is not None and len(hic_files) > 1
        
        # Prepare parameters
        if use_multiple_hics:
            hic_files_list = hic_files
            hic_cmaps_list = hic_cmaps if hic_cmaps else [hic_cmap] * len(hic_files)
            hic_labels_list = hic_labels if hic_labels else [f'Sample {i+1}' for i in range(len(hic_files))]
            flip_list = flip_vertical_list if flip_vertical_list else [False] * len(hic_files)
            tri_list = triangle_ratio_list if triangle_ratio_list else [1.0] * len(hic_files)
            insul_list = insulation_paths if insulation_paths else [None] * len(hic_files)
            win_list = window_size_list if window_size_list else [100000] * len(hic_files)
            bound_cmap_list = boundary_cmap_list if boundary_cmap_list else ['Blues_r'] * len(hic_files)
            bound_alpha_list = boundary_alpha_list if boundary_alpha_list else [0.6] * len(hic_files)
            loops_list = loops_paths if loops_paths else [None] * len(hic_files)
            loop_color_list = loop_color_list if loop_color_list else ['blue'] * len(hic_files)
            loop_alpha_list = loop_alpha_list if loop_alpha_list else [0.6] * len(hic_files)
            # 修复：多 sample 分支默认值也改 50(2026-06-13 修, 跟 164 行单 sample 分支保持一致)
            loop_size_list = loop_size_list if loop_size_list else [50] * len(hic_files)
        else:
            hic_files_list = [hic_file] if hic_files is None else hic_files
            hic_cmaps_list = [hic_cmap]
            hic_labels_list = hic_labels if hic_labels else [None]
            # 使用 flip_vertical_list 的值，而不是硬编码 [False]
            flip_list = [flip_vertical_list[0]] if flip_vertical_list else [False]
            tri_list = [triangle_ratio_list[0] if triangle_ratio_list else 1.0]
            # 修复：使用传入的 insulation_paths，而不是强制设为 None
            insul_list = insulation_paths if insulation_paths else [None]
            win_list = window_size_list if window_size_list else [100000]
            bound_cmap_list = boundary_cmap_list if boundary_cmap_list else ['Blues_r']
            bound_alpha_list = boundary_alpha_list if boundary_alpha_list else [0.6]
            # 修复：单样本模式也要支持 loops_paths
            loops_list = loops_paths if loops_paths else [None]
            loop_color_list = loop_color_list if loop_color_list else ['blue']
            loop_alpha_list = loop_alpha_list if loop_alpha_list else [0.6]
            # 修复：loop_size 默认值改为 50（原来 1.0 太小看不见）
            loop_size_list = loop_size_list if loop_size_list else [50]
        
        # Setup plot style
        setup_plot_style()
        
        # Create figure
        cm_to_inch = 1 / 2.54
        figsize = (layout.total_width_cm * cm_to_inch, layout.total_height_cm * cm_to_inch)
        fig = plt.figure(figsize=figsize)
        
        # Load and prepare Hi-C matrices
        matrices = []
        for hfile in hic_files_list:
            reader = read_cooler(hfile, resolution=resolution)
            matrix = reader.fetch(region.chrom, region.start, region.end, balance=balance)
            matrix = np.nan_to_num(matrix, nan=0.0)
            matrices.append(matrix)
        
        # Calculate global color range
        # 如果外部指定了 vmin/vmax，则使用外部值；否则自动计算
        all_values = np.concatenate([m[np.nonzero(m)] for m in matrices if m is not None and m.size > 0])
        if len(all_values) > 0:
            auto_vmin, auto_vmax = get_matrix_range(all_values)
            # 外部值优先
            vmin = kwargs.get('vmin', auto_vmin)
            vmax = kwargs.get('vmax', auto_vmax)
        else:
            vmin = kwargs.get('vmin', 0.001)
            vmax = kwargs.get('vmax', 1)
        
        # Setup normalization
        if hic_color_scale == 'log':
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)
        
        # Create Hi-C axes and plot heatmaps
        # 反转 axes 以匹配用户预期的 hics 列表顺序（hics[0] 在上方）
        # 但保持 matrices、flip_list 等的原始顺序
        hic_axes_reversed = layout.hic_axes_list[::-1] if len(layout.hic_axes_list) > 1 else layout.hic_axes_list
        
        hic_axes_list = []
        for i in range(len(hic_files_list)):
            if i < len(hic_axes_reversed):
                ax = fig.add_axes(hic_axes_reversed[i])
                hic_axes_list.append(ax)
                
                # Plot 45-degree heatmap
                # matrices[i] 保持原始顺序（与 hics 列表对应）
                _plot_rotated_heatmap(
                    ax,
                    matrices[i],
                    region.start,
                    region.end,
                    resolution,
                    norm=norm,
                    cmap=hic_cmaps_list[i] if i < len(hic_cmaps_list) else hic_cmap,
                    vmin=vmin,
                    vmax=vmax,
                    flip_vertical=flip_list[i] if i < len(flip_list) else False,
                    triangle_ratio=tri_list[i] if i < len(tri_list) else 1.0
                )
                
                # Add sample label
                if hic_labels_list[i]:
                    ax.text(-0.02, 0.5, hic_labels_list[i], ha='right', va='center',
                           fontsize=5, rotation=90, transform=ax.transAxes)
                
                # Plot TAD boundaries if insulation file is provided
                if i < len(insul_list) and insul_list[i] is not None:
                    _plot_tad_boundaries_on_heatmap(
                        ax=ax,
                        insulation_path=insul_list[i],
                        window_size=win_list[i] if i < len(win_list) else 100000,
                        matrix_start=region.start,
                        matrix_end=region.end,
                        resolution=resolution,
                        matrix_shape=matrices[i].shape,
                        boundary_cmap=bound_cmap_list[i] if i < len(bound_cmap_list) else 'Blues_r',
                        boundary_alpha=bound_alpha_list[i] if i < len(bound_alpha_list) else 0.6,
                        flip_vertical=flip_list[i] if i < len(flip_list) else False,
                        triangle_ratio=tri_list[i] if i < len(tri_list) else 1.0,
                        chrom=region.chrom  # 添加 chromosome 参数
                    )
                
                # Add loops annotations if loops file is provided
                if i < len(loops_list) and loops_list[i] is not None:
                    _add_loops_to_heatmap_45deg(
                        ax=ax,
                        loops_path=loops_list[i],
                        matrix=matrices[i],
                        start=region.start,
                        resolution=resolution,
                        chrom=region.chrom,
                        loop_color=loop_color_list[i] if i < len(loop_color_list) else 'blue',
                        loop_alpha=loop_alpha_list[i] if i < len(loop_alpha_list) else 0.6,
                        loop_size=loop_size_list[i] if i < len(loop_size_list) else 5,
                        flip_vertical=flip_list[i] if i < len(flip_list) else False
                    )
        
        # Add colorbar with horizontal label (for 45-degree rotated heatmap)
        # FIX: vmin/vmax might be None (external not specified + auto failed on empty matrix)
        _vmin = vmin if vmin is not None else 0.001
        _vmax = vmax if vmax is not None else 1.0
        if len(hic_axes_list) > 0:
            ax = hic_axes_list[0]
            cbar_left = (layout.left_margin_cm + layout.core_width_cm + 0.2) / layout.total_width_cm
            setup_colorbar(
                fig=fig,
                sc=ax.collections[0] if ax.collections else None,
                vmin=_vmin,
                vmax=_vmax,
                balance=balance,
                colorbar_left=cbar_left,
                colorbar_bottom=0.75,
                colorbar_width=0.2 / layout.total_width_cm,
                colorbar_height=0.15,
                label_fontsize=5,
                tick_fontsize=5,
                label_orientation='horizontal'  # 水平 label（适合 45 度旋转热图）
            )
        
        # Plot tracks (T-7.2 + T-7.3 + T-7.4)
        tracks_kwargs = kwargs.get('tracks_kwargs', None)
        HeatmapTracks._plot_tracks(
            fig=fig,
            layout=layout,
            region=region,
            track_files=track_files,
            track_colors=track_colors,
            track_names=track_names,
            min_value=min_value,
            max_value=max_value,    # T-7.4 新加(跟 min_value 对称,不反转,跟现有 min_value 一致)
            tracks_kwargs=tracks_kwargs,
        )

        # Add coordinate labels
        _add_coordinate_labels(fig, layout.coordinate_label_area, region)

        # Save output
        output_path = Path(output)
        fig.savefig(f"{output_path}.png", dpi=dpi)
        fig.savefig(f"{output_path}.svg")
        plt.close()

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    @staticmethod
    def _plot_tracks(
        fig,
        layout,
        region,
        track_files,
        track_colors=None,
        track_names=None,
        min_value=None,
        max_value=None,    # T-7.4 新加(跟 min_value 对称)
        tracks_kwargs=None,
    ):
        """
        绘制 genomic tracks 到 fig 上的 layout.tracks_axes 指定位置。

        """
        import os as _os
        from cfizz.api.integrated.tracks.simple import (
            create_track,
            format_y_axis_value,
            GenomeRange as TrackGenomeRange,
        )

        n_tracks = len(track_files)
        if n_tracks == 0:
            return

        for i, track_file in enumerate(track_files):
            if i >= len(layout.tracks_axes):
                break

            # T-7.6: 反转 layout.tracks_axes 引用,让 user tracks[0] 画到顶部(贴 Hi-C)
            # 因为 layout.tracks_axes[0] 在底部,axes[::-1][0] = axes[N-1] 在顶部
            ax = fig.add_axes(layout.tracks_axes[::-1][i])

            # 提取 per-track kwargs
            tk = {}
            if tracks_kwargs and i < len(tracks_kwargs):
                tk = dict(tracks_kwargs[i])

            # 提取 styling(从外层列表)
            color = track_colors[i] if track_colors and i < len(track_colors) else '#333333'
            name = track_names[i] if track_names and i < len(track_names) else None

            # 应用 color/name 到 tk(传给 create_track)
            tk['color'] = color
            if name is not None:
                tk['name'] = name

            # T-7.4: 统一提取 min_value / max_value(外层列表优先,tk 兜底)
            y_min = min_value[i] if min_value and i < len(min_value) else tk.get('min_value', None)
            y_max = max_value[i] if max_value and i < len(max_value) else tk.get('max_value', None)

            # 注入到 tk,create_track 会把它传给 TrackConfig.min_value/max_value
            if y_max is not None:
                tk['max_value'] = y_max
            if y_min is not None:
                tk['min_value'] = y_min

            # 创建 track
            try:
                track = create_track(track_file, **tk)
            except Exception as e:
                ax.text(
                    0.5, 0.5,
                    f'? {_os.path.basename(track_file)}: {e}',
                    ha='center', va='center',
                    transform=ax.transAxes, fontsize=5,
                )
                continue

            # 转 region 到 TrackGenomeRange(可能类型不一样)
            track_region = TrackGenomeRange(region.chrom, region.start, region.end)

            # 绘制
            try:
                track.plot(ax, track_region)
            except Exception as e:
                ax.text(
                    0.5, 0.5,
                    f'! plot failed: {e}',
                    ha='center', va='center',
                    transform=ax.transAxes, fontsize=5,
                )
                continue

            # title(右侧)
            if track.config.name:
                ax.text(
                    1.02, 0.5, str(track.config.name),
                    fontsize=tk.get('fontsize', 5),
                    ha='left', va='center',
                    transform=ax.transAxes,
                )

            # 样式
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)

            is_gtf_track = (
                hasattr(track, 'track_type')
                and track.track_type in ['gtf', 'bed']
            )

            if is_gtf_track:
                # GTF/BED 不需要 y 轴标签
                ax.spines['left'].set_visible(False)
                ax.set_yticks([])
                ax.set_ylabel('')
            else:
                ax.spines['left'].set_visible(True)
                ylim = ax.get_ylim()
                y_min_curr, y_max_curr = ylim

                if 'min_value' in tk and tk['min_value'] is not None:
                    y_min_use = tk['min_value']
                else:
                    y_min_use = max(y_min_curr, 0.0)

                # T-7.7: 分情况加 headroom
                # - 用户外部指定 max_value:忠实于用户设的值,不加 headroom
                # - 用户没指定 max_value:加 20% headroom 避免数据截断
                if 'max_value' in tk and tk['max_value'] is not None:
                    new_min = y_min_use
                    new_max = y_max_curr
                else:
                    y_range = y_max_curr - y_min_use
                    y_padding = y_range * 0.2
                    new_min = y_min_use
                    new_max = y_max_curr + y_padding

                ax.set_ylim(new_min, new_max)
                ax.set_yticks([new_min, new_max])
                ax.tick_params(
                    axis='y', which='major', length=3, labelsize=5,
                    left=True, labelleft=False,
                )

                # 手动加 y 轴标签(black text)
                y_offset = (new_max - new_min) * 0.2
                label_y_min = new_min + y_offset
                label_y_max = new_max - y_offset

                ax.text(
                    -0.05, label_y_min, format_y_axis_value(new_min),
                    ha='right', va='center', fontsize=5, color='black',
                    transform=ax.get_yaxis_transform(),
                )
                ax.text(
                    -0.05, label_y_max, format_y_axis_value(new_max),
                    ha='right', va='center', fontsize=5, color='black',
                    transform=ax.get_yaxis_transform(),
                )

            # 隐藏 x 轴(由底部 coordinate_label_area 统一画)
            ax.set_xticks([])
            ax.set_xlabel('')


def _plot_rotated_heatmap(
    ax,
    matrix: np.ndarray,
    start: int,
    end: int,
    resolution: int,
    norm=None,
    cmap='Reds',
    vmin=None,
    vmax=None,
    flip_vertical=False,
    triangle_ratio=1.0
):
    """
    Plot 45-degree rotated Hi-C heatmap.
    
    """
    # Handle colormap
    if isinstance(cmap, str):
        cmap = plt.colormaps.get_cmap(cmap)
    
    n = matrix.shape[0]
    import itertools
    
    # Setup coordinates using actual genomic positions
    start_pos_vector = [start + resolution * i for i in range(n + 1)]
    t = np.array([[1, 0.5], [-1, 0.5]])
    
    # Generate rotated coordinates
    if flip_vertical:
        # Lower triangle
        matrix_a = np.dot(
            np.array([
                (i[1], i[0])
                for i in itertools.product(start_pos_vector, start_pos_vector)
            ]),
            t
        )
    else:
        # Upper triangle (flipped)
        matrix_a = np.dot(
            np.array([
                (i[1], i[0])
                for i in itertools.product(start_pos_vector[::-1], start_pos_vector)
            ]),
            t
        )
    
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    
    # Process matrix
    processed_matrix = matrix.copy()
    mask = ~np.isnan(processed_matrix)
    if vmin is not None and vmax is not None:
        processed_matrix[mask] = np.clip(processed_matrix[mask], vmin, vmax)
    
    # Plot
    if flip_vertical:
        im = ax.pcolormesh(x, y, processed_matrix, norm=norm, cmap=cmap)
    else:
        im = ax.pcolormesh(x, y, np.flipud(processed_matrix), norm=norm, cmap=cmap)
    
    # Style axes
    ax.set_aspect(0.5)
    
    if flip_vertical:
        ax.set_ylim(-(end - start) * triangle_ratio, 0)
    else:
        ax.set_ylim(0, (end - start) * triangle_ratio)
    
    ax.set_xlim(start, end)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    im.set_rasterized(True)


def _add_coordinate_labels(fig, label_area, region):
    """Add genomic coordinate labels with formatted coordinates (e.g., 5M instead of 5,000,000)."""
    ax = fig.add_axes(label_area)
    ax.text(0, 0.5, print_coordinate(region.start), ha='left', va='center', fontsize=5, transform=ax.transAxes)
    ax.text(1, 0.5, print_coordinate(region.end), ha='right', va='center', fontsize=5, transform=ax.transAxes)
    ax.text(0.5, 0.5, region.chrom, ha='center', va='center', fontsize=5, transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _mark_boundaries_from_insulation(
    insulation_data: pd.DataFrame,
    matrix_data: np.ndarray,
    matrix_start: int,
    resolution: int
) -> np.ndarray:
    """
    Create a boundary matrix from insulation data.
    
    
    Key difference: When multiple boundaries exist, marks the regions BETWEEN
    boundaries (TAD blocks), not just the boundary lines. This creates proper
    TAD boundary visualization.
    
    Parameters
    ----------
    insulation_data : pd.DataFrame
        DataFrame with columns: chrom, start, end, insulation_score, is_boundary
    matrix_data : np.ndarray
        Contact matrix
    matrix_start : int
        Start position of the matrix region
    resolution : int
        Bin resolution
        
    Returns
    -------
    boundary_matrix : np.ndarray
        Matrix with boundary markers (1 for boundary regions, NaN elsewhere)
    """
    boundary_matrix = np.full_like(matrix_data, np.nan, dtype=float)
    
    if 'is_boundary' not in insulation_data.columns:
        return boundary_matrix
    
    boundaries = insulation_data[insulation_data['is_boundary']]
    if len(boundaries) == 0:
        return boundary_matrix
    
    # Convert boundary positions to bin indices
    boundary_indices = []
    for _, boundary in boundaries.iterrows():
        idx = int((boundary['start'] - matrix_start) // resolution)
        if 0 <= idx < matrix_data.shape[0]:
            boundary_indices.append(idx)
    
    # Mark boundaries based on number of boundary indices
    if len(boundary_indices) == 1:
        # Single boundary: mark the row and column
        i = boundary_indices[0]
        boundary_matrix[i, :] = 1
        boundary_matrix[:, i] = 1
    else:
        # Multiple boundaries: mark the regions BETWEEN adjacent boundaries
        # This creates proper TAD block visualization
        for k in range(len(boundary_indices) - 1):
            i = boundary_indices[k]
            j = boundary_indices[k + 1]
            # Mark the boundary lines
            boundary_matrix[i, i:j+1] = 1
            boundary_matrix[j, i:j+1] = 1
            boundary_matrix[i:j+1, i] = 1
            boundary_matrix[i:j+1, j] = 1
    
    return boundary_matrix


def _plot_tad_boundaries_on_heatmap(
    ax,
    insulation_path: str,
    window_size: int,
    matrix_start: int,
    matrix_end: int,
    resolution: int,
    matrix_shape: tuple,
    boundary_cmap: str = 'Blues_r',
    boundary_alpha: float = 0.6,
    flip_vertical: bool = False,
    triangle_ratio: float = 0.5,
    chrom: str = None
):
    """
    Plot TAD boundaries on a 45-degree rotated heatmap.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    insulation_path : str
        Path to insulation TSV file
    window_size : int
        Window size in bp
    matrix_start : int
        Start position of matrix region
    matrix_end : int
        End position of matrix region
    resolution : int
        Bin resolution
    matrix_shape : tuple
        Shape of the matrix (n_bins, n_bins)
    boundary_cmap : str
        Colormap for boundaries
    boundary_alpha : float
        Transparency of boundary markers
    flip_vertical : bool
        Whether heatmap is flipped vertically
    triangle_ratio : float
        Triangle ratio for layout
    chrom : str
        Chromosome name for filtering
    """
    try:
        # Read insulation data
        # 关键修复：必须传递 chrom 参数以确保只读取目标染色体的数据
        ins_data = read_insulation_scores(
            file_path=insulation_path,
            windows=window_size,
            chrom=chrom,  # 添加 chromosome 过滤
            start=matrix_start,
            end=matrix_end
        )
        
        if ins_data is None or len(ins_data) == 0:
            return
        
        # Create boundary matrix
        n = matrix_shape[0]
        mock_matrix = np.zeros((n, n))  # Dummy matrix for boundary calculation
        boundary_matrix = _mark_boundaries_from_insulation(
            ins_data, mock_matrix, matrix_start, resolution
        )
        
        # Handle colormap
        if isinstance(boundary_cmap, str):
            cmap = plt.colormaps.get_cmap(boundary_cmap)
        
        import itertools
        
        # Setup coordinates
        start_pos_vector = [matrix_start + resolution * i for i in range(n + 1)]
        t = np.array([[1, 0.5], [-1, 0.5]])
        
        # Generate rotated coordinates (same as _plot_rotated_heatmap)
        if flip_vertical:
            matrix_a = np.dot(
                np.array([
                    (i[1], i[0])
                    for i in itertools.product(start_pos_vector, start_pos_vector)
                ]),
                t
            )
        else:
            matrix_a = np.dot(
                np.array([
                    (i[1], i[0])
                    for i in itertools.product(start_pos_vector[::-1], start_pos_vector)
                ]),
                t
            )
        
        x = matrix_a[:, 1].reshape(n + 1, n + 1)
        y = matrix_a[:, 0].reshape(n + 1, n + 1)
        
        # Process boundary matrix
        processed_matrix = boundary_matrix.copy()
        mask = ~np.isnan(processed_matrix)
        processed_matrix[mask] = np.clip(processed_matrix[mask], 0, 1)
        
        # Plot boundary overlay
        if flip_vertical:
            ax.pcolormesh(x, y, processed_matrix, cmap=cmap, alpha=boundary_alpha, 
                         vmin=0, vmax=1, rasterized=True)
        else:
            ax.pcolormesh(x, y, np.flipud(processed_matrix), cmap=cmap, alpha=boundary_alpha,
                         vmin=0, vmax=1, rasterized=True)
        
    except Exception as e:
        print(f"Warning: Could not plot TAD boundaries: {e}")


def _add_loops_to_heatmap_45deg(
    ax,
    loops_path: str,
    matrix: np.ndarray,
    start: int,
    resolution: int,
    chrom: str,
    loop_color: str = 'blue',
    loop_alpha: float = 0.6,
    loop_size: float = 5,
    flip_vertical: bool = False
):
    """
    Add loops annotations to a 45-degree rotated Hi-C heatmap.

    This function uses the same coordinate transformation as _pcolormesh_45deg:
    - Rotation matrix: t = [[1, 0.5], [-1, 0.5]]
    - Transforms genomic coordinates (i, j) to rotated coordinates (x, y)
    - Key difference: coordinates are passed in REVERSED order to the rotation matrix

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    loops_path : str
        Path to loops BEDPE/TSV file
    matrix : np.ndarray
        Hi-C contact matrix
    start : int
        Region start position (bp)
    resolution : int
        Bin resolution (bp)
    chrom : str
        Chromosome name
    loop_color : str
        Color for loop markers
    loop_alpha : float
        Transparency of loop markers
    loop_size : float
        Size of loop markers (scatter s parameter)
    flip_vertical : bool
        Whether heatmap is flipped vertically
    """
    try:
        # Import read_loops function
        from cfizz.io.loops import read_loops

        # Read loops data
        loops_data = read_loops(
            file_path=loops_path,
            chrom=chrom,
            start=start,
            end=start + matrix.shape[0] * resolution
        )

        if loops_data is None or len(loops_data) == 0:
            print(f"No loops found in region {chrom}:{start:,}")
            return

        print(f"\n🔍 Adding {len(loops_data)} loops to 45-deg heatmap...")

        # Rotation matrix for 45-degree transformation (same as heatmap)
        rotation_matrix = np.array([[1, 0.5], [-1, 0.5]])
        region_start = start
        matrix_size = matrix.shape[0]

        # Process each loop
        loops_added = 0
        loops_skipped = 0
        for idx, loop in loops_data.iterrows():
            # Get loop anchors (use genomic coordinates, NOT bin indices!)
            anchor1_genome = (loop['start1'] + loop['end1']) / 2
            anchor2_genome = (loop['start2'] + loop['end2']) / 2

            # Check if anchors are within matrix bounds (using genomic coordinates)
            if not (region_start <= anchor1_genome < region_start + matrix_size * resolution and
                    region_start <= anchor2_genome < region_start + matrix_size * resolution):
                loops_skipped += 1
                continue

            # Add offset to convert center point to cell edge coordinates
            # Heatmap uses cell edges, loop uses center point
            offset = resolution / 2

            if flip_vertical:
                # 下三角热图 - v5 综合方案(王 scan + v4 bin-center)
                # 关键改进:
                #   1. 王 scan:在 anchor BEDPE 区间内扫描找最强信号 bin
                #      抗 BEDPE 坐标精度误差(pyHICCUPS 输出可能 1-2 bin 区间)
                #   2. v4 bin-center:用 heatmap 自己的 x_mesh, y_mesh 算 bin 视觉中心
                #      精准对齐菱形 bin 中央,不靠数学推导
                
                import itertools
                n = matrix.shape[0]
                start_pos_vector = [start + resolution * i for i in range(n + 1)]
                t = np.array([[1, 0.5], [-1, 0.5]])
                
                # 下三角 - 跟 _plot_rotated_heatmap L326-332 一致(不翻转)
                matrix_a = np.dot(
                    np.array([
                        (i[1], i[0])
                        for i in itertools.product(start_pos_vector, start_pos_vector)
                    ]),
                    t
                )
                x_mesh = matrix_a[:, 1].reshape(n + 1, n + 1)
                y_mesh = matrix_a[:, 0].reshape(n + 1, n + 1)
                
                # 王 scan:在 anchor BEDPE 区间找最强信号 bin
                s_l = range(int(anchor1_genome // resolution),
                            int(np.ceil((anchor1_genome + offset) / float(resolution))) + 1)
                e_l = range(int(anchor2_genome // resolution),
                            int(np.ceil((anchor2_genome + offset) / float(resolution))) + 1)
                
                si, ei, max_signal = None, None, -np.inf
                region_start_bin = int(start // resolution)
                for i in s_l:
                    for j in e_l:
                        st = i - region_start_bin
                        et = j - region_start_bin
                        if (0 <= st < matrix.shape[0]) and (0 <= et < matrix.shape[0]):
                            if matrix[st, et] > max_signal:
                                si, ei = st, et
                                max_signal = matrix[st, et]
                
                if si is None or ei is None:
                    loops_skipped += 1
                    continue
                
                # 下三角:heatmap (k1, k2) = (si, ei)
                k1, k2 = si, ei
                
                # bin 视觉中心 = 对角 2 个 cell edge 顶点平均
                x1 = (x_mesh[k1, k2] + x_mesh[k1+1, k2+1]) / 2
                y1 = (y_mesh[k1, k2] + y_mesh[k1+1, k2+1]) / 2
                
                # 对称点(transpose)
                k1_t, k2_t = ei, si
                x2 = (x_mesh[k1_t, k2_t] + x_mesh[k1_t+1, k2_t+1]) / 2
                y2 = (y_mesh[k1_t, k2_t] + y_mesh[k1_t+1, k2_t+1]) / 2
            else:
                # 上三角热图 - v5 综合方案(王 scan + v4 bin-center)
                # 关键改进:
                #   1. 王 scan:在 anchor BEDPE 区间内扫描找最强信号 bin
                #      抗 BEDPE 坐标精度误差(pyHICCUPS 输出可能 1-2 bin 区间)
                #   2. v4 bin-center:用 heatmap 自己的 x_mesh, y_mesh 算 bin 视觉中心
                #      精准对齐菱形 bin 中央,不靠数学推导
                
                import itertools
                n = matrix.shape[0]
                start_pos_vector = [start + resolution * i for i in range(n + 1)]
                t = np.array([[1, 0.5], [-1, 0.5]])
                
                # 上三角 - 跟 _plot_rotated_heatmap L335-341 一致(start_pos_vector[::-1] 翻转 x)
                matrix_a = np.dot(
                    np.array([
                        (i[1], i[0])
                        for i in itertools.product(start_pos_vector[::-1], start_pos_vector)
                    ]),
                    t
                )
                x_mesh = matrix_a[:, 1].reshape(n + 1, n + 1)
                y_mesh = matrix_a[:, 0].reshape(n + 1, n + 1)
                
                # 王 scan:在 anchor BEDPE 区间找最强信号 bin
                s_l = range(int(anchor1_genome // resolution),
                            int(np.ceil((anchor1_genome + offset) / float(resolution))) + 1)
                e_l = range(int(anchor2_genome // resolution),
                            int(np.ceil((anchor2_genome + offset) / float(resolution))) + 1)
                
                si, ei, max_signal = None, None, -np.inf
                region_start_bin = int(start // resolution)
                for i in s_l:
                    for j in e_l:
                        st = i - region_start_bin
                        et = j - region_start_bin
                        if (0 <= st < matrix.shape[0]) and (0 <= et < matrix.shape[0]):
                            if matrix[st, et] > max_signal:
                                si, ei = st, et
                                max_signal = matrix[st, et]
                
                if si is None or ei is None:
                    loops_skipped += 1
                    continue
                
                # 上三角:heatmap (k1, k2) = (n-1-si, ei)
                k1, k2 = n - 1 - si, ei
                
                # bin 视觉中心 = 对角 2 个 cell edge 顶点平均
                x1 = (x_mesh[k1, k2] + x_mesh[k1+1, k2+1]) / 2
                y1 = (y_mesh[k1, k2] + y_mesh[k1+1, k2+1]) / 2
                
                # 对称点(transpose:交换 si 和 ei)
                k1_t, k2_t = n - 1 - ei, si
                x2 = (x_mesh[k1_t, k2_t] + x_mesh[k1_t+1, k2_t+1]) / 2
                y2 = (y_mesh[k1_t, k2_t] + y_mesh[k1_t+1, k2_t+1]) / 2

            # Plot the two circular markers
            # Marker at anchor position
            ax.scatter(x1, y1, s=loop_size, c='none', marker='o',
                      edgecolors=loop_color, alpha=loop_alpha, linewidths=0.8,
                      zorder=10)  # Ensure loops are on top

            # Marker at symmetric position
            ax.scatter(x2, y2, s=loop_size, c='none', marker='o',
                      edgecolors=loop_color, alpha=loop_alpha, linewidths=0.8,
                      zorder=10)  # Ensure loops are on top

            loops_added += 1

        print(f"✓ Added {loops_added} loops, skipped {loops_skipped}")

    except Exception as e:
        print(f"Warning: Could not add loops: {e}")


class GenomeRange:
    """Simple genomic range class."""
    
    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        self.start = start
        self.end = end
    
    @property
    def size(self) -> int:
        return self.end - self.start
    
    def __repr__(self):
        return f"GenomeRange('{self.chrom}', {self.start:,}, {self.end:,})"


def quick_plot_integrated(
    hic_file: Optional[str] = None,
    hic: Optional[dict] = None,
    hics: Optional[List[dict]] = None,
    track_files: Optional[List[str]] = None,
    tracks: Optional[List[dict]] = None,
    region: Optional[GenomeRange] = None,
    output: Optional[str] = None,
    n_tracks: int = 2,
    width_cm: float = 5.0,
    gap_cm: float = 0.1,
    left_margin_cm: float = 1.0,
    right_margin_cm: float = 2.0,
    dpi: int = 300,
    triangle_ratio: float = 0.5,
    hic_cmap: str = 'Reds',
    hic_color_scale: str = 'linear',
    balance: bool = True,
    resolution: int = 10000,
    track_colors: Optional[List[str]] = None,
    track_names: Optional[List[str]] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    # TAD boundary parameters
    insulation_paths: Optional[List[str]] = None,
    window_size_list: Optional[List[int]] = None,
    boundary_cmap_list: Optional[List[str]] = None,
    boundary_alpha_list: Optional[List[float]] = None,
    # Loops annotation parameters
    loops_paths: Optional[List[str]] = None,
    loop_color_list: Optional[List[str]] = None,
    loop_alpha_list: Optional[List[float]] = None,
    loop_size_list: Optional[List[float]] = None,
    **kwargs
) -> None:
    """
    Quick plot function for integrated Hi-C heatmap and tracks.
    
    This is a convenience wrapper for common use cases.
    
    Parameters
    ----------
    hic_file : str
        Hi-C file path
    hic : dict
        Hi-C configuration dictionary
    hics : list
        List of Hi-C configuration dictionaries
    track_files : list
        Track file paths
    tracks : list
        Track configuration dictionaries
    region : GenomeRange
        Genomic region
    output : str
        Output file path
    n_tracks : int
        Number of tracks
    width_cm : float
        Figure width in cm
    gap_cm : float
        Gap between elements
    left_margin_cm : float
        Left margin in cm
    right_margin_cm : float
        Right margin in cm
    dpi : int
        Resolution
    triangle_ratio : float
        Heatmap height ratio
    hic_cmap : str
        Colormap
    hic_color_scale : str
        Color scale
    balance : bool
        Balance matrix
    resolution : int
        Resolution in bp
    track_colors : list
        Track colors
    track_names : list
        Track names
    vmin : float, optional
        Minimum value for color scale (auto-calculated if not provided)
    vmax : float, optional
        Maximum value for color scale (auto-calculated if not provided)
    """
    if region is None:
        raise ValueError("region is required")
    if output is None:
        raise ValueError("output is required")
    
    # Process Hi-C input
    hic_files = None
    hic_cmaps = None
    hic_labels = None
    flip_vertical_list = None
    
    if hics is not None:
        hic_files = [h.get('file') for h in hics]
        hic_cmaps = [h.get('cmap', hic_cmap) for h in hics]
        hic_labels = [h.get('name') for h in hics]
        flip_vertical_list = [h.get('flip_vertical', False) for h in hics]
        # ⚠️ 修复: 提取 loops_paths
        loops_paths = [h.get('loops_path') for h in hics]
        loop_color_list = [h.get('loop_color', 'blue') for h in hics]
        loop_alpha_list = [h.get('loop_alpha', 0.6) for h in hics]
        loop_size_list = [h.get('loop_size', 10) for h in hics]
    elif hic is not None:
        hic_files = [hic.get('file')]
        hic_cmaps = [hic.get('cmap', hic_cmap)]
        hic_labels = [hic.get('name')]
        flip_vertical_list = [hic.get('flip_vertical', False)]
        # ⚠️ 修复: 提取 loops_path
        loops_paths = [hic.get('loops_path')]
        loop_color_list = [hic.get('loop_color', 'blue')]
        loop_alpha_list = [hic.get('loop_alpha', 0.6)]
        loop_size_list = [hic.get('loop_size', 50)]
    elif hic_file is not None:
        hic_files = [hic_file]
        hic_cmaps = [hic_cmap]
        hic_labels = [None]
        flip_vertical_list = [False]
    else:
        raise ValueError("hic_file, hic, or hics is required")
    
    # Process track input
    if tracks is not None:
        track_files_out = [t.get('file') for t in tracks]
        if track_colors is None:
            track_colors = [t.get('color') for t in tracks]
        if track_names is None:
            track_names = [t.get('name') for t in tracks]
        n_tracks = len(tracks)
    else:
        track_files_out = track_files if track_files else []
    
    # Calculate layout
    # Handle case where n_tracks is 0
    track_heights = [1.0] * n_tracks if n_tracks > 0 else []
    layout = calculate_integrated_layout(
        n_tracks=n_tracks,
        track_heights_cm=track_heights,
        width_cm=width_cm,
        gap_cm=gap_cm,
        triangle_ratio=triangle_ratio,
        left_margin_cm=left_margin_cm,
        right_margin_cm=right_margin_cm,
        n_hics=len(hic_files),
        hic_gap_cm=0.1
    )
    
    # Plot
    HeatmapTracks.plot(
        hic_file=hic_files[0],
        track_files=track_files_out,
        region=region,
        layout=layout,
        output=output,
        dpi=dpi,
        hic_cmap=hic_cmap,
        hic_color_scale=hic_color_scale,
        balance=balance,
        resolution=resolution,
        hic_files=hic_files,
        hic_cmaps=hic_cmaps,
        hic_labels=hic_labels,
        flip_vertical_list=flip_vertical_list,
        track_colors=track_colors,
        track_names=track_names,
        vmin=vmin,  # 外部指定的最小值
        vmax=vmax,  # 外部指定的最大值
        # TAD boundary parameters
        insulation_paths=insulation_paths,
        window_size_list=window_size_list,
        boundary_cmap_list=boundary_cmap_list,
        boundary_alpha_list=boundary_alpha_list,
        # Loops annotation parameters
        loops_paths=loops_paths,
        loop_color_list=loop_color_list,
        loop_alpha_list=loop_alpha_list,
        loop_size_list=loop_size_list,
        **kwargs
    )
