"""
Visualization module for cfizz package.

Functions for creating publication-quality Hi-C visualizations.
"""

from .heatmap import (
    plot_heatmap, 
    plot_multi_heatmap,
    plot_single_heatmap,
    plot_oe_heatmap,
    get_matrix_range,
    generate_multi_heatmap,       # 多 sample 热图一条龙(T-3.1)
    generate_heatmap,             # 单 sample 一条龙(包装,P3.1)
)

from .compartment import (
    calculate_chromosome_data,     # 算 eig+O/E 一条龙
    plot_heatmap_with_e1,           # 单 sample 拼图
    plot_compartment,                # 单 sample 一条龙
    slice_compartment_region,        # 多 sample 并行 + 区域切片(T-3.4 改名)
    plot_multi_compartment,          # 多 sample 拼图(T-3.2)
    generate_multi_compartment,      # 多 sample 一条龙(T-3.2)
    plot_eigenvector,
)

from .tad import (
    plot_insulation_track,
    plot_tad_heatmap_with_boundaries,
    get_tad_boundary_positions,
    create_tad_summary,
)

# TAD 热图（含 TAD boundaries 标注）从 heatmap_tad_ext 导入
from .heatmap_tad_ext import (
    plot_heatmap_with_tad_boundaries,
    pcolormesh_45deg,
    add_tad_boundaries_to_heatmap,
)

# TAD pileup 从 pileup 导入
from .pileup import (
    plot_multi_tad_boundary_pileup,
)

from .loop import (
    plot_multi_heatmap_with_loops,  # 多 sample loop + 热图
    plot_loops, 
    plot_loop_comparison
)

from .pileup import (
    plot_pileup, 
    plot_multi_pileup, 
)

from .tracks import (
    plot_tracks, 
    add_tracks_to_heatmap
)

from .layout import (
    setup_plot_style,
    setup_ratio_colormap,
    calculate_heatmap_layout,
    calculate_heatmap_with_tracks_layout,
    calculate_rotated_heatmap_layout,
    setup_axes,
    add_coordinate_labels,
    add_rotated_coordinate_labels,
    setup_colorbar,
    setup_horizontal_colorbar,
    log2_and_mask,
    log10_and_mask,
    mask_diagonal,
    save_figure,
    generate_output_filename,
    save_figure_multi_format,
)

__all__ = [
    # Heatmap
    "plot_heatmap",
    "plot_multi_heatmap",
    "plot_single_heatmap",
    "plot_oe_heatmap",
    "get_matrix_range",
    "generate_multi_heatmap",
    "generate_heatmap",
    # Compartment (T-3.2: 新增 5 个函数)
    "calculate_chromosome_data",
    "plot_heatmap_with_e1",
    "plot_compartment",
    "slice_compartment_region",
    "plot_multi_compartment",
    "generate_multi_compartment",
    "plot_eigenvector",
    # TAD
    "plot_insulation_track",
    "plot_tad_heatmap_with_boundaries",
    "get_tad_boundary_positions",
    "create_tad_summary",
    # TAD 热图扩展
    "plot_heatmap_with_tad_boundaries",
    "pcolormesh_45deg",
    "add_tad_boundaries_to_heatmap",
    # TAD pileup
    "plot_multi_tad_boundary_pileup",
    # Loop
    "plot_loops",
    "plot_loop_comparison",
    "plot_multi_heatmap_with_loops",
    # Pileup
    "plot_pileup",
    "plot_multi_pileup",
    # Tracks
    "plot_tracks",
    "add_tracks_to_heatmap",
    "setup_plot_style",
    "setup_ratio_colormap",
    "calculate_heatmap_layout",
    "calculate_heatmap_with_tracks_layout",
    "calculate_rotated_heatmap_layout",
    "setup_axes",
    "add_coordinate_labels",
    "add_rotated_coordinate_labels",
    "setup_colorbar",
    "setup_horizontal_colorbar",
    "log2_and_mask",
    "log10_and_mask",
    "mask_diagonal",
    "save_figure",
    "generate_output_filename",
    "save_figure_multi_format",
]
