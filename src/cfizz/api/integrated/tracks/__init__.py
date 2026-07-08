"""
cfizz tracks module.

Provides simple track drawing (BigWig / GTF / BED), ported from
the upstream reference implementation for visual style alignment.

T-7.3 port: only the minimal subset needed by HeatmapTracks._plot_tracks
is exposed here. Standalone plotting functions
(plot_bw_tracks / plot_gtf_tracks / plot_bed_tracks / plot_mixed_tracks /
quick_plot / etc.) are intentionally NOT re-exported.
"""

from cfizz.api.integrated.tracks.simple import (
    TrackConfig,
    SimpleTrack,
    format_y_axis_value,
    create_track,
    GenomeRange,
)

__all__ = [
    "TrackConfig",
    "SimpleTrack",
    "format_y_axis_value",
    "create_track",
    "GenomeRange",
]