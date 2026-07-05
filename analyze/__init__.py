"""
Analysis module for cfizz package.

Functions for Hi-C data analysis:
- oe: Observed/Expected normalization
- compartment: A/B compartment analysis
- tad: TAD identification and analysis
- loop: Chromatin loop detection
- distance: Distance decay analysis
"""

from .oe import (
    calculate_oe_matrix,
    compute_expected,
    quick_oe,
    calculate_decay,
    log2_and_mask
)
from .compartment import (
    assign_compartments,
    calculate_compartment_strength,
    # Diff analysis
    load_compartment_data,
    calculate_compartment_diff,
    plot_compartment_scatter,
    analyze_single_comparison as analyze_compartment_diff,
)
from .tad import (
    calculate_insulation,
    identify_boundaries,
    call_tads,
    calculate_tad_statistics,
    # Diff analysis
    calculate_boundary_context,
    generate_dual_direction_pairing,
    generate_final_classification,
    plot_tad_stacked_bar,
    analyze_single_comparison_window as analyze_tad_diff,
)
from .loop import (
    call_loops_hiccups,
    calculate_loop_score,
    filter_loops_by_distance,
    # Diff analysis
    load_loops_data,
    find_matching_loops,
    plot_loops_stacked_bar,
    analyze_single_comparison_loops as analyze_loop_diff,
)
from .distance import (
    calculate_distance_decay,
    fit_power_law,
    fit_exponential,
    calculate_decay_rate,
    calculate_decay_profile,
    compare_decay_rates
)

__all__ = [
    # O/E
    "calculate_oe_matrix",
    "compute_expected",
    "quick_oe",
    "calculate_decay",
    "log2_and_mask",
    # Compartment
    "assign_compartments",
    "calculate_compartment_strength",
    # TAD
    "calculate_insulation",
    "identify_boundaries",
    "call_tads",
    "calculate_tad_statistics",
    # Loop
    "call_loops_hiccups",
    "calculate_loop_score",
    "filter_loops_by_distance",
    # Distance
    "calculate_distance_decay",
    "fit_power_law",
    "fit_exponential",
    "calculate_decay_rate",
    "calculate_decay_profile",
    "compare_decay_rates",
]
