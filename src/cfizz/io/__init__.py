"""
IO module for cfizz package.

Functions for reading and writing Hi-C data in various formats:
- cooler: .cool/.mcool files
- npy: NumPy array files
- npz: NumPy compressed archive files
- insulation: Insulation score data
- loops: Chromatin loop calls
- tad: TAD coordinates
"""

from .cooler import read_cooler, CoolerReader, list_resolutions
from .npy import read_npy, save_npy, NpyReader
from .npz import read_npz, save_npz, read_npz_matrix, NpzReader
from .insulation import (
    calculate_insulation_score,
    compute_insulation_from_cooler,
    compute_insulation_and_save,
    read_insulation_from_file,
    get_available_chroms,
    format_resolution,
)
from .loops import read_loops_bedpe, write_loops_bedpe, LoopSet, loops_to_anchor_pairs
from .tad import read_tads_bed, write_tads_bed, TADSet, tads_to_boundaries
from .paths import (
    Stage, ComputeFeature, VizFeature, AggregateFeature,
    StagePath, FeaturePath, SubstepRegistry, PathBuilder,
    make_filename, extract_qualifiers,
    compartment, tad, loop,
)

__all__ = [
    # Cooler
    "read_cooler",
    "CoolerReader",
    "list_resolutions",
    # NPY
    "read_npy",
    "save_npy",
    "NpyReader",
    # NPZ
    "read_npz",
    "save_npz",
    "read_npz_matrix",
    "NpzReader",
    # Insulation
    "calculate_insulation_score",
    "compute_insulation_from_cooler",
    "compute_insulation_and_save",
    "read_insulation_from_file",
    "get_available_chroms",
    "format_resolution",
    # Loops
    "read_loops_bedpe",
    "write_loops_bedpe",
    "LoopSet",
    "loops_to_anchor_pairs",
    # TAD
    "read_tads_bed",
    "write_tads_bed",
    "TADSet",
    "tads_to_boundaries",
    # Paths
    "Stage",
    "ComputeFeature",
    "VizFeature",
    "AggregateFeature",
    "StagePath",
    "FeaturePath",
    "SubstepRegistry",
    "PathBuilder",
    "make_filename",
    "extract_qualifiers",
]
