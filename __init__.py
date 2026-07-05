"""
cfizz - Chromosome Flow Integration with Zone-Terminal Visualization.

CFIZZ stands for Chromosome Flow Integration with Zone-Terminal Visualization.

Its core function is end-to-end integration and differential analysis of
chromatin conformation data; the zone-terminal visualization module resides
at the terminal step of the entire pipeline, serving as the last-mile
auxiliary component for generating publication-ready figures.

Name breakdown (CFIZZ, five capital letters):
    C  = Chromosome      (data domain: Hi-C chromatin contact maps)
    F  = Flow            (end-to-end pipeline, one-liner reachability)
    I  = Integration     (multi-omics: Hi-C + BigWig + GTF + BED)
    Z  = Zone-Terminal   (Zonal overlays + Terminal pipeline step)
                          - Zonal:    Adapts to A/B compartments, TADs, loops as overlays
                          - Terminal: Last-mile module that converts analytical results
                                     into publication-ready figures
    Z  = VisualiZation   (the V/Z is folded into the trailing Z)
    (cfizz spells the concept as a single five-letter token, like YAML/JSON,
    without internal separators.)

A unified package for end-to-end Hi-C differential analysis and
zone-terminal visualization.
"""

__version__ = "0.1.0"
__author__ = "cfizz developers"

# Core imports
from . import io
from . import analyze
from . import viz
from . import api
from . import utils

# Version check
def version():
    """Return the package version."""
    return __version__

__all__ = [
    "io",
    "analyze", 
    "viz",
    "api",
    "utils",
    "version",
]
