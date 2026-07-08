"""
Configuration settings for cfizz package.
"""

# Default resolutions for Hi-C data (in bp)
DEFAULT_RESOLUTIONS = [10000, 25000, 50000, 100000, 250000, 500000, 1000000]

# Default chromosome names
DEFAULT_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# Color schemes
DEFAULT_CMAP = "Reds"
DIFF_CMAP = "RdBu_r"
SADDLE_CMAP = "coolwarm"

# Plot settings
DEFAULT_DPI = 300
PUBLICATION_DPI = 600

# Parallel processing
DEFAULT_N_JOBS = 8
MAX_N_JOBS = 32

# File extensions
SUPPORTED_COOL_EXTENSIONS = [".cool", ".mcool", ".cool.gz"]
SUPPORTED_MATRIX_EXTENSIONS = [".npy", ".npz"]
SUPPORTED_TRACK_EXTENSIONS = [".bw", ".bigwig", ".bed", ".gtf", ".gff"]

# Cache settings
CACHE_DIR = "~/.cfizz/cache"
