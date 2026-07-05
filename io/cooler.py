"""
Cooler file reader for Hi-C data.

This module provides functions to read .cool and .mcool files using the cooler library.
"""

import cooler
import numpy as np
from pathlib import Path
from typing import Optional, List, Tuple, Union, Any


class CoolerReader:
    """
    Reader for cooler format files.
    
    Provides a high-level interface for reading Hi-C contact matrices
    from .cool and .mcool files.
    
    Examples
    --------
    >>> reader = CoolerReader("sample.mcool")
    >>> reader.resolutions  # [10000, 25000, 50000, 100000, ...]
    >>> matrix = reader.fetch("chr1", start=1000000, end=2000000)
    >>> reader.info()  # File metadata
    """
    
    def __init__(self, cool_path: str, resolution: Optional[int] = None):
        """
        Initialize CoolerReader.
        
        Parameters
        ----------
        cool_path : str
            Path to .cool or .mcool file
        resolution : int, optional
            Resolution to use for .mcool files. If None, uses the first resolution.
        """
        self.path = Path(cool_path)
        
        if not self.path.exists():
            raise FileNotFoundError(f"File not found: {cool_path}")
        
        # Open the cooler file
        if resolution:
            self.uri = f'{str(self.path)}::resolutions/{resolution}'
            self.clr = cooler.Cooler(self.uri)
        else:
            self.clr = cooler.Cooler(str(self.path))
            self.uri = self.clr.uri
        
        self._resolution = resolution
        
    @property
    def is_mcool(self) -> bool:
        """Check if this is a multi-resolution .mcool file."""
        return '.mcool' in str(self.path)
    
    @property
    def resolutions(self) -> List[int]:
        """Get list of available resolutions."""
        if self.is_mcool:
            return list(self.clr.store['resolutions'][:])
        else:
            return [self.clr.binsize]
    
    @property
    def chromosomes(self) -> List[str]:
        """Get list of chromosome names."""
        return list(self.clr.chromnames)
    
    @property
    def chromsizes(self) -> dict:
        """Get chromosome sizes."""
        return dict(self.clr.chromsizes)
    
    @property
    def binsize(self) -> int:
        """Get current bin size (resolution)."""
        return self.clr.binsize
    
    @property
    def n_bins(self) -> int:
        """Get total number of bins."""
        return len(self.clr.bins()[:])
    
    def fetch(
        self, 
        chrom: str, 
        start: Optional[int] = None, 
        end: Optional[int] = None,
        balance: bool = True,
        sparse: bool = False
    ) -> np.ndarray:
        """
        Fetch matrix for a region.
        
        Parameters
        ----------
        chrom : str
            Chromosome name (e.g., "chr1", "1")
        start : int, optional
            Start position (bp)
        end : int, optional
            End position (bp)
        balance : bool
            Whether to return balanced matrix
        sparse : bool
            Whether to return sparse matrix
            
        Returns
        -------
        matrix : np.ndarray or sparse matrix
            Contact matrix
            
        Examples
        --------
        >>> # Fetch whole chromosome
        >>> matrix = reader.fetch("chr1")
        >>> # Fetch region
        >>> matrix = reader.fetch("chr1", start=1000000, end=2000000)
        """
        # Normalize chromosome name
        chrom = self._normalize_chrom(chrom)
        
        if start is not None and end is not None:
            region = (chrom, start, end)
        else:
            region = chrom
        
        if sparse:
            return self.clr.matrix(balance=balance, sparse=True).fetch(region)
        else:
            return self.clr.matrix(balance=balance, sparse=False).fetch(region)
    
    def fetch_pixels(
        self,
        chrom: str,
        start: Optional[int] = None,
        end: Optional[int] = None
    ) -> Any:
        """
        Fetch pixel (contact) data for a region.
        
        Parameters
        ----------
        chrom : str
            Chromosome name
        start : int, optional
            Start position
        end : int, optional
            End position
            
        Returns
        -------
        pixels : iterator
            Pixel data as DataFrame-compatible iterator
        """
        chrom = self._normalize_chrom(chrom)
        
        if start is not None and end is not None:
            region = (chrom, start, end)
        else:
            region = chrom
        
        return self.clr.pixels().fetch(region)
    
    def get_chrom_length(self, chrom: str) -> int:
        """Get the length of a chromosome."""
        chrom = self._normalize_chrom(chrom)
        return self.clr.chromsizes[chrom]
    
    def info(self) -> dict:
        """
        Get cooler file metadata.
        
        Returns
        -------
        info : dict
            Dictionary containing file metadata
        """
        return self.clr.info
    
    def _normalize_chrom(self, chrom: str) -> str:
        """
        Normalize chromosome name to match cooler file.
        
        Parameters
        ----------
        chrom : str
            Chromosome name (e.g., "chr1", "1", "CHR1")
            
        Returns
        -------
        chrom : str
            Normalized chromosome name
        """
        chrom = chrom.strip()
        chrom_names = self.clr.chromnames
        
        # If already in chromnames, return as is
        if chrom in chrom_names:
            return chrom
        
        # Try adding/removing "chr" prefix
        if chrom.startswith('chr'):
            chrom_no_chr = chrom[3:]
            if chrom_no_chr in chrom_names:
                return chrom_no_chr
        else:
            chrom_with_chr = f'chr{chrom}'
            if chrom_with_chr in chrom_names:
                return chrom_with_chr
        
        raise ValueError(f"Unknown chromosome: {chrom}. Available: {chrom_names}")
    
    def __repr__(self) -> str:
        return f"CoolerReader(path='{self.path}', binsize={self.binsize}, chromosomes={len(self.chromosomes)})"
    
    def __str__(self) -> str:
        return f"""CoolerReader
  Path: {self.path}
  URI: {self.uri}
  Resolution: {self.binsize:,} bp
  Chromosomes: {len(self.chromosomes)}
  Total bins: {self.n_bins:,}
  Is mcool: {self.is_mcool}
"""


def read_cooler(cool_path: str, resolution: Optional[int] = None) -> CoolerReader:
    """
    Read a cooler file and return a CoolerReader.
    
    Parameters
    ----------
    cool_path : str
        Path to .cool or .mcool file
    resolution : int, optional
        Resolution to use for .mcool files
        
    Returns
    -------
    reader : CoolerReader
        Cooler reader instance
        
    Examples
    --------
    >>> reader = read_cooler("sample.mcool", resolution=10000)
    >>> matrix = reader.fetch("chr1", start=1000000, end=2000000)
    """
    return CoolerReader(cool_path, resolution)


def extract_contact_matrix(
    clr: Union[CoolerReader, cooler.Cooler],
    chrom: str,
    balance: bool = False
) -> Tuple[np.ndarray, int]:
    """
    Extract contact matrix for a chromosome.
    
    Parameters
    ----------
    clr : CoolerReader or cooler.Cooler
        Cooler object
    chrom : str
        Chromosome name
        
    balance : bool
        Whether to return balanced matrix
        
    Returns
    -------
    matrix : np.ndarray
        Contact matrix
    chrom_length : int
        Chromosome length
    """
    if isinstance(clr, CoolerReader):
        cooler_obj = clr.clr
        chrom = clr._normalize_chrom(chrom)
    else:
        cooler_obj = clr
    
    chrom_length = cooler_obj.chromsizes[chrom]
    M = cooler_obj.matrix(balance=balance, sparse=False).fetch((chrom, 0, chrom_length))
    
    # Replace NaN with 0
    M[np.isnan(M)] = 0
    
    return M, chrom_length


def list_resolutions(mcool_path: str) -> List[int]:
    """
    List all available resolutions in an .mcool file.
    
    Parameters
    ----------
    mcool_path : str
        Path to .mcool file
        
    Returns
    -------
    resolutions : list
        List of available resolutions
    """
    clr = cooler.Cooler(mcool_path)
    return list(clr.store['resolutions'][:])


# Convenience function for quick access
def quick_fetch(
    cool_path: str,
    chrom: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    resolution: int = 10000,
    balance: bool = True
) -> np.ndarray:
    """
    Quick function to fetch a Hi-C matrix region.
    
    This is a convenience function for common use cases.
    
    Parameters
    ----------
    cool_path : str
        Path to .cool or .mcool file
    chrom : str
        Chromosome name
    start : int, optional
        Start position
    end : int, optional
        End position
    resolution : int
        Resolution for .mcool files
    balance : bool
        Whether to return balanced matrix
        
    Returns
    -------
    matrix : np.ndarray
        Contact matrix
        
    Examples
    --------
    >>> matrix = quick_fetch("sample.mcool", "chr1", 1000000, 2000000)
    """
    reader = CoolerReader(cool_path, resolution=resolution)
    return reader.fetch(chrom, start, end, balance=balance)
