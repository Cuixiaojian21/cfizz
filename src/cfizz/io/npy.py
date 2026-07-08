"""
NPY file reader for Hi-C matrix data.

This module provides functions to read and write .npy format matrix files.
"""

import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Dict, Any


class NpyReader:
    """
    Reader for .npy format matrix files.
    
    Provides a simple interface for reading Hi-C contact matrices
    stored in NumPy's .npy format.
    
    Examples
    --------
    >>> reader = NpyReader("matrix.npy")
    >>> matrix = reader.read()
    >>> print(f"Matrix shape: {matrix.shape}")
    """
    
    def __init__(self, npy_path: str):
        """
        Initialize NpyReader.
        
        Parameters
        ----------
        npy_path : str
            Path to .npy file
        """
        self.path = Path(npy_path)
        
        if not self.path.exists():
            raise FileNotFoundError(f"File not found: {npy_path}")
    
    def read(self) -> np.ndarray:
        """
        Read the matrix from file.
        
        Returns
        -------
        matrix : np.ndarray
            Contact matrix
        """
        return np.load(str(self.path))
    
    @property
    def shape(self) -> Tuple[int, int]:
        """Get matrix shape without loading data."""
        return np.load(str(self.path), mmap_mode='r').shape
    
    def __repr__(self) -> str:
        return f"NpyReader(path='{self.path}')"


def read_npy(npy_path: str) -> np.ndarray:
    """
    Read a .npy file.
    
    Parameters
    ----------
    npy_path : str
        Path to .npy file
        
    Returns
    -------
    matrix : np.ndarray
        Matrix data
        
    Examples
    --------
    >>> matrix = read_npy("contact_matrix.npy")
    >>> print(f"Matrix shape: {matrix.shape}")
    """
    return np.load(npy_path)


def save_npy(matrix: np.ndarray, npy_path: str) -> None:
    """
    Save a matrix to .npy file.
    
    Parameters
    ----------
    matrix : np.ndarray
        Matrix to save
    npy_path : str
        Output path for .npy file
        
    Examples
    --------
    >>> save_npy(matrix, "output_matrix.npy")
    """
    np.save(npy_path, matrix)


def read_npy_with_metadata(
    npy_path: str,
    metadata_path: Optional[str] = None
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Read NPY matrix with associated metadata.
    
    Parameters
    ----------
    npy_path : str
        Path to .npy file
    metadata_path : str, optional
        Path to metadata JSON file
        
    Returns
    -------
    matrix : np.ndarray
        Matrix data
    metadata : dict
        Metadata dictionary
    """
    matrix = np.load(npy_path)
    
    metadata = {}
    if metadata_path and Path(metadata_path).exists():
        import json
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
    
    return matrix, metadata


def save_npy_with_metadata(
    matrix: np.ndarray,
    npy_path: str,
    metadata: Optional[Dict[str, Any]] = None,
    metadata_path: Optional[str] = None
) -> None:
    """
    Save NPY matrix with associated metadata.
    
    Parameters
    ----------
    matrix : np.ndarray
        Matrix to save
    npy_path : str
        Output path for .npy file
    metadata : dict, optional
        Metadata dictionary
    metadata_path : str, optional
        Output path for metadata JSON file
    """
    np.save(npy_path, matrix)
    
    if metadata and metadata_path:
        import json
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
