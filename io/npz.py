"""
NPZ file reader for Hi-C matrix data.

This module provides functions to read and write .npz format files,
which can store multiple NumPy arrays in a single file.
"""

import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List


class NpzReader:
    """
    Reader for .npz format matrix files.
    
    NPZ format can store multiple matrices (e.g., observed, expected)
    in a single file.
    
    Examples
    --------
    >>> reader = NpzReader("matrix.npz")
    >>> obs = reader.get_matrix("observed")
    >>> exp = reader.get_matrix("expected")
    """
    
    def __init__(self, npz_path: str):
        """
        Initialize NpzReader.
        
        Parameters
        ----------
        npz_path : str
            Path to .npz file
        """
        self.path = Path(npz_path)
        
        if not self.path.exists():
            raise FileNotFoundError(f"File not found: {npz_path}")
        
        self._data = np.load(str(self.path))
        self._keys = list(self._data.keys())
    
    @property
    def keys(self) -> List[str]:
        """Get list of array names in the file."""
        return self._keys
    
    def get_matrix(self, name: str) -> np.ndarray:
        """
        Get a matrix by name.
        
        Parameters
        ----------
        name : str
            Name of the array
            
        Returns
        -------
        matrix : np.ndarray
        """
        if name not in self._data:
            raise KeyError(f"Array '{name}' not found. Available: {self._keys}")
        return self._data[name]
    
    def read(self) -> Dict[str, np.ndarray]:
        """
        Read all matrices from file.
        
        Returns
        -------
        matrices : dict
            Dictionary of all arrays
        """
        return {key: self._data[key] for key in self._keys}
    
    def close(self):
        """Close the NPZ file."""
        self._data.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def __repr__(self) -> str:
        return f"NpzReader(path='{self.path}', keys={self._keys})"


def read_npz(npz_path: str) -> Dict[str, np.ndarray]:
    """
    Read a .npz file and return all matrices.
    
    Parameters
    ----------
    npz_path : str
        Path to .npz file
        
    Returns
    -------
    matrices : dict
        Dictionary of all arrays
        
    Examples
    --------
    >>> matrices = read_npz("hic_data.npz")
    >>> print(f"Available: {list(matrices.keys())}")
    """
    data = np.load(npz_path)
    result = {key: data[key] for key in data.keys()}
    data.close()
    return result


def read_npz_matrix(npz_path: str, matrix_name: str = "matrix") -> np.ndarray:
    """
    Read a specific matrix from .npz file.
    
    Parameters
    ----------
    npz_path : str
        Path to .npz file
    matrix_name : str
        Name of the matrix to read
        
    Returns
    -------
    matrix : np.ndarray
    """
    with np.load(npz_path) as data:
        if matrix_name not in data:
            available = list(data.keys())
            raise KeyError(f"Matrix '{matrix_name}' not found. Available: {available}")
        return data[matrix_name]


def save_npz(
    matrices: Dict[str, np.ndarray],
    npz_path: str,
    compressed: bool = True
) -> None:
    """
    Save multiple matrices to .npz file.
    
    Parameters
    ----------
    matrices : dict
        Dictionary of matrices to save
    npz_path : str
        Output path for .npz file
    compressed : bool
        Whether to use compressed format
        
    Examples
    --------
    >>> matrices = {
    ...     "observed": obs_matrix,
    ...     "expected": exp_matrix
    ... }
    >>> save_npz(matrices, "output.npz")
    """
    if compressed:
        np.savez_compressed(npz_path, **matrices)
    else:
        np.savez(npz_path, **matrices)


def save_single_npz(
    matrix: np.ndarray,
    npz_path: str,
    matrix_name: str = "matrix"
) -> None:
    """
    Save a single matrix to .npz file.
    
    Parameters
    ----------
    matrix : np.ndarray
        Matrix to save
    npz_path : str
        Output path
    matrix_name : str
        Name for the matrix in the file
    """
    np.savez_compressed(npz_path, **{matrix_name: matrix})


def read_npz_with_metadata(
    npz_path: str,
    metadata_path: Optional[str] = None
) -> Tuple[Dict[str, np.ndarray], Dict[str, Any]]:
    """
    Read NPZ matrices with associated metadata.
    
    Parameters
    ----------
    npz_path : str
        Path to .npz file
    metadata_path : str, optional
        Path to metadata JSON file
        
    Returns
    -------
    matrices : dict
        Matrix data
    metadata : dict
        Metadata dictionary
    """
    matrices = read_npz(npz_path)
    
    metadata = {}
    if metadata_path and Path(metadata_path).exists():
        import json
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
    
    return matrices, metadata


def save_npz_with_metadata(
    matrices: Dict[str, np.ndarray],
    npz_path: str,
    metadata: Optional[Dict[str, Any]] = None,
    metadata_path: Optional[str] = None
) -> None:
    """
    Save NPZ matrices with associated metadata.
    
    Parameters
    ----------
    matrices : dict
        Dictionary of matrices to save
    npz_path : str
        Output path for .npz file
    metadata : dict, optional
        Metadata dictionary
    metadata_path : str, optional
        Output path for metadata JSON file
    """
    np.savez_compressed(npz_path, **matrices)
    
    if metadata and metadata_path:
        import json
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
