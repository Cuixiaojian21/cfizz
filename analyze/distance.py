"""
Distance decay analysis module.

This module provides functions for analyzing the distance-dependent
decay of chromatin interactions in Hi-C data.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Any
from scipy.optimize import curve_fit
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)


def calculate_distance_decay(
    matrix: np.ndarray,
    resolution: int,
    max_distance: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate distance decay curve for Hi-C contacts.
    
    The distance decay follows a power law: P(s) ∝ s^(-α)
    where s is the genomic distance and α is typically around 1.
    
    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    resolution : int
        Bin resolution (bp)
    max_distance : int, optional
        Maximum distance to include (bp)
        
    Returns
    -------
    distances : np.ndarray
        Distance values (bp)
    contacts : np.ndarray
        Average contact frequency at each distance
        
    Examples
    --------
    >>> from cfizz.io import read_cooler
    >>> reader = read_cooler("sample.mcool")
    >>> matrix = reader.fetch("chr1", 0, 50000000)
    >>> distances, contacts = calculate_distance_decay(matrix, 100000)
    """
    n = matrix.shape[0]
    
    if max_distance is None:
        max_distance = n * resolution
    
    max_distance_bins = min(max_distance // resolution, n - 1)
    
    # Calculate average contact frequency at each distance
    distances = []
    contacts = []
    
    for d in range(1, max_distance_bins):
        values = []
        for i in range(n - d):
            values.append(matrix[i, i + d])
        
        if values:
            distances.append(d * resolution)
            contacts.append(np.mean(values))
    
    return np.array(distances), np.array(contacts)


def fit_power_law(
    distances: np.ndarray,
    contacts: np.ndarray
) -> Tuple[float, float, float]:
    """
    Fit power law to distance decay data.
    
    Model: P(s) = A * s^(-α)
    
    Parameters
    ----------
    distances : np.ndarray
        Distance values
    contacts : np.ndarray
        Contact frequencies
        
    Returns
    -------
    A : float
        Prefactor
    alpha : float
        Power law exponent
    R2 : float
        R-squared value for fit
        
    Examples
    --------
    >>> dist, contacts = calculate_distance_decay(matrix, 100000)
    >>> A, alpha, R2 = fit_power_law(dist, contacts)
    >>> print(f"Power law exponent: {alpha:.3f}")
    """
    # Log-transform for linear fitting
    log_dist = np.log10(distances)
    log_contacts = np.log10(contacts)
    
    # Linear regression: log(P) = log(A) - α * log(s)
    n = len(log_dist)
    sum_x = np.sum(log_dist)
    sum_y = np.sum(log_contacts)
    sum_xy = np.sum(log_dist * log_contacts)
    sum_x2 = np.sum(log_dist ** 2)
    
    alpha = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x ** 2)
    log_A = (sum_y - alpha * sum_x) / n
    A = 10 ** log_A
    
    # Calculate R-squared
    y_mean = np.mean(log_contacts)
    SS_tot = np.sum((log_contacts - y_mean) ** 2)
    SS_res = np.sum((log_contacts - log_A + alpha * log_dist) ** 2)
    R2 = 1 - SS_res / SS_tot if SS_tot > 0 else 0
    
    return A, alpha, R2


def fit_exponential(
    distances: np.ndarray,
    contacts: np.ndarray
) -> Tuple[float, float, float]:
    """
    Fit exponential decay to distance decay data.
    
    Model: P(s) = A * exp(-s/λ)
    
    Parameters
    ----------
    distances : np.ndarray
        Distance values
    contacts : np.ndarray
        Contact frequencies
        
    Returns
    -------
    A : float
        Prefactor
    lambda_val : float
        Decay length
    R2 : float
        R-squared value for fit
    """
    def exp_func(s, A, lam):
        return A * np.exp(-s / lam)
    
    try:
        popt, _ = curve_fit(exp_func, distances, contacts, p0=[1, 1000000], maxfev=5000)
        A, lambda_val = popt
        
        # Calculate R-squared
        y_mean = np.mean(contacts)
        SS_tot = np.sum((contacts - y_mean) ** 2)
        predicted = exp_func(distances, A, lambda_val)
        SS_res = np.sum((contacts - predicted) ** 2)
        R2 = 1 - SS_res / SS_tot if SS_tot > 0 else 0
        
        return A, lambda_val, R2
    except:
        return 0, 0, 0


def calculate_decay_rate(
    contacts: np.ndarray,
    distances: np.ndarray,
    window_size: int = 10
) -> float:
    """
    Calculate local decay rate.
    
    Parameters
    ----------
    contacts : np.ndarray
        Contact frequencies
    distances : np.ndarray
        Distance values
    window_size : int
        Window size for local calculation
        
    Returns
    -------
    decay_rate : float
        Local decay rate (exponent)
    """
    if len(contacts) < window_size:
        return 0
    
    # Calculate decay rate in sliding windows
    rates = []
    
    for i in range(len(contacts) - window_size):
        window_dist = distances[i:i+window_size]
        window_contacts = contacts[i:i+window_size]
        
        try:
            _, alpha, _ = fit_power_law(window_dist, window_contacts)
            rates.append(alpha)
        except:
            pass
    
    if rates:
        return np.median(rates)
    return 0


def calculate_decay_profile(
    matrix: np.ndarray,
    resolution: int,
    bin_size: int = 1000000,
    max_distance: int = 20000000
) -> pd.DataFrame:
    """
    Calculate comprehensive decay profile.
    
    Parameters
    ----------
    matrix : np.ndarray
        Hi-C contact matrix
    resolution : int
        Bin resolution (bp)
    bin_size : int
        Bin size for output (bp)
    max_distance : int
        Maximum distance (bp)
        
    Returns
    -------
    profile : pd.DataFrame
        DataFrame with decay statistics
    """
    distances, contacts = calculate_distance_decay(matrix, resolution, max_distance)
    
    if len(distances) == 0:
        return pd.DataFrame(columns=['distance', 'contact', 'log_distance', 'log_contact'])
    
    # Bin by distance
    bin_edges = np.arange(0, max_distance + bin_size, bin_size)
    binned_distances = []
    binned_contacts = []
    
    for i in range(len(bin_edges) - 1):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i + 1])
        if np.any(mask):
            binned_distances.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            binned_contacts.append(np.mean(contacts[mask]))
    
    profile = pd.DataFrame({
        'distance': binned_distances,
        'contact': binned_contacts,
        'log_distance': np.log10(binned_distances),
        'log_contact': np.log10(binned_contacts)
    })
    
    # Fit power law
    if len(binned_distances) > 2:
        A, alpha, R2 = fit_power_law(
            np.array(binned_distances),
            np.array(binned_contacts)
        )
        profile['power_law_A'] = A
        profile['power_law_alpha'] = alpha
        profile['power_law_R2'] = R2
    
    return profile


def compare_decay_rates(
    matrix1: np.ndarray,
    matrix2: np.ndarray,
    resolution: int,
    max_distance: int = 20000000
) -> Dict[str, Any]:
    """
    Compare distance decay rates between two matrices.
    
    Parameters
    ----------
    matrix1 : np.ndarray
        First contact matrix
    matrix2 : np.ndarray
        Second contact matrix
    resolution : int
        Bin resolution (bp)
    max_distance : int
        Maximum distance (bp)
        
    Returns
    -------
    comparison : dict
        Dictionary with comparison metrics
    """
    dist1, contacts1 = calculate_distance_decay(matrix1, resolution, max_distance)
    dist2, contacts2 = calculate_distance_decay(matrix2, resolution, max_distance)
    
    _, alpha1, R2_1 = fit_power_law(dist1, contacts1)
    _, alpha2, R2_2 = fit_power_law(dist2, contacts2)
    
    return {
        'alpha1': alpha1,
        'alpha2': alpha2,
        'alpha_diff': alpha2 - alpha1,
        'R2_1': R2_1,
        'R2_2': R2_2
    }
