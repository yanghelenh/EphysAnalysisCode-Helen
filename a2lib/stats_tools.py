"""
Statistics tools
"""
import numpy as np


def mad(x):
    """Calculate the median absolute deviation (MAD), a robust statistic

    Parameters
    ----------
    x : np.array
        Input to calculate the MAD on

    Returns
    -------
    x_mad : float
        median absolute deviation
    """
    return np.median(np.abs(x - np.median(x, axis=0)), axis=0) * 1.4826


def mad_standardize(x):
    """ Standardize to the median absolute deviation

    Parameters
    ----------
    x : np.array
        Input to median subtract and divide by the MAD

    Returns
    -------
    x_mad : np.array
        MAD standardized array
    """
    return (x - np.median(x)) / mad(x)


def conditional_zscore(x, inds_or_bool):
    """
    Parameters
    ----------
    x : np.array
        Array to be standardized

    inds_or_bool : bool, np.array
        Index or boolean over which to calculate standard deviation (for standardization of x)

    Returns
    -------
    x_zscored : np.array
        Standardized version of x

    """
    return x / np.std(x[inds_or_bool] - np.mean(x[inds_or_bool]))

