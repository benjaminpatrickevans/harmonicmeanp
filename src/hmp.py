"""
This file intends to provide a simple python implementation of harmonic mean p-values.

The original r file is available on cran as: 'harmonicmeanp'.
Wilson, D. J. (2019). The harmonic mean p-value for combining dependent tests. Proceedings of the National Academy of Sciences, 116(4), 1195-1200.

All implementations (and errors) here are my own, with no affiliation to Wilson.
"""

import numpy as np
from scipy import integrate


def stat(p_values, weights=None):
    """
    Uses the harmonic mean to combine p-values. This is the raw harmonic mean p_value calculation,
    but using hmp() should be preferred which transforms to an asymptotically exact test.

    For the original reference, see:
       Wilson, D. J. (2019). The harmonic mean p-value for combining dependent tests. Proceedings of the National Academy of Sciences, 116(4), 1195-1200.
    
    """
    L = len(p_values)

    if weights is None:
        weights = np.asarray([1. / L] * L)  # Assume uniform

    assert len(weights) == len(p_values)
    assert np.all(weights >= 0)
    assert np.isclose(np.sum(weights), 1)  # Should sum to 1

    p_values = np.asarray(p_values)
    return 1 / np.sum(weights / p_values)


def upper_bound(p_values, weights=None):
    """
    A worst case upper bound can be estimated by noting the harmonic mean
    is never more anti-conservative than a factor of e*ln(L).
    where L=|p_values|. 

    See: Vovk, V., & Wang, R. (2020). Combining p-values via averaging. Biometrika, 107(4), 791-808.
    """
    L = len(p_values)
    assert L >= 2  # Only defined with multiple values
    harmonic_mean = stat(p_values, weights)

    # Note: as k->∞, the e term can be dropped. Todo: Implement?
    return harmonic_mean * np.e * np.log(L)


def _landau_density(x, mu, sigma):
    """
    Computes the density of the Landau distribution. Note: This could be improved for efficiency,
    for example, using pylandau. Here I have used a naive implementation to keep dependencies small.
    """
    fn = lambda t:  np.exp(-t * ((x-mu)/sigma) - (2/np.pi) * t * np.log(t)) * np.sin(2 * t)
    return 1 / (mu * sigma) * integrate.quad(fn, 0, np.inf)[0]  # Note: may need to increase upper bound for sub intervals ('limit' keyword for quad)


def hmp(p_values, weights=None):
    """
    Uses the harmonic mean to combine p-values. This is an asymptotically exact transformation of HMP,
    using the landau distribution.

    For the original reference, see:
       Wilson, D. J. (2019). The harmonic mean p-value for combining dependent tests. Proceedings of the National Academy of Sciences, 116(4), 1195-1200.
    """
    harmonic_mean = stat(p_values, weights)

    L = len(p_values)
    mu = np.log(L) + 0.874
    sigma = np.pi / 2

    landau = lambda x: _landau_density(x, mu, sigma)
    return integrate.quad(landau, 1/harmonic_mean, np.inf)[0]
