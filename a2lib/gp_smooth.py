"""
Use the pymc3 library to perform gaussian process smoothing on a timeseries.

Separated out due to potential environment configuration difficulties.
"""
import pymc3 as pm
from theano import shared
from pymc3.distributions.timeseries import GaussianRandomWalk


def gp_smooth(var, smooth_ratio=0.5, prior_std=10000, progressbar=False):
    """Smooth a timeseries using a gaussian random walk, taking the MAP at each timepoint

    Parameters
    ----------
    var : np.array
        Variable to smooth

    smooth_ratio : float (optional)
        Ratio to assign between noise and walk. 0.5 works well in most cases

    prior_std : float (optional)
        The prior on random walk standard deviation, should be very large.

    progressbar : bool (optional)
        Flag to display progressbar

    Returns
    -------
    var_smooth : np.array
        GP smoothed timeseries

    """
    model = pm.Model()
    with model:
        alpha = shared(smooth_ratio)
        mu = pm.Normal("mu", sd=prior_std)
        mu_0 = pm.Normal("mu_0", sd=1/prior_std)
        tau = pm.HalfCauchy("tau", 1/prior_std)
        z = GaussianRandomWalk("z", mu=mu+mu_0, tau=tau / (1.0 - alpha), shape=var.shape)
        obs = pm.Normal("obs", mu=z, tau=tau/alpha, observed=var)

        # Get the MAP estimates
        res = pm.find_MAP(vars=[z], progressbar=progressbar)
        var_smooth = res['z']  # MAP estimate
    return var_smooth
