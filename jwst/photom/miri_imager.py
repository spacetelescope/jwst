import logging
import numpy as np

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Routine to find a time-dependent correction to the PHOTOM value.

def time_corr_photom(param, t):
    """
    Short Summary
    --------------
    Time dependent PHOTOM function.

    The model parameters are amplitude, tau, t0. t0 is the reference day 
    from which the time-dependent parameters were derived. This function will return 
    a correction to apply to the PHOTOM value at a given MJD.

    Parameters
    ----------
    param : numpy array
        Set of parameters for the PHOTOM value
    t : int
        Modified Julian Day (MJD) of the observation

    Returns
    -------
    corr: float
        The time-dependent correction to the photmjsr term.
    """
    
    amplitude, tau, t0 = param["amplitude"], param["tau"], param["t0"]
    corr = amplitude * np.exp(-(t - t0)/tau)

    return corr
