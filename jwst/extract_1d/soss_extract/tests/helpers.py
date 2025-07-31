"""Test data and functions useful for testing."""

PWCPOS = 245.85932900002442
DATA_SHAPE = (25, 200)
WAVE_BNDS_O1 = [2.8, 0.8]
WAVE_BNDS_O2 = [1.4, 0.5]
WAVE_BNDS_GRID = [0.7, 2.7]
ORDER1_SCALING = 20.0
ORDER2_SCALING = 2.0
TRACE_END_IDX = [DATA_SHAPE[1], 180]
SPECTRAL_SLOPE = 2

__all__ = ["f_lam"]


def f_lam(wl, m=SPECTRAL_SLOPE, b=0):
    """
    Return a linear model of flux as function of wavelength.

    Parameters
    ----------
    wl : array[float]
        Wavelength
    m : float
        Slope of linear model
    b : float
        Intercept of linear model

    Returns
    -------
    array[float]
        Flux
    """
    return m * wl + b
