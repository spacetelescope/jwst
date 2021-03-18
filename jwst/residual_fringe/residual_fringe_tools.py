import os
import gc
#
import numpy as np
import math
import numpy.polynomial.polynomial as poly
from scipy.signal import savgol_filter
import scipy.signal as signal
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.interpolate import pchip
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp
#
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.modeling import models
from astropy.timeseries import LombScargle
#from specutils import Spectrum1D
#from specutils.manipulation import FluxConservingResampler, SplineInterpolatedResampler
from jwst import datamodels
#from BayesicFitting import SplinesModel, Fitter, SineModel, LevenbergMarquardtFitter, \
#    PolynomialModel, RobustShell, SineSplineModel, SineAmpModel
from numpy.linalg.linalg import LinAlgError
#

# hard coded parameters, have been selected based on testing but can be changed
NUM_KNOTS = 80  # number of knots for bkg model if no other info provided


def find_nearest(array, value):
    """ Utility function to find the index of pixel with value in 'array' nearest to 'value'.

    Used by det_pixel_trace function
    """
    idx = (np.abs(array-value)).argmin()
    return idx


def slice_info(slice_map, c):


    slice_inventory = np.unique(slice_map)
    slices_in_band = slice_inventory[np.where((slice_inventory >= 100 * c)
                                              & (slice_inventory < 100 * (c + 1)))]
    print('slices in band',slices_in_band.shape)

    slice_x_ranges = np.zeros((slices_in_band.shape[0], 3))
    all_slice_masks = np.zeros((slices_in_band.shape[0], slice_map.shape[0], slice_map.shape[1]))
    for n, s in enumerate(slices_in_band):
        # create a mask of the slice
        pixels = np.where(slice_map == s)
        slice = np.zeros(slice_map.shape)
        print('slice shape',slice.shape)
        slice[pixels] = 1

        # add this to the all_slice_mask array
        all_slice_masks[n] = slice

        # get the indices at the start and end of the slice
        collapsed_slice = np.sum(slice, axis=0)
        indices = np.where(collapsed_slice[:-1] != collapsed_slice[1:])[0]
        slice_x_ranges[n, 0], slice_x_ranges[n, 1], slice_x_ranges[n, 2] = s, np.amin(indices), np.amax(indices) + 1

        print('slice_x_ranges',n,s,slice_x_ranges[n,0],slice_x_ranges[n,1],slice_x_ranges[n,2])

    print(np.amin(slice_x_ranges[:,1]), np.amax(slice_x_ranges[:,2]))
    xrange_channel = np.zeros(2)
    xrange_channel[0] = np.amin(slice_x_ranges[:,1])
    xrange_channel[1] = np.amax(slice_x_ranges[:,2])
    
    result = (slices_in_band, xrange_channel, all_slice_masks)
    return result
