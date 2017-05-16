"""
Create a median image from the singly drizzled images.

NOTE: This version is simplified from astrodrizzle's version in the
        following ways:
        - type of combination: fixed to 'median'
        - 'minmed' not implemented as an option
        - does not use buffers to try to minimize memory usage

:Authors: Warren Hack

:License:

"""
from __future__ import (division, print_function, unicode_literals,
    absolute_import)

import numpy as np

from stsci.image import numcombine
from stsci.imagestats import ImageStats

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_median(drizzle_groups_sci, drizzle_groups_wht, **pars):
    # start by interpreting input parameters
    nlow = pars.get('nlow', 0)
    nhigh = pars.get('nhigh', 0)
    high_threshold = pars.get('hthresh', None)
    low_threshold = pars.get('lthresh', None)
    nsigma = pars.get('nsigma', '4 3')
    maskpt = pars.get('maskpt', 0.7)

    # Perform additional interpretation of some parameters
    nsigma1 = float(nsigma.split()[0])
    nsigma2 = float(nsigma.split()[1])

    if high_threshold and high_threshold < 0:
        high_threshold = None
    if low_threshold and low_threshold < 0:
        low_threshold = None

    weight_mask_list = []

    for weight_arr in drizzle_groups_wht:
        # Initialize an output mask array to ones
        # This array will be reused for every output weight image
        weight_mask = np.zeros(weight_arr.shape, dtype=np.uint8)
        try:
            tmp_mean_value = ImageStats(weight_arr, lower=1e-8,
                fields="mean", nclip=0).mean
        except ValueError:
            tmp_mean_value = 0.0
        wht_mean = tmp_mean_value * maskpt
        # 0 means good, 1 means bad here...
        np.putmask(weight_mask, np.less(weight_arr, wht_mean), 1)
        weight_mask_list.append(weight_mask)

    # Create the combined array object using the numcombine task
    result = numcombine.numCombine(drizzle_groups_sci,
        numarrayMaskList=weight_mask_list, combinationType="median", nlow=nlow,
        nhigh=nhigh, upper=high_threshold, lower=low_threshold)
    median_array = result.combArrObj

    return median_array
