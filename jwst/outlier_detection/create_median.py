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
import numpy as np

from stsci.image import numcombine
from stsci.imagestats import ImageStats


def do_median(drizzle_groups_sci, drizzle_groups_wht, **pars):
    # start by interpreting input parameters
    nlow = pars.get('nlow',0)
    nhigh = pars.get('nhigh',0)
    high_threshold = pars.get('hthresh',None)
    low_threshold = pars.get('lthresh',None)
    nsigma = pars.get('nsigma','4 3')
    maskpt = pars.get('maskpt',0.7)
    
    # Perform additional interpretation of some parameters
    sigmaSplit=nsigma.split()
    nsigma1 = float(sigmaSplit[0])
    nsigma2 = float(sigmaSplit[1])

    if high_threshold is not None and (high_threshold.strip() == "" or high_threshold < 0): 
        high_threshold = None
    if low_threshold is not None and (low_threshold.strip() == "" or low_threshold < 0): 
        low_threshold = None
    
    if high_threshold is not None: high_threshold = float(high_threshold)
    if low_threshold is not None: low_threshold = float(low_threshold)

    
    _weight_mask_list = []
    
    for weight_arr in drizzle_groups_wht:
        # Initialize an output mask array to ones
        # This array will be reused for every output weight image
        _weight_mask = np.zeros(weight_arr.shape,dtype=np.uint8)
        try:
            tmp_mean_value = ImageStats(weight_arr, lower=1e-8,
                fields="mean", nclip=0).mean
        except ValueError:
            tmp_mean_value = 0.0
        _wht_mean = tmp_mean_value * maskpt
        # 0 means good, 1 means bad here...
        np.putmask(_weight_mask, np.less(weight_arr,_wht_mean), 1)
        #_weight_mask.info()
        _weight_mask_list.append(_weight_mask)

    # Create the combined array object using the numcombine task
    result = numcombine.numCombine(drizzle_groups_sci,
                            numarrayMaskList=_weight_mask_list,
                            combinationType="median",
                            nlow=nlow,
                            nhigh=nhigh,
                            upper=high_threshold,
                            lower=low_threshold
                        )
    median_array = result.combArrObj
    
    del _weight_mask_list
    
    return median_array
    
