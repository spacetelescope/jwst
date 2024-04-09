import pytest
import numpy as np
from jwst.wfss_contam.disperse import dispersed_pixel
from jwst.wfss_contam.tests.test_observations import grism_wcs, segmentation_map, direct_image


def test_oversample_same_result(grism_wcs, segmentation_map):
    '''
    Coverage for bug where wavelength oversampling led to double-counted fluxes

    note: segmentation_map fixture needs to be able to find module-scoped direct_image
    fixture, so it must be imported here
    '''

    # manual input of input params set the same as test_observations.py
    x0 = 300.5
    y0 = 300.5
    order = 1
    width = 1.0
    height = 1.0
    lams = [2.0]
    flxs = [1.0]
    ID = 0
    naxis = (300, 500)
    sens_waves = np.linspace(1.708, 2.28, 100)
    wmin, wmax = np.min(sens_waves), np.max(sens_waves)
    sens_resp = np.ones(100)
    seg_wcs = segmentation_map.meta.wcs
    0, (300, 500), 2, False, 
    xoffset = 2200
    yoffset = 1000


    xs, ys, areas, lams_out, counts_1, ID = dispersed_pixel(
                    x0, y0, width, height, lams, flxs, order, wmin, wmax,
                    sens_waves, sens_resp, seg_wcs, grism_wcs, ID, naxis,
                    oversample_factor=1, extrapolate_sed=False, xoffset=xoffset,
                    yoffset=yoffset)
    
    xs, ys, areas, lams_out, counts_3, ID = dispersed_pixel(
                x0, y0, width, height, lams, flxs, order, wmin, wmax,
                sens_waves, sens_resp, seg_wcs, grism_wcs, ID, naxis,
                oversample_factor=3, extrapolate_sed=False, xoffset=xoffset,
                yoffset=yoffset)

    assert np.isclose(np.sum(counts_1), np.sum(counts_3), rtol=1e-2)