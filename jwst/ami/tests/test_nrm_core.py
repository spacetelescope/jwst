"""
Unit tests for nrm_core and oifits modules.

Since nrm_core makes oifits objects, we also test oifits here.
"""

import pytest
import numpy as np
from scipy.signal import convolve
from numpy.testing import assert_allclose
import stdatamodels.jwst.datamodels as dm

from jwst.ami.nrm_core import FringeFitter
from jwst.ami.instrument_data import NIRISS
from jwst.ami.bp_fix import filtwl_d


def test_fringe_fitter(example_model, nrm_model, bandpass, nrm_psf):
    """
    Test generating a fringe fitter object, then using it to fit example data.

    This touches all the methods of FringeFitter, as well as creating an OIModel instance.
    """

    filt = example_model.meta.instrument.filter
    niriss = NIRISS(
        filt,
        nrm_model,
        bandpass,
    )

    # Need data to be convolved with the PSF
    sci_data = example_model.data[0]
    sci_data_conv = convolve(sci_data, nrm_psf, mode="same")
    example_model.data[0] = sci_data_conv
    example_model.data[1] = sci_data_conv

    # Create a FringeFitter instance with all the default options., then fit the fringes
    fitter = FringeFitter(niriss)
    output_model, output_model_multi, lgfit = fitter.fit_fringes_all(example_model)

    # output_model is an oifits model.
    assert isinstance(output_model, dm.AmiOIModel)
    assert isinstance(output_model_multi, dm.AmiOIModel)
    for i, oimodel in enumerate([output_model, output_model_multi]):
        assert len(oimodel.vis) == 21  # number of baselines: 7 holes choose 2
        assert len(oimodel.vis2) == 21
        assert len(oimodel.t3) == 35  # number of triples: 7 holes choose 3
        vis_amp = oimodel.vis["VISAMP"]
        vis_amp_err = oimodel.vis["VISAMPERR"]
        vis_phase = oimodel.vis["VISPHI"]
        vis_phase_err = oimodel.vis["VISPHIERR"]
        vis2 = oimodel.vis2["VIS2DATA"]
        vis2_err = oimodel.vis2["VIS2ERR"]
        t3_amp = oimodel.t3["T3AMP"]
        t3_amp_err = oimodel.t3["T3AMPERR"]
        t3_phase_err = oimodel.t3["T3PHIERR"]

        for arr in [
            vis_amp,
            vis_amp_err,
            vis_phase,
            vis_phase_err,
            vis2,
            vis2_err,
            t3_amp,
            t3_amp_err,
            t3_phase_err,
        ]:
            assert arr.dtype == np.float64
            if i == 0:
                assert arr.ndim == 1  # the mean data
            if i == 1:
                assert arr.ndim == 2
                assert arr.shape[1] == example_model.data.shape[0]  # one value per integration

        eff_wave, eff_bandwidth = oimodel.wavelength[0]
        assert np.isclose(eff_wave, filtwl_d[filt])

    # lgfit is an AmiLgFitModel
    assert isinstance(lgfit, dm.AmiLgFitModel)
    for att_str in [
        "centered_image",
        "norm_centered_image",
        "fit_image",
        "norm_fit_image",
        "resid_image",
        "norm_resid_image",
    ]:
        assert hasattr(lgfit, att_str)
        im = getattr(lgfit, att_str)
        assert isinstance(im, np.ndarray)
        assert im.ndim == 3

    coeffs = lgfit.solns_table["coeffs"]
    # identical because input data are identical in both planes
    # Why is the shape hard-coded to 44?
    assert coeffs.shape == (example_model.data.shape[0], 44)
    assert np.allclose(coeffs[0], coeffs[1])
