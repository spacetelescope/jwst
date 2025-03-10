"""Unit tests for nrm_core module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from jwst.ami.nrm_core import FringeFitter
from jwst.ami.instrument_data import NIRISS


def test_fringe_fitter(example_model, nrm_model, bandpass):
    
    filt = example_model.meta.instrument.filter 
    niriss = NIRISS(filt, nrm_model, bandpass=bandpass,)

    # Create a FringeFitter instance with all the default options., then fit the fringes
    fitter = FringeFitter(niriss)
    output_model, output_model_multi, lgfit = fitter.fit_fringes_all(example_model)

    # print(output_model, output_model_multi, lgfit)

    # # output_model is an oifits model.
    # print(output_model.array)
    # print(output_model.target)
    # print(output_model.vis)
    # print(output_model.vis2)
    # print(output_model.t3)
    # print(output_model.wavelength)

    # # output_model_multi is an oifits model
    # print(output_model_multi.array)
    # print(output_model_multi.target)
    # print(output_model_multi.vis)
    # print(output_model_multi.vis2)
    # print(output_model_multi.t3)
    # print(output_model_multi.wavelength)

    # lgfit is an lgfit model

    # test all the initialization options