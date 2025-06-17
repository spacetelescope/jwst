"""Regression tests for NIRSpec IFU"""

import warnings

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.flatfield import FlatFieldStep
from jwst.flatfield.flat_field import nirspec_ifu
from jwst.pathloss import PathLossStep
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

# Define artifactory source and truth
INPUT_PATH = "nirspec/ifu"
TRUTH_PATH = "truth/test_nirspec_ifu"


@pytest.mark.bigdata
def test_nirspec_ifu_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test using predefined interpolated flat"""
    basename = "jw01251004001_03107_00001_nrs1"
    output_file = f"{basename}_flat_from_user_model.fits"
    with dm.open(rtdata.get_data(f"nirspec/ifu/{basename}_assign_wcs.fits")) as data:
        with dm.open(
            rtdata.get_data(f"nirspec/ifu/{basename}_interpolatedflat.fits")
        ) as user_supplied_flat:
            # Call the flat field function directly with a user flat
            nirspec_ifu(data, None, None, None, None, user_supplied_flat=user_supplied_flat)

    rtdata.output = output_file
    data.save(rtdata.output)

    rtdata.get_truth(TRUTH_PATH + "/" + output_file)
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    basename = "jw01251004001_03107_00001_nrs1"
    output_file = f"{basename}_flat_from_user_file.fits"

    data = rtdata.get_data(f"nirspec/ifu/{basename}_assign_wcs.fits")
    user_supplied_flat = rtdata.get_data(f"nirspec/ifu/{basename}_interpolatedflat.fits")

    # Call the step with a user flat
    data_flat_fielded = FlatFieldStep.call(
        data, user_supplied_flat=user_supplied_flat, save_results=False
    )

    rtdata.output = output_file
    data_flat_fielded.save(rtdata.output)
    del data_flat_fielded

    rtdata.get_truth(TRUTH_PATH + "/" + output_file)
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.slow
@pytest.mark.bigdata
def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    with dm.open(
        rtdata.get_data("nirspec/ifu/jw01251004001_03107_00001_nrs1_assign_wcs.fits")
    ) as data:
        flatted = FlatFieldStep.call(data)
        unflatted = FlatFieldStep.call(flatted, inverse=True)

    # flat fielding may set some new NaN values - ignore these in test
    is_nan = np.isnan(unflatted.data)
    assert np.allclose(data.data[~is_nan], unflatted.data[~is_nan]), "Inversion failed"

    # make sure NaNs are only at do_not_use pixels
    assert np.all(unflatted.dq[is_nan] & dm.dqflags.pixel["DO_NOT_USE"])

    # make sure NaNs at science pixels are also NaN in error and var arrays
    assert np.all(np.isnan(unflatted.err[is_nan]))
    assert np.all(np.isnan(unflatted.var_poisson[is_nan]))
    assert np.all(np.isnan(unflatted.var_rnoise[is_nan]))
    assert np.all(np.isnan(unflatted.var_flat[is_nan]))


@pytest.mark.slow
@pytest.mark.bigdata
def test_pathloss_corrpars(rtdata):
    """Test PathLossStep using correction_pars"""
    with dm.open(
        rtdata.get_data("nirspec/ifu/jw01251004001_03107_00001_nrs1_flat_field.fits")
    ) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.use_correction_pars = True
        corrected_corrpars = pls.run(data)

    assert np.allclose(corrected.data, corrected_corrpars.data, equal_nan=True)


@pytest.mark.slow
@pytest.mark.bigdata
def test_pathloss_inverse(rtdata):
    """Test PathLossStep using correction_pars"""
    with dm.open(
        rtdata.get_data("nirspec/ifu/jw01251004001_03107_00001_nrs1_flat_field.fits")
    ) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.inverse = True
        corrected_inverse = pls.run(corrected)
        non_nan = ~np.isnan(corrected_inverse.data)

    assert np.allclose(corrected.data[non_nan], corrected_inverse.data[non_nan])


@pytest.mark.slow
@pytest.mark.bigdata
def test_pathloss_source_type(rtdata):
    """Test PathLossStep forcing source type"""
    with dm.open(
        rtdata.get_data("nirspec/ifu/jw01251004001_03107_00001_nrs1_flat_field.fits")
    ) as data:
        pls = PathLossStep()
        pls.source_type = "extended"
        pls.run(data)

    assert np.allclose(
        pls.correction_pars.data, pls.correction_pars.pathloss_uniform, equal_nan=True
    )
