import os
import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import stdatamodels.jwst.datamodels as dm

from jwst.flatfield import FlatFieldStep
from jwst.lib.suffix import replace_suffix
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_tso_spec2_pipeline(jail, rtdata_module, request):
    """Run the calwebb_spec2 pipeline performed on NIRSpec
        fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit)
    """

    rtdata = rtdata_module

    # Get the input exposure
    rtdata.get_data('nirspec/tso/jw02420001001_04101_00001-seg001_nrs1_rateints.fits')

    # Run the calwebb_spec2 pipeline;
    args = ["calwebb_spec2", rtdata.input,
            "--steps.assign_wcs.save_results=True",
            "--steps.extract_2d.save_results=True",
            "--steps.wavecorr.save_results=True",
            "--steps.flat_field.save_results=True",
            "--steps.photom.save_results=True"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ['assign_wcs', 'extract_2d', 'wavecorr', 'flat_field', 'photom', 'calints', 'x1dints'])
def test_nirspec_brightobj_spec2(run_tso_spec2_pipeline, fitsdiff_default_kwargs, suffix):
    """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec
        fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit).
    """
    rtdata = run_tso_spec2_pipeline
    output = replace_suffix(
        os.path.splitext(os.path.basename(rtdata.input))[0], suffix) + '.fits'
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_brightobj_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    data = rtdata.get_data('nirspec/tso/nrs2_wavecorr.fits')
    user_supplied_flat = rtdata.get_data('nirspec/tso/nrs2_interpolatedflat.fits')

    data_flat_fielded = FlatFieldStep.call(data, user_supplied_flat=user_supplied_flat)
    rtdata.output = 'flat_fielded_step_user_supplied.fits'
    data_flat_fielded.write(rtdata.output)

    rtdata.get_truth('truth/test_nirspec_brightobj_spec2/flat_fielded_step_user_supplied.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_bots_interp_flat(rtdata, fitsdiff_default_kwargs):
    """Test the interpolated flat for a NRS BOTS exposure"""
    data = rtdata.get_data('nirspec/tso/jw93056001001_short_nrs1_wavecorr.fits')

    FlatFieldStep.call(data, save_interpolated_flat=True)
    rtdata.output = 'jw93056001001_short_nrs1_wavecorr_interpolatedflat.fits'

    rtdata.get_truth('truth/test_nirspec_brightobj_spec2/jw93056001001_short_nrs1_wavecorr_interpolatedflat.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    with dm.open(rtdata.get_data('nirspec/tso/nrs2_wavecorr.fits')) as data:
        flatted = FlatFieldStep.call(data)
        unflatted = FlatFieldStep.call(flatted, inverse=True)

    # flat fielding may set some new NaN values - ignore these in test
    is_nan = np.isnan(unflatted.data)
    assert np.allclose(data.data[~is_nan], unflatted.data[~is_nan]), 'Inversion failed'

    # make sure NaNs are only at do_not_use pixels
    assert np.all(unflatted.dq[is_nan] & dm.dqflags.pixel['DO_NOT_USE'])
