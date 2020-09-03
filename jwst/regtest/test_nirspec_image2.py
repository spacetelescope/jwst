import pytest

from astropy.io.fits.diff import FITSDiff
import numpy as np

import jwst.datamodels as dm
from jwst.flatfield import FlatFieldStep
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

@pytest.mark.bigdata
def test_nirspec_image2(_jail, rtdata, fitsdiff_default_kwargs):
    rtdata.get_data("nirspec/imaging/jw84600010001_02102_00001_nrs2_rate.fits")

    collect_pipeline_cfgs("config")
    args = ["config/calwebb_image2.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw84600010001_02102_00001_nrs2_cal.fits"

    rtdata.get_truth("truth/test_nirspec_image2/jw84600010001_02102_00001_nrs2_cal.fits")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    data = rtdata.get_data('nirspec/imaging/usf_assign_wcs.fits')
    user_supplied_flat = rtdata.get_data('nirspec/imaging/usf_flat.fits')

    data_flat_fielded = FlatFieldStep.call(data, user_supplied_flat=user_supplied_flat)
    rtdata.output = 'flat_fielded_step_user_supplied.fits'
    data_flat_fielded.write(rtdata.output)

    rtdata.get_truth('truth/test_nirspec_image2/flat_fielded_step_user_supplied.fits')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    data = dm.open(rtdata.get_data('nirspec/imaging/usf_assign_wcs.fits'))

    flatted = FlatFieldStep.call(data)
    unflatted = FlatFieldStep.call(flatted, inverse=True)

    assert np.allclose(data.data, unflatted.data), 'Inversion failed'


@pytest.mark.bigdata
def test_correction_pars(rtdata, fitsdiff_default_kwargs):
    """Test use of correction parameters"""
    data = dm.open(rtdata.get_data('nirspec/imaging/usf_assign_wcs.fits'))

    # First use of FlatFieldStep will store the correction.
    # The next use will use that correction
    step = FlatFieldStep()
    flatted = step.run(data)
    assert step.correction_pars['flat'] is not None

    step.use_correction_pars = True
    reflatted = step.run(data)

    assert np.allclose(flatted.data,reflatted.data), 'Re-run with correction parameters failed'
