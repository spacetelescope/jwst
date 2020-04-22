import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.lib.set_telescope_pointing import update_mt_kwds
from jwst import datamodels

from jwst.stpipe import Step


@pytest.mark.bigdata
def test_nircam_image_moving_target(rtdata, fitsdiff_default_kwargs):
    """Test resampled i2d of moving target exposures for NIRCam imaging"""
    collect_pipeline_cfgs("config")
    rtdata.get_asn("nircam/image/mt_asn.json")
    rtdata.output = "mt_assoc_i2d.fits"
    args = ["config/calwebb_image3.cfg", rtdata.input]
    Step.from_cmdline(args)
    rtdata.get_truth("truth/test_nircam_mtimage/mt_assoc_i2d.fits")

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

@pytest.mark.bigdata
def test_nircam_image_moving_target_kwds_oot(rtdata, fitsdiff_default_kwargs):
    """Test that if the requested time does not overlap the times in the moving
       target table the input data are not altered."""

    # Get the input level-1b file
    rtdata.get_data("nircam/image/jw00634_nrcblong_mttest_tnotinrange_uncal.fits")

    # The add_mt_kwds function overwrites its input, so output = input
    rtdata.output = rtdata.input

    with datamodels.open(rtdata.output) as model:
        update_mt_kwds(model)
        # since the model is updated in place we need to write out the update
        model.write(rtdata.input)

    rtdata.get_truth("truth/test_nircam_mtimage/jw00634_nrcblong_mttest_tnotinrange_uncal.fits")

    # Compare the results and the truth
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_no_moving_target_table_kwds(rtdata, fitsdiff_default_kwargs):
    """Test that if the input data does not have a moving target table the
       input data are returned without the MT_RA & MT_DEC kwywords."""

    # Get the input level-1b file
    rtdata.get_data("nircam/image/jw00634_nrcblong_no_mtt_uncal.fits")

    # The add_mt_kwds function overwrites its input, so output = input
    rtdata.output = rtdata.input

    with datamodels.open(rtdata.output) as model:
        update_mt_kwds(model)
        # since the model is updated in place we need to write out the update
        model.write(rtdata.input)

    rtdata.get_truth("truth/test_nircam_mtimage/jw00634_nrcblong_no_mtt_uncal.fits")

    # Compare the results and the truth
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

@pytest.mark.bigdata
def test_nircam_moving_target_table_kwds(rtdata, fitsdiff_default_kwargs):
    """Test that if the input data does not have a moving target table the
       input data are returned without the MT_RA & MT_DEC kwywords."""

    # Get the input level-1b file
    rtdata.get_data("nircam/image/jw00634_nrcblong_mttest_uncal.fits")

    # The add_mt_kwds function overwrites its input, so output = input
    rtdata.output = rtdata.input

    with datamodels.open(rtdata.input) as model:
        update_mt_kwds(model)
        # since the model is updated in place we need to write out the update
        model.write(rtdata.input)

    rtdata.get_truth("truth/test_nircam_mtimage/jw00634_nrcblong_mttest_uncal.fits")

    # Compare the results and the truth
    fitsdiff_default_kwargs['ignore_keywords'].append('HISTORY')
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
