import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.lib.set_telescope_pointing import update_mt_kwds
from jwst import datamodels
from jwst.stpipe import Step


@pytest.mark.bigdata
def test_nircam_image_moving_target_i2d(rtdata, fitsdiff_default_kwargs):
    """Test resampled i2d of moving target exposures for NIRCam imaging"""
    rtdata.get_asn("nircam/image/mt_asn.json")
    rtdata.output = "mt_assoc_i2d.fits"
    args = ["calwebb_image3", rtdata.input]
    Step.from_cmdline(args)
    rtdata.get_truth("truth/test_nircam_mtimage/mt_assoc_i2d.fits")

    fitsdiff_default_kwargs["atol"] = 1e-5
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("input_file",
                         [
                             'jw00634_nrcblong_mttest_tnotinrange_uncal.fits',
                             'jw00634_nrcblong_no_mtt_uncal.fits',
                             'jw00634_nrcblong_mttest_uncal.fits',
                         ],
                         ids=["midpt_not_in_mt_table_range", "no_mt_table", "with_mt_table"]
                         )
@pytest.mark.bigdata
def test_nircam_image_moving_target_kwds(input_file, rtdata, fitsdiff_default_kwargs):
    """Tests for moving target table nkeyword additions"""

    # Get the input file
    rtdata.get_data(f"nircam/image/{input_file}")

    # The add_mt_kwds function overwrites its input, so output = input
    rtdata.output = rtdata.input

    with datamodels.open(rtdata.output) as model:
        update_mt_kwds(model)
        # since the model is updated in place we need to write out the update
        model.write(rtdata.input)

    rtdata.get_truth(f"truth/test_nircam_mtimage/{input_file}")

    # Compare the results and the truth
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
