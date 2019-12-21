import os

import pytest
from astropy.io.fits.diff import FITSDiff
from astropy.table import Table, setdiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipelines(jail, rtdata_module):
    """Run stages 1-3 pipelines on NIRCam imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("nircam/image/jw42424001001_01101_00001_nrca5_uncal.fits")

    collect_pipeline_cfgs("config")

    # Run detector1 pipeline only on one of the _uncal files
    args = ["config/calwebb_detector1.cfg", rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.jump.save_results=True",
        "--steps.jump.rejection_threshold=50.0",
        ]
    Step.from_cmdline(args)

    # And run image2 pipeline on the produced _rate file
    rtdata.input = "jw42424001001_01101_00001_nrca5_rate.fits"
    args = ["config/calwebb_image2.cfg", rtdata.input,
        "--steps.assign_wcs.save_results=True",
        ]
    Step.from_cmdline(args)

    # Grab rest of _rate files for the asn and run image2 pipeline on each to
    # produce fresh _cal files for the image3 pipeline.  We won't check these
    # or look at intermediate products
    rate_files = [
    "nircam/image/jw42424001001_01101_00001_nrcb5_rate.fits",
    "nircam/image/jw42424001001_01101_00002_nrca5_rate.fits",
    "nircam/image/jw42424001001_01101_00002_nrcb5_rate.fits",
    "nircam/image/jw42424001001_01101_00003_nrca5_rate.fits",
    "nircam/image/jw42424001001_01101_00003_nrcb5_rate.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["config/calwebb_image2.cfg", rtdata.input,
            "--steps.resample.skip=True"]
        Step.from_cmdline(args)

    # Get the level3 assocation json file (though not its members) and run
    # image3 pipeline on all _cal files listed in association
    rtdata.get_data("nircam/image/jw42424-o002_20191220t214154_image3_001_asn.json")
    args = ["config/calwebb_image3.cfg", rtdata.input,
        # Comment out following lines, as the dataset is currently broken
        # "--steps.tweakreg.save_results=True",
        # "--steps.skymatch.save_results=True",
        ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["dq_init", "saturation", "superbias",
    "refpix", "linearity", "trapsfilled", "dark_current", "jump", "rate",
    "assign_wcs", "cal", "i2d"])
def test_nircam_image_stages12(run_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 and image2 pipelines performed on NIRCam data."""
    rtdata = run_pipelines
    rtdata.input = "jw42424001001_01101_00001_nrca5_uncal.fits"
    output = "jw42424001001_01101_00001_nrca5_" + suffix + ".fits"
    rtdata.output = output
    assert os.path.exists(rtdata.output)

    rtdata.get_truth("truth/test_nircam_image_stages/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_image_stage3_i2d(run_pipelines, fitsdiff_default_kwargs):
    rtdata = run_pipelines
    rtdata.input = "jw42424-o002_20191220t214154_image3_001_asn.json"
    rtdata.output = "jw42424-o002_t001_nircam_clear-f444w_i2d.fits"
    rtdata.get_truth("truth/test_nircam_image_stages/jw42424-o002_t001_nircam_clear-f444w_i2d.fits")

    fitsdiff_default_kwargs['ignore_fields'] = ['date', 'filename']
    fitsdiff_default_kwargs['ignore_keywords'] += ['naxis1', 'tform*']
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_nircam_image_stage3_catalog(run_pipelines):
    rtdata = run_pipelines
    rtdata.input = "jw42424-o002_20191220t214154_image3_001_asn.json"
    rtdata.output = "jw42424-o002_t001_nircam_clear-f444w_cat.ecsv"
    rtdata.get_truth("truth/test_nircam_image_stages/jw42424-o002_t001_nircam_clear-f444w_cat.ecsv")

    t = Table.read(rtdata.output)
    tt = Table.read(rtdata.truth)

    # Compare the first 3 columns only, as the RA/DEC columns cannot be sorted
    # and thus setdiff cannot work on the whole table
    table = Table([t[col] for col in ['id', 'xcentroid', 'ycentroid']])
    table_truth = Table([tt[col] for col in ['id', 'xcentroid', 'ycentroid']])

    # setdiff returns a table of length zero if there is no difference
    assert len(setdiff(table, table_truth)) == 0
