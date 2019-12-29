import pytest
from astropy.io.fits.diff import FITSDiff
from astropy.table import Table, setdiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.fixture(scope="module")
def run_pipelines(jail, rtdata_module):
    """Run stage 1-3 pipelines on MIRI imaging data."""
    rtdata = rtdata_module
    rtdata.get_data("miri/image/det_image_1_MIRIMAGE_F770Wexp1_5stars_uncal.fits")

    collect_pipeline_cfgs("config")

    # Run detector1 pipeline only on one of the _uncal files
    args = ["config/calwebb_detector1.cfg", rtdata.input,
        "--save_calibrated_ramp=True",
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.rscd.save_results=True",
        "--steps.lastframe.save_results=True",
        "--steps.firstframe.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.jump.rejection_threshold=10.0",
        ]
    Step.from_cmdline(args)

    # Now run image2 pipeline on the _rate file, saving intermediate products
    rtdata.input = 'det_image_1_MIRIMAGE_F770Wexp1_5stars_rate.fits'
    args = ["config/calwebb_image2.cfg", rtdata.input,
        "--steps.assign_wcs.save_results=True",
        "--steps.flat_field.save_results=True"
        ]
    Step.from_cmdline(args)

    # Grab rest of _rate files for the asn and run image2 pipeline on each to
    # produce fresh _cal files for the image3 pipeline.  We won't check these
    # or look at intermediate products, and skip resample (don't need i2d image)
    rate_files = [
    "miri/image/det_image_1_MIRIMAGE_F770Wexp2_5stars_rate.fits",
    "miri/image/det_image_2_MIRIMAGE_F770Wexp1_5stars_rate.fits",
    "miri/image/det_image_2_MIRIMAGE_F770Wexp2_5stars_rate.fits",
    ]
    for rate_file in rate_files:
        rtdata.get_data(rate_file)
        args = ["config/calwebb_image2.cfg", rtdata.input,
            "--steps.resample.skip=True"]
        Step.from_cmdline(args)

    # Get the level3 assocation json file (though not its members) and run
    # image3 pipeline on all _cal files listed in association
    rtdata.get_data("miri/image/det_dithered_5stars_image3_asn.json")
    args = ["config/calwebb_image3.cfg", rtdata.input,
        # Set some unique param values needed for these data
        "--steps.tweakreg.snr_threshold=200",
        "--steps.tweakreg.use2dhist=False",
        "--steps.source_catalog.snr_threshold=20",
        ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["a3001_crf", "assign_wcs", "cal", "dark_current",
    "dq_init", "firstframe", "flat_field", "i2d", "lastframe", "linearity", "ramp",
    "rate", "rateints", "refpix", "rscd", "saturation"])
def test_miri_image_stages12(run_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 and image2 pipelines performed on MIRI data."""
    rtdata = run_pipelines
    rtdata.input = "det_image_1_MIRIMAGE_F770Wexp1_5stars_uncal.fits"
    output = "det_image_1_MIRIMAGE_F770Wexp1_5stars_" + suffix + ".fits"
    rtdata.output = output
    #assert os.path.exists(rtdata.output)

    rtdata.get_truth("truth/test_miri_image_stages/" + output)

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_image_stage3_i2d(run_pipelines, fitsdiff_default_kwargs):
    rtdata = run_pipelines
    rtdata.input = "det_dithered_5stars_image3_asn.json"
    rtdata.output = "det_dithered_5stars_f770w_i2d.fits"
    rtdata.get_truth("truth/test_miri_image_stages/det_dithered_5stars_f770w_i2d.fits")

    fitsdiff_default_kwargs['ignore_fields'] = ['date', 'filename']
    fitsdiff_default_kwargs['ignore_keywords'] += ['naxis1', 'tform*']
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_miri_image_stage3_catalog(run_pipelines):
    rtdata = run_pipelines
    rtdata.input = "det_dithered_5stars_image3_asn.json"
    rtdata.output = "det_dithered_5stars_f770w_cat.ecsv"
    rtdata.get_truth("truth/test_miri_image_stages/det_dithered_5stars_f770w_cat.ecsv")

    t = Table.read(rtdata.output)
    tt = Table.read(rtdata.truth)

    # Compare the first 3 columns only, as the RA/DEC columns cannot be sorted
    # and thus setdiff cannot work on the whole table
    table = Table([t[col] for col in ['id', 'xcentroid', 'ycentroid']])
    table_truth = Table([tt[col] for col in ['id', 'xcentroid', 'ycentroid']])

    # setdiff returns a table of length zero if there is no difference
    assert len(setdiff(table, table_truth)) == 0
