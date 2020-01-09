import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


DATASET_ID = "det_image_MIRIMAGE_P750Lexp1"


@pytest.fixture(scope="module")
def run_tso1_pipeline(jail, rtdata_module):
    """Run the calwebb_tso1 pipeline on a MIRI LRS slitless exposure."""
    rtdata = rtdata_module
    collect_pipeline_cfgs("config")
    rtdata.get_data(f"miri/lrs/{DATASET_ID}.fits")

    args = [
        "config/calwebb_tso1.cfg",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.rscd.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.jump.rejection_threshold=25",
        "--steps.jump.save_results=True"
    ]
    Step.from_cmdline(args)
    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("step_suffix", ['rate', 'rateints', 'linearity',
                                         'rscd',
                                         'dq_init',
                                     'saturation', 'dark_current', 'refpix',
                                         'jump'])
def test_miri_lrs_slitless_tso1(run_tso1_pipeline, fitsdiff_default_kwargs, step_suffix):
    """
    Regression test of calwebb_detector1 pipeline performed on MIRI data.
    """
    rtdata = run_tso1_pipeline
    output_filename = f"{DATASET_ID}_{step_suffix}.fits"
    rtdata.output = output_filename

    rtdata.get_truth(f"truth/test_miri_lrs_slitless_tso1/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
