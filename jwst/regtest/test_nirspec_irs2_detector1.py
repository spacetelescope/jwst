""" Test for the detector1 pipeline using NIRSpec data in IRS2 mode. This takes
    an uncal file and generates the stage 1 FITS files (rate) along with the
    intermediate products."""

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.stpipe import Step
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs


@pytest.fixture(scope="module")
def run_detector1pipeline(rtdata_module, jail):
    """Run calwebb_detector1 pipeline on NIRSpec data with IRS2 readout mode."""
    rtdata = rtdata_module
    rtdata.get_data("nirspec/irs2/jw0010010_11010_nrs1_chimera_uncal.fits")

    collect_pipeline_cfgs("config")
    Step.from_cmdline([
        "config/calwebb_detector1.cfg",
        rtdata.input,
        "--steps.dq_init.save_results=True",
        "--steps.saturation.save_results=True",
        "--steps.superbias.save_results=True",
        "--steps.refpix.save_results=True",
        "--steps.rscd.save_results=True",
        "--steps.linearity.save_results=True",
        "--steps.dark_current.save_results=True",
        "--steps.jump.save_results=True",
        "--steps.ramp_fit.save_results=True",
        "--steps.gain_scale.save_results=True",
        "--steps.jump.rejection_threshold=200",
    ])


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ['dq_init', 'saturation', 'superbias',
    'refpix', 'linearity', 'dark_current', 'jump', '0_ramp_fit', 'gain_scale',
    'rate'])
def test_nirspec_irs2_detector1(run_detector1pipeline, rtdata_module,
    fitsdiff_default_kwargs, suffix):
    """
    Regression test of calwebb_detector1 pipeline performed on NIRSpec IRS2 data.
    """
    rtdata = rtdata_module

    output_filename = f"jw0010010_11010_nrs1_chimera_{suffix}.fits"
    rtdata.output = output_filename
    rtdata.get_truth(f"truth/test_nirspec_irs2_detector1/{output_filename}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
