""" Test for the detector1 pipeline using NIRISS image mode, starting with
    an uncal file. Results from all intermediate steps, including
    charge_migration, are saved for comparisons with truth files.
"""

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst import datamodels
from jwst.stpipe import Step
from jwst.tweakreg import TweakRegStep


@pytest.fixture(scope="module")
def run_detector1(rtdata_module):
    """Run calwebb_detector1 pipeline on NIRISS imaging data."""
    rtdata = rtdata_module

    rtdata.get_data("niriss/imaging/jw01094001002_02107_00001_nis_uncal.fits")

    # Run detector1 pipeline on an _uncal files
    args = ["calwebb_detector1", rtdata.input,
            "--steps.persistence.save_trapsfilled=False",
            "--steps.dq_init.save_results=True",
            "--steps.saturation.save_results=True",
            "--steps.superbias.save_results=True",
            "--steps.refpix.save_results=True",
            "--steps.linearity.save_results=True",
            "--steps.dark_current.save_results=True",
            "--steps.charge_migration.skip=False",
            "--steps.charge_migration.save_results=True",
            "--steps.jump.save_results=True"
            ]


    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_detector1_with_clean_flicker_noise(rtdata_module):
    """Run detector1 pipeline on NIRISS imaging data with noise cleaning."""
    rtdata_module.get_data("niriss/imaging/jw01094001002_02107_00001_nis_uncal.fits")

    # Run detector1 pipeline only on one of the _uncal files.
    # Run optional clean_flicker_noise step,
    # saving extra outputs and masking to science regions
    args = ["jwst.pipeline.Detector1Pipeline", rtdata_module.input,
            "--output_file=jw01094001002_02107_00001_nis_cfn",
            "--save_calibrated_ramp=True",
            "--steps.clean_flicker_noise.skip=False",
            "--steps.clean_flicker_noise.save_results=True",
            "--steps.clean_flicker_noise.save_mask=True",
            "--steps.clean_flicker_noise.save_background=True",
            "--steps.clean_flicker_noise.save_noise=True",
            ]
    Step.from_cmdline(args)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["dq_init", "saturation", "superbias",
                                    "refpix", "linearity", "dark_current",
                                    "charge_migration", "jump", "rate", "rateints"])
def test_niriss_image_detector1(run_detector1, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Test detector1 pipeline for NIRISS imaging data with noise cleaning.
    """
    truth_dir = 'test_niriss_image'
    _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix, truth_dir)


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", ["cfn_clean_flicker_noise", "mask",
                                    "flicker_bkg", "flicker_noise",
                                    "cfn_ramp", "cfn_rate", "cfn_rateints"])
def test_niriss_image_detector1_with_clean_flicker_noise(
        run_detector1_with_clean_flicker_noise, rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 pipeline performed on NIRISS imaging data.
    """
    truth_dir = 'test_niriss_image_clean_flicker_noise'
    _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix, truth_dir)


@pytest.mark.bigdata
def test_niriss_tweakreg_no_sources(rtdata, fitsdiff_default_kwargs):
    """Make sure tweakreg is skipped when sources are not found.
    """

    rtdata.input = "niriss/imaging/jw01537-o003_20240406t164421_image3_00004_asn.json"
    rtdata.get_asn("niriss/imaging/jw01537-o003_20240406t164421_image3_00004_asn.json")

    args = [
        "jwst.tweakreg.TweakRegStep",
        rtdata.input,
        "--abs_refcat='GAIADR3'",
        "--save_results=True",
    ]

    # run the test from the command line:
    result = Step.from_cmdline(args)

    # Check the status of the step is set correctly in the files.
    mc = datamodels.ModelContainer(rtdata.input)

    for model in mc:
        assert model.meta.cal_step.tweakreg != 'SKIPPED'

    result = TweakRegStep.call(mc)
    with result:
        for model in result:
            assert model.meta.cal_step.tweakreg == 'SKIPPED'
            result.shelve(model, modify=False)


def _assert_is_same(rtdata_module, fitsdiff_default_kwargs, suffix, truth_dir):
    """Assertion helper for the above tests"""
    rtdata = rtdata_module
    rtdata.input = "jw01094001002_02107_00001_nis_uncal.fits"
    output = f"jw01094001002_02107_00001_nis_{suffix}.fits"
    rtdata.output = output

    rtdata.get_truth(f"truth/{truth_dir}/{output}")

    # Set tolerances so the crf, rscd and rateints file comparisons work across
    # architectures
    fitsdiff_default_kwargs["rtol"] = 1e-4
    fitsdiff_default_kwargs["atol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
