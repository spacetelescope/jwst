import os
import shutil
from pathlib import Path

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.tests.helpers import make_nirspec_ifu_rate_model
from jwst.stpipe import Step
from jwst.targ_centroid.tests.helpers import (
    make_empty_lrs_model,
    make_slit_data,
    make_ta_association,
)
from jwst.tests.helpers import _help_pytest_warns

INPUT_FILE = "test_rate.fits"
INPUT_FILE_2 = "test2_rate.fits"
INPUT_ASN = "test_asn.json"
OUTPUT_FILE = "custom_name.fits"
OUTPUT_FILE_ASN = "custom_name_asn.fits"  # cannot reuse because everything runs in same cwd
LOGFILE = "run_asn.log"


@pytest.fixture(scope="module")
def make_test_rate_file(tmp_cwd_module):
    """
    Make and save a test rate file in the temporary working directory
    Partially copied from test_background.py
    """
    image = make_nirspec_ifu_rate_model()
    with image as dm:
        dm.save(INPUT_FILE)


@pytest.fixture(scope="module")
def make_test_association(make_test_rate_file):
    shutil.copy(INPUT_FILE, INPUT_FILE_2)
    os.system(f"asn_from_list -o {INPUT_ASN} -r DMSLevel2bBase {INPUT_FILE} {INPUT_FILE_2}")


@pytest.fixture(scope="module", params=[OUTPUT_FILE])
def run_spec2_pipeline(make_test_rate_file, request):
    """
    Run pipeline, saving one intermediate step
    and skipping most of the rest for runtime
    """
    args = [
        "calwebb_spec2",
        INPUT_FILE,
        "--steps.badpix_selfcal.skip=true",
        "--steps.msa_flagging.skip=true",
        "--steps.clean_flicker_noise.skip=true",
        "--steps.flat_field.skip=true",
        "--steps.pathloss.skip=true",
        "--steps.photom.skip=true",
        "--steps.pixel_replace.skip=true",
        "--steps.cube_build.save_results=true",
        "--steps.extract_1d.skip=true",
        f"--output_file={request.param}",
    ]

    Step.from_cmdline(args)


@pytest.fixture(scope="module", params=[OUTPUT_FILE_ASN])
def run_spec2_pipeline_asn(make_test_association, request):
    """
    Two-product association passed in. This should trigger a warning
    and the output_file parameter should be ignored.
    """
    # save warnings to logfile so can be checked later
    args = [
        "calwebb_spec2",
        INPUT_ASN,
        "--log-level=INFO",
        f"--log-file={LOGFILE}",
        "--steps.badpix_selfcal.skip=true",
        "--steps.msa_flagging.skip=true",
        "--steps.clean_flicker_noise.skip=true",
        "--steps.flat_field.skip=true",
        "--steps.pathloss.skip=true",
        "--steps.photom.skip=true",
        "--steps.pixel_replace.skip=true",
        "--steps.cube_build.save_results=true",
        "--steps.extract_1d.skip=true",
        f"--output_file={request.param}",
    ]

    Step.from_cmdline(args)


def test_output_file_rename(run_spec2_pipeline):
    """
    Covers a bug where the output_file parameter was not being
    respected in calls to Spec2Pipeline.
    _s3d file expected from cube_build save_results=true
    _cal file expected from pipeline finishing
    """
    assert os.path.exists(INPUT_FILE)
    custom_stem = OUTPUT_FILE.split(".")[0]
    for extension in ["s3d", "cal"]:
        assert os.path.exists(f"{custom_stem}_{extension}.fits")


def test_output_file_norename_asn(run_spec2_pipeline_asn):
    """
    Ensure output_file parameter is ignored, with warning,
    when multiple products are in the same association.
    """
    # ensure tmp_cwd_module is successfully keeping all files in cwd
    assert os.path.exists(INPUT_ASN)
    assert os.path.exists(INPUT_FILE)
    assert os.path.exists(INPUT_FILE_2)

    custom_stem = OUTPUT_FILE_ASN.split(".")[0]
    input_stem = INPUT_FILE.split("_")[0]
    input_stem_2 = INPUT_FILE_2.split("_")[0]

    # ensure default filenames were written, and not the custom one
    for extension in ["s3d", "cal"]:
        assert not os.path.exists(f"{custom_stem}_{extension}.fits")
        assert os.path.exists(f"{input_stem}_{extension}.fits")
        assert os.path.exists(f"{input_stem_2}_{extension}.fits")

    # ensure warning goes to log file
    assert os.path.exists(LOGFILE)
    with open(LOGFILE, "r") as f:
        log_content = f.read()

    assert (
        "Multiple products in input association. Output file name will be ignored." in log_content
    )


def test_filenotfounderror_raised(capsys):
    # Verify the failure is in the traceback message
    with pytest.raises(RuntimeError, match="FileNotFoundError"):
        Spec2Pipeline().run("file_does_not_exist.fits")

    # Verify the failure is printed to stderr
    captured = capsys.readouterr()
    assert "FileNotFoundError" in captured.err


def test_bsub_deprecated(make_test_rate_file):
    """
    Ensure that the deprecated save_bsub parameter raises a
    DeprecationWarning when set to True.
    """
    args = [
        "calwebb_spec2",
        INPUT_FILE,
        "--save_bsub=true",
        "--steps.badpix_selfcal.skip=true",
        "--steps.msa_flagging.skip=true",
        "--steps.clean_flicker_noise.skip=true",
        "--steps.flat_field.skip=true",
        "--steps.pathloss.skip=true",
        "--steps.photom.skip=true",
        "--steps.pixel_replace.skip=true",
        "--steps.extract_1d.skip=true",
    ]

    with (
        _help_pytest_warns(),
        pytest.warns(DeprecationWarning, match="The --save_bsub parameter is deprecated"),
    ):
        Step.from_cmdline(args)


@pytest.mark.parametrize("use_asn", [True, False])
def test_targ_centroid_logic(use_asn, tmp_cwd):
    sci_model = make_empty_lrs_model("MIR_LRS-FIXEDSLIT")
    sci_model = datamodels.ImageModel(sci_model)  # assign_wcs can't operate on SlitModel

    if use_asn:
        ta_model = make_slit_data(offset=(0, 0))
        step_input = make_ta_association(sci_model, ta_model)
    else:
        step_input = sci_model

    all_steps = Spec2Pipeline.step_defs.keys()
    do_steps = {step_name: {"skip": True} for step_name in all_steps}
    do_steps["assign_wcs"] = {"skip": False}
    do_steps["targ_centroid"] = {"skip": False}
    result = Spec2Pipeline.call(
        step_input,
        steps=do_steps,
    )
    status = result[0].meta.cal_step.targ_centroid
    if use_asn:
        assert status == "COMPLETE"
    else:
        # TA centering image not provided
        assert status == "SKIPPED"


def test_spec2_nirspec_ifu_model(make_test_rate_file):
    input_model = datamodels.open(INPUT_FILE)
    input_model_copy = input_model.copy()

    # Skip all steps but the first for speed
    all_steps = Spec2Pipeline.step_defs.keys()
    do_steps = {step_name: {"skip": True} for step_name in all_steps}
    do_steps["assign_wcs"] = {"skip": False}
    Spec2Pipeline.call(input_model, save_results=True, steps=do_steps)

    # Check for expected output
    assert Path("test_cal.fits").exists()

    # Input model is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step._instance == input_model_copy.meta.cal_step._instance

    input_model.close()
    input_model_copy.close()
