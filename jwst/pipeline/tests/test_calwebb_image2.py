import os
import shutil
from pathlib import Path

import numpy as np
import pytest
from stdatamodels.jwst import datamodels
from stpipe import _log

from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.pipeline import Image2Pipeline
from jwst.pipeline.tests.helpers import make_nircam_rate_model
from jwst.stpipe import Step

INPUT_FILE = "test_rate.fits"
INPUT_FILE_2 = "test2_rate.fits"
OUTPUT_FILE = "custom_name.fits"
OUTPUT_FILE_ASN = "custom_name_asn.fits"  # cannot reuse because everything runs in same cwd
LOGFILE = "run_asn.log"


@pytest.fixture(scope="module")
def make_rate_file(tmp_cwd_module):
    """
    Make and save a rate file in the temporary working directory
    Partially copied from test_background.py
    """
    image = make_nircam_rate_model()
    with image as dm:
        dm.save(INPUT_FILE)


@pytest.fixture(scope="module")
def make_association(make_rate_file):
    shutil.copy(INPUT_FILE, INPUT_FILE_2)
    return asn_from_list([INPUT_FILE, INPUT_FILE_2], rule=DMSLevel2bBase)


@pytest.fixture(scope="module", params=[OUTPUT_FILE])
def run_image2_pipeline_file(make_rate_file, request):
    """
    Run pipeline, skipping most steps
    """
    args = [
        "calwebb_image2",
        INPUT_FILE,
        "--steps.flat_field.skip=true",
        "--steps.photom.skip=true",
        "--steps.resample.skip=true",
        f"--output_file={request.param}",
    ]

    Step.from_cmdline(args)


@pytest.fixture(scope="module")
def run_image2_pipeline_asn(make_association):
    """
    Two-product association passed in. This should trigger a warning
    and the output_file parameter should be ignored.
    """
    # save warnings to logfile so can be checked later
    log_cfg = _log.load_configuration(log_level="INFO", log_file=LOGFILE)
    log_cfg.set_recording_formatter(Image2Pipeline._log_records_formatter)
    with log_cfg.context(Image2Pipeline.get_stpipe_loggers()):
        Image2Pipeline.call(
            make_association,
            save_results=True,
            steps={
                "flat_field": {"skip": True},
                "photom": {"skip": True},
                "resample": {"skip": True},
            },
            output_file=OUTPUT_FILE_ASN,
        )


def test_output_file_rename_file(run_image2_pipeline_file):
    """
    Covers a bug where the output_file parameter was not being
    respected in calls to Image2Pipeline.
    """
    assert os.path.exists(INPUT_FILE)  # ensures tmp_cwd_module is working
    custom_stem = OUTPUT_FILE.split(".")[0]
    for extension in ["cal"]:
        assert os.path.exists(f"{custom_stem}_{extension}.fits")


def test_output_file_norename_asn(run_image2_pipeline_asn):
    """
    Ensure output_file parameter is ignored, with warning,
    when multiple products are in the same association.
    """
    # ensure tmp_cwd_module is successfully keeping all files in cwd
    assert os.path.exists(INPUT_FILE)
    assert os.path.exists(INPUT_FILE_2)

    custom_stem = OUTPUT_FILE_ASN.split(".")[0]
    input_stem = INPUT_FILE.split("_")[0]
    input_stem_2 = INPUT_FILE_2.split("_")[0]

    # ensure default filenames were written, and not the custom one
    for extension in ["cal"]:
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


def test_image2_nircam_model(make_rate_file):
    input_model = datamodels.open(INPUT_FILE)
    input_model_copy = input_model.copy()
    Image2Pipeline.call(input_model, save_results=True, steps={"resample": {"skip": True}})

    # Check for expected output
    assert Path("test_cal.fits").exists()

    # Input model is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step.instance == input_model_copy.meta.cal_step.instance

    input_model.close()
    input_model_copy.close()
