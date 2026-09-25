import shutil

import pytest
from stpipe import _log

from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline import Image3Pipeline
from jwst.pipeline.tests.helpers import make_nircam_image_cal_model
from jwst.stpipe import Step

INPUT_FILE = "mock_cal.fits"
INPUT_FILE_2 = "mock2_cal.fits"
OUTPUT_PRODUCT = "custom_name"
LOGFILE = "run_asn.log"


@pytest.fixture(scope="module")
def make_mock_cal_model():
    """Make a mock cal model."""
    return make_nircam_image_cal_model()


@pytest.fixture(scope="module")
def make_mock_cal_file(tmp_cwd_module, make_mock_cal_model):
    """Make and save a mock cal file in the temporary working directory."""
    with make_mock_cal_model as dm:
        dm.save(INPUT_FILE)


@pytest.fixture(scope="module")
def make_mock_association(make_mock_cal_file):
    shutil.copy(INPUT_FILE, INPUT_FILE_2)
    return asn_from_list(
        [INPUT_FILE, INPUT_FILE_2], product_name=OUTPUT_PRODUCT, rule=DMS_Level3_Base
    )


@pytest.mark.parametrize("in_memory", [True, False])
def test_run_image3_pipeline(make_mock_association, in_memory):
    """
    Two-product association passed in, run pipeline, skipping most steps
    """
    # save warnings to logfile so can be checked later
    log_cfg = _log.load_configuration(log_level="INFO", log_file=LOGFILE)
    log_cfg.set_recording_formatter(Image3Pipeline._log_records_formatter)
    with log_cfg.context(Image3Pipeline.get_stpipe_loggers()):
        Image3Pipeline.call(
            make_mock_association,
            steps={
                "tweakreg": {"skip": True},
                "skymatch": {"skip": True},
                "outlier_detection": {"skip": True},
                "resample": {"skip": True},
                "source_catalog": {"skip": True},
            },
            in_memory=in_memory,
        )

    _is_run_complete(LOGFILE)


def test_run_image3_single_file(make_mock_cal_file):
    args = [
        "calwebb_image3",
        INPUT_FILE,
        "--log-level=INFO",
        f"--log-file={LOGFILE}",
        "--steps.tweakreg.skip=true",
        "--steps.skymatch.skip=true",
        "--steps.outlier_detection.skip=true",
        "--steps.resample.skip=true",
        "--steps.source_catalog.skip=true",
    ]

    Step.from_cmdline(args)
    _is_run_complete(LOGFILE)


def test_run_image3_single_model(caplog, make_mock_cal_model):
    model_copy = make_mock_cal_model.copy()
    all_steps = Image3Pipeline.step_defs.keys()
    do_steps = {step_name: {"skip": True} for step_name in all_steps}
    Image3Pipeline.call(model_copy, steps=do_steps)

    # The pipeline ran, but all steps were skipped
    assert "Step Image3Pipeline done" in caplog.text
    for step in do_steps:
        if step in ["assign_mtwcs", "source_catalog"]:
            # These steps are not even attempted for the test data
            assert getattr(model_copy.meta.cal_step, step) is None
        else:
            assert getattr(model_copy.meta.cal_step, step) == "SKIPPED"

        # Input is not modified
        assert getattr(make_mock_cal_model.meta.cal_step, step) is None


def _is_run_complete(logfile):
    """
    Check that the pipeline runs to completion
    """
    msg = "Step Image3Pipeline done"
    with open(LOGFILE, "r") as f:
        log = f.read()
    assert msg in log
