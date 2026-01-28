import os
import shutil

import pytest

from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel  # type: ignore[attr-defined]
from jwst.pipeline import Image3Pipeline
from jwst.stpipe import Step

INPUT_FILE = "mock_cal.fits"
INPUT_FILE_2 = "mock2_cal.fits"
INPUT_ASN = "mock_asn.json"
OUTPUT_PRODUCT = "custom_name"
LOGFILE = "run_asn.log"


@pytest.fixture(scope="module")
def make_mock_cal_model():
    """
    Make a mock cal model.

    Partially copied from test_calwebb_image2.py
    """

    image = ImageModel((2048, 2048))
    image.data[:, :] = 1
    image.meta.instrument.name = "NIRCAM"
    image.meta.instrument.filter = "F210M"
    image.meta.instrument.pupil = "CLEAR"
    image.meta.exposure.type = "NRC_IMAGE"
    image.meta.observation.date = "2024-02-27"
    image.meta.observation.time = "13:37:18.548"
    image.meta.date = "2024-02-27T13:37:18.548"
    image.meta.subarray.xstart = 1
    image.meta.subarray.ystart = 1

    image.meta.subarray.xsize = image.data.shape[-1]
    image.meta.subarray.ysize = image.data.shape[-2]

    image.meta.instrument.channel = "SHORT"
    image.meta.instrument.module = "A"
    image.meta.instrument.detector = "NRCA1"

    # bare minimum wcs info to get assign_wcs step to pass
    image.meta.wcsinfo.crpix1 = 693.5
    image.meta.wcsinfo.crpix2 = 512.5
    image.meta.wcsinfo.v2_ref = -453.37849
    image.meta.wcsinfo.v3_ref = -373.810549
    image.meta.wcsinfo.roll_ref = 272.3237653262276
    image.meta.wcsinfo.ra_ref = 80.54724018120017
    image.meta.wcsinfo.dec_ref = -69.5081101864959

    image = AssignWcsStep.call(image)

    return image


@pytest.fixture(scope="module")
def make_mock_cal_file(tmp_cwd_module, make_mock_cal_model):
    """
    Make and save a mock cal file in the temporary working directory.

    Partially copied from test_calwebb_image2.py
    """
    with make_mock_cal_model as dm:
        dm.save(INPUT_FILE)


@pytest.fixture(scope="module")
def make_mock_association(make_mock_cal_file):
    shutil.copy(INPUT_FILE, INPUT_FILE_2)
    os.system(
        f"asn_from_list -o {INPUT_ASN} --product-name {OUTPUT_PRODUCT} -r DMS_Level3_Base {INPUT_FILE} {INPUT_FILE_2}"
    )


@pytest.mark.parametrize("in_memory", [True, False])
def test_run_image3_pipeline(make_mock_association, in_memory):
    """
    Two-product association passed in, run pipeline, skipping most steps
    """
    # save warnings to logfile so can be checked later
    args = [
        "calwebb_image3",
        INPUT_ASN,
        "--log-level=INFO",
        f"--log-file={LOGFILE}",
        "--steps.tweakreg.skip=true",
        "--steps.skymatch.skip=true",
        "--steps.outlier_detection.skip=true",
        "--steps.resample.skip=true",
        "--steps.source_catalog.skip=true",
        f"--in_memory={str(in_memory)}",
    ]

    Step.from_cmdline(args)

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
