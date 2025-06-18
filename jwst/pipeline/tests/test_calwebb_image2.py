import os
import shutil

import pytest

from jwst.stpipe import Step
from jwst.datamodels import ImageModel  # type: ignore[attr-defined]


INPUT_FILE = "dummy_rate.fits"
INPUT_FILE_2 = "dummy2_rate.fits"
INPUT_ASN = "dummy_asn.json"
OUTPUT_FILE = "custom_name.fits"
OUTPUT_FILE_ASN = "custom_name_asn.fits"  # cannot reuse because everything runs in same cwd
LOGFILE = "run_asn.log"
LOGCFG = "test_logs.cfg"


@pytest.fixture(scope="module")
def make_dummy_rate_file(tmp_cwd_module):
    """
    Make and save a dummy rate file in the temporary working directory
    Partially copied from test_background.py
    """

    image = ImageModel((2048, 2048))
    image.data[:, :] = 1
    image.meta.instrument.name = "NIRCAM"
    image.meta.instrument.filter = "CLEAR"
    image.meta.exposure.type = "NRC_IMAGE"
    image.meta.observation.date = "2019-02-27"
    image.meta.observation.time = "13:37:18.548"
    image.meta.date = "2019-02-27T13:37:18.548"
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

    with image as dm:
        dm.save(INPUT_FILE)


@pytest.fixture(scope="module")
def make_dummy_association(make_dummy_rate_file):
    shutil.copy(INPUT_FILE, INPUT_FILE_2)
    os.system(f"asn_from_list -o {INPUT_ASN} -r DMSLevel2bBase {INPUT_FILE} {INPUT_FILE_2}")


@pytest.fixture(scope="module", params=[OUTPUT_FILE])
def run_image2_pipeline_file(make_dummy_rate_file, request):
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


@pytest.fixture(scope="module", params=[OUTPUT_FILE_ASN])
def run_image2_pipeline_asn(make_dummy_association, request):
    """
    Two-product association passed in. This should trigger a warning
    and the output_file parameter should be ignored.
    """
    # save warnings to logfile so can be checked later
    logcfg_content = f"[*] \n \
        level = INFO \n \
        handler = file:{LOGFILE}"
    with open(LOGCFG, "w") as f:
        f.write(logcfg_content)

    args = [
        "calwebb_image2",
        INPUT_ASN,
        f"--logcfg={LOGCFG}",
        "--steps.flat_field.skip=true",
        "--steps.photom.skip=true",
        "--steps.resample.skip=true",
        f"--output_file={request.param}",
    ]

    Step.from_cmdline(args)


def test_output_file_rename_file(run_image2_pipeline_file):
    """
    Covers a bug where the output_file parameter was not being
    respected in calls to Image2Pipeline.
    """
    assert os.path.exists(INPUT_FILE)  # ensures tmp_cwd_module is working
    custom_stem = OUTPUT_FILE.split(".")[0]
    for extension in ["cal"]:
        assert os.path.exists(f"{custom_stem}_{extension}.fits")


@pytest.mark.filterwarnings("ignore::ResourceWarning")
def test_output_file_norename_asn(run_image2_pipeline_asn):
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
