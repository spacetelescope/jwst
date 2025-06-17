import pytest
import os
import shutil
from jwst.stpipe import Step
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import ImageModel  # type: ignore[attr-defined]


INPUT_FILE = "dummy_cal.fits"
INPUT_FILE_2 = "dummy2_cal.fits"
INPUT_ASN = "dummy_asn.json"
OUTPUT_PRODUCT = "custom_name"
LOGFILE = "run_asn.log"
LOGCFG = "test_logs.cfg"
LOGCFG_CONTENT = f"[*] \n \
        handler = file:{LOGFILE}"


@pytest.fixture(scope="module")
def make_dummy_cal_file(tmp_cwd_module):
    """
    Make and save a dummy cal file in the temporary working directory
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

    with image as dm:
        dm.save(INPUT_FILE)


@pytest.fixture(scope="module")
def make_dummy_association(make_dummy_cal_file):
    shutil.copy(INPUT_FILE, INPUT_FILE_2)
    os.system(
        f"asn_from_list -o {INPUT_ASN} --product-name {OUTPUT_PRODUCT} -r DMS_Level3_Base {INPUT_FILE} {INPUT_FILE_2}"
    )


@pytest.mark.filterwarnings("ignore::ResourceWarning")  # in_memory=False
@pytest.mark.parametrize("in_memory", [True, False])
def test_run_image3_pipeline(make_dummy_association, in_memory):
    """
    Two-product association passed in, run pipeline, skipping most steps
    """
    # save warnings to logfile so can be checked later
    with open(LOGCFG, "w") as f:
        f.write(LOGCFG_CONTENT)

    args = [
        "calwebb_image3",
        INPUT_ASN,
        f"--logcfg={LOGCFG}",
        "--steps.tweakreg.skip=true",
        "--steps.skymatch.skip=true",
        "--steps.outlier_detection.skip=true",
        "--steps.resample.skip=true",
        "--steps.source_catalog.skip=true",
        f"--in_memory={str(in_memory)}",
    ]

    Step.from_cmdline(args)

    _is_run_complete(LOGFILE)


@pytest.mark.filterwarnings("ignore::ResourceWarning")
def test_run_image3_single_file(make_dummy_cal_file):
    with open(LOGCFG, "w") as f:
        f.write(LOGCFG_CONTENT)

    args = [
        "calwebb_image3",
        INPUT_FILE,
        f"--logcfg={LOGCFG}",
        "--steps.tweakreg.skip=true",
        "--steps.skymatch.skip=true",
        "--steps.outlier_detection.skip=true",
        "--steps.resample.skip=true",
        "--steps.source_catalog.skip=true",
    ]

    Step.from_cmdline(args)
    _is_run_complete(LOGFILE)


def _is_run_complete(logfile):
    """
    Check that the pipeline runs to completion
    """
    msg = "Step Image3Pipeline done"
    with open(LOGFILE, "r") as f:
        log = f.read()
    assert msg in log
