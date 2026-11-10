"""Regression tests for FGS Guidestar in ID and FINEGUIDE modes"""

import pytest

from jwst.regtest import regtestdata as rt
from jwst.stpipe import Step

INPUT_DATA_PATH = "rtdata/fgs/level1b"
INPUT_DATA = {
    "fgs_acq1": rt.RTFile(
        file_name="jw01029001001_gs-acq1_2022142180746_uncal.fits",
        path=INPUT_DATA_PATH,
    ),
    "fgs_fineguide": rt.RTFile(
        file_name="jw01029001001_gs-fg_2022142181502_uncal.fits",
        path=INPUT_DATA_PATH,
    ),
    "fgs_id-image": rt.RTFile(
        file_name="jw01029001001_gs-id_1_image_uncal.fits",
        path=INPUT_DATA_PATH,
    ),
    "fgs_id-stack": rt.RTFile(
        file_name="jw01029001001_gs-id_1_stacked_uncal.fits",
        path=INPUT_DATA_PATH,
    ),
}
GUIDER_SUFFIXES = ["cal", "dq_init", "guider_cds"]


@pytest.fixture(
    scope="module",
    params=[rtobj.full_path for rtobj in INPUT_DATA.values()],
    ids=[k for k in INPUT_DATA.keys()],
)
def run_guider_pipelines(rtdata_module, request):
    """Run pipeline for guider data"""
    rtdata = rtdata_module
    rtdata.get_data(request.param)

    args = [
        "calwebb_guider",
        rtdata.input,
        "--steps.dq_init.save_results=true",
        "--steps.guider_cds.save_results=true",
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.bigdata
@pytest.mark.parametrize("suffix", GUIDER_SUFFIXES, ids=GUIDER_SUFFIXES)
def test_fgs_guider(run_guider_pipelines, fitsdiff_default_kwargs, suffix):
    """Regression for FGS Guider data"""
    rt.is_like_truth(
        run_guider_pipelines,
        fitsdiff_default_kwargs,
        suffix,
        "rtdata/truth/test_fgs_guider",
        is_suffix=True,
    )
