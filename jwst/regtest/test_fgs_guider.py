"""Regression tests for FGS Guidestar in ID and FINEGUIDE modes"""

import pytest

from jwst.regtest import regtestdata as rt
from jwst.stpipe import Step

EXP_TYPES = ["fgs_acq1", "fgs_fineguide", "fgs_id-image", "fgs_id-stack"]
FILE_ROOTS = [
    "jw01029001001_gs-acq1_2022142180746",
    "jw01029001001_gs-fg_2022142181502",
    "jw01029001001_gs-id_1_image",
    "jw01029001001_gs-id_1_stacked",
]
GUIDER_SUFFIXES = ["cal", "dq_init", "guider_cds"]


@pytest.fixture(scope="module", params=FILE_ROOTS, ids=FILE_ROOTS)
def run_guider_pipelines(rtdata_module, request):
    """Run pipeline for guider data"""
    rtdata = rtdata_module
    rtdata.get_data("fgs/level1b/" + request.param + "_uncal.fits")

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
        "truth/test_fgs_guider",
        is_suffix=True,
    )
