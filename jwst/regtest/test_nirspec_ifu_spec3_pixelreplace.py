"""Regression tests for NIRSpec IFU"""

import warnings

import pytest

from jwst.regtest import regtestdata as rt

TRUTH_PATH_PIXEL_REPLACE = "truth/test_nirspec_ifu_pixel_replace"


@pytest.fixture(scope="module")
def run_spec3_multi_pixel_replace(rtdata_module):
    """Run Spec3Pipeline with pixel replacement"""
    rtdata = rtdata_module

    step_params = {
        "input_path": "nirspec/ifu/jw01249-o005_20230622t074431_rp_spec3_00001_asn.json",
        "step": "calwebb_spec3",
        "args": {
            "--steps.pixel_replace.skip=false",
            "--steps.pixel_replace.save_results=true",
            "--steps.pixel_replace.algorithm=mingrad",
            "--steps.cube_build.save_results=true",
            "--steps.extract_1d.save_results=true",
        },
    }
    rtdata = rt.run_step_from_dict(rtdata, **step_params)
    return rtdata


@pytest.mark.slow
@pytest.mark.bigdata
@pytest.mark.parametrize(
    "output",
    [
        "jw01249005001_03101_00001_nrs1_o005_pixel_replace.fits",
        "jw01249005001_03101_00001_nrs2_o005_pixel_replace.fits",
        "jw01249005001_03101_00002_nrs1_o005_pixel_replace.fits",
        "jw01249005001_03101_00002_nrs2_o005_pixel_replace.fits",
        "jw01249005001_03101_00003_nrs1_o005_pixel_replace.fits",
        "jw01249005001_03101_00003_nrs2_o005_pixel_replace.fits",
        "jw01249005001_03101_00004_nrs1_o005_pixel_replace.fits",
        "jw01249005001_03101_00004_nrs2_o005_pixel_replace.fits",
        "jw01249-o005_t001_nirspec_rp_test_g395h-f290lp_s3d.fits",
        "jw01249-o005_t001_nirspec_rp_test_g395h-f290lp_x1d.fits",
    ],
)
def test_spec3_multi_pixel_replace(run_spec3_multi_pixel_replace, fitsdiff_default_kwargs, output):
    """Regression test matching output files"""
    rt.is_like_truth(
        run_spec3_multi_pixel_replace,
        fitsdiff_default_kwargs,
        output,
        truth_path="truth/test_nirspec_ifu_pixel_replace",
        is_suffix=False,
    )
