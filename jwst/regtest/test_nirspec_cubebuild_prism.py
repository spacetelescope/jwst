"""Test CubeBuildStep on NIRSpec Prims data"""

import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step


@pytest.mark.bigdata
def test_cube_build_nirspec_prism_linear(rtdata, fitsdiff_default_kwargs):
    """Run cube_build on prism data and produce a linear wavelength"""
    input_file = "jw01208062001_07101_00001_nrs1_cal.fits"
    rtdata.get_data(f"nirspec/ifu/{input_file}")

    args = ["jwst.cube_build.CubeBuildStep", input_file]
    Step.from_cmdline(args)

    output = input_file.replace("cal", "prism-clear_s3d")
    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_cubebuild_prism/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.bigdata
def test_cube_build_nirspec_prism_nonlinear(rtdata, fitsdiff_default_kwargs):
    """Run cube_build on prism data and produce a non-linear wavelength"""
    input_file = "jw01208062001_07101_00001_nrs1_cal.fits"
    rtdata.get_data(f"nirspec/ifu/{input_file}")

    output_name = "jw01208062001_07101_00001_nrs1_nonlinear"
    args = [
        "jwst.cube_build.CubeBuildStep",
        input_file,
        "--output_type=multi",
        "--output_file=" + output_name,
    ]
    Step.from_cmdline(args)
    output = output_name + "_prism-clear_s3d.fits"

    rtdata.output = output

    rtdata.get_truth(f"truth/test_nirspec_cubebuild_prism/{output}")

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
