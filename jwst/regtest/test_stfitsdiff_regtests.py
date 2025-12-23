"""Regression tests for STFitsDiff"""

import warnings

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.regtest.st_fitsdiff import STFITSDiffBeta as STFITSDiff
from jwst.regtest.test_stfitsdiff import report_to_list
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


def astropy_fitsdiff(file, truth_file, fitsdiff_default_kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        diff = FITSDiff(file, truth_file, **fitsdiff_default_kwargs)
        apresult = diff.identical
        # Turn the report into a list of strings for easier comparison
        apreport = report_to_list(diff.report())
        return apresult, apreport


def get_stfitsdiff_reports(file, truth_file, fitsdiff_default_kwargs):
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    stdiff = STFITSDiff(file, truth_file, **fitsdiff_default_kwargs)
    result = stdiff.identical
    # The function report_to_list returns two lists when report_pixel_loc_diffs
    # is True. The first list is the ST ad hoc report and the second is the
    # report that  matches astropy's, with pixel or table location differences.
    # This report is the one we want to compare to validate stfitsdiff.
    _, report = report_to_list(stdiff.report(), report_pixel_loc_diffs=True)
    return result, report


@pytest.mark.parametrize(
    "suffix",
    ["s3d"],
)
def test_nirspec_ifu(rtdata_module, suffix, fitsdiff_default_kwargs):
    rtdata = rtdata_module
    rtdata.get_data("nirspec/ifu/jw01251004001_03107_00001_nrs1_rate.fits")
    args = [
        "calwebb_spec2",
        rtdata.input,
    ]
    Step.from_cmdline(args)

    output = "jw01251004001_03107_00001_nrs1_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_nirspec_ifu/" + output)
    apresult, apreport = astropy_fitsdiff(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)
    result, report = get_stfitsdiff_reports(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)

    assert result == apresult
    assert report == apreport


@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize(
    "source_id,slit_name",
    [
        ("s000000024", "s1600a1"),
    ],
)
def test_nirspec_fs_spec3(rtdata_module, fitsdiff_default_kwargs, suffix, source_id, slit_name):
    """Test spec3 pipeline on a set of NIRSpec FS exposures."""
    rtdata = rtdata_module

    # Get the ASN file and input exposures
    rtdata.get_asn("nirspec/fs/jw01309-o022_spec3_regtest_asn.json")

    # Run the calwebb_spec3 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec3",
        rtdata.input,
        "--steps.outlier_detection.save_results=true",
        "--steps.resample_spec.save_results=true",
        "--steps.extract_1d.save_results=true",
    ]
    Step.from_cmdline(args)

    output = f"jw01309-o022_{source_id}_nirspec_f290lp-g395h-{slit_name}-allslits_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_fs_spec3/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-2
        fitsdiff_default_kwargs["atol"] = 2e-4

    # Compare the results
    apresult, apreport = astropy_fitsdiff(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)
    result, report = get_stfitsdiff_reports(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)

    assert result == apresult
    assert report == apreport


@pytest.mark.parametrize(
    "suffix",
    ["esec"],
)
def test_nis_wfss_spec2(rtdata_module, fitsdiff_default_kwargs, suffix):
    """Regression test for calwebb_spec2 applied to NIRISS WFSS data"""
    rtdata = rtdata_module
    spec2_asns = [
        "niriss/wfss/jw01324-o001_20220629t171902_spec2_001_asn.json",
        "niriss/wfss/jw01324-o001_20220629t171902_spec2_002_asn.json",
        "niriss/wfss/jw01324-o001_20220629t171902_spec2_005_asn.json",
        "niriss/wfss/jw01324-o001_20220629t171902_spec2_007_asn.json",
    ]
    rtdata.get_asn(spec2_asns[0])
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--save_wfss_esec=true",
    ]
    Step.from_cmdline(args)

    # Run the remaining exposures without doing comparisons, just so that
    # fresh results are available for level-3 processing
    for asn in spec2_asns[1:]:
        rtdata.get_asn(asn)
        args = ["calwebb_spec2", rtdata.input, "--steps.extract_2d.wfss_nbright=10"]
        Step.from_cmdline(args)

    rtdata.input = "jw01324001001_03101_00001_nis_rate.fits"
    output = "jw01324001001_03101_00001_nis_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_niriss_wfss/{output}")

    # Compare the results
    apresult, apreport = astropy_fitsdiff(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)
    result, report = get_stfitsdiff_reports(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)

    assert result == apresult
    assert report == apreport


def test_miri_mrs_extract1d_nominal(rtdata, fitsdiff_default_kwargs):
    """Test running extract_1d on an s3d cube containing a point source"""
    # input s3d are created using the same data that was used in test_miri_mrs_spec3_ifushort:
    # run calwebb_spec3 on miri/mrs/jw01024_ifushort_mediumlong_spec3_00001_asn.json to create
    # the input data for this test.

    rtdata.get_data("miri/mrs/jw01024-c1000_t002_miri_ch2-mediumlong_s3d.fits")

    args = ["jwst.extract_1d.Extract1dStep", rtdata.input]
    Step.from_cmdline(args)
    rtdata.output = "jw01024-c1000_t002_miri_ch2-mediumlong_extract1dstep.fits"

    # Get the truth file
    rtdata.get_truth(
        "truth/test_miri_mrs_extract1d/jw01024-c1000_t002_miri_ch2-mediumlong_extract1dstep.fits"
    )

    # Compare the results
    apresult, apreport = astropy_fitsdiff(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)
    result, report = get_stfitsdiff_reports(rtdata.output, rtdata.truth, fitsdiff_default_kwargs)

    assert result == apresult
    assert report == apreport


def test_one_nan_in_table(rtdata, fitsdiff_default_kwargs):
    """Test running STFITSDiff on a table with one NaN difference."""

    rtdata.get_data("nirspec/mos/nrs1_one_nan_x1d.fits")
    rtdata.get_truth("nirspec/mos/nrs1_no_nans_x1d.fits")

    apresult, apreport = astropy_fitsdiff(rtdata.input, rtdata.truth, fitsdiff_default_kwargs)
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    result, report = get_stfitsdiff_reports(rtdata.input, rtdata.truth, fitsdiff_default_kwargs)
    assert result == apresult
    assert report == apreport

    diff = STFITSDiff(rtdata.input, rtdata.truth, **fitsdiff_default_kwargs)
    assert "    1 failed the (atol, rtol) test" in diff.report()
    assert "    Found 1 different table data element(s). " in diff.report()
    assert "    WAVELENGTH    f8         0     nan      nan     nan" in diff.report()
