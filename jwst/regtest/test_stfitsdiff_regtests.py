"""Regression tests for STFitsDiff"""

import pytest

from stdatamodels.jwst import datamodels
from jwst.stpipe import Step
from jwst.regtest.st_fitsdiff import STFITSDiffBeta as STFITSDiff
from jwst.regtest.test_stfitsdiff import report_to_list
from jwst.flatfield.flat_field import nirspec_ifu
from astropy.io.fits.diff import FITSDiff


# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


def test_nirspec_ifu_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test using predefined interpolated flat"""
    basename = "jw01251004001_03107_00001_nrs1"
    output_file = f"{basename}_flat_from_user_model.fits"
    with datamodels.open(rtdata.get_data(f"nirspec/ifu/{basename}_assign_wcs.fits")) as data:
        with datamodels.open(
            rtdata.get_data(f"nirspec/ifu/{basename}_interpolatedflat.fits")
        ) as user_supplied_flat:
            # Call the flat field function directly with a user flat
            nirspec_ifu(data, None, None, None, None, user_supplied_flat=user_supplied_flat)
    rtdata.output = output_file
    data.save(rtdata.output)

    rtdata.get_truth("truth/test_nirspec_ifu/" + output_file)
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    apresult, apreport = diff.identical, report_to_list(diff.report())

    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    stdiff = STFITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    result = stdiff.identical
    _, report = report_to_list(stdiff.report(), report_pixel_loc_diffs=True)

    assert result == apresult
    assert report == apreport


@pytest.mark.parametrize("suffix", ["cal", "crf", "s2d", "x1d"])
@pytest.mark.parametrize(
    "source_id",
    [
        "b000000030",
        "b000000031",
        "s000004385",
        "s000007380",
        "v000000048",
        "v000000049",
        "v000000053",
        "v000000056",
    ],
)
def test_nirspec_mos_spec3(
    rtdata_module, resource_tracker, suffix, source_id, fitsdiff_default_kwargs
):
    """Check results of calwebb_spec3"""
    rtdata = rtdata_module
    rtdata.get_asn("nirspec/mos/jw01345-o066_20230831t181155_spec3_00002_asn.json")

    # Run the calwebb_spec3 pipeline on the association
    args = ["calwebb_spec3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)

    output = f"jw01345-o066_{source_id}_nirspec_f170lp-g235m_{suffix}.fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nirspec_mos_spec3/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "s2d":
        fitsdiff_default_kwargs["rtol"] = 1e-4
        fitsdiff_default_kwargs["atol"] = 1e-5

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    apresult, apreport = diff.identical, report_to_list(diff.report())

    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    stdiff = STFITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    result = stdiff.identical
    _, report = report_to_list(stdiff.report(), report_pixel_loc_diffs=True)

    assert result == apresult
    assert report == apreport


@pytest.mark.parametrize(
    "suffix",
    ["assign_wcs", "esec", "x1d"],
)
def test_nis_wfss_spec2(rtdata_module, resource_tracker, fitsdiff_default_kwargs, suffix):
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
        "--steps.assign_wcs.save_results=true",
        "--steps.extract_1d.save_results=true",
        "--save_wfss_esec=true",
        rtdata.input,
    ]
    with resource_tracker.track():
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

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    apresult, apreport = diff.identical, report_to_list(diff.report())

    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    stdiff = STFITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    result = stdiff.identical
    _, pixreport = report_to_list(stdiff.report(), report_pixel_loc_diffs=True)

    assert result == apresult
    assert pixreport == apreport


def test_nircam_tsgrism_stage3_x1dints(rtdata_module, resource_tracker, fitsdiff_default_kwargs):
    rtdata = rtdata_module
    # Run spec2 pipeline on the _rateints file, saving intermediate products
    rtdata.get_data("nircam/tsgrism/jw01366002001_04103_00001-seg001_nrcalong_rateints.fits")
    args = [
        "calwebb_spec2",
        rtdata.input,
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)
    # Run spec3
    rtdata.get_data("nircam/tsgrism/jw01366-o002_20230107t004627_tso3_00001_asn.json")
    args = ["calwebb_tso3", rtdata.input]
    with resource_tracker.track():
        Step.from_cmdline(args)

    # Compare
    rtdata.input = "jw01366-o002_20230107t004627_tso3_00001_asn.json"
    rtdata.output = "jw01366-o002_t001_nircam_f322w2-grismr-subgrism256_x1dints.fits"
    rtdata.get_truth(
        "truth/test_nircam_tsgrism_stages/jw01366-o002_t001_nircam_f322w2-grismr-subgrism256_x1dints.fits"
    )
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    apresult, apreport = diff.identical, report_to_list(diff.report())

    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    stdiff = STFITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    result = stdiff.identical
    _, pixreport = report_to_list(stdiff.report(), report_pixel_loc_diffs=True)

    assert result == apresult
    assert pixreport == apreport


@pytest.mark.parametrize(
    "suffix",
    [
        "rate",
        "cal",
        "i2d",
    ],
)
def test_nircam_image_stages12(rtdata_module, resource_tracker, fitsdiff_default_kwargs, suffix):
    """Regression test of detector1 and image2 pipelines performed on NIRCam data."""
    rtdata = rtdata_module
    rtdata.get_data("nircam/image/jw01345001001_10201_00001_nrca3_uncal.fits")
    # Run detector1 pipeline only on one of the _uncal files
    args = [
        "calwebb_detector1",
        rtdata.input,
        "--output_file=jw01345001001_10201_00001_nrca3_likely",
        "--steps.ramp_fit.algorithm=LIKELY",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)
    # Run image2
    rtdata.input = "jw01538046001_03105_00001_nrcalong_rate.fits"
    args = [
        "calwebb_image2",
        rtdata.input,
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    # Compare
    rtdata.input = "jw01538046001_03105_00001_nrcalong_uncal.fits"
    output = "jw01538046001_03105_00001_nrcalong_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth(f"truth/test_nircam_image_stages/{output}")

    # Adjust tolerance for machine precision with float32 drizzle code
    if suffix == "i2d":
        fitsdiff_default_kwargs["rtol"] = 5e-5
        fitsdiff_default_kwargs["atol"] = 1e-4

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    apresult, apreport = diff.identical, report_to_list(diff.report())

    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    stdiff = STFITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    result = stdiff.identical
    _, pixreport = report_to_list(stdiff.report(), report_pixel_loc_diffs=True)

    assert result == apresult
    assert pixreport == apreport


@pytest.mark.parametrize(
    "suffix",
    [
        "cal",
        "s2d",
        "x1d",
    ],
)
def test_miri_lrs_slit_spec2(rtdata_module, resource_tracker, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on MIRI
    LRS fixedslit data using along-slit-nod pattern for
    background subtraction."""
    rtdata = rtdata_module
    # Get the spec2 ASN and its members
    rtdata.get_asn("miri/lrs/jw01530-o005_20221202t204827_spec2_00001_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)
    output = f"jw01530005001_03103_00001_mirimage_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(f"truth/test_miri_lrs_slit_spec2/{output}")

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    apresult, apreport = diff.identical, report_to_list(diff.report())

    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    stdiff = STFITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    result = stdiff.identical
    _, pixreport = report_to_list(stdiff.report(), report_pixel_loc_diffs=True)

    assert result == apresult
    assert pixreport == apreport
