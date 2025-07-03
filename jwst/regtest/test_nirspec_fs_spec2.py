import os
import warnings

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff

from jwst.lib.suffix import replace_suffix
from jwst.pathloss import PathLossStep
from jwst.stpipe import Step

file_roots = [
    "jw02072-o002_20221206t143745_spec2_00001_asn.json",
    "jw01245-o002_20230107t223023_spec2_00001_asn.json",
    "jw01309-o022_20230113t025924_spec2_00001_asn.json",
]

asn_memberdict = {
    "jw01309-o022_20230113t025924_spec2_00001_asn.json": [
        "jw01309022001_04102_00004_nrs2_rate.fits",
        "jw01309022001_04102_00002_nrs2_rate.fits",
        "jw01309022001_04102_00001_nrs2_rate.fits",
    ],
    "jw01245-o002_20230107t223023_spec2_00001_asn.json": [
        "jw01245002001_04102_00002_nrs1_rate.fits",
        "jw01245002001_04102_00001_nrs1_rate.fits",
    ],
    "jw02072-o002_20221206t143745_spec2_00001_asn.json": [
        "jw02072002001_05101_00001_nrs1_rate.fits",
        "jw02072002001_05101_00002_nrs1_rate.fits",
        "jw02072002001_05101_00003_nrs1_rate.fits",
    ],
}
# ids = ["fullframe", "S400A1-subarray", "ALLSLITS-subarray"]

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module", params=file_roots)  # ids=ids)
def run_pipeline(rtdata_module, request, resource_tracker):
    """Run the calwebb_spec2 pipeline on NIRSpec Fixed-Slit exposures.
    We currently test the following types of inputs:
      1) Full-frame exposure (all slits will be extracted)
      2) ALLSLITS subarray exposure (all slits will be extracted)
      3) S400A1 subarray exposure (1 slit extracted)"""

    rtdata = rtdata_module
    for filename in asn_memberdict[request.param]:
        rtdata.get_data("nirspec/fs/" + filename)
    # Get the input exposure
    rtdata.get_data("nirspec/fs/" + request.param)

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=true",
        "--steps.extract_2d.save_results=true",
        "--steps.wavecorr.save_results=true",
        "--steps.srctype.save_results=true",
        "--steps.flat_field.save_results=true",
        "--steps.pathloss.save_results=true",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def run_pipeline_nsclean(rtdata_module):
    """Run the calwebb_spec2 pipeline on NIRSpec Fixed-Slit exposures with nsclean."""

    rtdata = rtdata_module

    filename = "jw01245-o002_20230107t223023_spec2_00001_asn.json"
    rtdata.get_asn("nirspec/fs/" + filename)

    # Run the calwebb_spec2 pipeline with nsclean
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--output_file=jw01245002001_04102_00002_nrs1_nsc",
        "--steps.nsclean.skip=False",
        "--steps.nsclean.save_results=true",
    ]
    Step.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def run_pipeline_pixel_replace(rtdata_module):
    """Run the calwebb_spec2 pipeline on NIRSpec Fixed-Slit exposures including pixel replacement.
    The run_pipeline fixture saves the output for step before pixel_replace, so no need to retest.
    We currently test the following types of inputs:
      Full-frame exposure (all slits will be extracted)
    """
    rtdata = rtdata_module
    rtdata.get_asn("nirspec/fs/jw01309_prtest_spec2_00001_asn.json")

    # Run the calwebb_spec2 pipeline; save results from intermediate steps
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.pixel_replace.save_results=true",
        "--steps.pixel_replace.skip=false",
    ]
    Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_spec2(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "assign_wcs",
        "extract_2d",
        "wavecorr",
        "flat_field",
        "pathloss",
        "srctype",
        "cal",
        "s2d",
        "x1d",
    ],
)
def test_nirspec_fs_spec2(run_pipeline, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline on a
    NIRSpec FS exposures."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline
    output = (
        replace_suffix(
            os.path.splitext(asn_memberdict[os.path.basename(rtdata.input)][0])[0], suffix
        )
        + ".fits"
    )
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_fs_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("suffix", ["nsclean", "cal", "s2d", "x1d"])
def test_nirspec_fs_nsclean(run_pipeline_nsclean, fitsdiff_default_kwargs, suffix):
    """Regression test of the calwebb_spec2 pipeline with nsclean."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline_nsclean
    basename = "jw01245002001_04102_00002_nrs1_nsc"
    output = f"{basename}_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_fs_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize(
    "output",
    [
        "jw013090_prtest_04102_00004_nrs2_pixel_replace.fits",
        "jw013090_prtest_04102_00004_nrs2_cal.fits",
        "jw013090_prtest_04102_00004_nrs2_s2d.fits",
        "jw013090_prtest_04102_00004_nrs2_x1d.fits",
    ],
)
def test_nirspec_fs_spec2_pixel_replace(
    run_pipeline_pixel_replace, fitsdiff_default_kwargs, output
):
    """Regression test of the calwebb_spec2 pipeline on a
    NIRSpec FS exposures."""

    # Run the pipeline and retrieve outputs
    rtdata = run_pipeline_pixel_replace
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_fs_spec2_pixel_replace", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_pathloss_corrpars(rtdata):
    """Test PathLossStep using correction_pars"""
    basename = "jw02072002001_05101_00001_nrs1_flatfieldstep"
    with dm.open(rtdata.get_data(f"nirspec/fs/{basename}.fits")) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.use_correction_pars = True
        corrected_corrpars = pls.run(data)

    bad_slits = []
    for idx, slits in enumerate(zip(corrected.slits, corrected_corrpars.slits)):
        corrected_slit, corrected_corrpars_slit = slits
        if not np.allclose(corrected_slit.data, corrected_corrpars_slit.data, equal_nan=True):
            bad_slits.append(idx)
    assert not bad_slits, f"correction_pars failed for slits {bad_slits}"


def test_pathloss_inverse(rtdata):
    """Test PathLossStep using inversion"""
    basename = "jw02072002001_05101_00001_nrs1_flatfieldstep"
    with dm.open(rtdata.get_data(f"nirspec/fs/{basename}.fits")) as data:
        pls = PathLossStep()
        corrected = pls.run(data)

        pls.inverse = True
        corrected_inverse = pls.run(corrected)

        bad_slits = []
        for idx, slits in enumerate(zip(data.slits, corrected_inverse.slits)):
            data_slit, corrected_inverse_slit = slits
            non_nan = ~np.isnan(corrected_inverse_slit.data)
            if not np.allclose(data_slit.data[non_nan], corrected_inverse_slit.data[non_nan]):
                bad_slits.append(idx)

    assert not bad_slits, f"Inversion failed for slits {bad_slits}"


def test_pathloss_source_type(rtdata):
    """Test PathLossStep forcing source type"""
    basename = "jw02072002001_05101_00001_nrs1_flatfieldstep"
    with dm.open(rtdata.get_data(f"nirspec/fs/{basename}.fits")) as data:
        pls = PathLossStep()
        pls.source_type = "extended"
        pls.run(data)

    bad_slits = []
    for idx, slit in enumerate(pls.correction_pars.slits):
        if slit:
            if not np.allclose(slit.data, slit.pathloss_uniform, equal_nan=True):
                bad_slits.append(idx)
    assert not bad_slits, f"Force to uniform failed for slits {bad_slits}"


def test_nirspec_fs_rateints_spec2(rtdata_module):
    """Run the calwebb_spec2 pipeline on a NIRSpec Fixed-Slit _rateints exposure.
    This is a test that the pipeline completes when processing this
    multi-integration input.
    """
    rtdata = rtdata_module

    rtdata.get_data("nirspec/fs/jw01128004001_03102_00001_nrs1_rateints.fits")

    # Run the spec2 pipeline on a (3D) _rateints file
    args = ["calwebb_spec2", rtdata.input]
    Step.from_cmdline(args)
