import os

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.flatfield import FlatFieldStep
from jwst.lib.suffix import replace_suffix
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_tso_spec2_pipeline(rtdata_module, resource_tracker):
    """Run the calwebb_spec2 pipeline performed on NIRSpec
    fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit)
    """

    rtdata = rtdata_module

    # Get the input exposure
    # Input data is from jw02420001001_04101_00001-seg001_nrs1_rateints.fits,
    # modified to truncate the data to the first 100 integrations, for
    # faster processing.
    rtdata.get_data("nirspec/tso/jw02420001001_04101_00001-first100_nrs1_rateints.fits")

    # Run the calwebb_spec2 pipeline;
    args = [
        "calwebb_spec2",
        rtdata.input,
        "--steps.assign_wcs.save_results=True",
        "--steps.extract_2d.save_results=True",
        "--steps.wavecorr.save_results=True",
        "--steps.flat_field.save_results=True",
        "--steps.flat_field.save_interpolated_flat=True",
        "--steps.photom.save_results=True",
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


@pytest.fixture(scope="module")
def run_tso3_pipeline(rtdata_module, resource_tracker):
    """Run the calwebb_tso3 pipeline on NRS_BRIGHTOBJ (S1600A1 slit) data."""

    rtdata = rtdata_module

    # Get the input ASN and data
    # Input data is from jw01118005001_04101_00001-seg001_nrs1_calints.fits
    # and jw01118005001_04101_00001-seg001_nrs2_calints.fits,
    # modified to truncate the data to the first 20 integrations, for
    # faster processing.
    rtdata.get_asn("nirspec/tso/jw01118005001_04101_tso3_asn.json")

    # Run the calwebb_tso3 pipeline
    args = ["calwebb_tso3", rtdata.input]

    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_spec2(log_tracked_resources, run_tso_spec2_pipeline):
    log_tracked_resources()


def test_log_tracked_resources_tso3(log_tracked_resources, run_tso3_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize(
    "suffix",
    [
        "assign_wcs",
        "extract_2d",
        "wavecorr",
        "flat_field",
        "photom",
        "calints",
        "x1dints",
        "interpolatedflat",
    ],
)
def test_nirspec_brightobj_spec2(run_tso_spec2_pipeline, fitsdiff_default_kwargs, suffix):
    """
    Regression test of calwebb_spec2 pipeline performed on NIRSpec
    fixed-slit data that uses the NRS_BRIGHTOBJ mode (S1600A1 slit).
    """
    rtdata = run_tso_spec2_pipeline
    basename = "jw02420001001_04101_00001-first100_nrs1"
    output = f"{basename}_{suffix}.fits"
    rtdata.output = output

    # Get the truth files
    rtdata.get_truth(os.path.join("truth/test_nirspec_brightobj_spec2", output))

    # Compare the results
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_tso3_x1dints(run_tso3_pipeline, fitsdiff_default_kwargs):
    rtdata = run_tso3_pipeline
    rtdata.output = "jw01118005001_04101_00001-first20_x1dints.fits"
    rtdata.get_truth(
        "truth/test_nirspec_brightobj_tso3/jw01118005001_04101_00001-first20_x1dints.fits"
    )

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_tso3_whtlt(run_tso3_pipeline, diff_astropy_tables):
    rtdata = run_tso3_pipeline
    rtdata.output = "jw01118005001_04101_00001-first20_whtlt.ecsv"
    rtdata.get_truth(
        "truth/test_nirspec_brightobj_tso3/jw01118005001_04101_00001-first20_whtlt.ecsv"
    )

    assert diff_astropy_tables(rtdata.output, rtdata.truth)


def test_flat_field_step_user_supplied_flat(rtdata, fitsdiff_default_kwargs):
    """Test providing a user-supplied flat field to the FlatFieldStep"""
    basename = "jw02420001001_04101_00001-first100_nrs1"
    output_file = f"{basename}_flat_from_user_file.fits"

    data = rtdata.get_data(f"nirspec/tso/{basename}_wavecorr.fits")
    user_supplied_flat = rtdata.get_data(f"nirspec/tso/{basename}_user_flat.fits")

    data_flat_fielded = FlatFieldStep.call(
        data, user_supplied_flat=user_supplied_flat, save_results=False
    )
    rtdata.output = output_file
    data_flat_fielded.save(rtdata.output)

    rtdata.get_truth(f"truth/test_nirspec_brightobj_spec2/{output_file}")
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_ff_inv(rtdata, fitsdiff_default_kwargs):
    """Test flat field inversion"""
    basename = "jw02420001001_04101_00001-first100_nrs1"
    with dm.open(rtdata.get_data(f"nirspec/tso/{basename}_wavecorr.fits")) as data:
        flatted = FlatFieldStep.call(data)
        unflatted = FlatFieldStep.call(flatted, inverse=True)

    # flat fielding may set some new NaN values - ignore these in test
    is_nan = np.isnan(unflatted.data)
    assert np.allclose(data.data[~is_nan], unflatted.data[~is_nan]), "Inversion failed"

    # make sure NaNs are only at do_not_use pixels
    assert np.all(unflatted.dq[is_nan] & dm.dqflags.pixel["DO_NOT_USE"])
