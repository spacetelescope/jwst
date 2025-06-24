"""Tests for STFitsDiff"""

import numpy as np
import pytest

from stdatamodels.jwst import datamodels
from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff


@pytest.fixture(scope="module")
def mock_rampfiles(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("stfitsdiff_ramps")
    truth = tmp_dir / "truth_ramp.fits"
    keyword_mod = tmp_dir / "keyword_mod_ramp.fits"
    sci_mod = tmp_dir / "sci_mod_ramp.fits"
    nan_in_sci = tmp_dir / "nan_in_sci_ramp.fits"

    nints = 3
    ngroups = 10
    nrows = 1032
    ncols = 1024
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint8)
    with datamodels.RampModel(data=data, pixeldq=pixdq, groupdq=gdq) as rampmodel:
        rampmodel.meta.instrument.name = "MIRI"
        rampmodel.meta.instrument.detector = "MIRIMAGE"
        rampmodel.meta.instrument.filter = "F480M"
        rampmodel.meta.observation.date = "2015-10-13"
        rampmodel.meta.exposure.type = "MIR_IMAGE"
        rampmodel.meta.subarray.name = "FULL"
        rampmodel.meta.instrument.gwa_tilt = 37.0610
        # Save the 'truth' file
        rampmodel.save(truth)

        # Change a keyword value only and save
        rampmodel.meta.observation.date = "2025-10-13"
        rampmodel.meta.instrument.gwa_tilt = 37.0510
        rampmodel.save(keyword_mod)
        # return to truth value
        rampmodel.meta.observation.date = "2015-10-13"
        rampmodel.meta.instrument.gwa_tilt = 37.0610

        # Change the data in a couple of places in the SCI extension and save
        rampmodel.data[1, 5, 100, 100] = 100.0
        rampmodel.data[1, 5, 101, 101] = 101.0
        rampmodel.pixeldq[100, 100] = 1.0
        rampmodel.pixeldq[101, 101] = 1.0
        rampmodel.save(sci_mod)
        # return to truth value
        rampmodel.data[1, 5, 100, 100] = 0.0
        rampmodel.data[1, 5, 101, 101] = 0.0
        rampmodel.pixeldq[100, 100] = 0.0
        rampmodel.pixeldq[101, 101] = 0.0

        # Set first group to nans and save
        rampmodel.data[0, 0, ...] = np.nan
        rampmodel.save(nan_in_sci)

    return truth, keyword_mod, sci_mod, nan_in_sci


def report_to_list(report, from_line=11):
    rsplit = report.split(sep="\n")
    rsplit = [rs.strip() for rs in rsplit if rs]
    # Remove the lines of fits diff version and comparison filenames, as well
    # HDUs, keywords and columns not compared, and the maximum number of
    # different data values to be reported and the abs and rel tolerances
    return rsplit[from_line:]


def test_identical(mock_rampfiles):
    truth = mock_rampfiles[0]
    diff = FITSDiff(truth, truth)
    assert diff.identical, diff.report()


def test_keyword_change(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    keyword_mod = mock_rampfiles[1]
    diff = FITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword DATE-OBS has different values:",
        "a> 2025-10-13",
        "?   ^",
        "b> 2015-10-13",
        "?   ^",
        "Keyword FILENAME has different values:",
        "a> keyword_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Keyword GWA_TILT has different values:",
        "a> 37.051",
        "?     ^",
        "b> 37.061",
        "?     ^",
    ]

    extension_tolerances = {"primary": {"rtol": 1, "atol": 2}}
    fitsdiff_default_kwargs["extension_tolerances"] = extension_tolerances
    diff2 = FITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    result2 = diff2.identical
    report2 = report_to_list(diff2.report(), from_line=10)
    expected_report2 = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1.0e+00, Absolute tolerance: 2.0e+00",
        "Headers contain differences:",
        "Keyword DATE-OBS has different values:",
        "a> 2025-10-13",
        "?   ^",
        "b> 2015-10-13",
        "?   ^",
        "Keyword FILENAME has different values:",
        "a> keyword_mod_ramp.fits",
        "b> truth_ramp.fits",
    ]

    assert result is False
    assert report == expected_report
    assert result2 is False
    assert report2 == expected_report2


def test_sci_change(mock_rampfiles, fitsdiff_default_kwargs):
    truth_file = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    diff = FITSDiff(sci_mod, truth_file, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Values in a and b",
        "Quantity     a        b",
        "---------- -------- --------",
        "zeros 31703038 31703040",
        "nan        0        0",
        "no-nan 31703040 31703040",
        "min_value        0        0",
        "max_value      101        0",
        "mean_value 6.34e-06        0",
        "Difference stats: abs(a - b)",
        "Quantity abs_diff no_rel_stats_available",
        "-------- -------- ----------------------",
        "max      101                    nan",
        "min        0                    nan",
        "mean 6.34e-06                    nan",
        "std_dev  0.02524                    nan",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff%",
        "--------- ---------",
        "0.1 6.309e-06",
        "0.01 6.309e-06",
        "0.001 6.309e-06",
        "0.0001 6.309e-06",
        "1e-05 6.309e-06",
        "1e-06 6.309e-06",
        "1e-07 6.309e-06",
        "0.0 6.309e-06",
        "Extension HDU 2 (PIXELDQ, 1):",
        "Values in a and b",
        "Quantity      a        b",
        "---------- --------- -------",
        "zeros   1056766 1056768",
        "nan         0       0",
        "no-nan   1056768 1056768",
        "min_value         0       0",
        "max_value         1       0",
        "mean_value 1.893e-06       0",
        "Difference stats: abs(a - b)",
        "Quantity  abs_diff no_rel_stats_available",
        "-------- --------- ----------------------",
        "max         1                    nan",
        "min         0                    nan",
        "mean 1.893e-06                    nan",
        "std_dev  0.001376                    nan",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff%",
        "--------- ---------",
        "0.1 0.0001893",
        "0.01 0.0001893",
        "0.001 0.0001893",
        "0.0001 0.0001893",
        "1e-05 0.0001893",
        "1e-06 0.0001893",
        "1e-07 0.0001893",
        "0.0 0.0001893",
    ]

    assert result is False
    assert report == expected_report


def test_nan_in_sci(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    nan_in_sci = mock_rampfiles[3]
    diff = FITSDiff(nan_in_sci, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> nan_in_sci_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Values in a and b",
        "Quantity     a        b",
        "---------- -------- --------",
        "zeros 30646272 31703040",
        "nan  1056768        0",
        "no-nan 30646272 31703040",
        "min_value        0        0",
        "max_value        0        0",
        "mean_value        0        0",
        "Difference stats: abs(a - b)",
        "Quantity abs_diff no_rel_stats_available",
        "-------- -------- ----------------------",
        "max        0                    nan",
        "min        0                    nan",
        "mean        0                    nan",
        "std_dev        0                    nan",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff%",
        "--------- ---------",
        "0.1         0",
        "0.01         0",
        "0.001         0",
        "0.0001         0",
        "1e-05         0",
        "1e-06         0",
        "1e-07         0",
        "0.0         0",
    ]

    assert result is False
    assert report == expected_report


def test_change_tols(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Inflate all tolerances so only file names are different
    fitsdiff_default_kwargs["rtol"] = 1e4
    fitsdiff_default_kwargs["atol"] = 1e5
    diff = FITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result1 = diff.identical
    report1 = report_to_list(diff.report())
    # The expected result is False
    # The report should look like this
    expected_report1 = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
    ]
    # Only fail the PIXELDQ extension
    extension_tolerances = {
        "sci": {"rtol": 1e2, "atol": 1e3},
        "pixeldq": {"rtol": 1e-3, "atol": 1e-5},
        "default": {"rtol": 1e5, "atol": 1e7},
    }
    fitsdiff_default_kwargs["extension_tolerances"] = extension_tolerances
    diff2 = FITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result2 = diff.identical
    report2 = report_to_list(diff2.report(), from_line=10)
    # The expected result is False
    # The report should look like this
    expected_report2 = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1.0e+05, Absolute tolerance: 1.0e+07",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 2 (PIXELDQ, 1):",
        "Relative tolerance: 1.0e-03, Absolute tolerance: 1.0e-05",
        "Values in a and b",
        "Quantity      a        b",
        "---------- --------- -------",
        "zeros   1056766 1056768",
        "nan         0       0",
        "no-nan   1056768 1056768",
        "min_value         0       0",
        "max_value         1       0",
        "mean_value 1.893e-06       0",
        "Difference stats: abs(a - b)",
        "Quantity  abs_diff no_rel_stats_available",
        "-------- --------- ----------------------",
        "max         1                    nan",
        "min         0                    nan",
        "mean 1.893e-06                    nan",
        "std_dev  0.001376                    nan",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff%",
        "--------- ---------",
        "0.1 0.0001893",
        "0.01 0.0001893",
        "0.001 0.0001893",
        "0.0001 0.0001893",
        "1e-05 0.0001893",
        "1e-06 0.0001893",
        "1e-07 0.0001893",
        "0.0 0.0001893",
    ]

    assert result1 is False
    assert report1 == expected_report1
    assert result2 is False
    assert report2 == expected_report2


@pytest.fixture(scope="module")
def mock_table(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("stfitsdiff_table")
    truth = tmp_dir / "truth_x1d.fits"
    keyword_mod = tmp_dir / "keyword_mod_x1d.fits"
    data_mod = tmp_dir / "data_mod_x1d.fits"
    nan_in_data = tmp_dir / "nan_in_data_x1d.fits"
    nan_column = tmp_dir / "nan_column_x1d.fits"

    with datamodels.SpecModel() as spec:
        spec.spectral_order = 2
        spec.meta.soss_extract1d.type = "OBSERVATION"
        spec.meta.soss_extract1d.factor = np.nan
        spec.spec_table = np.zeros((100,), dtype=datamodels.SpecModel().spec_table.dtype)
        spec.spec_table["WAVELENGTH"] = np.arange(100) * 0.1
        spec.spec_table["FLUX"] = np.ones(100)
        spec.spec_table["DQ"] = np.ones(100, dtype=int)
        spec.save(truth)

        # Change data in a couple of places and save
        spec.spec_table["FLUX"][20] = 100.0
        spec.spec_table["FLUX"][50] = 1e-5
        spec.save(data_mod)
        # return to truth value
        spec.spec_table["FLUX"][20] = 1.0
        spec.spec_table["FLUX"][50] = 1.0

        # Set a few values to zeros and a few to nan
        spec.spec_table["FLUX"][10] = 0.0
        spec.spec_table["FLUX"][15] = 0.0
        spec.spec_table["FLUX"][20] = np.nan
        spec.spec_table["FLUX"][25] = np.nan
        spec.save(nan_in_data)
        # return to truth value
        spec.spec_table["FLUX"][10] = 1.0
        spec.spec_table["FLUX"][15] = 1.0
        spec.spec_table["FLUX"][20] = 1.0
        spec.spec_table["FLUX"][25] = 1.0

        # Set a whole column to nan
        spec.spec_table["WAVELENGTH"][:] = np.nan
        spec.save(nan_column)
        # return to truth
        spec.spec_table["WAVELENGTH"] = np.arange(100) * 0.1

        # Add a keyword and save
        spec.meta.cal_step.flat_field = "SKIPPED"
        spec.save(keyword_mod)

    return truth, keyword_mod, data_mod, nan_in_data, nan_column


def test_identical_tables(mock_table):
    truth = mock_table[0]
    diff = FITSDiff(truth, truth)
    assert diff.identical, diff.report()


def test_table_keyword_change(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    keyword_mod = mock_table[1]
    diff = FITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Headers have different number of cards:",
        "a: 10",
        "b: 8",
        "Extra keyword 'S_FLAT' in a: 'SKIPPED'",
        "Inconsistent duplicates of keyword ''      :",
        "Occurs 2 time(s) in a, 1 times in (b)",
        "Keyword FILENAME has different values:",
        "a> keyword_mod_x1d.fits",
        "b> truth_x1d.fits",
    ]

    assert result is False
    assert report == expected_report


def test_table_data_mod(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    data_mod = mock_table[2]
    diff = FITSDiff(data_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> data_mod_x1d.fits",
        "b> truth_x1d.fits",
        "Extension HDU 1 (EXTRACT1D, 1):",
        "Found 2 different table data element(s). Reporting percentages above respective tolerances:",
        "- absolute .... 0.1111%",
        "- relative .... 0.1111%",
        "Values in a and b",
        "col_name zeros_a zeros_b nan_a nan_b no-nan_a no-nan_b max_a max_b min_a min_b mean_a mean_b",
        "-------- ------- ------- ----- ----- -------- -------- ----- ----- ----- ----- ------ ------",
        "FLUX       0       0     0     0      100      100   100     1 1e-05     1   1.98      1",
        "Difference stats: abs(a - b)",
        "col_name dtype abs_diffs abs_max abs_mean abs_std rel_diffs rel_max rel_mean rel_std",
        "-------- ----- --------- ------- -------- ------- --------- ------- -------- -------",
        "FLUX    f8         2      99        1    9.85         2      99       50      49",
    ]

    assert result is False
    assert report == expected_report


def test_table_nan_in_data(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    nan_in_data = mock_table[3]
    diff = FITSDiff(nan_in_data, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> nan_in_data_x1d.fits",
        "b> truth_x1d.fits",
        "Extension HDU 1 (EXTRACT1D, 1):",
        "Found 100 different table data element(s). Reporting percentages above respective tolerances:",
        "- absolute .... 5.556%",
        "* Unable to calculate relative differences and stats due to data types",
        "Values in a and b",
        "col_name zeros_a zeros_b nan_a nan_b no-nan_a no-nan_b max_a max_b min_a min_b mean_a mean_b",
        "-------- ------- ------- ----- ----- -------- -------- ----- ----- ----- ----- ------ ------",
        "FLUX       2       0     2     0       98      100     1     1     0     1 0.9796      1",
        "Difference stats: abs(a - b)",
        "col_name dtype abs_diffs abs_max abs_mean abs_std rel_diffs rel_max rel_mean rel_std",
        "-------- ----- --------- ------- -------- ------- --------- ------- -------- -------",
        "FLUX    f8         0       0        0       0         0       0        0       0",
    ]

    assert result is False
    assert report == expected_report


def test_table_nan_column(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    nan_column = mock_table[4]
    diff = FITSDiff(nan_column, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> nan_column_x1d.fits",
        "b> truth_x1d.fits",
        "Extension HDU 1 (EXTRACT1D, 1):",
        "Found 100 different table data element(s). Reporting percentages above respective tolerances:",
        "- absolute .... 5.556%",
        "* Unable to calculate relative differences and stats due to data types",
        "Values in a and b",
        "col_name  zeros_a zeros_b nan_a nan_b no-nan_a no-nan_b max_a max_b min_a min_b mean_a mean_b",
        "---------- ------- ------- ----- ----- -------- -------- ----- ----- ----- ----- ------ ------",
        "WAVELENGTH       0       1   100     0        0      100   nan   9.9   nan     0    nan   4.95",
        "Difference stats: abs(a - b)",
        "col_name  dtype abs_diffs abs_max abs_mean abs_std rel_diffs rel_max rel_mean rel_std",
        "---------- ----- --------- ------- -------- ------- --------- ------- -------- -------",
        "WAVELENGTH    f8         0       0        0       0         0       0        0       0",
    ]

    assert result is False
    assert report == expected_report
