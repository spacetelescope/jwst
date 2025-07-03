"""Tests for STFitsDiff"""

import numpy as np
import pytest

from stdatamodels.jwst import datamodels
from jwst.regtest.st_fitsdiff import STFITSDiff
from astropy.io.fits.diff import FITSDiff
from astropy.io import fits


@pytest.fixture(scope="module")
def mock_rampfiles(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("stfitsdiff_ramps")
    truth = tmp_dir / "truth_ramp.fits"
    keyword_mod = tmp_dir / "keyword_mod_ramp.fits"
    sci_mod = tmp_dir / "sci_mod_ramp.fits"
    nan_in_sci = tmp_dir / "nan_in_sci_ramp.fits"
    diff_exts = tmp_dir / "ext_removed_ramp.fits"
    diff_dim = tmp_dir / "diff_dim.fits"

    nints = 2
    ngroups = 2
    nrows = 3
    ncols = 3
    data = np.ones(shape=(nints, ngroups, nrows, ncols))
    data[0, 0, ...] += np.arange(3) * 0.001
    data[1, 1, ...] += np.arange(3) * 0.0001
    pixdq = np.ones(shape=(nrows, ncols), dtype=int)
    gdq = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=int)
    with datamodels.open(data=data) as model:
        model.meta.observation.date = "2015-10-13"
        model.meta.instrument.gwa_tilt = 37.0610
        model.save(diff_exts)

    with datamodels.RampModel(data=data, pixeldq=pixdq, groupdq=gdq) as rampmodel:
        rampmodel.meta.observation.date = "2015-10-13"
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
        rampmodel.data[1, 1, 1, 0] = 1.5
        rampmodel.data[1, 1, 1, 2] = 0.5
        rampmodel.data[1, 1, 1, 2] = 0.01
        rampmodel.data[1, 1, 2, 0] = 1.0005
        rampmodel.data[1, 1, 2, 1] = 0.0005
        rampmodel.data[1, 1, 2, 2] = 1.0000005
        rampmodel.save(sci_mod)
        # return to truth value
        rampmodel.data = data

        # To test nans and zeros
        rampmodel.data[0, 0, ...] = np.nan
        rampmodel.data[1, 1, ...] = 0.0
        rampmodel.save(nan_in_sci)
        # return to truth value
        rampmodel.data = data

        # Test different shapes
        rampmodel.data = np.ones(shape=(2, 2, 4, 4))
        rampmodel.save(diff_dim)

    return truth, keyword_mod, sci_mod, nan_in_sci, diff_exts, diff_dim, tmp_dir


def report_to_list(report, from_line=11, report_pixel_loc_diffs=False):
    rsplit = report.split(sep="\n")
    rsplit = [rs.strip() for rs in rsplit if rs]
    # Remove the lines of fits diff version and comparison filenames, as well
    # HDUs, keywords and columns not compared, and the maximum number of
    # different data values to be reported and the abs and rel tolerances
    if not report_pixel_loc_diffs:
        return rsplit[from_line:]
    else:
        # Match the astropy report
        # Remove the legend '* Pixel indices below are 1-based.'
        report = [line for line in rsplit if "Pixel indices below are 1-based." not in line]
        primary_diffs, pixidx = None, None
        for idx, line in enumerate(report):
            if "Values" in line or "Found" in line:
                if primary_diffs is None:
                    primary_diffs = idx
            if "differs" in line or "differ:" in line and pixidx is None:
                pixidx = idx
                break
        streport = report[from_line:pixidx]
        # If primary_diffs is still None, means that no further comparison
        # was made and so no stats were calculated
        if primary_diffs is None:
            pixelreport = report[from_line:pixidx]
        else:
            pixelreport = report[from_line:primary_diffs]
        pixelreport.append("Data contains differences:")
        pixelreport.extend(report[pixidx:])
        return streport, pixelreport


def test_identical(mock_rampfiles):
    truth = mock_rampfiles[0]
    apdiff = FITSDiff(truth, truth)
    diff = STFITSDiff(truth, truth)
    assert apdiff.identical, apdiff.report()
    assert diff.identical, diff.report()


def test_filename(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    tmp_dir = mock_rampfiles[6]
    hdu = fits.open(truth)
    del hdu[0].header["FILENAME"]
    test1 = tmp_dir / "test1.fits"
    hdu.writeto(test1)
    hdu.close()
    diff = STFITSDiff(test1, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report(), from_line=11)
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Headers have different number of cards:",
        "a: 10",
        "b: 11",
        "Extra keyword 'FILENAME' in b: 'truth_ramp.fits'",
    ]
    assert result is False
    assert report == expected_report

    diff = STFITSDiff(truth, test1, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report(), from_line=11)
    expected_report2 = [
        "Primary HDU:",
        "Headers contain differences:",
        "Headers have different number of cards:",
        "a: 11",
        "b: 10",
        "Extra keyword 'FILENAME' in a: 'truth_ramp.fits'",
    ]
    assert result is False
    assert report == expected_report2


def test_diff_exts(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    diff_exts = mock_rampfiles[4]
    apdiff = FITSDiff(diff_exts, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    diff = STFITSDiff(diff_exts, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    asptropy_expected_report = [
        "Files contain different numbers of HDUs:",
        "a: 2",
        "b: 5",
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword DATAMODL has different values:",
        "a> JwstDataModel",
        "b> RampModel",
        "Keyword FILENAME has different values:",
        "a> ext_removed_ramp.fits",
        "b> truth_ramp.fits",
    ]
    expected_report = [
        "Files contain different numbers of HDUs:",
        "a: 2, ['ASDF', 'PRIMARY']",
        "b: 5, ['ASDF', 'GROUPDQ', 'PIXELDQ', 'PRIMARY', 'SCI']",
        "Common HDUs: ['ASDF', 'PRIMARY']",
        "Missing HDUs: ['GROUPDQ', 'PIXELDQ', 'SCI']",
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword DATAMODL has different values:",
        "a> JwstDataModel",
        "b> RampModel",
        "Keyword FILENAME has different values:",
        "a> ext_removed_ramp.fits",
        "b> truth_ramp.fits",
    ]

    assert result == apresult
    assert apreport == asptropy_expected_report
    assert report == expected_report


def test_keyword_change(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    keyword_mod = mock_rampfiles[1]
    apdiff = FITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    diff = STFITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())

    # change the tolerance for the primary header
    extension_tolerances = {"primary": {"rtol": 1, "atol": 2}}
    fitsdiff_default_kwargs["extension_tolerances"] = extension_tolerances
    diff2 = STFITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    result2 = diff2.identical
    report2 = report_to_list(diff2.report(), from_line=10)
    expected_report2 = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1, Absolute tolerance: 2",
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

    assert result == apresult
    assert report == apreport
    assert result2 is False
    assert report2 == expected_report2


def test_ignore(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Ignoring these keywords (in addition to 'CAL_VER', 'CAL_VCS', 'CRDS_VER',
    # 'CRDS_CTX', 'NAXIS1', 'TFORM*'), files should be identical
    fitsdiff_default_kwargs["ignore_keywords"].extend(["FILENAME", "DATE"])
    # Ignoring these extensions but with indices (indices only work for STFITSDIFF)
    fitsdiff_default_kwargs["ignore_hdus"].extend([1, 2])
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()

    # Populate ignore_hdu_patterns, files should be identical
    fitsdiff_default_kwargs["ignore_hdus"].extend(["S*"])
    apdiff = FITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    assert apdiff.identical, apdiff.report()
    assert diff.identical, diff.report()

    # Ignoring these extensions, files should be identical
    fitsdiff_default_kwargs["ignore_hdus"] = ["ASDF", "PIXELDQ", "SCI"]
    apdiff = FITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    assert apdiff.identical, apdiff.report()
    assert diff.identical, diff.report()


def test_array_diffs(mock_rampfiles, fitsdiff_default_kwargs):
    truth_file = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    apdiff = FITSDiff(sci_mod, truth_file, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    # The report from astropy crashes due to dividing by 0
    diff = STFITSDiff(sci_mod, truth_file, **fitsdiff_default_kwargs)
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
        "Quantity   a        b",
        "-------- ------ ---------",
        "zeros      0         0",
        "nans      0         0",
        "no-nans     36        36",
        "min 0.0005         1",
        "max    1.5     1.002",
        "mean 0.9589         1",
        "std_dev 0.2454 0.0005875",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max   0.9996   0.9995",
        "min        0        0",
        "mean  0.06918  0.06917",
        "std_dev   0.2391    0.239",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.1     8.333     8.333",
        "0.01     8.333     8.333",
        "0.001     8.333     8.333",
        "0.0001     13.89     13.89",
        "1e-05     13.89     13.89",
        "1e-06     13.89     13.89",
        "1e-07     13.89     13.89",
        "0.0     13.89     13.89",
    ]

    assert result == apresult
    assert report == expected_report

    tmp_dir = mock_rampfiles[6]
    gdq = tmp_dir / "groupdq.fits"
    fitsdiff_default_kwargs["numdiffs"] = -1
    with datamodels.RampModel(truth_file) as model:
        model.groupdq[1, 1, 1, 0] = 2
        model.groupdq[1, 1, 1, 2] = 4
        model.groupdq[1, 1, 1, 2] = 8
        model.save(gdq)
    diff = STFITSDiff(gdq, truth_file, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> groupdq.fits",
        "b> truth_ramp.fits",
        "Extension HDU 3 (GROUPDQ, 1):",
        "Values in a and b",
        "Quantity   a    b",
        "-------- ----- ---",
        "zeros     0   0",
        "nans     0   0",
        "no-nans    36  36",
        "min     1   1",
        "max     8   1",
        "mean 1.222   1",
        "std_dev 1.157   0",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max      255      255",
        "min        0        0",
        "mean       14       14",
        "std_dev    57.73    57.73",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.1     5.556     5.556",
        "0.01     5.556     5.556",
        "0.001     5.556     5.556",
        "0.0001     5.556     5.556",
        "1e-05     5.556     5.556",
        "1e-06     5.556     5.556",
        "1e-07     5.556     5.556",
        "0.0     5.556     5.556",
    ]
    assert result is False
    assert report == expected_report

    truth3d = tmp_dir / "truth3d.fits"
    test3d = tmp_dir / "test3d.fits"
    with datamodels.CubeModel((4, 4, 4)) as model:
        model.save(truth3d)
        model.data[1, 1, 1] = 4.5
        model.data[1, 1, 2] = 8.00023
        model.save(test3d)
    diff = STFITSDiff(test3d, truth3d, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> test3d.fits",
        "?  ^^",
        "b> truth3d.fits",
        "?  ^^ +",
        "Extension HDU 1 (SCI, 1):",
        "Values in a and b",
        "Quantity   a     b",
        "-------- ------ ---",
        "zeros     62  64",
        "nans      0   0",
        "no-nans     64  64",
        "min      0   0",
        "max      8   0",
        "mean 0.1953   0",
        "std_dev  1.131   0",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max        8      nan",
        "min        0      nan",
        "mean   0.1953      nan",
        "std_dev    1.131      nan",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff%",
        "--------- ---------",
        "0.1     3.125",
        "0.01     3.125",
        "0.001     3.125",
        "0.0001     3.125",
        "1e-05     3.125",
        "1e-06     3.125",
        "1e-07     3.125",
        "0.0     3.125",
    ]
    assert result is False
    assert report == expected_report

    truth2d = tmp_dir / "truth2d.fits"
    test2d = tmp_dir / "test2d.fits"
    with datamodels.ImageModel((4, 4)) as model:
        model.save(truth2d)
        model.data[1, 1] = 4.5
        model.data[1, 2] = 8.00023
        model.save(test2d)
    diff = STFITSDiff(test2d, truth2d, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> test2d.fits",
        "?  ^^",
        "b> truth2d.fits",
        "?  ^^ +",
        "Extension HDU 1 (SCI, 1):",
        "Values in a and b",
        "Quantity   a     b",
        "-------- ------ ---",
        "zeros     14  16",
        "nans      0   0",
        "no-nans     16  16",
        "min      0   0",
        "max      8   0",
        "mean 0.7813   0",
        "std_dev  2.158   0",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max        8      nan",
        "min        0      nan",
        "mean   0.7813      nan",
        "std_dev    2.158      nan",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff%",
        "--------- ---------",
        "0.1      12.5",
        "0.01      12.5",
        "0.001      12.5",
        "0.0001      12.5",
        "1e-05      12.5",
        "1e-06      12.5",
        "1e-07      12.5",
        "0.0      12.5",
    ]
    assert result is False
    assert report == expected_report


def test_nan_in_sci(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    nan_in_sci = mock_rampfiles[3]
    apdiff = FITSDiff(nan_in_sci, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(nan_in_sci, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report, pixelreport = report_to_list(diff.report(), report_pixel_loc_diffs=True)
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
        "Quantity   a        b",
        "-------- ------ ---------",
        "zeros      9         0",
        "nans      9         0",
        "no-nans     27        36",
        "min      0         1",
        "max      1     1.002",
        "mean 0.6667         1",
        "std_dev 0.4714 0.0005875",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max        1        1",
        "min        0        0",
        "mean   0.3334   0.3333",
        "std_dev   0.4715   0.4714",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.1        25        25",
        "0.01        25        25",
        "0.001        25        25",
        "0.0001        25        25",
        "1e-05        25        25",
        "1e-06        25        25",
        "1e-07        25        25",
        "0.0        25        25",
    ]

    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report

    tmp_dir = mock_rampfiles[6]
    all_sci_nan = tmp_dir / "all_sci_nan.fits"
    hdu = fits.open(truth)
    shape = np.shape(hdu[1].data)
    hdu[1].data = np.ones(shape) * np.nan
    hdu.writeto(all_sci_nan)
    hdu.close()
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = False
    diff = STFITSDiff(all_sci_nan, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report(), from_line=11)
    expected_report = [
        "Extension HDU 1 (SCI, 1):",
        "Headers contain differences:",
        "Keyword BITPIX   has different values:",
        "a> -64",
        "b> -32",
        "Values in a and b",
        "Quantity  a      b",
        "-------- --- ---------",
        "zeros   0         0",
        "nans  36         0",
        "no-nans   0        36",
        "min   -         1",
        "max   -     1.002",
        "mean   -         1",
        "std_dev   - 0.0005875",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff",
        "-------- --------",
        "max      nan",
        "min      nan",
        "mean      nan",
        "std_dev      nan",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.0       100       100",
    ]
    assert result is False
    assert report == expected_report


def test_change_tols(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Only change the abs tolerance
    fitsdiff_default_kwargs["extension_tolerances"] = {"sci": {"atol": 0.01}}
    diff1 = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result1 = diff1.identical
    report1 = report_to_list(diff1.report(), from_line=10)
    # The expected result is False
    # The report should look like this
    expected_report1 = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1e-05, Absolute tolerance: 1e-07",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Relative tolerance: 1e-05, Absolute tolerance: 0.01",
        "Values in a and b",
        "Quantity   a        b",
        "-------- ------ ---------",
        "zeros      0         0",
        "nans      0         0",
        "no-nans     36        36",
        "min 0.0005         1",
        "max    1.5     1.002",
        "mean 0.9589         1",
        "std_dev 0.2454 0.0005875",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max   0.9996   0.9995",
        "min        0        0",
        "mean  0.06918  0.06917",
        "std_dev   0.2391    0.239",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.1     8.333     8.333",
        "0.01     8.333     8.333",
        "0.001     8.333     8.333",
        "0.0001     8.333     13.89",
        "1e-05     8.333     13.89",
        "1e-06     8.333     13.89",
        "1e-07     8.333     13.89",
        "0.0     8.333     13.89",
    ]

    # Only change the rtol tolerance
    fitsdiff_default_kwargs["extension_tolerances"] = {1: {"rtol": 0.001}}
    diff2 = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result2 = diff2.identical
    report2 = report_to_list(diff2.report(), from_line=10)
    # The expected result is False
    # The report should look like this
    expected_report2 = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1e-05, Absolute tolerance: 1e-07",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Relative tolerance: 0.001, Absolute tolerance: 1e-07",
        "Values in a and b",
        "Quantity   a        b",
        "-------- ------ ---------",
        "zeros      0         0",
        "nans      0         0",
        "no-nans     36        36",
        "min 0.0005         1",
        "max    1.5     1.002",
        "mean 0.9589         1",
        "std_dev 0.2454 0.0005875",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max   0.9996   0.9995",
        "min        0        0",
        "mean  0.06918  0.06917",
        "std_dev   0.2391    0.239",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.1     8.333     8.333",
        "0.01     8.333     8.333",
        "0.001     8.333     8.333",
        "0.0001     13.89     8.333",
        "1e-05     13.89     8.333",
        "1e-06     13.89     8.333",
        "1e-07     13.89     8.333",
        "0.0     13.89     8.333",
    ]

    # Fail the SCI extension a little less bad
    extension_tolerances = {
        "sci": {"rtol": 1e-3, "atol": 1e-5},
        "pixeldq": {"rtol": 1, "atol": 2},
        "headers": {"rtol": 1e5, "atol": 1e6},
        "default": {"rtol": 1e3, "atol": 1e4},
    }
    fitsdiff_default_kwargs["extension_tolerances"] = extension_tolerances
    diff3 = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result3 = diff3.identical
    report3 = report_to_list(diff3.report(), from_line=10)
    # The expected result is False
    # The report should look like this
    expected_report3 = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1000, Absolute tolerance: 1e+04",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Relative tolerance: 0.001, Absolute tolerance: 1e-05",
        "Values in a and b",
        "Quantity   a        b",
        "-------- ------ ---------",
        "zeros      0         0",
        "nans      0         0",
        "no-nans     36        36",
        "min 0.0005         1",
        "max    1.5     1.002",
        "mean 0.9589         1",
        "std_dev 0.2454 0.0005875",
        "Difference stats: abs(b - a)",
        "Quantity abs_diff rel_diff",
        "-------- -------- --------",
        "max   0.9996   0.9995",
        "min        0        0",
        "mean  0.06918  0.06917",
        "std_dev   0.2391    0.239",
        "Percentages of difference above (tolerance + threshold)",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.1     8.333     8.333",
        "0.01     8.333     8.333",
        "0.001     8.333     8.333",
        "0.0001     13.89     8.333",
        "1e-05     13.89     8.333",
        "1e-06     13.89     8.333",
        "1e-07     13.89     8.333",
        "0.0     13.89     8.333",
    ]

    # Inflate all tolerances so only file names are different
    fitsdiff_default_kwargs["extension_tolerances"] = None
    fitsdiff_default_kwargs["rtol"] = 1e4
    fitsdiff_default_kwargs["atol"] = 1e5
    diff4 = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result4 = diff4.identical
    report4 = report_to_list(diff4.report())
    # The expected result is False but only for the file name diffs
    # The report should look like this
    expected_report4 = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
    ]
    assert result1 is False
    assert report1 == expected_report1
    assert result2 is False
    assert report2 == expected_report2
    assert result3 is False
    assert report3 == expected_report3
    assert result4 is False
    assert report4 == expected_report4


def test_diff_dim(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    diff_dim = mock_rampfiles[5]
    apdiff = FITSDiff(diff_dim, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(diff_dim, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    _, pixelreport = report_to_list(diff.report(), report_pixel_loc_diffs=True)
    assert result == apresult
    assert pixelreport == apreport


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
    apdiff = FITSDiff(truth, truth)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    diff = STFITSDiff(truth, truth)
    result = diff.identical
    report = report_to_list(diff.report())
    assert result == apresult
    assert report == apreport


def test_table_keyword_change(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    keyword_mod = mock_table[1]
    apdiff = FITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    diff = STFITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
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

    assert result == apresult
    assert report == apreport
    assert report == expected_report


def test_table_data_mod(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    data_mod = mock_table[2]
    apdiff = FITSDiff(data_mod, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(data_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report, pixelreport = report_to_list(diff.report(), report_pixel_loc_diffs=True)

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
        "- relative .... 0%",
        "Values in a and b",
        "col_name zeros_a zeros_b nan_a nan_b no-nan_a no-nan_b max_a max_b min_a min_b mean_a mean_b",
        "-------- ------- ------- ----- ----- -------- -------- ----- ----- ----- ----- ------ ------",
        "FLUX       0       0     0     0      100      100   100     1 1e-05     1   1.98      1",
        "Difference stats: abs(b - a)",
        "col_name dtype abs_diffs abs_max abs_mean abs_std rel_diffs rel_max rel_mean rel_std",
        "-------- ----- --------- ------- -------- ------- --------- ------- -------- -------",
        "FLUX    f8         2      99        1    9.85         2      99       50      49",
    ]

    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report


def test_table_nan_in_data(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    nan_in_data = mock_table[3]
    apdiff = FITSDiff(nan_in_data, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(nan_in_data, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report, pixelreport = report_to_list(diff.report(), report_pixel_loc_diffs=True)

    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> nan_in_data_x1d.fits",
        "b> truth_x1d.fits",
        "Extension HDU 1 (EXTRACT1D, 1):",
        "Found 4 different table data element(s). Reporting percentages above respective tolerances:",
        "- absolute .... 0.2222%",
        "- relative .... 0%",
        "Values in a and b",
        "col_name zeros_a zeros_b nan_a nan_b no-nan_a no-nan_b max_a max_b min_a min_b mean_a mean_b",
        "-------- ------- ------- ----- ----- -------- -------- ----- ----- ----- ----- ------ ------",
        "FLUX       2       0     2     0       98      100     1     1     0     1 0.9796      1",
        "Difference stats: abs(b - a)",
        "col_name dtype abs_diffs abs_max abs_mean abs_std rel_diffs rel_max rel_mean rel_std",
        "-------- ----- --------- ------- -------- ------- --------- ------- -------- -------",
        "FLUX    f8         0       0        0       0         0       0        0       0",
    ]

    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report


def test_table_nan_column(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    nan_column = mock_table[4]
    apdiff = FITSDiff(nan_column, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report())
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(nan_column, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report, pixelreport = report_to_list(diff.report(), report_pixel_loc_diffs=True)

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
        "- relative .... 0%",
        "Values in a and b",
        "col_name  zeros_a zeros_b nan_a nan_b no-nan_a no-nan_b max_a max_b min_a min_b mean_a mean_b",
        "---------- ------- ------- ----- ----- -------- -------- ----- ----- ----- ----- ------ ------",
        "WAVELENGTH       0       1   100     0        0      100   nan   9.9   nan     0    nan   4.95",
        "Difference stats: abs(b - a)",
        "col_name  dtype abs_diffs abs_max abs_mean abs_std rel_diffs rel_max rel_mean rel_std",
        "---------- ----- --------- ------- -------- ------- --------- ------- -------- -------",
        "WAVELENGTH    f8         0       0        0       0         0       0        0       0",
    ]

    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report
