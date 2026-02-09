"""Tests for STFitsDiff."""

import numpy as np
import pytest
from astropy.io import fits
from astropy.io.fits.diff import FITSDiff
from stdatamodels.jwst import datamodels

from jwst.regtest.st_fitsdiff import STFITSDiffBeta as STFITSDiff


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
    with datamodels.JwstDataModel(data=data) as model:
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
    """
    Turn the fitsdiff or stfitsdiff report into a comparable list of strings.

    Parameters
    ----------
    report : str
        Full string report.
    from_line : int
        The index of the line at which to start the comparison (this is to skip the lines
        that contain the astropy version and directory path).
    report_pixel_loc_diffs : bool
        If True, it is expected that the pixel or column differences are part of the
        report (at the end), as well as the ST ad hoc report (at the beginning).

    Returns
    -------
    report or streport: list
        This is either the no differences report, or the ST ad hoc report.
    pixelreport : list
        Report containing the pixel or column location of the differences, also
        only returned when report_pixel_loc_diffs is True.
    """
    rsplit = report.split(sep="\n")
    rsplit = [rs.strip() for rs in rsplit if len(rs) > 0 and "STScI Custom FITSDiff" not in rs]
    # Remove the lines of fits diff version and comparison filenames, as well
    # HDUs, keywords and columns not compared, and the maximum number of
    # different data values to be reported and the abs and rel tolerances
    report = rsplit[from_line:]
    # Remove the max absolute and max relative to pass old astropy version tests
    no_diffs = False
    new_report = []
    stidx, pixidx = False, False
    for line in report:
        if "No differences found" in line:
            no_diffs = True
            break
        elif "Maximum relative difference" in line:
            continue
        elif "Maximum absolute difference" in line:
            continue
        new_report.append(line)
    if no_diffs:
        if report_pixel_loc_diffs:
            return report, report
        else:
            return report
    report = new_report
    if not report_pixel_loc_diffs:
        return report
    else:
        # Match the astropy report
        streport, pixelreport = [], []
        for idx, line in enumerate(report):
            # Ignore the ST added legends
            if "These values are calculated" in line:
                continue
            elif "Pixel indices below are 1-based." in line:
                continue
            elif "Found" in line and "table data element(s)" in line:
                continue
            elif "(atol, rtol)" in line:
                continue
            elif "Primary" in line or "Extension" in line:
                stidx = True
                pixidx = True
            elif "Values in" in line or "Reporting percentages" in line:
                stidx = True
                pixidx = False
            elif "differs" in line or "differ:" in line or "Extra column" in line:
                stidx = False
                pixidx = True
            if stidx:
                streport.append(line)
            if pixidx:
                pixelreport.append(line)

        return streport, pixelreport


def test_identical(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    apdiff = FITSDiff(truth, truth)
    diff = STFITSDiff(truth, truth, **fitsdiff_default_kwargs)
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
    report = report_to_list(diff.report())
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


def test_filename_reversed(mock_rampfiles, fitsdiff_default_kwargs):
    # Reverse the order of files to compare
    truth = mock_rampfiles[0]
    tmp_dir = mock_rampfiles[6]
    diff = STFITSDiff(truth, tmp_dir / "test1.fits", **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
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
    apreport = report_to_list(apdiff.report(), from_line=11)
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
        "Files contain different HDUs:",
        "a: 2 HDUs: ['ASDF', 'PRIMARY']",
        "b: 5 HDUs: ['ASDF', 'GROUPDQ', 'PIXELDQ', 'PRIMARY', 'SCI']",
        "Common HDUs: ['ASDF', 'PRIMARY']",
        "HDUs in a but not b: []",
        "HDUs in b but not a: ['GROUPDQ', 'PIXELDQ', 'SCI']",
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


def test_rm_sci(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    tmp_dir = mock_rampfiles[6]
    hdu = fits.open(truth)
    del hdu[1]
    rmsci = tmp_dir / "rm_sci.fits"
    hdu.writeto(rmsci)
    hdu.close()
    diff = STFITSDiff(rmsci, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    expected_report = [
        "Files contain different HDUs:",
        "a: 4 HDUs: ['ASDF', 'GROUPDQ', 'PIXELDQ', 'PRIMARY']",
        "b: 5 HDUs: ['ASDF', 'GROUPDQ', 'PIXELDQ', 'PRIMARY', 'SCI']",
        "Common HDUs: ['ASDF', 'GROUPDQ', 'PIXELDQ', 'PRIMARY']",
        "HDUs in a but not b: []",
        "HDUs in b but not a: ['SCI']",
        "No differences found between common HDUs.",
    ]
    assert result is False
    assert report == expected_report


def test_keyword_change(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    keyword_mod = mock_rampfiles[1]
    apdiff = FITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report(), from_line=11)
    diff = STFITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    assert result == apresult
    assert report == apreport


def test_keyword_change_tol(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    keyword_mod = mock_rampfiles[1]
    # change the tolerance for the primary header
    extension_tolerances = {"primary": {"rtol": 1, "atol": 2}}
    fitsdiff_default_kwargs["extension_tolerances"] = extension_tolerances
    diff = STFITSDiff(keyword_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report(), from_line=10)
    expected_report = [
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
    assert result is False
    assert report == expected_report


def test_ignore_kw(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Ignoring these keywords (in addition to 'CAL_VER', 'CAL_VCS', 'CRDS_VER',
    # 'CRDS_CTX', 'NAXIS1', 'TFORM*'), files should be identical
    fitsdiff_default_kwargs["ignore_keywords"].extend(["FILENAME", "DATE"])
    # Ignoring these extensions but with indices (indices only work for STFITSDIFF)
    fitsdiff_default_kwargs["ignore_hdus"].extend([1, 2])
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


def test_ignore_wild_card(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Populate ignore_hdu_patterns, files should be identical
    fitsdiff_default_kwargs["ignore_keywords"].extend(["FILENAME", "DATE"])
    fitsdiff_default_kwargs["ignore_hdus"].extend(["S*"])
    apdiff = FITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    assert apdiff.identical, apdiff.report()
    assert diff.identical, diff.report()


def test_ignore_hdu(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Ignoring these extensions, files should be identical
    fitsdiff_default_kwargs["ignore_keywords"].extend(["FILENAME", "DATE"])
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
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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


def test_array4d_diffs(mock_rampfiles, fitsdiff_default_kwargs):
    truth_file = mock_rampfiles[0]
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
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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


def test_array3d_diffs(mock_rampfiles, fitsdiff_default_kwargs):
    tmp_dir = mock_rampfiles[6]
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
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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


def test_array2d_diffs(mock_rampfiles, fitsdiff_default_kwargs):
    tmp_dir = mock_rampfiles[6]
    truth2d = tmp_dir / "truth2d.fits"
    test2d = tmp_dir / "test2d.fits"
    with datamodels.ImageModel((4, 4)) as model:
        model.save(truth2d)
        model.data[1, 1] = 4.5
        model.data[1, 2] = 8.00023
        model.save(test2d)
    fitsdiff_default_kwargs["numdiffs"] = -1
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
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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
    fitsdiff_default_kwargs["numdiffs"] = 5
    apdiff = FITSDiff(nan_in_sci, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report(), from_line=11)
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
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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
    assert "Pixel indices below are 1-based." in diff.report()
    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report


def test_allnan_sci(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    tmp_dir = mock_rampfiles[6]
    all_sci_nan = tmp_dir / "all_sci_nan.fits"
    hdu = fits.open(truth)
    shape = np.shape(hdu[1].data)
    hdu[1].data = np.ones(shape) * np.nan
    hdu.writeto(all_sci_nan)
    hdu.close()
    diff = STFITSDiff(all_sci_nan, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    expected_report = [
        "Extension HDU 1 (SCI, 1):",
        "Headers contain differences:",
        "Keyword BITPIX   has different values:",
        "a> -64",
        "b> -32",
        "Data contains differences:",
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
        "Percentages of difference above threshold",
        "threshold abs_diff% rel_diff%",
        "--------- --------- ---------",
        "0.0       100       100",
    ]
    assert result is False
    assert report == expected_report


def test_change_sci_atol(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Only change the abs tolerance
    fitsdiff_default_kwargs["extension_tolerances"] = {"sci": {"atol": 0.01}}
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report(), from_line=10)
    # The expected result is False
    # The report should look like this
    expected_report = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1e-05, Absolute tolerance: 1e-07",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Relative tolerance: 1e-05, Absolute tolerance: 0.01",
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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
    assert result is False
    assert report == expected_report


def test_change_sci_rtol(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Only change the rtol tolerance
    fitsdiff_default_kwargs["extension_tolerances"] = {1: {"rtol": 0.001}}
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report(), from_line=10)
    # The expected result is False
    # The report should look like this
    expected_report = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1e-05, Absolute tolerance: 1e-07",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Relative tolerance: 0.001, Absolute tolerance: 1e-07",
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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
    assert result is False
    assert report == expected_report


def test_change_all_tols(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Fail the SCI extension a little less bad
    extension_tolerances = {
        "sci": {"rtol": 1e-3, "atol": 1e-5},
        "pixeldq": {"rtol": 1, "atol": 2},
        "headers": {"rtol": 1e5, "atol": 1e6},
        "default": {"rtol": 1e3, "atol": 1e4},
    }
    fitsdiff_default_kwargs["extension_tolerances"] = extension_tolerances
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report(), from_line=10)
    # The expected result is False
    # The report should look like this
    expected_report = [
        "Extension HDU 0 (PRIMARY, 1):",
        "Relative tolerance: 1000, Absolute tolerance: 1e+04",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
        "Extension HDU 1 (SCI, 1):",
        "Relative tolerance: 0.001, Absolute tolerance: 1e-05",
        "Data contains differences:",
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
        "Percentages of difference above threshold",
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
    assert result is False
    assert report == expected_report


def test_change_tols(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    sci_mod = mock_rampfiles[2]
    # Inflate all tolerances so only file names are different
    fitsdiff_default_kwargs["extension_tolerances"] = None
    fitsdiff_default_kwargs["rtol"] = 1e4
    fitsdiff_default_kwargs["atol"] = 1e5
    diff = STFITSDiff(sci_mod, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    # The expected result is False but only for the file name diffs
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> sci_mod_ramp.fits",
        "b> truth_ramp.fits",
    ]
    assert result is False
    assert report == expected_report


def test_diff_dim(mock_rampfiles, fitsdiff_default_kwargs):
    truth = mock_rampfiles[0]
    diff_dim = mock_rampfiles[5]
    apdiff = FITSDiff(diff_dim, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report(), from_line=11)
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
    diff_column = tmp_dir / "diff_column_x1d.fits"

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

        # Different number of columns
        spec.spec_table = np.zeros((105,), dtype=datamodels.SpecModel().spec_table.dtype)
        spec.spec_table["WAVELENGTH"] = np.arange(105) * 0.1
        spec.spec_table["FLUX"] = np.ones(105)
        spec.spec_table["DQ"] = np.ones(105, dtype=int)
        spec.save(diff_column)
        # return to truth
        spec.spec_table = np.zeros((100,), dtype=datamodels.SpecModel().spec_table.dtype)
        spec.spec_table["WAVELENGTH"] = np.arange(100) * 0.1
        spec.spec_table["FLUX"] = np.ones(100)
        spec.spec_table["DQ"] = np.ones(100, dtype=int)

        # Add a keyword and save
        spec.meta.cal_step.flat_field = "SKIPPED"
        spec.save(keyword_mod)

    testing_tables = [
        truth,
        keyword_mod,
        data_mod,
        nan_in_data,
        nan_column,
        diff_column,
        tmp_dir,
    ]

    return testing_tables


def test_identical_tables(mock_table):
    truth = mock_table[0]
    apdiff = FITSDiff(truth, truth)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report(), from_line=11)
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
    apreport = report_to_list(apdiff.report(), from_line=11)
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
    apreport = report_to_list(apdiff.report(), from_line=11)
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
        "Data contains differences:",
        "Values in a and b",
        "col_name zeros_a_b nan_a_b no-nan_a_b       max_a_b             min_a_b             mean_a_b",
        "-------- --------- ------- ---------- ------------------- ------------------- -------------------",
        "FLUX       0 0     0 0    100 100       100         1     1e-05         1      1.98         1",
        "Difference stats for non-NaN diffs that fail the [atol, rtol] test: abs(b - a)",
        "col_name dtype rel_diffs rel_max rel_mean rel_std",
        "-------- ----- --------- ------- -------- -------",
        "FLUX    f8         2      99       50      49",
        "Columns ['BACKGROUND', 'BKGD_ERROR', 'BKGD_VAR_FLAT', 'BKGD_VAR_POISSON', "
        "'BKGD_VAR_RNOISE', 'DQ', 'FLUX_ERROR', 'FLUX_VAR_FLAT', "
        "'FLUX_VAR_POISSON', 'FLUX_VAR_RNOISE', 'NPIXELS', 'SB_ERROR', "
        "'SB_VAR_FLAT', 'SB_VAR_POISSON', 'SB_VAR_RNOISE', 'SURF_BRIGHT', "
        "'WAVELENGTH'] are identical",
    ]
    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report


def test_table_nan_in_data(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    nan_in_data = mock_table[3]
    apdiff = FITSDiff(nan_in_data, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report(), from_line=11)
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
        "Data contains differences:",
        "Values in a and b",
        "col_name zeros_a_b nan_a_b no-nan_a_b       max_a_b             min_a_b             mean_a_b",
        "-------- --------- ------- ---------- ------------------- ------------------- -------------------",
        "FLUX       2 0     2 0     98 100         1         1         0         1    0.9796         1",
        "Difference stats for non-NaN diffs that fail the [atol, rtol] test: abs(b - a)",
        "col_name dtype rel_diffs rel_max rel_mean rel_std",
        "-------- ----- --------- ------- -------- -------",
        "FLUX    f8         2       1        1       0",
        "Columns ['BACKGROUND', 'BKGD_ERROR', 'BKGD_VAR_FLAT', 'BKGD_VAR_POISSON', "
        "'BKGD_VAR_RNOISE', 'DQ', 'FLUX_ERROR', 'FLUX_VAR_FLAT', "
        "'FLUX_VAR_POISSON', 'FLUX_VAR_RNOISE', 'NPIXELS', 'SB_ERROR', "
        "'SB_VAR_FLAT', 'SB_VAR_POISSON', 'SB_VAR_RNOISE', 'SURF_BRIGHT', "
        "'WAVELENGTH'] are identical",
    ]
    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report


def test_table_nan_column(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    nan_column = mock_table[4]
    apdiff = FITSDiff(nan_column, truth, **fitsdiff_default_kwargs)
    apresult = apdiff.identical
    apreport = report_to_list(apdiff.report(), from_line=11)
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
        "Data contains differences:",
        "Values in a and b",
        "col_name  zeros_a_b nan_a_b no-nan_a_b       max_a_b             min_a_b             mean_a_b",
        "---------- --------- ------- ---------- ------------------- ------------------- -------------------",
        "WAVELENGTH       0 1   100 0      0 100       nan       9.9       nan         0       nan      4.95",
        "Difference stats for non-NaN diffs that fail the [atol, rtol] test: abs(b - a)",
        "col_name  dtype rel_diffs rel_max rel_mean rel_std",
        "---------- ----- --------- ------- -------- -------",
        "WAVELENGTH    f8         0     nan      nan     nan",
        "Columns ['BACKGROUND', 'BKGD_ERROR', 'BKGD_VAR_FLAT', 'BKGD_VAR_POISSON', "
        "'BKGD_VAR_RNOISE', 'DQ', 'FLUX', 'FLUX_ERROR', 'FLUX_VAR_FLAT', "
        "'FLUX_VAR_POISSON', 'FLUX_VAR_RNOISE', 'NPIXELS', 'SB_ERROR', "
        "'SB_VAR_FLAT', 'SB_VAR_POISSON', 'SB_VAR_RNOISE', 'SURF_BRIGHT'] are "
        "identical",
    ]
    assert result == apresult
    assert pixelreport == apreport
    assert report == expected_report


def test_table_diff_column(mock_table, fitsdiff_default_kwargs):
    truth = mock_table[0]
    diff_column = mock_table[5]
    diff = STFITSDiff(diff_column, truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    # The expected result is False
    # The report should look like this
    expected_report = [
        "Primary HDU:",
        "Headers contain differences:",
        "Keyword FILENAME has different values:",
        "a> diff_column_x1d.fits",
        "b> truth_x1d.fits",
        "Extension HDU 1 (EXTRACT1D, 1):",
        "Headers contain differences:",
        "Keyword NAXIS2   has different values:",
        "a> 105",
        "b> 100",
        "Data contains differences:",
        "Table rows differ:",
        "a: 105",
        "b: 100",
        "No further data comparison performed.",
    ]
    assert result is False
    assert report == expected_report


def test_table_pq_coltype(mock_table, fitsdiff_default_kwargs):
    tmp_dir = mock_table[6]
    diff_coltype_truth = tmp_dir / "diff_coltype_truth.fits"
    diff_coltype = tmp_dir / "diff_coltype.fits"

    arr1 = [[0, 1], [2, 3, 4]]
    arr2 = [[0.0, 1.0], [2.0, 3.0, np.nan, np.nan]]
    c1 = fits.Column(name="col_1", array=arr1, format="PI(2)")
    c2 = fits.Column(name="col_2", array=arr2, format="QD(2)")
    tab = fits.BinTableHDU.from_columns([c1, c2], name="test")
    outfile = fits.HDUList()
    outfile.append(tab)
    outfile.writeto(diff_coltype_truth)

    arr1 = [[0, 11], [2, 3, 4]]
    arr2 = [[0.0, 1.0], [np.nan, 23.0, 4.0, np.nan]]
    c1 = fits.Column(name="col_1", array=arr1, format="PI(2)")
    c2 = fits.Column(name="col_2", array=arr2, format="QD(2)")
    tab = fits.BinTableHDU.from_columns([c1, c2], name="test")
    outfile = fits.HDUList()
    outfile.append(tab)
    outfile.writeto(diff_coltype)

    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(diff_coltype, diff_coltype_truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    # The expected result is False
    # The report should look like this
    expected_report = [
        "Extension HDU 1 (TEST, 1):",
        "Data contains differences:",
        "Found 2 different table data element(s).",
        "2 failed the (atol, rtol) test",
        "Values in a and b",
        "col_name zeros_a_b nan_a_b no-nan_a_b    max_a_b       min_a_b       mean_a_b",
        "-------- --------- ------- ---------- ------------- ------------- -------------",
        "col_1       1 1     0 0        5 5     11      4      0      0      4      2",
        "col_2       1 1     2 2        4 4     23      3      0      0      7    1.5",
        "Difference stats for non-NaN diffs that fail the [atol, rtol] test: abs(b - a)",
        "col_name dtype  rel_diffs rel_max rel_mean rel_std",
        "-------- ------ --------- ------- -------- -------",
        "col_1 object         1      10       10       0",
        "col_2 object         1      20       20       0",
        "* Pixel indices below are 1-based.",
        "Column col_1 data differs in row 0:",
        "at [1]:",
        "a> 11",
        "b> 1",
        "Column col_2 data differs in row 1:",
        "at [0]:",
        "a> nan",
        "b> 2.0",
        "at [1]:",
        "a> 23.0",
        "? -",
        "b> 3.0",
        "at [2]:",
        "a> 4.0",
        "b> nan",
        "2 different table data element(s) found (50.00% different).",
    ]
    assert result is False
    assert report == expected_report
    # Make sure pixels differences aren't reported when report_pixel_loc_diffs isn't True
    del fitsdiff_default_kwargs["report_pixel_loc_diffs"]
    diff = STFITSDiff(diff_coltype, diff_coltype_truth, **fitsdiff_default_kwargs)
    assert "data differs in row" not in report


def test_hdulists_different_extnames(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    data_array = np.arange((1000)).reshape((10, 10, 10))
    a.append(fits.ImageHDU(data=data_array, name="EXTENSION1"))
    b.append(fits.ImageHDU(data=data_array, name="EXTENSION2"))
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "a: <HDUList object at " in report
    assert "b: <HDUList object at " in report
    assert "HDUs in a but not b: ['EXTENSION1']" in report
    assert "HDUs in b but not a: ['EXTENSION2']" in report
    assert "No differences found between common HDUs." in report


def test_hdulists_with_ignores(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    data_array = np.arange((1000)).reshape((10, 10, 10))
    a.append(fits.ImageHDU(data=data_array, name="EXTENSION1"))
    b.append(fits.ImageHDU(data=data_array, name="EXTENSION2"))
    fitsdiff_default_kwargs["ignore_hdus"].extend(["EXT*"])
    fitsdiff_default_kwargs["ignore_comments"] = ["EXTNAME"]
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "HDU(s) not to be compared:\n  EXT*\n" in report
    assert "Keyword(s) whose comments are not to be compared:\n  EXTNAME\n" in report


def test_hdus_different_array_sizes(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    data_array = np.arange((1000)).reshape((10, 10, 10))
    a.append(fits.ImageHDU(data=data_array, name="SCI"))
    b.append(fits.ImageHDU(data=data_array[:-1, :-1], name="SCI"))
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "Data dimensions differ:" in report


def test_hdus_different_headers(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    data_array = np.arange((1000)).reshape((10, 10, 10))
    hdu1 = fits.ImageHDU(data=data_array)
    hdu2 = fits.ImageHDU(data=data_array + 1.0)
    hdu2.header["EXPTIME"] = 10.0
    hdu1.header["XTENSION"] = "SPECTRUM1"
    hdu2.header["XTENSION"] = "SPECTRUM2"
    a.append(hdu1)
    b.append(hdu2)
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "Extension HDU 1:" in report


def test_hdus_1d_arrays(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    hdu1 = fits.hdu.base.NonstandardExtHDU(name="SCI")
    hdu2 = fits.hdu.base.NonstandardExtHDU(name="SCI")
    hdu1.header["XTENSION"] = "Spectrum1"
    hdu2.header["XTENSION"] = "Spectrum1"
    hdu2.header["EXTLEVEL"] = 2
    hdu1.data = np.arange(100.0)
    hdu2.data = np.arange(100.0) + 1.0
    a.append(hdu1)
    b.append(hdu2)
    hdu3 = fits.hdu.base.NonstandardExtHDU(name="ERR")
    hdu4 = fits.hdu.base.NonstandardExtHDU(name="ERR")
    hdu3.header["XTENSION"] = "Spectrum2"
    hdu4.header["XTENSION"] = "Spectrum2"
    hdu4.header["EXTLEVEL"] = 3
    hdu3.data = np.arange(100.0)
    hdu4.data = np.arange(99.0)
    a.append(hdu3)
    b.append(hdu4)
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "Extension HDU 1 (SCI, 1)" in report
    # Earlier version didn't print pixel differences even when report_pixel_loc_diffs was set to True
    assert "Data differs at" in report


def test_hdus_image_zeros_rpl(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    data_array = np.zeros((6, 6), dtype=np.float32)
    a.append(fits.ImageHDU(data=data_array, name="ZERO"))
    b.append(fits.ImageHDU(data=data_array, name="ZERO"))
    data_a = data_array.copy()
    data_b = data_array.copy()
    data_a[2:] = 1.0
    data_b[4:] = 1.0
    a.append(fits.ImageHDU(data=data_a, name="SCI"))
    b.append(fits.ImageHDU(data=data_b, name="SCI"))
    data_a2 = data_array.copy()
    data_b2 = data_array.copy()
    data_a2[4:] = 1.0
    data_b2[2:] = 1.0
    a.append(fits.ImageHDU(data=data_a2, name="ERR"))
    b.append(fits.ImageHDU(data=data_b2, name="ERR"))
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "0.1     33.33         0" in report
    assert "0.1     33.33     33.33" in report


def test_hdu_tables_diff_col_numbers():
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    cf = fits.Column(name="FLUX", format="E", unit="erg /s /cm**2 /Angstrom")
    cd_a = fits.ColDefs([cw, cf])
    cd_b = fits.ColDefs([cw])
    table_a = fits.BinTableHDU.from_columns(cd_a, nrows=1)
    table_b = fits.BinTableHDU.from_columns(cd_b, nrows=1)
    a.append(table_a)
    b.append(table_b)
    diff = STFITSDiff(a, b)
    report = diff.report()
    assert "Tables have different number of columns:" in report


def test_hdu_tables_ignore_all_fields(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    cf = fits.Column(name="FLUX", format="E", unit="erg /s /cm**2 /Angstrom")
    cd_a = fits.ColDefs([cw, cf])
    cd_b = fits.ColDefs([cw])
    table_a = fits.BinTableHDU.from_columns(cd_a, nrows=100)
    table_b = fits.BinTableHDU.from_columns(cd_b, nrows=100)
    table_a.data["WAVELENGTH"] = np.arange(100.0)
    table_b.data["WAVELENGTH"] = np.arange(100.0) + 1.0
    a.append(table_a)
    b.append(table_b)
    fitsdiff_default_kwargs["ignore_fields"].extend("*")
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "Values in a and b" not in report


def test_hdus_tables_ignored_columns(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    cf = fits.Column(name="FLUX", format="E", unit="erg /s /cm**2 /Angstrom")
    cd_a = fits.ColDefs([cw, cf])
    cd_b = fits.ColDefs([cw, cf])
    table_a = fits.BinTableHDU.from_columns(cd_a, nrows=100)
    table_b = fits.BinTableHDU.from_columns(cd_b, nrows=100)
    table_a.data["WAVELENGTH"] = np.arange(100.0)
    table_b.data["WAVELENGTH"] = np.arange(100.0) + 1.0
    a.append(table_a)
    b.append(table_b)
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    fitsdiff_default_kwargs["ignore_fields"].extend(["WAVELENGTH"])
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "No differences found." in report


def test_hdus_tables_extra_cols_in_b(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    cf = fits.Column(name="FLUX", format="E", unit="erg /s /cm**2 /Angstrom")
    ce = fits.Column(name="ERROR", format="E")
    cd_a = fits.ColDefs([cw, cf])
    cd_b = fits.ColDefs([cw, cf, ce])
    table_a = fits.BinTableHDU.from_columns(cd_a, nrows=100)
    table_b = fits.BinTableHDU.from_columns(cd_b, nrows=100)
    table_a.data["WAVELENGTH"] = np.arange(100.0)
    table_b.data["WAVELENGTH"] = np.arange(100.0) + 1.0
    a.append(table_a)
    b.append(table_b)
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "Extra column ERROR of format E in b" in report


def test_hdus_tables_no_rows(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    cf = fits.Column(name="FLUX", format="E", unit="erg /s /cm**2 /Angstrom")
    ce = fits.Column(name="ERROR", format="E")
    cd_a = fits.ColDefs([cw, cf, ce])
    cd_b = fits.ColDefs([cw, cf, ce])
    table_a = fits.BinTableHDU.from_columns(cd_a)
    table_b = fits.BinTableHDU.from_columns(cd_b)
    a.append(table_a)
    b.append(table_b)
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "No differences found." in report


def test_hdus_tables_lowercase_column_names(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    clcw = fits.Column(name="wavelength", format="E", unit="m")
    cf = fits.Column(name="FLUX", format="E", unit="erg /s /cm**2 /Angstrom")
    ce = fits.Column(name="ERROR", format="E")
    cd_a = fits.ColDefs([clcw, cf, ce])
    cd_b = fits.ColDefs([cw, cf, ce])
    table_a = fits.BinTableHDU.from_columns(cd_a, nrows=100)
    table_b = fits.BinTableHDU.from_columns(cd_b, nrows=100)
    table_a.data["WAVELENGTH"] = np.arange(100.0)
    table_b.data["wavelength"] = np.arange(100.0) + 1.0
    a.append(table_a)
    b.append(table_b)
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "Found 100 different table data element(s)" in report


def test_hdus_tables_misc(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    cf = fits.Column(name="FLUX", format="E", unit="erg /s /cm**2 /Angstrom")
    ce = fits.Column(name="ERROR", format="E")
    ci = fits.Column(name="INDEX", format="J")
    cd_a = fits.ColDefs([cw, cf, ce, ci])
    cd_b = fits.ColDefs([cw, cf, ce, ci])
    table_a = fits.BinTableHDU.from_columns(cd_a, nrows=100)
    table_b = fits.BinTableHDU.from_columns(cd_b, nrows=100)
    table_a.data["WAVELENGTH"] = np.arange(100.0)
    table_b.data["WAVELENGTH"] = np.arange(100.0)
    table_b.data["WAVELENGTH"][:5] = np.arange(5) + 1.0
    table_b.data["WAVELENGTH"][40] = np.nan
    table_a.data["INDEX"] = np.arange(100)
    table_b.data["INDEX"] = np.arange(100)
    table_b.data["INDEX"][:5] = np.arange(5) + 1
    table_a.data["ERROR"] = np.arange(100.0)
    table_b.data["ERROR"] = np.arange(100.0)
    table_b.data["ERROR"][90:95] = np.nan
    table_a.data["FLUX"] = np.arange(100.0)
    table_b.data["FLUX"] = np.arange(100.0)
    table_b.data["FLUX"][10:20] *= 1.0000001
    a.append(table_a)
    b.append(table_b)
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    fitsdiff_default_kwargs["numdiffs"] = -1
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = report_to_list(diff.report())
    expected_report = [
        "Extension HDU 1:",
        "Data contains differences:",
        "Found 26 different table data element(s).",
        "16 failed the (atol, rtol) test",
        "Values in a and b",
        "col_name  zeros_a_b nan_a_b no-nan_a_b       max_a_b             min_a_b             mean_a_b",
        "---------- --------- ------- ---------- ------------------- ------------------- -------------------",
        "ERROR       1 1     0 5     100 95        99        99         0         0      49.5     47.26",
        "INDEX       1 0     0 0    100 100        99        99         0         1      49.5     49.55",
        "WAVELENGTH       1 0     0 1     100 99        99        99         0         1      49.5     49.65",
        "Difference stats for non-NaN diffs that fail the [atol, rtol] test: abs(b - a)",
        "col_name   dtype  rel_diffs rel_max rel_mean rel_std",
        "---------- ------- --------- ------- -------- -------",
        "ERROR float32         0     nan      nan     nan",
        "INDEX   int32         5       1        1       0",
        "WAVELENGTH float32         5       1        1       0",
        "* Pixel indices below are 1-based.",
        "Column ERROR data differs in row 90:",
        "a> 90.0",
        "b> nan",
        "Column ERROR data differs in row 91:",
        "a> 91.0",
        "b> nan",
        "Column ERROR data differs in row 92:",
        "a> 92.0",
        "b> nan",
        "Column ERROR data differs in row 93:",
        "a> 93.0",
        "b> nan",
        "Column ERROR data differs in row 94:",
        "a> 94.0",
        "b> nan",
        "Column INDEX data differs in row 0:",
        "a> 0",
        "b> 1",
        "Column INDEX data differs in row 1:",
        "a> 1",
        "b> 2",
        "Column INDEX data differs in row 2:",
        "a> 2",
        "b> 3",
        "Column INDEX data differs in row 3:",
        "a> 3",
        "b> 4",
        "Column INDEX data differs in row 4:",
        "a> 4",
        "b> 5",
        "Column WAVELENGTH data differs in row 0:",
        "a> 0.0",
        "b> 1.0",
        "Column WAVELENGTH data differs in row 1:",
        "a> 1.0",
        "b> 2.0",
        "Column WAVELENGTH data differs in row 2:",
        "a> 2.0",
        "b> 3.0",
        "Column WAVELENGTH data differs in row 3:",
        "a> 3.0",
        "b> 4.0",
        "Column WAVELENGTH data differs in row 4:",
        "a> 4.0",
        "b> 5.0",
        "Column WAVELENGTH data differs in row 40:",
        "a> 40.0",
        "b> nan",
        "...",
        "16 different table data element(s) found (4.00% different).",
    ]
    assert report == expected_report


def test_hdus_tables_non_numeric(fitsdiff_default_kwargs):
    a = fits.HDUList()
    b = fits.HDUList()
    a.append(fits.PrimaryHDU())
    b.append(fits.PrimaryHDU())
    cw = fits.Column(name="WAVELENGTH", format="E", unit="Angstrom")
    ci = fits.Column(name="INDEX", format="J")
    ct = fits.Column(name="NAME", format="10A")
    cidentical = fits.Column(name="IDENTICAL", format="10A")
    cd_a = fits.ColDefs([cw, ci, ct, cidentical])
    cd_b = fits.ColDefs([cw, ci, ct, cidentical])
    table_a = fits.BinTableHDU.from_columns(cd_a, nrows=10)
    table_b = fits.BinTableHDU.from_columns(cd_b, nrows=10)
    table_a.data["WAVELENGTH"] = np.arange(10.0)
    table_b.data["WAVELENGTH"] = np.arange(10.0)
    table_b.data["WAVELENGTH"][6:] = 4.0 - np.arange(4.0)
    table_a.data["INDEX"] = np.arange(10)
    table_b.data["INDEX"] = np.arange(10)
    table_b.data["INDEX"][6:] = 4 - np.arange(4)
    table_a.data["NAME"] = np.array(
        ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"]
    )
    table_b.data["NAME"] = np.array(
        ["Zero", "One", "Two", "Three", "Four", "Five", "Four", "Three", "Two", "One"]
    )
    table_a.data["IDENTICAL"] = np.array(
        ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"]
    )
    table_b.data["IDENTICAL"] = np.array(
        ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"]
    )
    a.append(table_a)
    b.append(table_b)
    diff = STFITSDiff(a, b)
    report = diff.report()
    assert "Column NAME has 4 different non-numeric entries" in report
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(a, b, **fitsdiff_default_kwargs)
    report = diff.report()
    assert "12 different table data element(s) found (30.00% different)." in report
    assert "Column ['IDENTICAL'] is identical" in report


def test_table_pq_different_array_sizes(mock_table, fitsdiff_default_kwargs):
    tmp_dir = mock_table[6]
    diff_arraysizes_truth = tmp_dir / "diff_arraysizes_truth.fits"
    diff_arraysizes = tmp_dir / "diff_arraysizes.fits"

    arr1 = [[0, 1], [2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4, 5]]
    arr2 = [[0.0, 1.0], [2.0, 3.0, np.nan, np.nan]]
    arr3 = [[0, 1], [2, 3, 4], [1, 2, 3, 4]]
    c1 = fits.Column(name="col_1", array=arr1, format="PI(2)")
    c2 = fits.Column(name="col_2", array=arr2, format="QD(2)")
    c3 = fits.Column(name="col_3", array=arr3, format="PI(2)")
    tab = fits.BinTableHDU.from_columns([c1, c2, c3], name="test")
    outfile = fits.HDUList()
    outfile.append(tab)
    outfile.writeto(diff_arraysizes_truth)

    arr1 = [[0, 11], [2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]]
    arr2 = [[0.0, 1.0], [np.nan, 23.0, 4.0, np.nan]]
    c1 = fits.Column(name="col_1", array=arr1, format="PI(2)")
    c2 = fits.Column(name="col_2", array=arr2, format="QD(2)")
    c3 = fits.Column(name="col_3", array=arr3, format="PI(2)")
    tab = fits.BinTableHDU.from_columns([c1, c2, c3], name="test")
    outfile = fits.HDUList()
    outfile.append(tab)
    outfile.writeto(diff_arraysizes)
    fitsdiff_default_kwargs["report_pixel_loc_diffs"] = True
    diff = STFITSDiff(diff_arraysizes, diff_arraysizes_truth, **fitsdiff_default_kwargs)
    result = diff.identical
    report = report_to_list(diff.report())
    assert result is False
    assert "Extra column col_1 of format PI(4) in a" in report
    assert "Extra column col_1 of format PI(5) in b" in report
    assert "Column ['col_3'] is identical" in report
