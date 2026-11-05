"""Test the RTFile class."""

import pytest

from jwst.regtest.regtestdata import RTFile


def test_good_fits_input_rtfile():
    rtobj = RTFile(
        file_name="jw_mock1_mod_uncal.fits",
        path="rtdata/fgs/image3",
        from_mast=True,
        callback_function="N/A",
    )
    rtobj.validate_file_name()
    assert rtobj.root_name == "jw_mock1_mod"
    assert rtobj.full_path == "rtdata/fgs/image3/jw_mock1_mod_uncal.fits"


def test_good_asn_input_rtfile():
    rtobj = RTFile(
        file_name="jw_mock_mod_asn.json",
        path="rtdata/fgs/image3",
        from_mast=False,
        callback_function="fake_function",
        asn_files=["jw_mock_file2.fits", "jw_mock_file3.fits"],
        asn_files_from_mast=True,
    )
    rtobj.validate_file_name()
    assert rtobj.full_path == "rtdata/fgs/image3/jw_mock_mod_asn.json"


def test_bad_fits_input_rtfile():
    # Test bad place of 'mod' in file name
    rtobj = RTFile(
        file_name="jw_mod_mock1_uncal.fits",
        path="rtdata/fgs/image3",
        from_mast=True,
        callback_function="fake_function",
    )
    pytest.raises(ValueError, lambda: rtobj.validate_file_name())


def test_bad_asn_input_rtfile():
    # Test no 'mod' in file name
    rtobj = RTFile(
        file_name="jw_mock_asn.json",
        path="rtdata/fgs/image3",
        from_mast=True,
        callback_function="fake_function",
    )
    pytest.raises(ValueError, lambda: rtobj.validate_file_name())
    pytest.raises(ValueError, lambda: rtobj.root_name)
