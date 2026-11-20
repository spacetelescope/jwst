import os
from glob import glob

import astropy.units as u
import numpy as np
import pytest
from astropy.table import Table

from jwst.regtest.regtestdata import RTFile, text_diff

INPUT_DATA = {
    "mock_asn": RTFile(
        file_name="my_asn.json",
        path="rtdata/infrastructure/test_regtestdata",
        from_mast=False,
        asn_files=["file1_rate.fits", "file2_rate.fits", "file3_rate.fits"],
        asn_files_from_mast=False,
        callback_function="mk_empty_files",
        comment="This code was used to create the files in the association.",
    )
}


def mk_empty_files():
    """Create empty FITS files for the test."""
    for i in range(3):
        file = "file" + str(i + 1) + "_rate.fits"
        with open(file, "w") as f:
            f.write("")


@pytest.mark.bigdata
def test_regtestdata_get_data(rtdata, tmp_cwd):
    rate_file_path = INPUT_DATA["mock_asn"].path + "/" + INPUT_DATA["mock_asn"].asn_files[0]
    rtdata.get_data(rate_file_path)
    rate_file = INPUT_DATA["mock_asn"].asn_files[0]
    rtdata.output = rate_file

    assert rtdata.input == str(tmp_cwd / rate_file)


@pytest.mark.bigdata
def test_regtestdata_get_truth(rtdata, tmp_cwd):
    rate_file_path = INPUT_DATA["mock_asn"].path + "/" + INPUT_DATA["mock_asn"].asn_files[0]
    rtdata.get_truth(rate_file_path)
    rate_file = INPUT_DATA["mock_asn"].asn_files[0]
    rtdata.output = rate_file

    assert rtdata.truth == str(tmp_cwd / "truth" / rate_file)


@pytest.mark.bigdata
def test_regtestdata_get_asn(rtdata):
    rtdata.get_asn(INPUT_DATA["mock_asn"])
    files = glob("*.fits")
    rtdata.output = INPUT_DATA["mock_asn"].asn_files[0]

    assert os.path.isfile(INPUT_DATA["mock_asn"].file_name)
    assert len(files) == 3


def test_fitsdiff_defaults(fitsdiff_default_kwargs):
    assert "ASDF" in fitsdiff_default_kwargs["ignore_hdus"]


@pytest.fixture
def two_tables(tmp_path):
    """Return identical astropy tables written to 2 .ecsv file paths"""
    path1 = tmp_path / "catalog1.ecsv"
    path2 = tmp_path / "catalog2.ecsv"
    a = np.array([1, 4, 5], dtype=float)
    b = [2.0, 5.0, 8.5]
    c = ["x", "y", "z"]
    d = [10, 20, 30] * u.m / u.s
    names = ("a", "b", "c", "d")
    meta = {"name": "first table"}
    t = Table([a, b, c, d], names=names, meta=meta)
    t.write(path1, format="ascii.ecsv")
    t.write(path2, format="ascii.ecsv")

    return path1, path2


def test_diff_astropy_tables_same(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_length(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1.add_row([7, 5.0, "q", 40 * u.m / u.s])
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="Row count"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_columns(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    del t1["d"]
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="Column names"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_allclose(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1["a"] += 2e-5
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="Not equal to tolerance"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_dtype(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1["a"] = np.array([1, 4, 6], dtype=int)
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="dtype does not match"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_all_equal(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1["a"] = np.array([1, 4, 6], dtype=int)
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    t2 = Table.read(path2)
    t2["a"] = np.array([1, 4, 7], dtype=int)
    t2.write(path2, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="values do not match"):
        assert diff_astropy_tables(path1, path2)


def test_text_diff(tmp_path):
    path1 = tmp_path / "test1.txt"
    path2 = tmp_path / "test2.txt"
    with open(path1, "w") as text_file:
        print("foo", file=text_file)
    with open(path2, "w") as text_file:
        print("foo", file=text_file)

    assert text_diff(path1, path2)

    with open(path2, "w") as text_file:
        print("bar", file=text_file)

    with pytest.raises(AssertionError):
        assert text_diff(path1, path2)
