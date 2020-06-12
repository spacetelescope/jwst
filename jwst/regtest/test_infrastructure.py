from glob import glob
import os

import pytest
from astropy.table import Table
import astropy.units as u
import numpy as np
from .regtestdata import text_diff


@pytest.mark.bigdata
def test_regtestdata_get_data(rtdata):
    rtdata.get_data("infrastructure/test_regtestdata/file1_rate.fits")
    rtdata.output = "file1_cal.fits"

    assert rtdata.input == os.path.join(os.getcwd(), "file1_rate.fits")


@pytest.mark.bigdata
def test_regtestdata_get_truth(rtdata):
    rtdata.get_truth("infrastructure/test_regtestdata/file1_rate.fits")
    rtdata.output = "file1_rate.fits"

    assert rtdata.truth == os.path.join(os.getcwd(), "truth", "file1_rate.fits")


@pytest.mark.bigdata
def test_regtestdata_get_asn(rtdata):
    rtdata.get_asn("infrastructure/test_regtestdata/my_asn.json")
    files = glob("*.fits")
    rtdata.output = "file1_rate.fits"

    assert os.path.isfile("my_asn.json")
    assert len(files) == 3


def test_fitsdiff_defaults(fitsdiff_default_kwargs):
    assert 'ASDF' in fitsdiff_default_kwargs['ignore_hdus']


@pytest.fixture
def two_tables(tmpdir):
    """Return identical astropy tables written to 2 .ecsv file paths"""
    path1 = str(tmpdir.join("catalog1.ecsv"))
    path2 = str(tmpdir.join("catalog2.ecsv"))
    a = np.array([1, 4, 5], dtype=np.float)
    b = [2.0, 5.0, 8.5]
    c = ['x', 'y', 'z']
    d = [10, 20, 30] * u.m / u.s
    names = ('a', 'b', 'c', 'd')
    meta = {'name': 'first table'}
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
    t1.add_row([7, 5.0, 'q', 40 * u.m / u.s])
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="Row count"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_columns(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    del t1['d']
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="Column names"):
        assert diff_astropy_tables(path1, path2)


@pytest.mark.xfail(reason='table meta comparison currently deactivated')
def test_diff_astropy_tables_meta(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1.meta = {"name": "not the first table"}
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="Metadata does not match"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_allclose(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1['a'] += 2e-5
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="Not equal to tolerance"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_dtype(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1['a'] = np.array([1, 4, 6], dtype=np.int)
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="dtype does not match"):
        assert diff_astropy_tables(path1, path2)


def test_diff_astropy_tables_all_equal(diff_astropy_tables, two_tables):
    path1, path2 = two_tables

    t1 = Table.read(path1)
    t1['a'] = np.array([1, 4, 6], dtype=np.int)
    t1.write(path1, overwrite=True, format="ascii.ecsv")

    t2 = Table.read(path2)
    t2['a'] = np.array([1, 4, 7], dtype=np.int)
    t2.write(path2, overwrite=True, format="ascii.ecsv")

    with pytest.raises(AssertionError, match="values do not match"):
        assert diff_astropy_tables(path1, path2)


def test_text_diff(tmpdir):
    path1 = str(tmpdir.join("test1.txt"))
    path2 = str(tmpdir.join("test2.txt"))
    with open(path1, "w") as text_file:
        print("foo", file=text_file)
    with open(path2, "w") as text_file:
        print("foo", file=text_file)

    assert text_diff(path1, path2)

    with open(path2, "w") as text_file:
        print("bar", file=text_file)

    with pytest.raises(AssertionError):
        assert text_diff(path1, path2)
