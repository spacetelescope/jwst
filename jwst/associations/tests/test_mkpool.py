"""
Test mkpool
"""

import os
from glob import glob

import pytest
from astropy.io import fits
from astropy.utils.data import get_pkg_data_path

from jwst.associations import AssociationPool
from jwst.associations.mkpool import from_cmdline, mkpool, NON_HEADER_COLS

# Optional column settings
OPT_COLS = [
    ("asn_candidate", [("a3001", "coron")]),
    ("dms_note", "a note from dms"),
    ("is_imprt", "t"),
    ("pntgtype", "target_acquisition"),
]

# Required column names
REQUIRED_PARAMS = set(("program", "filename"))
NON_HEADER_PARAMS = set(NON_HEADER_COLS.keys())


def test_hdu(exposures):
    hdus = []
    for exposure in exposures:
        with fits.open(exposure) as hdulist:
            hdus.append(hdulist[0])
    pool = mkpool(hdus)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert NON_HEADER_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)


def test_hdulist(exposures):
    hduls = [fits.open(exposure) for exposure in exposures]
    pool = mkpool(hduls)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert NON_HEADER_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)
    [h.close() for h in hduls]


def test_mkpool(exposures):
    pool = mkpool(exposures)
    assert isinstance(pool, AssociationPool)
    assert REQUIRED_PARAMS.issubset(pool.colnames)
    assert NON_HEADER_PARAMS.issubset(pool.colnames)
    assert len(pool) == len(exposures)
    filenames = [filename for filename in pool["filename"]]
    assert set(exposures) == set(filenames)


@pytest.mark.parametrize("opt_cols", OPT_COLS)
def test_opt_cols(mkpool_with_args, opt_cols):
    """Ensure that optional arguments are properly used"""
    _test_opt_cols(mkpool_with_args, opt_cols)


@pytest.mark.parametrize("opt_cols", OPT_COLS)
def test_opt_cols_cmdline(mkpool_cmdline, opt_cols):
    """Ensure that the command line with optional arguments are properly used"""
    _test_opt_cols(mkpool_cmdline, opt_cols)


def test_nonheader_cols(mkpool_with_args):
    """Ensure that non-header defaults are properly set."""
    _test_nonheader_cols(mkpool_with_args)


def test_nonheader_cols_cmdline(mkpool_cmdline):
    """Ensure that non-header defaults are properly set with the command line."""
    _test_nonheader_cols(mkpool_cmdline)


# ####################
# Fixtures & Utilities
# ####################
@pytest.fixture(
    scope="module",
    params=[
        "data/exposures_nopsf",
        "data/exposures",
    ],
)
def exposures(request):
    exposure_path = get_pkg_data_path(request.param, package="jwst.associations.tests")
    return sorted(glob(os.path.join(exposure_path, "*.fits")))


@pytest.fixture(scope="module")
def mkpool_cmdline(exposures):
    """Create a pool with optional arguments from the commandline"""
    args = ["pool.csv"]
    for column, value in OPT_COLS:
        args.append(f"--{column.replace('_', '-')}")
        args.append(f"{value}")
    for exposure in exposures:
        args.append(exposure)

    mkpool_args = from_cmdline(args)
    pool = mkpool(**mkpool_args)
    return pool


@pytest.fixture(scope="module")
def mkpool_with_args(exposures):
    """Create a pool with all optional arguments specified"""
    kargs = {column: value for column, value in OPT_COLS}
    pool = mkpool(exposures, **kargs)

    return pool


def _test_opt_cols(mkpool_with_args, opt_cols):
    """Ensure that optional arguments are properly used"""
    column, expected = opt_cols
    if column == "asn_candidate":
        expected = "[" + str(("o001", "observation")) + ", " + str(expected)[1:]
    assert mkpool_with_args[0][column] == expected


def _test_nonheader_cols(mkpool_with_args):
    """Ensure that non-header defaults are properly set"""
    opt_dict = dict(OPT_COLS)
    for column, expected in NON_HEADER_COLS.items():
        if column in opt_dict:
            expected = opt_dict[column]
        if column == "asn_candidate":
            expected = "[" + str(("o001", "observation")) + ", " + str(expected)[1:]
        assert mkpool_with_args[0][column] == expected
