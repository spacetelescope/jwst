import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.photom import time_dependence as td

# Some useful constants to reuse
NROW = 10
MIDTIME = 59700.0
TIMEDIFF = np.arange(1, NROW + 1)
T0 = MIDTIME - TIMEDIFF


def create_photom():
    """
    Create a photom model with timecoeff extensions.

    Returns
    -------
    ftab : `~jwst.datamodels.MirImgPhotomModel`
        An open data model for a MIRI image photom reference file.
    """
    # Filler values for filter/subarray - they would actually vary by row.
    filter_list = ["F2550W"] * NROW
    subarray = ["SUB256"] * NROW
    nrows = len(filter_list)

    photmjsr = np.linspace(3.1, 3.1 + (nrows - 1.0) * 0.1, nrows)
    uncertainty = np.zeros(nrows, dtype=np.float32)
    dtype = np.dtype(
        [("filter", "S12"), ("subarray", "S15"), ("photmjsr", "<f4"), ("uncertainty", "<f4")]
    )
    reftab = np.array(
        list(zip(filter_list, subarray, photmjsr, uncertainty, strict=True)), dtype=dtype
    )

    # linear time coefficients
    lossperyear = np.full(nrows, 0.01)
    dtypec = np.dtype(
        [("filter", "S12"), ("subarray", "S15"), ("t0", "<f4"), ("lossperyear", "<f4")]
    )
    reftab_lin = np.array(
        list(zip(filter_list, subarray, T0, lossperyear, strict=True)), dtype=dtypec
    )

    # exponential time coefficients
    amp = np.linspace(2.1, 2.1 + (nrows - 1.0) * 0.1, nrows) / photmjsr
    tau = np.full(nrows, 145)
    const = np.full(nrows, 1.0)
    dtypec = np.dtype(
        [
            ("filter", "S12"),
            ("subarray", "S15"),
            ("t0", "<f4"),
            ("amplitude", "<f4"),
            ("tau", "<f4"),
            ("const", "<f4"),
        ]
    )
    reftab_exp = np.array(
        list(zip(filter_list, subarray, T0, amp, tau, const, strict=True)), dtype=dtypec
    )

    # power law time coefficients
    year1value = np.full(nrows, 0.01)
    tsoft = np.full(nrows, 1.0)
    alpha = np.full(nrows, -1.01)
    dtypec = np.dtype(
        [
            ("filter", "S12"),
            ("subarray", "S15"),
            ("t0", "<f4"),
            ("year1value", "<f4"),
            ("tsoft", "<f4"),
            ("alpha", "<f4"),
        ]
    )
    reftab_plaw = np.array(
        list(zip(filter_list, subarray, T0, year1value, tsoft, alpha, strict=True)), dtype=dtypec
    )

    ftab = datamodels.MirImgPhotomModel(
        phot_table=reftab,
        timecoeff_linear=reftab_lin,
        timecoeff_exponential=reftab_exp,
        timecoeff_powerlaw=reftab_plaw,
    )
    return ftab


def test_linear_correction():
    lossperyear = np.full(NROW, 1.0)
    correction = td.linear_correction(MIDTIME, T0, lossperyear)
    assert correction.shape == (NROW,)

    # Linear dropoff with yearly slope
    expected = 1 - TIMEDIFF / 365
    np.testing.assert_allclose(correction, expected)


@pytest.mark.parametrize("bounded", [True, False])
def test_linear_correction_bounded(bounded):
    lossperyear = np.full(NROW, -1.0)
    correction = td.linear_correction(MIDTIME, T0, lossperyear, bounded=bounded)
    assert correction.shape == (NROW,)

    if not bounded:
        # unrealistic linear gain to test bounding
        expected = 1 + TIMEDIFF / 365
        np.testing.assert_allclose(correction, expected)
    else:
        # bounding sets all >1 values to 1.0
        np.testing.assert_equal(correction, 1.0)


def test_exponential_correction():
    amplitude = np.full(NROW, 1.0)
    tau = np.full(NROW, 10.0)
    const = np.full(NROW, 0.0)
    correction = td.exponential_correction(MIDTIME, T0, amplitude, tau, const)
    assert correction.shape == (NROW,)

    # exponential dropoff
    expected = np.exp(-TIMEDIFF / 10)
    np.testing.assert_allclose(correction, expected)


@pytest.mark.parametrize("bounded", [True, False])
def test_exponential_correction_bounded(bounded):
    amplitude = np.full(NROW, 1.0)
    tau = np.full(NROW, -10.0)
    const = np.full(NROW, 0.0)
    correction = td.exponential_correction(MIDTIME, T0, amplitude, tau, const, bounded=bounded)
    assert correction.shape == (NROW,)

    if not bounded:
        # unrealistic exponential increase to test bounding
        expected = np.exp(TIMEDIFF / 10)
        np.testing.assert_allclose(correction, expected)
    else:
        # bounding sets all >1 values to 1.0
        np.testing.assert_equal(correction, 1.0)


def test_powerlaw_correction():
    tsoft = np.full(NROW, 1.0)
    alpha = np.full(NROW, -1.1)
    year1value = np.full(NROW, (365 + tsoft) ** alpha)
    correction = td.powerlaw_correction(MIDTIME, T0, year1value, tsoft, alpha)
    assert correction.shape == (NROW,)

    # powerlaw dropoff
    expected = (TIMEDIFF + tsoft) ** alpha
    np.testing.assert_allclose(correction, expected)


@pytest.mark.parametrize("bounded", [True, False])
def test_powerlaw_correction_bounded(bounded):
    tsoft = np.full(NROW, 1.0)
    alpha = np.full(NROW, 1.1)
    year1value = np.full(NROW, (365 + tsoft) ** alpha)
    correction = td.powerlaw_correction(MIDTIME, T0, year1value, tsoft, alpha, bounded=bounded)
    assert correction.shape == (NROW,)

    if not bounded:
        # unrealistic powerlaw increase to test bounding
        expected = (TIMEDIFF + tsoft) ** alpha
        np.testing.assert_allclose(correction, expected)
    else:
        # bounding sets all >1 values to 1.0
        np.testing.assert_equal(correction, 1.0)


@pytest.mark.parametrize("bounded", [True, False])
def test_get_correction_table(bounded):
    """Check specific expected output values given synthetic input."""
    ftab = create_photom()
    correction = td.get_correction_table(ftab, MIDTIME, bounded=bounded)
    if bounded:
        # Bounding happens in each individual correction, and at least one of
        # them is unreasonable, so composed correction is different for all values.
        expected = [
            0.999973,
            0.999945,
            0.957193,
            0.764027,
            0.635512,
            0.54387,
            0.475239,
            0.421925,
            0.379323,
            0.344501,
        ]
    else:
        # Correction is allowed to be >1 for each
        expected = [
            3.224766,
            2.147894,
            1.610667,
            1.288666,
            1.074063,
            0.920744,
            0.80569,
            0.71613,
            0.644408,
            0.585658,
        ]
    np.testing.assert_allclose(correction, expected, atol=1e-6)


def test_get_correction_table_composition():
    """Check that composed correction is the individual corrections multiplied."""
    ftab = create_photom()
    composed_correction = td.get_correction_table(ftab, MIDTIME)

    lin_cor = td.linear_correction(
        MIDTIME, ftab.timecoeff_linear["t0"], ftab.timecoeff_linear["lossperyear"]
    )
    exp_cor = td.exponential_correction(
        MIDTIME,
        ftab.timecoeff_exponential["t0"],
        ftab.timecoeff_exponential["amplitude"],
        ftab.timecoeff_exponential["tau"],
        ftab.timecoeff_exponential["const"],
    )
    plaw_cor = td.powerlaw_correction(
        MIDTIME,
        ftab.timecoeff_powerlaw["t0"],
        ftab.timecoeff_powerlaw["year1value"],
        ftab.timecoeff_powerlaw["tsoft"],
        ftab.timecoeff_powerlaw["alpha"],
    )
    np.testing.assert_allclose(composed_correction, lin_cor * exp_cor * plaw_cor)
