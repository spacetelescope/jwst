import numpy as np

from jwst.adaptive_trace_model.bspline import bspline_fit
from jwst.adaptive_trace_model.tests.helpers import profile_1d


def test_bspline_fit_flat_data():
    # Make some flat test data
    xvec = np.arange(100, dtype=float)
    yvec = np.full(100, 1.0)

    # Fit the data with as many breakpoints as possible, one iteration
    bspline = bspline_fit(xvec, yvec)

    # The fit spline should exactly reproduce the data
    np.testing.assert_allclose(bspline(xvec), yvec)


def test_bspline_fit_spikes():
    # Make some test data with spikes
    xvec = np.arange(100, dtype=float)
    yvec = np.zeros(100, dtype=float)
    yvec[5::10] = 1

    # Fit the data with as many breakpoints as possible, one iteration
    bspline = bspline_fit(xvec, yvec, nbkpts=90, wrapiter=1)

    # The fit spline with regular breakpoints
    # should reproduce the spikes pretty well, but not perfectly.
    np.testing.assert_allclose(bspline(xvec), yvec, atol=0.3)


def test_bspline_fit_smooth_profile():
    # Make some test data with a smooth profile
    xvec = np.arange(100, dtype=float)
    yvec = profile_1d(xvec)

    # Fit the data with default breakpoints, 1 iteration
    bspline = bspline_fit(xvec, yvec, wrapiter=1)

    # The fit spline with regular breakpoints
    # should reproduce the smooth profile very well
    np.testing.assert_allclose(bspline(xvec), yvec, rtol=1e-3)


def test_bspline_fit_spacing_gap():
    # Make some test data with one data point far from the others
    xvec = np.arange(100, dtype=float)
    xvec[-1] = 1000
    yvec = profile_1d(xvec)

    # The extra point should be rejected and fit should be good,
    # except for the last point
    bspline = bspline_fit(xvec, yvec, wrapiter=1)
    np.testing.assert_allclose(bspline(xvec[:-1]), yvec[:-1], rtol=1e-3)


def test_bspline_fit_too_many_gaps():
    xvec = np.arange(0, 1000, 50, dtype=float)
    yvec = profile_1d(xvec)

    # too many gaps compared to breakpoint spacing: spline fit returns None
    bspline = bspline_fit(xvec, yvec, wrapiter=1)
    assert bspline is None


def test_bspline_fit_too_few_points():
    # Vector is too small for cubic fit
    xvec = np.arange(10, dtype=float)
    yvec = np.zeros(10, dtype=float)
    bspline = bspline_fit(xvec, yvec, wrapiter=1)
    assert bspline is None


def test_bspline_fit_failed():
    # Vector contains bad data: fit will fail
    xvec = np.arange(100, dtype=float)
    yvec = np.zeros(100, dtype=float)
    yvec[10:50] = np.nan

    bspline = bspline_fit(xvec, yvec, wrapiter=1)
    assert bspline is None


def test_bspline_fit_outlier_rejection(caplog):
    caplog.set_level("DEBUG", "jwst")

    # Make test data with an outlier to reject
    xvec = np.arange(0, 100, dtype=float)
    yvec = profile_1d(xvec)
    yvec[10] = 10

    # Fit with two iterations
    bspline = bspline_fit(xvec, yvec, wrapiter=2, verbose=True)

    # The outlier (plus surrounding points) is rejected in iteration 1;
    # any further iterations would reject the profile peak
    assert "Iter 0 Rejected 3 Kept 97" in caplog.text

    # Fit is good except for the outlier pixel
    np.testing.assert_allclose(bspline(xvec[2:9]), yvec[2:9], rtol=1e-3)
    np.testing.assert_allclose(bspline(xvec[11:]), yvec[11:], rtol=1e-3)


def test_bspline_fit_reject_too_many(caplog):
    caplog.set_level("DEBUG", "jwst")

    # Make noisy test data with spikes to successively reject
    npt = 1000
    xvec = np.arange(0, npt, dtype=float)
    yvec = np.full(npt, 1.0)
    yvec[20::14] = 50
    yvec[3::25] = 40
    yvec[7::27] = 30
    yvec[10::29] = 20
    yvec[13::29] = 10

    # Fit with 10 iterations
    bspline = bspline_fit(xvec, yvec, wrapiter=10, verbose=True)

    # Stops before reaching the end of the iterations because
    # too many data points are rejected
    assert "Iter 9" not in caplog.text
    assert "too many points rejected" in caplog.text
    assert np.mean(bspline(xvec)) > 1.0
