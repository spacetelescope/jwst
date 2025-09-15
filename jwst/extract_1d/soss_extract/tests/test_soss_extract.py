from functools import partial

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import SossWaveGridModel, SpecModel

from jwst.extract_1d.soss_extract.soss_extract import (
    SHORT_CUTOFF,
    _build_null_spec_table,
    _compute_box_weights,
    _model_image,
)
from jwst.extract_1d.soss_extract.tests.helpers import DATA_SHAPE


@pytest.mark.parametrize("order", [1, 2, 3])
def test_build_null_spec_table(wave_grid, order):
    expected_cut = SHORT_CUTOFF[order - 1]
    spec = _build_null_spec_table(wave_grid, order)

    assert isinstance(spec, SpecModel)
    assert spec.spectral_order == order
    assert spec.spec_table.shape == (93,)
    assert np.all(np.isnan(spec.spec_table["FLUX"]))
    assert np.all(spec.spec_table["DQ"] == 1)


@pytest.fixture
def monkeypatch_setup(
    monkeypatch,
    wave_map,
    trace_profile,
    throughput,
    webb_kernels,
    trace1d,
):
    """Monkeypatch get_ref_file_args and get_trace_1d to return the miniature model detector"""

    def mock_get_ref_file_args(wave, trace, thru, kern, reffiles, orders_requested):
        """Return the arrays from conftest instead of querying CRDS"""
        if (orders_requested is None) or (orders_requested == [1, 2, 3]):
            return [wave, trace, thru, kern]
        elif orders_requested == [1, 2]:
            return [wave[:2], trace[:2], thru[:2], kern[:2]]
        else:
            raise ValueError(f"Unexpected orders_requested for mock: {orders_requested}")

    def mock_trace1d(trace, reffiles, order):
        """Return the traces from conftest instead of doing math that requires a full-sized detector"""
        return trace[int(order) - 1]

    monkeypatch.setattr(
        "jwst.extract_1d.soss_extract.soss_extract.get_ref_file_args",
        partial(mock_get_ref_file_args, wave_map, trace_profile, throughput, webb_kernels),
    )
    monkeypatch.setattr(
        "jwst.extract_1d.soss_extract.soss_extract._get_trace_1d", partial(mock_trace1d, trace1d)
    )


@pytest.mark.parametrize("order_list", [[1, 2], None])
def test_model_image(monkeypatch_setup, imagemodel, detector_mask, ref_files, order_list):
    scidata, scierr = imagemodel

    if order_list is None:
        orders_expected = [1, 2, 3]
    else:
        orders_expected = order_list
    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, _wavelengths = _compute_box_weights(
        ref_files, DATA_SHAPE, box_width, orders_requested=orders_expected
    )

    tracemodels, tikfac, logl, wave_grid, spec_list = _model_image(
        scidata,
        scierr,
        detector_mask,
        refmask,
        ref_files,
        box_weights,
        tikfac=None,
        threshold=1e-4,
        n_os=2,
        wave_grid=None,
        estimate=None,
        rtol=1e-3,
        max_grid_size=1000000,
        order_list=order_list,
    )

    # check output basics, types and shapes
    assert len(tracemodels) == len(orders_expected)
    for order in orders_expected:
        tm = tracemodels[f"Order {order}"]
        assert tm.dtype == np.float64
        assert tm.shape == DATA_SHAPE
        # should be some nans in the trace model but not all
        assert 0 < np.sum(np.isfinite(tm)) < tm.size
    for x in [tikfac, logl]:
        assert isinstance(x, float)
        assert np.isfinite(x)
    assert logl < 0
    assert wave_grid.dtype == np.float64
    for spec in spec_list:
        assert isinstance(spec, SpecModel)

    factors = np.array([getattr(spec.meta.soss_extract1d, "factor", np.nan) for spec in spec_list])
    chi2s = np.array([getattr(spec.meta.soss_extract1d, "chi2", np.nan) for spec in spec_list])
    orders = np.array([spec.spectral_order for spec in spec_list])
    colors = np.array([spec.meta.soss_extract1d.color_range for spec in spec_list])

    assert tikfac in factors

    # ensure outputs have the shapes we expect for each order and blue/red
    n_good = []
    for order in [1, 2]:
        for color in ["RED", "BLUE"]:
            good = (order == orders) & (color == colors)
            # check that there's at least one good spectrum for all valid order-color combinations
            if not np.any(good):
                assert order == 1
                assert color == "BLUE"
                continue
            n_good.append(np.sum(good))

            this_factors = factors[good]
            this_chi2s = chi2s[good]
            this_spec = np.array(spec_list)[good]
            nochi = np.isnan(this_chi2s)

            # _model_single_order is set up so that the final/best spectrum is last in the list
            # it lacks chi2 calculations
            assert np.sum(nochi) == 1
            assert np.where(nochi)[0][0] == len(this_chi2s) - 1

            # it represents the best tikhonov factor for that order-color combination
            # which is not necessarily the same as the top-level tikfac for the blue part of order 2
            # but it is the same for the red part of order 1 and the red part of order 2
            if color == "RED":
                assert this_factors[-1] == tikfac

            # check that the output spectra contain good data
            for spec in this_spec:
                spec = np.array([[s[0], s[1]] for s in spec.spec_table])
                assert np.sum(np.isfinite(spec)) == spec.size

    # check that all order-color combinations have the same number of spectra
    n_good = np.array(n_good)
    assert np.all(n_good >= 1)
    assert np.all(n_good - n_good[0] == 0)


def test_model_image_tikfac_specified(
    monkeypatch_setup,
    imagemodel,
    detector_mask,
    ref_files,
):
    """Ensure spec_list is a single-element list per order if tikfac is specified"""
    scidata, scierr = imagemodel

    order_list = [1, 2]
    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, wavelengths = _compute_box_weights(
        ref_files, DATA_SHAPE, box_width, orders_requested=order_list
    )

    tikfac_in = 1e-7
    tracemodels, tikfac, logl, wave_grid, spec_list = _model_image(
        scidata,
        scierr,
        detector_mask,
        refmask,
        ref_files,
        box_weights,
        tikfac=tikfac_in,
        threshold=1e-4,
        n_os=2,
        wave_grid=None,
        estimate=None,
        rtol=1e-3,
        max_grid_size=1000000,
        order_list=order_list,
    )
    # check that spec_list is a single-element list per order in this case
    assert len(spec_list) == 3
    assert tikfac == tikfac_in


def test_model_image_wavegrid_specified(
    monkeypatch_setup,
    imagemodel,
    detector_mask,
    ref_files,
):
    """Ensure wave_grid is used if specified.
    Also specify tikfac because it makes the code run faster to not have to re-derive it.

    Note the failure with SossWaveGridModel passed as input. What should be done about that?
    """
    scidata, scierr = imagemodel

    order_list = [1, 2]
    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, wavelengths = _compute_box_weights(
        ref_files, DATA_SHAPE, box_width, orders_requested=order_list
    )

    tikfac_in = 1e-7
    # test np.array input
    wave_grid_in = np.linspace(1.0, 2.5, 100)
    tracemodels, tikfac, logl, wave_grid, spec_list = _model_image(
        scidata,
        scierr,
        detector_mask,
        refmask,
        ref_files,
        box_weights,
        tikfac=tikfac_in,
        threshold=1e-4,
        n_os=2,
        wave_grid=wave_grid_in,
        estimate=None,
        rtol=1e-3,
        max_grid_size=1000000,
        order_list=order_list,
    )
    assert np.allclose(wave_grid, wave_grid_in)

    # test SossWaveGridModel input
    # the docs on main say this works, but I don't think it does even on main
    with pytest.raises(ValueError):
        wave_grid_in = SossWaveGridModel()
        wave_grid_in.wavegrid = np.linspace(1.0, 2.5, 100)
        tracemodels, tikfac, logl, wave_grid, spec_list = _model_image(
            scidata,
            scierr,
            detector_mask,
            refmask,
            ref_files,
            box_weights,
            tikfac=tikfac_in,
            threshold=1e-4,
            n_os=2,
            wave_grid=wave_grid_in,
            estimate=None,
            rtol=1e-3,
            max_grid_size=1000000,
            order_list=order_list,
        )
