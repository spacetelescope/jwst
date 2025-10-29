import numpy as np
import pytest
from stdatamodels.jwst.datamodels import SossWaveGridModel, SpecModel

from jwst.extract_1d.soss_extract.soss_extract import (
    SHORT_CUTOFF,
    DetectorModelOrder,
    Integration,
    _build_null_spec_table,
    _compute_box_weights,
    _process_one_integration,
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
):
    """Monkeypatch get_trace_1d to return the miniature model detector"""

    def mock_make_background_mask(*args, **kwargs):
        return np.zeros(DATA_SHAPE, dtype="bool")

    monkeypatch.setattr(
        "jwst.extract_1d.soss_extract.soss_extract.make_background_mask", mock_make_background_mask
    )
    monkeypatch.setattr("jwst.extract_1d.soss_extract.soss_extract.CUTOFFS", [200, 200, 200])


@pytest.fixture
def detector_models(
    wave_map,
    trace_profile,
    throughput,
    webb_kernels,
    trace1d,
):
    """Return the reference file arguments as DetectorModels"""
    detector_models = []
    for order in [1, 2, 3]:
        detector_model = DetectorModelOrder(
            spectral_order=order,
            wavemap=wave_map[order - 1],
            spectrace=trace1d[order - 1],
            specprofile=trace_profile[order - 1],
            throughput=throughput[order - 1],
            kernel=webb_kernels[order - 1],
            kernel_func=webb_kernels[order - 1],
            kernel_native=webb_kernels[order - 1],
            subarray="SUBSTRIP256",
        )
        detector_models.append(detector_model)
    return detector_models


@pytest.mark.parametrize("order_list", [[1, 2], [1, 2, 3]])
def test_model_image(monkeypatch_setup, imagemodel, detector_mask, detector_models, order_list):
    scidata, scierr = imagemodel

    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, _wavelengths = _compute_box_weights(
        detector_models, DATA_SHAPE, box_width, orders_requested=order_list
    )

    if len(order_list) == 3:
        extract_order3 = True
    else:
        extract_order3 = False

    integration = Integration(
        scidata,
        scierr,
        detector_mask,
        refmask,
        detector_models,
        box_weights,
        do_bkgsub=True,
        extract_order3=extract_order3,
    )

    estimate = integration.estim_flux_first_order()
    wave_grid = integration.make_decontamination_grid(
        estimate,
        rtol=1e-3,
        max_grid_size=1000000,
        n_os=2,
    )
    tracemodels, spec_list, tikfacs_out = integration.model_image(
        wave_grid,
        tikfacs_in=None,
        threshold=1e-4,
        estimate=estimate,
    )

    # check output basics, types and shapes
    assert len(tracemodels) == len(order_list)
    for order in order_list:
        tm = tracemodels[f"Order {order}"]
        assert tm.dtype == np.float64
        assert tm.shape == DATA_SHAPE
        # should be some nans in the trace model but not all
        assert 0 < np.sum(np.isfinite(tm)) < tm.size
    for q in [tikfacs_out["Order 1"], tikfacs_out["Order 2"]]:
        assert isinstance(q, float)
        assert np.isfinite(q)
    assert wave_grid.dtype == np.float64
    for spec in spec_list:
        assert isinstance(spec, SpecModel)

    factors = np.array([getattr(spec.meta.soss_extract1d, "factor", -1) for spec in spec_list])
    chi2s = np.array([getattr(spec.meta.soss_extract1d, "chi2", -1) for spec in spec_list])
    orders = np.array([spec.spectral_order for spec in spec_list])
    colors = np.array([spec.meta.soss_extract1d.color_range for spec in spec_list])

    assert tikfacs_out["Order 1"] in factors

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
            nochi = this_chi2s < 0  # we set this flag above if getattr(spec, "chi2") fails

            # _model_single_order is set up so that the final/best spectrum is last in the list
            # it lacks chi2 calculations
            assert np.sum(nochi) == 1
            assert np.where(nochi)[0][0] == len(this_chi2s) - 1

            # it represents the best tikhonov factor for that order-color combination
            # which is not necessarily the same as the top-level tikfac for the blue part of order 2
            # but it is the same for the red part of order 1 and the red part of order 2
            if color == "RED":
                assert this_factors[-1] == tikfacs_out["Order 1"]

            # check that the output spectra contain good data
            for spec in this_spec:
                spec = np.array([[s[0], s[1]] for s in spec.spec_table])
                assert np.sum(np.isfinite(spec[:, 0])) > 0

    # check that all order-color combinations have the same number of spectra
    n_good = np.array(n_good)
    assert np.all(n_good >= 1)
    assert np.all(n_good - n_good[0] == 0)


def test_model_image_tikfac_specified(
    monkeypatch_setup,
    imagemodel,
    detector_mask,
    detector_models,
):
    """Ensure spec_list is a single-element list per order if tikfac is specified"""
    scidata, scierr = imagemodel

    order_list = [1, 2]
    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, wavelengths = _compute_box_weights(
        detector_models, DATA_SHAPE, box_width, orders_requested=order_list
    )

    integration = Integration(
        scidata,
        scierr,
        detector_mask,
        refmask,
        detector_models,
        box_weights,
        do_bkgsub=False,
        extract_order3=False,
    )
    estimate = integration.estim_flux_first_order()
    wave_grid = integration.make_decontamination_grid(
        estimate,
        rtol=1e-3,
        max_grid_size=1000000,
        n_os=2,
    )

    tikfacs_in = {"Order 1": 1e-7, "Order 2": 1e-6}
    tracemodels, spec_list, tikfacs_out = integration.model_image(
        wave_grid,
        tikfacs_in=tikfacs_in,
        threshold=1e-4,
        estimate=estimate,
    )
    # check that spec_list is a single-element list per order in this case
    assert len(spec_list) == 3
    assert tikfacs_out["Order 1"] == tikfacs_in["Order 1"]
    assert tikfacs_out["Order 2"] == tikfacs_in["Order 2"]


def test_model_image_wavegrid_specified(
    monkeypatch_setup,
    imagemodel,
    detector_mask,
    detector_models,
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
        detector_models, DATA_SHAPE, box_width, orders_requested=order_list
    )

    integration = Integration(
        scidata,
        scierr,
        detector_mask,
        refmask,
        detector_models,
        box_weights,
        do_bkgsub=False,
        extract_order3=False,
    )
    estimate = integration.estim_flux_first_order()

    tikfacs_in = {"Order 1": 1e-7, "Order 2": 1e-6}
    # test np.array input
    wave_grid_in = np.linspace(1.0, 2.5, 100)
    tracemodels, spec_list, tikfacs_out = integration.model_image(
        tikfacs_in=tikfacs_in,
        threshold=1e-4,
        wave_grid=wave_grid_in,
        estimate=estimate,
    )

    # test SossWaveGridModel input
    # the docs on main say this works, but I don't think it does even on main
    wave_grid_in = SossWaveGridModel()
    wave_grid_in.wavegrid = np.linspace(1.0, 2.5, 100)
    with pytest.raises(ValueError):
        tracemodels, tikfacs_out, logl, wave_grid, spec_list = integration.model_image(
            tikfacs_in=tikfacs_in,
            threshold=1e-4,
            wave_grid=wave_grid_in,
            estimate=estimate,
        )


@pytest.mark.parametrize("tikfacs_in", [None, {"Order 1": 1e-7, "Order 2": 1e-6, "Order 3": 1e-5}])
@pytest.mark.parametrize("bad_pix", ["masking", "model"])
@pytest.mark.parametrize("generate_model", [False, True])
@pytest.mark.parametrize("intermediate", [False, True])
def test_process_one_integration(
    monkeypatch_setup,
    imagemodel,
    detector_mask,
    detector_models,
    tikfacs_in,
    bad_pix,
    generate_model,
    intermediate,
):
    """Smoke test for a bunch of configurations of _process_one_integration"""
    scidata, scierr = imagemodel
    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, wavelengths = _compute_box_weights(
        detector_models, DATA_SHAPE, box_width, orders_requested=[1, 2, 3]
    )

    soss_kwargs = {
        "subtract_background": True,
        "order_3": True,
        "bad_pix": bad_pix,
        "box_width": box_width,
        "estimate": None,
        "threshold": 1e-4,
        "rtol": 1e-3,
        "max_grid_size": 1000000,
        "n_os": 2,
        "model": intermediate,
    }

    tracemodels, spec_list, atoca_list, tikfacs_out, wave_grid = _process_one_integration(
        scidata,
        scierr,
        detector_mask,
        refmask,
        detector_models,
        box_weights,
        wavelengths,
        soss_kwargs,
        wave_grid=None,
        tikfacs_in=tikfacs_in,
        generate_model=generate_model,
        int_num=5,
    )
