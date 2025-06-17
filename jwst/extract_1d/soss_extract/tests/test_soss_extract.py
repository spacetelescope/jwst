from functools import partial
import pytest
import numpy as np
from stdatamodels.jwst.datamodels import SpecModel, SossWaveGridModel

from jwst.extract_1d.soss_extract.soss_extract import (
    _model_image,
    _compute_box_weights,
)
from .conftest import DATA_SHAPE


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

    def mock_get_ref_file_args(wave, trace, thru, kern, reffiles):
        """Return the arrays from conftest instead of querying CRDS"""
        return [wave, trace, thru, kern]

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


# slow because tests multiple Tikhonov factors
@pytest.mark.slow
def test_model_image(
    monkeypatch_setup,
    imagemodel,
    detector_mask,
    ref_files,
):
    scidata, scierr = imagemodel

    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, wavelengths = _compute_box_weights(ref_files, DATA_SHAPE, box_width)

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
    )

    # check output basics, types and shapes
    assert len(tracemodels) == 2
    for order in tracemodels:
        tm = tracemodels[order]
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

    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, wavelengths = _compute_box_weights(ref_files, DATA_SHAPE, box_width)

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

    refmask = np.zeros_like(detector_mask)
    box_width = 5.0
    box_weights, wavelengths = _compute_box_weights(ref_files, DATA_SHAPE, box_width)

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
        )
