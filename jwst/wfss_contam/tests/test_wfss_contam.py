import types

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.assign_wcs.tests.test_niriss import create_imaging_wcs
from jwst.wfss_contam.wfss_contam import (
    SlitOverlapError,
    UnmatchedSlitIDError,
    _apply_magnitude_limit,
    _build_simulated_image_from_slits,
    _cut_frame_to_match_slit,
    _find_matching_simul_slit,
    _validate_orders_against_reference,
    contam_corr,
    match_backplane_prefer_first,
)

VALID_ORDERS = [-1, 0, 1, 2, 3]


@pytest.fixture(scope="module")
def contam():
    return np.ones((10, 10)) * 0.1


@pytest.fixture(scope="module")
def slit0():
    slit = dm.SlitModel(data=np.ones((5, 3)))
    slit.xstart = 2
    slit.ystart = 3
    slit.xsize = 3
    slit.ysize = 5
    slit.meta.wcsinfo.spectral_order = 1
    slit.source_id = 1
    return slit


@pytest.fixture(scope="module")
def slit1():
    slit = dm.SlitModel(data=np.ones((4, 4)) * 0.5)
    slit.xstart = 3
    slit.ystart = 2
    slit.xsize = 4
    slit.ysize = 4
    return slit


@pytest.fixture(scope="module")
def slit2():
    slit = dm.SlitModel(data=np.ones((3, 5)) * 0.1)
    slit.xstart = 300
    slit.ystart = 200
    slit.xsize = 5
    slit.ysize = 3
    return slit


def test_find_matching_simul_slit(slit0):
    sids = [0, 1, 1]
    orders = [1, 1, 2]
    idx = _find_matching_simul_slit(slit0, sids, orders)
    assert idx == 1


def test_find_matching_simul_slit_no_match(slit0):
    sids = [0, 1, 1]
    orders = [1, 2, 2]
    with pytest.raises(UnmatchedSlitIDError):
        _find_matching_simul_slit(slit0, sids, orders)


def test_cut_frame_to_match_slit(slit0, contam):
    cut_contam = _cut_frame_to_match_slit(contam, slit0)
    assert cut_contam.shape == (5, 3)
    assert np.all(cut_contam == 0.1)


def test_match_backplane_prefer_first(slit0, slit1):
    slit1_final = match_backplane_prefer_first(slit0.copy(), slit1.copy())

    assert slit1_final.xstart == slit0.xstart
    assert slit1_final.ystart == slit0.ystart
    assert slit1_final.xsize == slit0.xsize
    assert slit1_final.ysize == slit0.ysize
    assert slit1_final.data.shape == slit0.data.shape
    assert np.count_nonzero(slit1_final.data) == 6


def test_match_backplane_prefer_first_yoffset():
    """
    Cover an indexing bug in match_backplane_prefer_first.

    When slit1 starts below slit0 (y1 > 0), the old code dropped the first y1 rows
    of data1 and placed the wrong rows in the output.
    """
    slit0 = dm.SlitModel(data=np.zeros((4, 4)))
    slit0.xstart, slit0.ystart = 0, 0  # origin at (0, 0)
    slit0.xsize, slit0.ysize = 4, 4

    # Give each row a distinct value so we can tell which rows end up where.
    data1 = np.array([[1.0] * 4, [2.0] * 4, [3.0] * 4, [4.0] * 4])
    wl1 = np.array([[0.1] * 4, [0.2] * 4, [0.3] * 4, [0.4] * 4])
    slit1 = dm.SlitModel(data=data1.copy())
    slit1.xstart, slit1.ystart = 0, 2  # starts 2 rows below slit0
    slit1.xsize, slit1.ysize = 4, 4

    slit1_out = match_backplane_prefer_first(slit0.copy(), slit1)

    # Overlap: slit1 rows 0-1 (values 1.0, 2.0) map to backplane rows 2-3.
    # Rows 0-1 of the backplane are outside slit1 so should be zero.
    assert slit1_out.data.shape == (4, 4)
    assert np.all(slit1_out.data[0:2, :] == 0.0)
    assert np.all(slit1_out.data[2, :] == 1.0)
    assert np.all(slit1_out.data[3, :] == 2.0)


def test_common_slit_prefer_expected_raise(slit0, slit2):
    with pytest.raises(SlitOverlapError):
        match_backplane_prefer_first(slit0.copy(), slit2.copy())


def test_apply_magnitude_limit(photom_ref_model, source_catalog):
    magnitude_limit = 23
    min_relresp_order1 = 1
    order = -1
    order_idx = np.where(photom_ref_model.phot_table["order"] == order)[0]
    sens_wave = photom_ref_model.phot_table["wavelength"][order_idx]
    sens_response = photom_ref_model.phot_table["relresponse"][order_idx]
    sources = _apply_magnitude_limit(
        order, source_catalog, sens_wave, sens_response, magnitude_limit, min_relresp_order1
    )
    assert sources == [17, 39, 51, 82, 93]


def test_apply_magnitude_limit_no_sources(photom_ref_model, source_catalog):
    magnitude_limit = 10
    min_relresp_order1 = 1
    order = -1
    order_idx = np.where(photom_ref_model.phot_table["order"] == order)[0]
    sens_wave = photom_ref_model.phot_table["wavelength"][order_idx]
    sens_response = photom_ref_model.phot_table["relresponse"][order_idx]
    sources = _apply_magnitude_limit(
        order, source_catalog, sens_wave, sens_response, magnitude_limit, min_relresp_order1
    )
    assert sources is None


def test_constrain_orders(log_watcher):
    # normal case
    orders = [1, 2]
    constrained_orders = _validate_orders_against_reference(orders, VALID_ORDERS)
    assert np.array_equal(constrained_orders, np.array(orders))

    # orders is None
    constrained_orders = _validate_orders_against_reference(None, VALID_ORDERS)
    assert np.array_equal(constrained_orders, np.array(VALID_ORDERS))


def test_constrain_orders_error_empty(log_watcher):
    # orders is empty
    orders = []
    watcher = log_watcher(
        "jwst.wfss_contam.wfss_contam",
        message="None of the requested spectral orders ",
        level="error",
    )
    constrained_orders = _validate_orders_against_reference(orders, VALID_ORDERS)
    assert not constrained_orders

    # none of the orders match spec_orders
    orders = [4, 5]
    watcher = log_watcher(
        "jwst.wfss_contam.wfss_contam",
        message="None of the requested spectral orders ",
        level="error",
    )
    constrained_orders = _validate_orders_against_reference(orders, VALID_ORDERS)
    assert not constrained_orders


def test_constrain_orders_warn_subset(log_watcher):
    # only a subset of orders match
    orders = [1, 4]
    watcher = log_watcher(
        "jwst.wfss_contam.wfss_contam", message="Skipping undefined orders", level="warning"
    )
    constrained_orders = _validate_orders_against_reference(orders, VALID_ORDERS)
    assert np.array_equal(constrained_orders, np.array([1]))
    watcher.assert_seen()


def test_build_simulated_image_from_slits():
    shape = (10, 10)
    simulated_slits = dm.MultiSlitModel()

    slit_a = dm.SlitModel(data=np.ones((3, 4)) * 2.0)
    slit_a.xstart = 2
    slit_a.ystart = 2
    slit_a.xsize = 4
    slit_a.ysize = 3

    slit_b = dm.SlitModel(data=np.ones((3, 4)) * 5.0)
    slit_b.xstart = 4
    slit_b.ystart = 3
    slit_b.xsize = 4
    slit_b.ysize = 3

    simulated_slits.slits.append(slit_a)
    simulated_slits.slits.append(slit_b)

    full_image = _build_simulated_image_from_slits(simulated_slits, shape)

    assert full_image.shape == shape
    assert np.all(full_image[1:4, 1:3] == 2.0)
    assert np.all(full_image[2:5, 5:7] == 5.0)
    assert np.all(full_image[2:4, 3:5] == 7.0)  # overlap region: values add
    # pixels not covered by either slit are zero
    covered = np.zeros(shape, dtype=bool)
    covered[1:5, 1:7] = True  # bounding box enclosing both slits
    assert np.all(full_image[~covered] == 0.0)


def test_build_simulated_image_from_slits_overflow():
    """Slit data extending beyond the frame boundary should be clipped without error."""
    shape = (5, 5)
    simulated_slits = dm.MultiSlitModel()

    slit = dm.SlitModel(data=np.ones((4, 4)) * 3.0)
    slit.xstart = 4
    slit.ystart = 4
    slit.xsize = 4
    slit.ysize = 4

    simulated_slits.slits.append(slit)

    full_image = _build_simulated_image_from_slits(simulated_slits, shape)

    assert full_image.shape == shape
    assert np.all(full_image[3:5, 3:5] == 3.0)
    # other pixels are zero
    assert np.all(full_image[:3, :] == 0.0)
    assert np.all(full_image[:, :3] == 0.0)


_ITER_NROWS, _ITER_NCOLS = 40, 80
_ITER_FRAME_SHAPE = (80, 90)
_ITER_OVERLAP = 20
_ITER_XA, _ITER_YA = 1, 1
_ITER_XB, _ITER_YB = 1, 21


def _iter_geometry():
    """Return the shared spectral geometry used by the iteration fixtures and test."""
    lam = np.linspace(1.0, 3.0, _ITER_NCOLS)
    flat = np.ones((_ITER_NROWS, _ITER_NCOLS))
    dlam = np.broadcast_to(lam - 2.0, (_ITER_NROWS, _ITER_NCOLS)).copy()
    tilt = flat * dlam

    true_A = 1.0 * flat + 0.6 * tilt
    true_B = 1.0 * flat - 0.8 * tilt
    return flat, tilt, true_A, true_B


@pytest.fixture
def two_source_input(tmp_cwd, grism_wcs):
    """
    MultiSlitModel with two slits whose grism traces partially overlap.

    Source A (ID=1) is placed at (xstart=1, ystart=1) and source B (ID=2) at
    (xstart=1, ystart=21) on the full-frame detector, giving a 20-row overlap.
    Observed slit data in each slit is contaminated by its neighbor.

    Direct image and segmentation map files corresponding to the same scenario
    are written to the ``tmp_cwd``. They must exist with proper file name and type
    but are trivial because the test monkeypatches ``Observation``.
    """
    _flat, _tilt, true_A, true_B = _iter_geometry()

    obs_A_data = true_A.copy()
    obs_A_data[_ITER_NROWS - _ITER_OVERLAP :, :] += true_B[:_ITER_OVERLAP, :]
    obs_B_data = true_B.copy()
    obs_B_data[:_ITER_OVERLAP, :] += true_A[_ITER_NROWS - _ITER_OVERLAP :, :]

    model = dm.MultiSlitModel()
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.filter = "GR150C"
    model.meta.instrument.pupil = "F200W"

    for sid, xstart, ystart, data in [
        (1, _ITER_XA, _ITER_YA, obs_A_data),
        (2, _ITER_XB, _ITER_YB, obs_B_data),
    ]:
        slit = dm.SlitModel()
        slit.source_id = sid
        slit.meta.wcs = grism_wcs
        slit.meta.wcsinfo.spectral_order = 1
        slit.xstart = xstart
        slit.ystart = ystart
        slit.data = data.astype(np.float32)
        slit.dq = np.zeros(data.shape, dtype=np.uint32)
        slit.xsize = data.shape[1]
        slit.ysize = data.shape[0]
        model.slits.append(slit)

    direct = dm.ImageModel(data=np.ones(_ITER_FRAME_SHAPE))
    direct.meta.wcs = create_imaging_wcs("F200W")
    direct.save("direct_image.fits")
    seg = dm.SegmentationMapModel(data=np.zeros(_ITER_FRAME_SHAPE, dtype=np.uint32))
    seg.save("seg_map.fits")
    direct.close()
    seg.close()

    model.meta.direct_image = "direct_image.fits"
    model.meta.segmentation_map = "seg_map.fits"

    yield model

    model.close()


def _make_slit(source_id, order, xstart, ystart, data, **fluxmodels):
    slit = dm.SlitModel()
    slit.source_id = source_id
    slit.meta.wcsinfo.spectral_order = order
    slit.xstart = xstart
    slit.ystart = ystart
    slit.data = data.astype(np.float32)
    slit.dq = np.zeros(data.shape, dtype=np.uint32)
    slit.xsize = data.shape[1]
    slit.ysize = data.shape[0]
    for k, v in fluxmodels.items():
        setattr(slit, k, v.astype(np.float32))
    return slit


def test_iteration_improves_contamination_correction(
    two_source_input, wavelengthrange_ref_model, photom_ref_model_niriss, monkeypatch
):
    """
    Test that iterative contamination correction is working.

    The fixture makes two sources (A and B) with partially overlapping grism traces.
    Source A has trace spanning detector rows 1-40.
    Source B has trace spanning rows 38-77 (a 3-row overlap).

    Source A has a linear spectrum: f_A = 1.0 + 0.6 * dlam
    Source B has a different linear spectrum: f_B = 1.0 - 0.8 * dlam

    In the overlapping rows each slit is contaminated by the other.
    A flat-spectrum simulation (``polyfit_degree=None``) leaves large spectral
    residuals in those rows.  With ``polyfit_degree=1`` the fitted shapes
    reduce the residual, and with ``n_iterations=3`` we converge to near
    exact recovery.
    """
    flat, tilt, true_A, true_B = _iter_geometry()

    def MockObservation():
        """Use SimpleNamespace to make a mock Observation object with pre-built simulated slits."""
        obs = types.SimpleNamespace()
        simul_A = _make_slit(1, 1, _ITER_XA, _ITER_YA, flat.copy(), fluxmodel_1=tilt.copy())
        simul_B = _make_slit(2, 1, _ITER_XB, _ITER_YB, flat.copy(), fluxmodel_1=tilt.copy())
        obs.simulated_slits = dm.MultiSlitModel()
        obs.simulated_slits.slits.extend([simul_A, simul_B])
        obs.simulated_image = np.zeros(_ITER_FRAME_SHAPE)
        obs.source_ids = {1, 2}
        # disperse() becomes a function that does nothing.
        # Normally it would set the attributes we pre-set above.
        obs.disperse_order = lambda *args, **kwargs: None
        return obs

    # Inject our fake Observation so no real dispersion is performed
    monkeypatch.setattr(
        "jwst.wfss_contam.wfss_contam.Observation",
        lambda *args, **kwargs: MockObservation(),
    )
    # Bypass the grism-WCS order-validity check
    monkeypatch.setattr(
        "jwst.wfss_contam.wfss_contam._validate_orders_against_transform",
        lambda _wcs, orders: orders,
    )

    def _run(n_iter, polyfit_degree):
        result, *_ = contam_corr(
            two_source_input.copy(),
            wavelengthrange_ref_model,
            photom_ref_model_niriss,
            "none",
            orders=[1],
            polyfit_degree=polyfit_degree,
            n_iterations=n_iter,
        )
        return [np.array(s.data) for s in result.slits]

    corrected_flat = _run(1, polyfit_degree=None)
    corrected_1iter = _run(1, polyfit_degree=1)
    corrected_3iter = _run(3, polyfit_degree=1)

    true_signals = [true_A, true_B]

    def _total_rms(corrected):
        return sum(
            np.sqrt(np.mean((c - t) ** 2)) for c, t in zip(corrected, true_signals, strict=True)
        )

    rms_flat = _total_rms(corrected_flat)
    rms_1iter = _total_rms(corrected_1iter)
    rms_3iter = _total_rms(corrected_3iter)

    assert rms_1iter < rms_flat
    assert rms_3iter < rms_1iter
