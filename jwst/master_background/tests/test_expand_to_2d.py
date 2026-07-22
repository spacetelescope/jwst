"""Test master_background.expand_to_2d."""

import gwcs
import numpy as np
import pytest
from astropy.modeling.models import Identity, Mapping, Planar2D
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.master_background import expand_to_2d


@pytest.fixture(scope="module")
def multislit_science_data():
    """
    Create multi-slit "science" data for testing.

    Yields
    ------
    input_model : `~jwst.datamodels.MultiSlitModel`
    """
    # Create a MultiSlitModel object.
    data_shape = (5, 9)
    data = np.zeros(data_shape, dtype=np.float32) + 10.0
    dq = np.zeros(data_shape, dtype=np.uint32)
    # One row of wavelengths.
    temp_wl = np.linspace(
        1.3, 4.8, num=data_shape[1], endpoint=True, retstep=False, dtype=np.float32
    )
    wavelength = np.zeros(data_shape, dtype=np.float32)
    # Shift the wavelength values from row to the next.
    dwl = 0.1
    for j in range(data_shape[0]):
        wavelength[j, :] = temp_wl + j * dwl
    wavelength = np.around(wavelength, 4)
    input_model = datamodels.MultiSlitModel()
    slit = datamodels.SlitModel(init=None, data=data, dq=dq, wavelength=wavelength)
    input_model.slits.append(slit)

    yield input_model
    input_model.close()


@pytest.fixture(scope="module")
def image_science_data():
    """
    Create "science" data for testing in ImageModel format.

    Yields
    ------
    input_model : `~jwst.datamodels.ImageModel`
    """
    data_shape = (9, 5)
    data = np.zeros(data_shape, dtype=np.float32) + 10.0
    dq = np.zeros(data_shape, dtype=np.uint32)
    input_model = datamodels.ImageModel(data=data, dq=dq)

    # Mock a wcs with a 2D linear plane for wavelengths
    w_start = 1.3
    w_stop = 4.8
    npts = data_shape[1]
    w_delta_col = (w_stop - w_start) / (npts - 1)
    w_delta_row = 0.1
    wave = Planar2D(slope_x=w_delta_col, slope_y=w_delta_row, intercept=w_start)
    det2world = Mapping((0, 1, 0, 1), n_inputs=2) | Identity(2) & wave
    input_frame = gwcs.Frame2D(name="detector")
    output_frame = gwcs.Frame2D(name="world")
    pipeline = [(input_frame, det2world), (output_frame, None)]
    wcs = gwcs.WCS(pipeline)

    input_model.meta.wcs = wcs

    yield input_model
    input_model.close()


@pytest.fixture(scope="module")
def multispec_background_data():
    """
    Create multi-spec "user-background" data for testing.

    Yields
    ------
    m_bkg_spec : `~jwst.datamodels.MultiSpecModel`
    """
    # This data type is used for creating a MultiSpecModel.
    spec_dtype = datamodels.SpecModel().get_dtype("spec_table")

    # m_bkg_spec doesn't have to be a MultiSpecModel, but that's an option.
    m_bkg_spec = datamodels.MultiSpecModel()
    wavelength = np.geomspace(1.5, 4.5, num=25, endpoint=True, dtype=np.float64)
    flux = np.zeros_like(wavelength)
    error = np.ones_like(wavelength)
    var_dummy = error.copy()
    surf_bright = np.linspace(13.0, 25.0, num=25, endpoint=True, retstep=False, dtype=np.float64)
    sb_error = np.ones_like(wavelength)
    dq = np.zeros(wavelength.shape, dtype=np.uint32)
    background = np.ones_like(wavelength)
    berror = np.ones_like(wavelength)
    # The npixels column should no longer be used.  Set it to a large value
    # to make it more obvious in case it actually is still used.
    npixels = np.zeros_like(wavelength) + 2000.0
    otab = np.array(
        list(
            zip(
                wavelength,
                flux,
                error,
                var_dummy,
                var_dummy,
                var_dummy,
                surf_bright,
                sb_error,
                var_dummy,
                var_dummy,
                var_dummy,
                dq,
                background,
                berror,
                var_dummy,
                var_dummy,
                var_dummy,
                npixels,
                strict=False,
            )
        ),
        dtype=spec_dtype,
    )
    spec = datamodels.SpecModel(spec_table=otab)
    m_bkg_spec.spec.append(spec)

    yield m_bkg_spec
    m_bkg_spec.close()


@pytest.fixture(scope="module")
def multispec_background_data_reversed():
    """
    Create reversed "user-background" data for testing.

    `expand_to_2d` uses `np.interp` for interpolation, and the wavelength
    array that is passed to `np.interp` must be increasing.  `expand_to_2d`
    is supposed to handle the case that the wavelengths are decreasing.

    Create data for checking that the results are the same even if the
    wavelength array in `m_bkg_spec` is reversed so that the values are
    decreasing, and the corresponding flux array is also reversed to retain
    the original (wavelength, flux) relation.

    Yields
    ------
    m_bkg_spec : `~jwst.datamodels.MultiSpecModel`
    """
    # This data type is used for creating a MultiSpecModel.
    spec_dtype = datamodels.SpecModel().get_dtype("spec_table")

    m_bkg_spec = datamodels.MultiSpecModel()
    wavelength = np.geomspace(1.5, 4.5, num=25, endpoint=True, dtype=np.float64)[::-1]
    flux = np.zeros_like(wavelength)
    error = np.ones_like(wavelength)
    var_dummy = error.copy()
    surf_bright = np.linspace(13.0, 25.0, num=25, endpoint=True, retstep=False, dtype=np.float64)[
        ::-1
    ]
    sb_error = np.ones_like(wavelength)
    dq = np.zeros(wavelength.shape, dtype=np.uint32)
    background = np.ones_like(wavelength)
    berror = np.ones_like(wavelength)
    npixels = np.ones_like(wavelength)
    otab = np.array(
        list(
            zip(
                wavelength,
                flux,
                error,
                var_dummy,
                var_dummy,
                var_dummy,
                surf_bright,
                sb_error,
                var_dummy,
                var_dummy,
                var_dummy,
                dq,
                background,
                berror,
                var_dummy,
                var_dummy,
                var_dummy,
                npixels,
                strict=False,
            )
        ),
        dtype=spec_dtype,
    )
    spec = datamodels.SpecModel(spec_table=otab)
    m_bkg_spec.spec.append(spec)

    yield m_bkg_spec
    m_bkg_spec.close()


@pytest.fixture(scope="module")
def combined_spec_background_data():
    """
    Create "user-background" data for testing in combined spec format.

    Yields
    ------
    m_bkg_spec : `~jwst.datamodels.CombinedSpecModel`
    """
    # This is the data type of an output table from combine_1d.
    spec_table_dtype = datamodels.CombinedSpecModel().get_dtype("spec_table")

    wavelength = np.geomspace(1.5, 4.5, num=25, endpoint=True, dtype=np.float64)
    flux = np.zeros_like(wavelength)
    error = np.ones_like(wavelength)
    surf_bright = np.linspace(13.0, 25.0, num=25, endpoint=True, retstep=False, dtype=np.float64)
    sb_error = np.ones_like(wavelength)
    dq = np.zeros(wavelength.shape, dtype=np.uint32)
    weight = np.ones_like(wavelength)
    n_input = np.ones_like(wavelength)  # yes, float64
    data = np.array(
        list(
            zip(wavelength, flux, error, surf_bright, sb_error, dq, weight, n_input, strict=False)
        ),
        dtype=spec_table_dtype,
    )
    m_bkg_spec = datamodels.CombinedSpecModel(spec_table=data)

    yield m_bkg_spec
    m_bkg_spec.close()


@pytest.fixture(scope="module")
def multispec_truth():
    """
    Create an array of comparison values, for testing multispec background input.

    Returns
    -------
    truth : ndarray, 2-D, float64
        An array to compare with the data in the output from `expand_to_2d`.
    """
    truth = np.array(
        [
            [0.0, 14.603571, 17.057364, 19.059263, 20.748844, 22.21304, 23.506573, 24.658548, 0.0],
            [0.0, 15.213889, 17.548525, 19.470137, 21.102207, 22.524107, 23.77871, 24.906878, 0.0],
            [13.0, 15.792757, 18.018988, 19.86396, 21.444334, 22.822334, 24.04857, 0.0, 0.0],
            [13.702182, 16.342768, 18.469252, 20.244974, 21.773643, 23.115166, 24.308529, 0.0, 0.0],
            [14.364901, 16.866348, 18.900745, 20.614536, 22.095966, 23.40005, 24.565424, 0.0, 0.0],
        ],
        dtype=np.float64,
    )

    return truth


@pytest.fixture(scope="module")
def combined_spec_truth():
    """Create an array of comparison values, for testing combined spec background input.

    Returns
    -------
    truth : ndarray, 2-D, float64
        An array to compare with the data in the output from `expand_to_2d`.
    """
    truth = np.array(
        [
            [0.0, 17.057364, 20.748844, 23.506573, 0.0],
            [0.0, 17.548523, 21.102207, 23.77871, 0.0],
            [13.0, 18.018988, 21.444334, 24.04857, 0.0],
            [13.702181, 18.469252, 21.773643, 24.308529, 0.0],
            [14.3649, 18.900745, 22.095966, 24.565424, 0.0],
            [14.9912815, 19.31606, 22.408163, 24.813755, 0.0],
            [15.5804825, 19.716778, 22.710499, 0.0, 0.0],
            [16.139992, 20.104378, 23.008335, 0.0, 0.0],
            [16.672644, 20.479303, 23.293217, 0.0, 0.0],
        ],
        dtype=np.float64,
    )

    return truth


def test_expand_to_2d_multispec(multislit_science_data, multispec_background_data, multispec_truth):
    input_data = multislit_science_data  # MultiSlitModel
    m_bkg_spec = multispec_background_data  # MultiSpecModel
    bkg = expand_to_2d.expand_to_2d(input_data, m_bkg_spec)

    assert np.allclose(bkg.slits[0].data, multispec_truth, rtol=1.0e-6)


def test_expand_to_2d_multispec_reversed(
    multislit_science_data, multispec_background_data_reversed, multispec_truth
):
    # Same input data as in the first test.
    input_data = multislit_science_data  # MultiSlitModel

    # Check that expand_to_2d works if the wavelength array is reversed.
    # The flux array is also reversed, so the results should be unchanged.
    m_bkg_spec = multispec_background_data_reversed
    bkg = expand_to_2d.expand_to_2d(input_data, m_bkg_spec)

    # Same truth array as in the first test.
    assert np.allclose(bkg.slits[0].data, multispec_truth, rtol=1.0e-6)


def test_expand_to_2d_multispec_ignore_extra(
    caplog, multislit_science_data, multispec_background_data, multispec_truth
):
    input_data = multislit_science_data  # MultiSlitModel
    m_bkg_spec = multispec_background_data  # MultiSpecModel

    # Add an extra spectrum
    m_bkg_spec.spec.append(m_bkg_spec.spec[0].instance.copy())

    # Expect the same results: only the first is used
    bkg = expand_to_2d.expand_to_2d(input_data, m_bkg_spec)
    assert np.allclose(bkg.slits[0].data, multispec_truth, rtol=1.0e-6)

    # Warning is logged
    assert "contains multiple spectra" in caplog.text


def test_expand_to_2d_combined_spec(
    image_science_data, combined_spec_background_data, combined_spec_truth
):
    input_data = image_science_data  # ImageModel
    m_bkg_spec = combined_spec_background_data  # CombinedSpecModel
    bkg = expand_to_2d.expand_to_2d(input_data, m_bkg_spec)

    assert np.allclose(bkg.data, combined_spec_truth, rtol=1.0e-6)


def test_expand_to_2d_container(multislit_science_data, multispec_background_data, multispec_truth):
    input_data = ModelContainer([multislit_science_data])  # Container of MultiSlitModel
    m_bkg_spec = multispec_background_data  # MultiSpecModel
    bkg = expand_to_2d.expand_to_2d(input_data, m_bkg_spec)

    # output is a container with a single MultiSlitModel
    assert isinstance(bkg, ModelContainer)
    assert len(bkg) == 1
    assert isinstance(bkg[0], datamodels.MultiSlitModel)

    # check against expected truth array
    assert np.allclose(bkg[0].slits[0].data, multispec_truth, rtol=1.0e-6)


def test_create_bkg_unsupported_model():
    # Unsupported model raises error
    input_model = datamodels.RampModel()
    with pytest.raises(TypeError, match="not supported"):
        expand_to_2d.create_bkg(input_model, None, None)


def test_bkg_for_multislit_missing_wl(multislit_science_data):
    input_data = multislit_science_data.copy()

    # Modify input to unset wavelengths
    input_data.slits[0].wavelength = None

    # Raises error
    with pytest.raises(RuntimeError, match="Can't determine wavelengths"):
        expand_to_2d.bkg_for_multislit(input_data, None, None)


@pytest.mark.parametrize("allow_mos", [True, False])
def test_bkg_for_multislit_allow_mos(
    caplog, multislit_science_data, multispec_background_data, allow_mos
):
    input_data = multislit_science_data.copy()
    input_data.meta.exposure.type = "NRS_MSASPEC"

    m_bkg_spec = multispec_background_data
    tab_wl = m_bkg_spec.spec[0].spec_table["WAVELENGTH"]
    tab_bg = m_bkg_spec.spec[0].spec_table["SURF_BRIGHT"]

    # Warns and sets all backgrounds and dq to 0 with allow_mos=False
    bg = expand_to_2d.bkg_for_multislit(input_data, tab_wl, tab_bg, allow_mos=allow_mos)
    if not allow_mos:
        assert "not supported for NIRSpec MOS spectra" in caplog.text
        for slit in bg.slits:
            assert np.allclose(slit.data, 0.0)
            assert np.all(slit.dq == 0)
    else:
        for slit in bg.slits:
            assert not np.allclose(slit.data, 0.0)


def test_bkg_for_multislit_nrs_fs_point(
    caplog, multislit_science_data, multispec_background_data, multispec_truth
):
    input_data = multislit_science_data.copy()
    input_data.meta.exposure.type = "NRS_FIXEDSLIT"
    input_data.slits[0].source_type = "POINT"
    shape = input_data.slits[0].data.shape
    input_data.slits[0].pathloss_point = np.full(shape, 2.0)
    input_data.slits[0].pathloss_uniform = np.full(shape, 3.0)
    input_data.slits[0].flatfield_point = np.full(shape, 4.0)
    input_data.slits[0].flatfield_uniform = np.full(shape, 5.0)
    input_data.slits[0].photom_point = np.full(shape, 6.0)
    input_data.slits[0].photom_uniform = np.full(shape, 7.0)

    m_bkg_spec = multispec_background_data
    tab_wl = m_bkg_spec.spec[0].spec_table["WAVELENGTH"]
    tab_bg = m_bkg_spec.spec[0].spec_table["SURF_BRIGHT"]

    bkg = expand_to_2d.bkg_for_multislit(input_data, tab_wl, tab_bg)

    # Background is corrected by ratio of point to uniform for 3 corrections.
    # Note: the pathloss and flat corrections are inverted.
    corrected = multispec_truth * (3 / 2) * (5 / 4) * (6 / 7)
    np.testing.assert_allclose(bkg.slits[0].data, corrected, rtol=1.0e-6)


def test_bkg_for_image_missing_wl(image_science_data):
    input_data = image_science_data.copy()

    # Modify input to unset wcs
    input_data.meta.wcs = None

    # Raises error
    with pytest.raises(RuntimeError, match="Can't determine wavelengths"):
        expand_to_2d.bkg_for_image(input_data, None, None)


def test_bkg_for_ifu_image_miri(
    image_science_data, combined_spec_background_data, combined_spec_truth
):
    input_data = image_science_data.copy()
    input_data.meta.instrument.name = "MIRI"

    m_bkg_spec = combined_spec_background_data
    tab_wl = m_bkg_spec.spec_table["WAVELENGTH"]
    tab_bg = m_bkg_spec.spec_table["SURF_BRIGHT"]

    bkg = expand_to_2d.bkg_for_ifu_image(input_data, tab_wl, tab_bg)
    assert np.allclose(bkg.data, combined_spec_truth, rtol=1.0e-6)


def test_bkg_for_ifu_unsupported_instrument(image_science_data):
    # Unsupported instrument raises error
    input_model = image_science_data.copy()
    input_model.meta.instrument.name = "NIRISS"
    with pytest.raises(TypeError, match="not supported"):
        expand_to_2d.bkg_for_ifu_image(input_model, None, None)
