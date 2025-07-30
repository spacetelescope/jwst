import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.io import fits
from astropy.table import Table

from jwst.extract_1d.tests import helpers


@pytest.fixture
def simple_wcs():
    """
    Mock a horizontal dispersion WCS with a simple callable function.

    Some other expected WCS attributes are also mocked with placeholder values:
       - bounding_box
       - get_transform
       - available_frames

    Returns
    -------
    callable
        A function that will return mock values for RA, Dec, wave,
        given x and y coordinates.
    """
    return helpers.simple_wcs_func()


def simple_wcs_transpose_func():
    """
    Mock a vertical dispersion WCS with a simple callable function.

    Some other expected WCS attributes are also mocked with placeholder values:
       - bounding_box
       - get_transform
       - backward_transform
       - available_frames

    Returns
    -------
    callable
        A function that will return mock values for RA, Dec, wave,
        given x and y coordinates.
    """
    return helpers.simple_wcs_transpose_func()


@pytest.fixture()
def simple_wcs_transpose():
    """
    Make simple_wcs_transpose_func available as a fixture.

    Returns
    -------
    callable
        A function that will return mock values for RA, Dec, wave,
        given x and y coordinates.
    """
    return simple_wcs_transpose_func()


@pytest.fixture()
def simple_wcs_ifu():
    """
    Mock an IFU WCS with a simple callable function.

    The bounding_box attribute is also mocked with a placeholder value.

    Returns
    -------
    callable
        A function that will return mock values for RA, Dec, wave,
        given x and y coordinates.
    """
    return helpers.simple_wcs_ifu_func()


@pytest.fixture
def mock_nirspec_fs_one_slit():
    """
    Mock one slit in NIRSpec FS mode.

    Yields
    ------
    SlitModel
        The mock model.
    """
    model = helpers.mock_nirspec_fs_one_slit_func()
    yield model
    model.close()


@pytest.fixture()
def mock_nirspec_mos():
    """
    Mock three slits in NIRSpec MOS mode.

    Yields
    ------
    MultiSlitModel
        The mock model.
    """
    model = helpers.mock_nirspec_mos_func()
    yield model
    model.close()


@pytest.fixture()
def mock_nirspec_bots():
    """
    Mock a single slit with 10 integrations in NIRSpec BOTS mode.

    Yields
    ------
    CubeModel
        The mock model.
    """
    model = helpers.mock_nirspec_bots_func()
    yield model
    model.close()


@pytest.fixture()
def mock_miri_lrs_fs():
    """
    Mock a spectral image in MIRI LRS FS mode.

    Yields
    ------
    ImageModel
        The mock model.
    """
    model = helpers.mock_miri_lrs_fs_func()
    yield model
    model.close()


@pytest.fixture()
def mock_miri_ifu():
    """
    Mock an IFU cube in MIRI MRS mode.

    Yields
    ------
    IFUCubeModel
        The mock model.
    """
    model = helpers.mock_miri_ifu_func()
    yield model
    model.close()


@pytest.fixture
def mock_niriss_wfss_l2():
    """
    Mock 3 slits in NIRISS WFSS mode, level 2 style.

    The slits correspond to a single exposure, with one slit per extracted source.

    Yields
    ------
    MultiSlitModel
        The mock model.
    """
    model = helpers.mock_nis_wfss_l2()
    yield model
    model.close()


@pytest.fixture()
def mock_niriss_wfss_l3():
    """
    Mock 3 slits in NIRISS WFSS mode, level 3 style.

    Here the container has one MultiSlitModel per source, and each model has one
    slit per exposure.

    Yields
    ------
    SourceModelContainer
        The mock model.
    """
    sources = helpers.mock_nis_wfss_l3()
    yield sources[0]
    for source in sources:
        source.close()


@pytest.fixture()
def mock_niriss_soss():
    """
    Mock a multi-integration cube with metadata for NIRISS SOSS mode.

    Yields
    ------
    CubeModel
        The mock model.
    """
    model = helpers.mock_niriss_soss_func()
    yield model
    model.close()


@pytest.fixture()
def mock_niriss_soss_256():
    """
    Make mock_niriss_soss_256_func available as a fixture.

    Yields
    ------
    CubeModel
        The mock model.
    """
    model = helpers.mock_niriss_soss_256_func()
    yield model
    model.close()


@pytest.fixture()
def mock_niriss_soss_96():
    """
    Make mock_niriss_soss_96_func available as a fixture.

    Yields
    ------
    CubeModel
        The mock model.
    """
    model = helpers.mock_niriss_soss_96_func()
    yield model
    model.close()


@pytest.fixture()
def mock_one_spec():
    """
    Mock one simple spectrum in a MultiSpecModel.

    Yields
    ------
    MultiSpecModel
        The mock model.
    """
    model = dm.MultiSpecModel()
    spec_model = helpers.make_spec_model()
    model.spec.append(spec_model)

    yield model
    model.close()


@pytest.fixture()
def mock_10_spec():
    """
    Mock 10 simple spectra in a MultiSpecModel.

    Yields
    ------
    MultiSpecModel
        The mock model.
    """
    model = dm.MultiSpecModel()

    for i in range(10):
        spec_model = helpers.make_spec_model(name=f"slit{i + 1}", value=i + 1)
        model.spec.append(spec_model)

    yield model
    model.close()


@pytest.fixture()
def mock_10_multi_int_spec():
    """
    Mock 10 simple spectra in a TSOMultiSpecModel.

    Yields
    ------
    TSOMultiSpecModel
        The mock model.
    """
    model = helpers.make_tso_spec_model(n_spectra=10)
    yield model
    model.close()


@pytest.fixture()
def mock_2_multi_int_spec():
    """
    Mock 2 simple spectra in a TSOMultiSpecModel.

    Used for generating spectra that do not match the int_times
    table in a 10-integration input.

    Yields
    ------
    TSOMultiSpecModel
        The mock model.
    """
    model = helpers.make_tso_spec_model(n_spectra=2)
    yield model
    model.close()


@pytest.fixture()
def miri_lrs_apcorr():
    """
    Mock a MIRI LRS aperture correction model.

    Yields
    ------
    MirLrsApcorrModel
        The mock model.
    """
    table = Table(
        {
            "subarray": ["FULL", "SLITLESSPRISM"],
            "wavelength": [[1, 2, 3], [1, 2, 3]],
            "nelem_wl": [3, 3],
            "size": [[1, 2, 3], [1, 2, 3]],
            "nelem_size": [3, 3],
            "apcorr": np.full((2, 3, 3), 0.5),
            "apcorr_err": np.full((2, 3, 3), 0.01),
        }
    )
    table = fits.table_to_hdu(table)
    table.header["EXTNAME"] = "APCORR"
    table.header["SIZEUNIT"] = "pixels"
    hdul = fits.HDUList([fits.PrimaryHDU(), table])

    apcorr_model = dm.MirLrsApcorrModel(hdul)
    yield apcorr_model
    apcorr_model.close()


@pytest.fixture()
def miri_lrs_apcorr_file(tmp_path, miri_lrs_apcorr):
    """
    Mock a MIRI LRS aperture correction reference file.

    Returns
    -------
    str
        Path to the reference file.
    """
    filename = str(tmp_path / "miri_lrs_apcorr.fits")
    miri_lrs_apcorr.save(filename)
    return filename


@pytest.fixture()
def nirspec_fs_apcorr():
    """
    Mock a NIRSpec FS aperture correction model.

    Yields
    ------
    NrsFsApcorrModel
        The mock model.
    """
    table = Table(
        {
            "filter": ["clear", "f290lp"],
            "grating": ["prism", "g395h"],
            "slit": ["S200A1", "S200A1"],
            "wavelength": [[1, 2, 3], [1, 2, 3]],
            "nelem_wl": [3, 3],
            "size": np.full((2, 3, 3), [0.3, 0.5, 1]),
            "nelem_size": [3, 3],
            "pixphase": [[0.0, 0.25, 0.5], [0.0, 0.25, 0.5]],
            "apcorr": np.full((2, 3, 3, 3), 0.5),
            "apcorr_err": np.full((2, 3, 3, 3), 0.01),
        }
    )
    table = fits.table_to_hdu(table)
    table.header["EXTNAME"] = "APCORR"
    table.header["SIZEUNIT"] = "pixels"
    hdul = fits.HDUList([fits.PrimaryHDU(), table])

    apcorr_model = dm.NrsFsApcorrModel(hdul)
    yield apcorr_model
    apcorr_model.close()


@pytest.fixture()
def nirspec_fs_apcorr_file(tmp_path, nirspec_fs_apcorr):
    """
    Mock a NIRSpec FS aperture correction reference file.

    Returns
    -------
    str
        Path to the reference file.
    """
    filename = str(tmp_path / "nirspec_fs_apcorr.fits")
    nirspec_fs_apcorr.save(filename)
    return filename


@pytest.fixture()
def psf_reference():
    """
    Mock a flat spectral PSF model.

    Yields
    ------
    SpecPsfModel
        The mock model.
    """
    psf_model = dm.SpecPsfModel()
    psf_model.data = np.ones((50, 50), dtype=float)
    psf_model.wave = np.linspace(0, 10, 50)
    psf_model.meta.psf.subpix = 1.0
    psf_model.meta.psf.center_col = 25
    psf_model.meta.psf.center_row = 25
    yield psf_model
    psf_model.close()


@pytest.fixture()
def psf_reference_file(tmp_path, psf_reference):
    """
    Mock a spectral PSF reference file.

    Returns
    -------
    str
        Path to the reference file.
    """
    filename = str(tmp_path / "psf_reference.fits")
    psf_reference.save(filename)
    return filename


@pytest.fixture()
def psf_reference_with_source():
    """
    Mock a spectral PSF model with a simple PSF structure.

    Yields
    ------
    SpecPsfModel
        The mock model.
    """
    psf_model = dm.SpecPsfModel()
    psf_model.data = np.full((50, 50), 1e-6)
    psf_model.data[:, 24:27] += 1.0

    psf_model.wave = np.linspace(0, 10, 50)
    psf_model.meta.psf.subpix = 1.0
    psf_model.meta.psf.center_col = 25
    psf_model.meta.psf.center_row = 25
    yield psf_model
    psf_model.close()


@pytest.fixture()
def psf_reference_file_with_source(tmp_path, psf_reference_with_source):
    """
    Mock a spectral PSF reference file with a simple PSF structure.

    Returns
    -------
    str
        Path to the reference file.
    """
    filename = str(tmp_path / "psf_reference_with_source.fits")
    psf_reference_with_source.save(filename)
    return filename


@pytest.fixture()
def simple_profile():
    """
    Create a simple profile with a central region marked for extraction.

    Returns
    -------
    ndarray of float
        The profile array, with shape (50, 50).
    """
    profile = np.zeros((50, 50), dtype=np.float32)
    profile[20:30, :] = 1.0
    return profile


@pytest.fixture()
def background_profile():
    """
    Create a simple profile with two regions marked for background.

    Returns
    -------
    ndarray of float
        The profile array, with shape (50, 50).
    """
    profile = np.zeros((50, 50), dtype=np.float32)
    profile[:10, :] = 1.0
    profile[40:, :] = 1.0
    return profile


@pytest.fixture()
def nod_profile():
    """
    Create a simple profile with an off-center region marked for a positive source.

    Returns
    -------
    ndarray of float
        The profile array, with shape (50, 50).
    """
    profile = np.zeros((50, 50), dtype=np.float32)
    profile[10:20, :] = 1.0 / 10
    return profile


@pytest.fixture()
def negative_nod_profile():
    """
    Create a simple profile with an off-center region marked for a negative source.

    Returns
    -------
    ndarray of float
        The profile array, with shape (50, 50).
    """
    profile = np.zeros((50, 50), dtype=np.float32)
    profile[30:40, :] = -1.0 / 10
    return profile
