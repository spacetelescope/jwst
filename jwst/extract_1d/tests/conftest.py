import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.io import fits
from astropy.table import Table

from jwst.assign_wcs.util import wcs_bbox_from_shape
from jwst.exp_to_source import multislit_to_container


@pytest.fixture()
def simple_wcs():
    shape = (50, 50)
    xcenter = shape[1] // 2.0

    def simple_wcs_function(x, y):
        """ Simple WCS for testing """
        crpix1 = xcenter
        crpix3 = 1.0
        cdelt1 = 0.1
        cdelt2 = 0.1
        cdelt3 = 0.01

        crval1 = 45.0
        crval2 = 45.0
        crval3 = 7.5

        wave = (x + 1 - crpix3) * cdelt3 + crval3
        ra = (x + 1 - crpix1) * cdelt1 + crval1
        dec = np.full_like(ra, crval2 + 1 * cdelt2)

        return ra, dec, wave

    # Add a bounding box
    simple_wcs_function.bounding_box = wcs_bbox_from_shape(shape)

    # Add a few expected attributes, so they can be monkeypatched as needed
    simple_wcs_function.get_transform = None
    simple_wcs_function.backward_transform = None
    simple_wcs_function.available_frames = []

    return simple_wcs_function


@pytest.fixture()
def simple_wcs_transpose():
    shape = (50, 50)
    ycenter = shape[0] // 2.0

    def simple_wcs_function(x, y):
        """ Simple WCS for testing """
        crpix2 = ycenter
        crpix3 = 1.0
        cdelt1 = 0.1
        cdelt2 = 0.1
        cdelt3 = -0.01

        crval1 = 45.0
        crval2 = 45.0
        crval3 = 7.5

        wave = (y + 1 - crpix3) * cdelt3 + crval3
        ra = (y + 1 - crpix2) * cdelt1 + crval1
        dec = np.full_like(ra, crval2 + 1 * cdelt2)

        return ra, dec, wave

    # Add a bounding box
    simple_wcs_function.bounding_box = wcs_bbox_from_shape(shape)

    # Add a few expected attributes, so they can be monkeypatched as needed
    simple_wcs_function.get_transform = None
    simple_wcs_function.backward_transform = None
    simple_wcs_function.available_frames = []

    return simple_wcs_function


@pytest.fixture()
def simple_wcs_ifu():
    shape = (10, 50, 50)
    xcenter = shape[1] // 2.0

    def simple_wcs_function(x, y, z):
        """ Simple WCS for testing """
        crpix1 = xcenter
        crpix3 = 1.0
        cdelt1 = 0.1
        cdelt2 = 0.1
        cdelt3 = 0.01

        crval1 = 45.0
        crval2 = 45.0
        crval3 = 7.5

        wave = (z + 1 - crpix3) * cdelt3 + crval3
        ra = (z + 1 - crpix1) * cdelt1 + crval1
        dec = np.full_like(ra, crval2 + 1 * cdelt2)

        return ra, dec, wave[::-1]

    simple_wcs_function.bounding_box = wcs_bbox_from_shape(shape)

    return simple_wcs_function


@pytest.fixture()
def mock_nirspec_fs_one_slit(simple_wcs):
    model = dm.SlitModel()
    model.meta.instrument.name = 'NIRSPEC'
    model.meta.instrument.detector = 'NRS1'
    model.meta.instrument.filter = 'F290LP'
    model.meta.instrument.grating = 'G395H'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.instrument.fixed_slit = 'S200A1'
    model.meta.exposure.nints = 1
    model.meta.exposure.type = 'NRS_FIXEDSLIT'
    model.meta.subarray.name = 'ALLSLITS'
    model.source_type = 'EXTENDED'

    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs

    model.data = np.arange(50 * 50, dtype=float).reshape((50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    yield model
    model.close()


@pytest.fixture()
def mock_nirspec_mos(mock_nirspec_fs_one_slit):
    model = dm.MultiSlitModel()
    model.meta.instrument.name = 'NIRSPEC'
    model.meta.instrument.detector = 'NRS1'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.type = 'NRS_MSASPEC'
    model.meta.exposure.nints = 1

    nslit = 3
    for i in range(nslit):
        slit = mock_nirspec_fs_one_slit.copy()
        slit.name = str(i + 1)
        model.slits.append(slit)

    yield model
    model.close()


@pytest.fixture()
def mock_nirspec_bots(simple_wcs):
    model = dm.CubeModel()
    model.meta.instrument.name = 'NIRSPEC'
    model.meta.instrument.detector = 'NRS1'
    model.meta.instrument.filter = 'F290LP'
    model.meta.instrument.grating = 'G395H'
    model.meta.observation.date = '2022-05-30'
    model.meta.observation.time = '01:03:16.369'
    model.meta.instrument.fixed_slit = 'S1600A1'
    model.meta.exposure.type = 'NRS_BRIGHTOBJ'
    model.meta.subarray.name = 'SUB2048'
    model.meta.exposure.nints = 10
    model.meta.visit.tsovisit = True

    model.name = 'S1600A1'
    model.meta.target.source_type = 'POINT'
    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs

    model.data = np.arange(10 * 50 * 50, dtype=float).reshape((10, 50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05

    # Add an int_times table
    integrations = [
        (1, 59729.04367729, 59729.04378181, 59729.04388632, 59729.04731706, 59729.04742158, 59729.04752609),
        (2, 59729.04389677, 59729.04400128, 59729.04410579, 59729.04753654, 59729.04764105, 59729.04774557),
        (3, 59729.04411625, 59729.04422076, 59729.04432527, 59729.04775602, 59729.04786053, 59729.04796504),
        (4, 59729.04433572, 59729.04444023, 59729.04454475, 59729.04797549, 59729.04808001, 59729.04818452),
        (5, 59729.0445552, 59729.04465971, 59729.04476422, 59729.04819497, 59729.04829948, 59729.048404),
        (6, 59729.04477467, 59729.04487918, 59729.0449837, 59729.04841445, 59729.04851896, 59729.04862347),
        (7, 59729.04499415, 59729.04509866, 59729.04520317, 59729.04863392, 59729.04873844, 59729.04884295),
        (8, 59729.04521362, 59729.04531813, 59729.04542265, 59729.0488534 , 59729.04895791, 59729.04906242),
        (9, 59729.0454331, 59729.04553761, 59729.04564212, 59729.04907288, 59729.04917739, 59729.0492819),
        (10, 59729.04565257, 59729.04575709, 59729.0458616, 59729.04929235, 59729.04939686, 59729.04950138),
    ]

    integration_table = np.array(integrations, dtype=[('integration_number', 'i4'),
                                                      ('int_start_MJD_UTC', 'f8'),
                                                      ('int_mid_MJD_UTC', 'f8'),
                                                      ('int_end_MJD_UTC', 'f8'),
                                                      ('int_start_BJD_TDB', 'f8'),
                                                      ('int_mid_BJD_TDB', 'f8'),
                                                      ('int_end_BJD_TDB', 'f8')])
    model.int_times = integration_table

    yield model
    model.close()


@pytest.fixture()
def mock_miri_lrs_fs(simple_wcs_transpose):
    model = dm.ImageModel()
    model.meta.instrument.name = 'MIRI'
    model.meta.instrument.detector = 'MIRIMAGE'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.nints = 1
    model.meta.exposure.type = 'MIR_LRS-FIXEDSLIT'
    model.meta.subarray.name = 'FULL'
    model.meta.target.source_type = 'EXTENDED'
    model.meta.dither.dithered_ra = 45.0
    model.meta.dither.dithered_ra = 45.0

    model.meta.wcsinfo.dispersion_direction = 2
    model.meta.wcs = simple_wcs_transpose

    model.data = np.arange(50 * 50, dtype=float).reshape((50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    yield model
    model.close()


@pytest.fixture()
def mock_miri_ifu(simple_wcs_ifu):
    model = dm.IFUCubeModel()
    model.meta.instrument.name = 'MIRI'
    model.meta.instrument.detector = 'MIRIFULONG'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.type = 'MIR_MRS'

    model.meta.wcsinfo.dispersion_direction = 2
    model.meta.photometry.pixelarea_steradians = 1.0
    model.meta.wcs = simple_wcs_ifu

    model.data = np.arange(10 * 50 * 50, dtype=float).reshape((10, 50, 50))
    model.var_poisson = model.data * 0.02
    model.var_rnoise = model.data * 0.02
    model.var_flat = model.data * 0.05
    model.weightmap = np.full_like(model.data, 1.0)
    yield model
    model.close()


@pytest.fixture()
def mock_niriss_wfss_l3(mock_nirspec_fs_one_slit):
    model = dm.MultiSlitModel()
    model.meta.instrument.name = 'NIRISS'
    model.meta.instrument.detector = 'NIS'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.type = 'NIS_WFSS'

    nslit = 3
    for i in range(nslit):
        slit = mock_nirspec_fs_one_slit.copy()
        slit.name = str(i + 1)
        slit.meta.exposure.type = 'NIS_WFSS'
        model.slits.append(slit)

    container = multislit_to_container([model])['0']

    yield container
    container.close()


@pytest.fixture()
def mock_niriss_soss(simple_wcs):
    model = dm.CubeModel()
    model.meta.instrument.name = 'NIRISS'
    model.meta.instrument.detector = 'NIS'
    model.meta.instrument.filter = 'CLEAR'
    model.meta.instrument.pupil_position = 245.79
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.exposure.type = 'NIS_SOSS'
    model.meta.exposure.nints = 3

    model.meta.target.source_type = 'POINT'
    model.meta.wcsinfo.dispersion_direction = 1
    model.meta.wcs = simple_wcs

    yield model
    model.close()


@pytest.fixture()
def mock_niriss_soss_256(mock_niriss_soss):
    model = mock_niriss_soss
    model.meta.subarray.name = 'SUBSTRIP256'

    shape = (3, 256, 2048)
    model.data = np.ones(shape, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.err = model.data * 0.02
    model.var_poisson = model.data * 0.001
    model.var_rnoise = model.data * 0.001
    model.var_flat = model.data * 0.001
    return model


@pytest.fixture()
def mock_niriss_soss_96(mock_niriss_soss):
    model = mock_niriss_soss
    model.meta.subarray.name = 'SUBSTRIP96'

    shape = (3, 96, 2048)
    model.data = np.ones(shape, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)
    model.err = model.data * 0.02
    model.var_poisson = model.data * 0.001
    model.var_rnoise = model.data * 0.001
    model.var_flat = model.data * 0.001
    return model


def make_spec_model(name='slit1', value=1.0):
    wavelength = np.arange(20, dtype=np.float32)
    flux = np.full(10, value)
    error = 0.05 * flux
    f_var_poisson = error ** 2
    f_var_rnoise = np.zeros_like(flux)
    f_var_flat = np.zeros_like(flux)
    surf_bright = flux / 10
    sb_error = error / 10
    sb_var_poisson = f_var_poisson / 10**2
    sb_var_rnoise = f_var_rnoise / 10**2
    sb_var_flat = f_var_flat / 10**2
    dq = np.zeros(20, dtype=np.uint32)
    background = np.zeros_like(flux)
    berror = np.zeros_like(flux)
    b_var_poisson = np.zeros_like(flux)
    b_var_rnoise = np.zeros_like(flux)
    b_var_flat = np.zeros_like(flux)
    npixels = np.full(20, 10)

    spec_dtype = dm.SpecModel().spec_table.dtype
    otab = np.array(
        list(
            zip(
                wavelength, flux, error, f_var_poisson, f_var_rnoise, f_var_flat,
                surf_bright, sb_error, sb_var_poisson, sb_var_rnoise, sb_var_flat,
                dq, background, berror, b_var_poisson, b_var_rnoise, b_var_flat,
                npixels
            ),
        ), dtype=spec_dtype
    )

    spec_model = dm.SpecModel(spec_table=otab)
    spec_model.name = name

    return spec_model


@pytest.fixture()
def mock_one_spec():
    model = dm.MultiSpecModel()
    spec_model = make_spec_model()
    model.spec.append(spec_model)

    yield model
    model.close()


@pytest.fixture()
def mock_10_spec(mock_one_spec):
    model = dm.MultiSpecModel()

    for i in range(10):
        spec_model = make_spec_model(name=f'slit{i + 1}', value=i + 1)
        model.spec.append(spec_model)

    yield model
    model.close()


@pytest.fixture()
def miri_lrs_apcorr():
    table = Table(
        {
            'subarray': ['FULL', 'SLITLESSPRISM'],
            'wavelength': [[1, 2, 3], [1, 2, 3]],
            'nelem_wl': [3, 3],
            'size': [[1, 2, 3], [1, 2, 3]],
            'nelem_size': [3, 3],
            'apcorr': np.full((2, 3, 3), 0.5),
            'apcorr_err': np.full((2, 3, 3), 0.01)
        }
    )
    table = fits.table_to_hdu(table)
    table.header['EXTNAME'] = 'APCORR'
    table.header['SIZEUNIT'] = 'pixels'
    hdul = fits.HDUList([fits.PrimaryHDU(), table])

    apcorr_model = dm.MirLrsApcorrModel(hdul)
    yield apcorr_model
    apcorr_model.close()


@pytest.fixture()
def miri_lrs_apcorr_file(tmp_path, miri_lrs_apcorr):
    filename = str(tmp_path / 'miri_lrs_apcorr.fits')
    miri_lrs_apcorr.save(filename)
    return filename


@pytest.fixture()
def nirspec_fs_apcorr():
    table = Table(
        {
            'filter': ['clear', 'f290lp'],
            'grating': ['prism', 'g395h'],
            'slit': ['S200A1', 'S200A1'],
            'wavelength': [[1, 2, 3], [1, 2, 3]],
            'nelem_wl': [3, 3],
            'size': np.full((2, 3, 3), [0.3, 0.5, 1]),
            'nelem_size': [3, 3],
            'pixphase': [[0.0, 0.25, 0.5], [0.0, 0.25, 0.5]],
            'apcorr': np.full((2, 3, 3, 3), 0.5),
            'apcorr_err': np.full((2, 3, 3, 3), 0.01)
        }
    )
    table = fits.table_to_hdu(table)
    table.header['EXTNAME'] = 'APCORR'
    table.header['SIZEUNIT'] = 'pixels'
    hdul = fits.HDUList([fits.PrimaryHDU(), table])

    apcorr_model = dm.NrsFsApcorrModel(hdul)
    yield apcorr_model
    apcorr_model.close()


@pytest.fixture()
def nirspec_fs_apcorr_file(tmp_path, nirspec_fs_apcorr):
    filename = str(tmp_path / 'nirspec_fs_apcorr.fits')
    nirspec_fs_apcorr.save(filename)
    return filename
