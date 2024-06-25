import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from stdatamodels.jwst.datamodels import ImageModel, CubeModel, MultiSlitModel, SlitModel

from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d.extract_2d_step import Extract2dStep


# WCS keywords, borrowed from NIRCam grism tests
WCS_KEYS = {'wcsaxes': 2, 'ra_ref': 53.1490299775, 'dec_ref': -27.8168745624,
            'v2_ref': 86.103458, 'v3_ref': -493.227512, 'roll_ref': 45.04234459270135,
            'v3i_yang': 0.0, 'vparity': -1,
            'crpix1': 1024.5, 'crpix2': 1024.5,
            'crval1': 53.1490299775, 'crval2': -27.8168745624,
            'cdelt1': 1.81661111111111e-05, 'cdelt2': 1.8303611111111e-05,
            'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
            'pc1_1': -0.707688557183348, 'pc1_2': 0.7065245261360363,
            'pc2_1': 0.7065245261360363, 'pc2_2': 1.75306861111111e-05,
            'cunit1': 'deg', 'cunit2': 'deg'}


def create_nirspec_hdul(detector='NRS1', grating='G395M', filter_name='F290LP',
                        exptype='NRS_MSASPEC', subarray='FULL', slit=None, nint=1,
                        wcskeys=None):
    if wcskeys is None:
        wcskeys = WCS_KEYS.copy()

    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['TELESCOP'] = 'JWST'
    phdu.header['INSTRUME'] = 'NIRSPEC'
    phdu.header['DETECTOR'] = detector
    phdu.header['FILTER'] = filter_name
    phdu.header['GRATING'] = grating
    phdu.header['PROGRAM'] = '01234'
    phdu.header['TIME-OBS'] = '8:59:37'
    phdu.header['DATE-OBS'] = '2023-01-05'
    phdu.header['EXP_TYPE'] = exptype
    phdu.header['PATT_NUM'] = 1
    phdu.header['SUBARRAY'] = subarray
    phdu.header['XOFFSET'] = 0.0
    phdu.header['YOFFSET'] = 0.0
    if subarray == 'SUBS200A1':
        phdu.header['SUBSIZE1'] = 2048
        phdu.header['SUBSIZE2'] = 64
        phdu.header['SUBSTRT1'] = 1
        phdu.header['SUBSTRT2'] = 1041
    elif subarray == 'SUB2048':
        phdu.header['SUBSIZE1'] = 2048
        phdu.header['SUBSIZE2'] = 32
        phdu.header['SUBSTRT1'] = 1
        phdu.header['SUBSTRT2'] = 946
    else:
        phdu.header['SUBSIZE1'] = 2048
        phdu.header['SUBSIZE2'] = 2048
        phdu.header['SUBSTRT1'] = 1
        phdu.header['SUBSTRT2'] = 1

    if exptype == 'NRS_MSASPEC':
        phdu.header['MSAMETID'] = 1
        phdu.header['MSAMETFL'] = 'test_msa_01.fits'

    if slit is not None:
        phdu.header['FXD_SLIT'] = slit
        phdu.header['APERNAME'] = f'NRS_{slit}_SLIT'

    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    scihdu.header.update(wcskeys)
    if nint > 1:
        scihdu.data = np.ones((nint, phdu.header['SUBSIZE2'], phdu.header['SUBSIZE1']))
    else:
        scihdu.data = np.ones((phdu.header['SUBSIZE2'], phdu.header['SUBSIZE1']))

    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def create_msa_hdul():
    # Two point sources, one in MSA, one fixed slit.
    # Source locations for the fixed slit are placeholders, not realistic.
    shutter_data = {
        'slitlet_id': [12, 12, 12, 100],
        'msa_metadata_id': [1, 1, 1, 1],
        'shutter_quadrant': [4, 4, 4, 0],
        'shutter_row': [251, 251, 251, 0],
        'shutter_column': [22, 23, 24, 0],
        'source_id': [1, 1, 1, 2],
        'background': ['Y', 'N', 'Y', 'N'],
        'shutter_state': ['OPEN', 'OPEN', 'OPEN', 'OPEN'],
        'estimated_source_in_shutter_x': [np.nan, 0.18283921, np.nan, 0.5],
        'estimated_source_in_shutter_y': [np.nan, 0.31907734, np.nan, 0.5],
        'dither_point_index': [1, 1, 1, 1],
        'primary_source': ['N', 'Y', 'N', 'Y'],
        'fixed_slit': ['NONE', 'NONE', 'NONE', 'S200A1']}

    source_data = {
        'program': [95065, 95065],
        'source_id': [1, 2],
        'source_name': ['95065_1', '95065_2'],
        'alias': ['2122', '2123'],
        'ra': [53.139904, 53.15],
        'dec': [-27.805002, -27.81],
        'preimage_id': ['95065001_000', '95065001_000'],
        'stellarity': [1.0, 1.0]}

    shutter_table = Table(shutter_data)
    source_table = Table(source_data)

    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(fits.ImageHDU())
    hdul.append(fits.table_to_hdu(shutter_table))
    hdul.append(fits.table_to_hdu(source_table))
    hdul[2].name = 'SHUTTER_INFO'
    hdul[3].name = 'SOURCE_INFO'

    return hdul


@pytest.fixture
def nirspec_msa_rate(tmp_path):
    hdul = create_nirspec_hdul()
    hdul[0].header['MSAMETFL'] = str(tmp_path / 'test_msa_01.fits')
    filename = str(tmp_path / 'test_nrs_msa_rate.fits')
    hdul.writeto(filename, overwrite=True)
    hdul.close()
    return filename


@pytest.fixture
def nirspec_fs_rate(tmp_path):
    hdul = create_nirspec_hdul(
        exptype='NRS_FIXEDSLIT', subarray='SUBS200A1', slit='S200A1')
    filename = str(tmp_path / 'test_nrs_fs_rate.fits')
    hdul.writeto(filename, overwrite=True)
    hdul.close()
    return filename


@pytest.fixture
def nirspec_bots_rateints(tmp_path):
    hdul = create_nirspec_hdul(
        exptype='NRS_BRIGHTOBJ', subarray='SUB2048', slit='S1600A1', nint=3)
    filename = str(tmp_path / 'test_nrs_bots_rateints.fits')
    hdul.writeto(filename, overwrite=True)
    hdul.close()
    return filename


@pytest.fixture
def nirspec_msa_metfl(tmp_path):
    hdul = create_msa_hdul()
    filename = str(tmp_path / 'test_msa_01.fits')
    hdul.writeto(filename, overwrite=True)
    hdul.close()
    return filename


def test_extract_2d_nirspec_msa_fs(nirspec_msa_rate, nirspec_msa_metfl):
    model = ImageModel(nirspec_msa_rate)
    result = AssignWcsStep.call(model)
    result = Extract2dStep.call(result)
    assert isinstance(result, MultiSlitModel)

    # there should be 2 slits extracted: one MSA, one FS
    assert len(result.slits) == 2

    # the MSA slit has an integer name, slitlet_id matches name
    assert result.slits[0].name == '12'
    assert result.slits[0].slitlet_id == 12
    assert result.slits[0].data.shape == (31, 1355)

    # the FS slit has a string name, slitlet_id matches shutter ID
    assert result.slits[1].name == 'S200A1'
    assert result.slits[1].slitlet_id == 0
    assert result.slits[1].data.shape == (45, 1254)

    model.close()
    result.close()


def test_extract_2d_nirspec_fs(nirspec_fs_rate):
    model = ImageModel(nirspec_fs_rate)
    model_wcs = AssignWcsStep.call(model)


    result = Extract2dStep.call(model_wcs)
    assert isinstance(result, MultiSlitModel)

    # there should be 1 slit extracted: FS, S200A1
    assert len(result.slits) == 1

    # the FS slit has a string name, slitlet_id matches shutter ID
    assert result.slits[0].name == 'S200A1'
    assert result.slits[0].slitlet_id == 0
    assert result.slits[0].data.shape == (45, 1254)

    # ensure x_offset, y_offset become zero when dither information is missing
    model_wcs.meta.dither = None
    result = Extract2dStep.call(model_wcs)
    assert result.slits[0].source_xpos == 0.0
    assert result.slits[0].source_ypos == 0.0
    model_wcs.meta.dither = {"x_offset": None, "y_offset": None}
    result = Extract2dStep.call(model_wcs)
    assert result.slits[0].source_xpos == 0.0
    assert result.slits[0].source_ypos == 0.0

    model.close()
    model_wcs.close()
    result.close()


def test_extract_2d_nirspec_bots(nirspec_bots_rateints):
    model = CubeModel(nirspec_bots_rateints)
    result = AssignWcsStep.call(model)
    result = Extract2dStep.call(result)

    # output is a single slit
    assert isinstance(result, SlitModel)

    # the BOTS slit has a string name, slitlet_id matches shutter ID
    assert result.name == 'S1600A1'
    assert result.data.shape == (3, 28, 1300)

    model.close()
    result.close()
