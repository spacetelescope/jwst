import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from stdatamodels.jwst.datamodels import ImageModel

from jwst.stpipe import query_step_status
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d.extract_2d_step import Extract2dStep
from jwst.master_background import MasterBackgroundMosStep
from jwst.master_background import nirspec_utils


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
    phdu.header['EFFEXPTM'] = 1.0
    phdu.header['DURATION'] = 1.0
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
    # Four shutters open on MSA, the first three are background slits only
    # with the fourth having a source
    shutter_data = {
        'slitlet_id': [12, 13, 14, 100],
        'msa_metadata_id': [1, 1, 1, 1],
        'shutter_quadrant': [4, 4, 4, 4],
        'shutter_row': [10, 10, 10, 10],
        'shutter_column': [22, 23, 24, 25],
        'source_id': [0, 0, 0, 2],
        'background': ['Y', 'Y', 'Y', 'N'],
        'shutter_state': ['OPEN', 'OPEN', 'OPEN', 'OPEN'],
        'estimated_source_in_shutter_x': [np.nan, np.nan, np.nan, 0.18283921],
        'estimated_source_in_shutter_y': [np.nan, np.nan, np.nan, 0.31907734],
        'dither_point_index': [1, 1, 1, 1],
        'primary_source': ['N', 'N', 'N', 'Y'],
        'fixed_slit': ['NONE', 'NONE', 'NONE', 'NONE']}

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
def nirspec_msa_metfl(tmp_path):
    hdul = create_msa_hdul()
    filename = str(tmp_path / 'test_msa_01.fits')
    hdul.writeto(filename, overwrite=True)
    hdul.close()
    return filename

@pytest.fixture
def nirspec_msa_extracted2d(nirspec_msa_rate, nirspec_msa_metfl):
    model = ImageModel(nirspec_msa_rate)
    model = AssignWcsStep.call(model)
    model = Extract2dStep.call(model)
    return model


def test_master_background_mos(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d

    result = MasterBackgroundMosStep.call(model)

    # Check that the master_background_mos step was run
    assert query_step_status(result, "master_background") == 'COMPLETE'

    finite = np.isfinite(result.slits[-1].data) & np.isfinite(model.slits[-1].data)
    sci_orig = model.slits[-1].data[finite]
    sci_bkgsub = result.slits[-1].data[finite]

    # Check that a background was subtracted from the science data
    assert not np.allclose(sci_orig, sci_bkgsub)

    model.close()
    result.close()


def test_create_background_from_multislit(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d

    # Insert a outliers into one of the background spectra
    nypix = len(model.slits[0].data)
    nxpix = len(model.slits[0].data)
    model.slits[0].data[nypix//2,nxpix//2-1:nxpix//2+1] = 10

    # First check that we can make a master background from the inputs

    # Check that with sigma_clip=None, the outlier is retained
    master_background, _ = nirspec_utils.create_background_from_multislit(
        model, sigma_clip=None)
    assert np.any(master_background.spec[0].spec_table['surf_bright'] > 1)

    # Confirm that using a median_filter will filter out the outlier
    master_background, _ = nirspec_utils.create_background_from_multislit(
        model, median_kernel=4)
    assert np.allclose(master_background.spec[0].spec_table['surf_bright'], 1)

    # Confirm that using a sigma clipping when combining background spectra
    # removes the outlier
    master_background, _ = nirspec_utils.create_background_from_multislit(
        model, sigma_clip=3)
    assert np.allclose(master_background.spec[0].spec_table['surf_bright'], 1)

    model.close()

def test_map_to_science_slits(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d

    master_background, _ = nirspec_utils.create_background_from_multislit(
        model)

    # Check that the master background is expanded to the shape of the input slits
    mb_multislit = nirspec_utils.map_to_science_slits(model, master_background)
    assert mb_multislit.slits[0].data.shape == model.slits[0].data.shape

    # background should be all ones, but won't be expanded to populate the whole
    # 2D array, so any non-zero background pixels should have a value of 1
    slit_data = mb_multislit.slits[0].data
    nonzero = slit_data != 0
    assert np.allclose(slit_data[nonzero], 1)

    model.close()

def test_apply_master_background(nirspec_msa_extracted2d):
    model = nirspec_msa_extracted2d

    master_background, _ = nirspec_utils.create_background_from_multislit(
    model)
    mb_multislit = nirspec_utils.map_to_science_slits(model, master_background)

    result = nirspec_utils.apply_master_background(model, mb_multislit, inverse=False)

    # where the background is applied to the science it should be 0 and elsewhere 1
    sci_data_orig = model.slits[-1].data
    sci_data_bkgsub = result.slits[-1].data
    diff = sci_data_orig - sci_data_bkgsub
    assert np.any(diff != 0)
    assert np.allclose(diff[diff != 0], 1)

    # Check inverse application
    result = nirspec_utils.apply_master_background(model, mb_multislit, inverse=True)

    # where the background is applied to the science it should be 0 and elsewhere 1
    sci_data_orig = model.slits[-1].data
    sci_data_bkgsub = result.slits[-1].data
    diff = sci_data_orig - sci_data_bkgsub

    # Background subtraction was inverted so the differences will be -1
    assert np.any(diff != 0)
    assert np.allclose(diff[diff != 0], -1)

    model.close()
    result.close()


