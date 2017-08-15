"""
Test suite for set_telescope_pointing

Notes
-----
This file has been named specifically so it is not
automatically found by py.test. This is because, to test,
a connection to the internal engineering service is needed,
which is generally not available.
"""

from __future__ import absolute_import

from backports.tempfile import TemporaryDirectory
import copy
import numpy as np
import os
import pytest
import requests_mock

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from .. import engdb_tools
from .. import set_telescope_pointing as stp

# Setup mock engineering service
GOOD_MNEMONIC = 'INRSI_GWA_Y_TILT_AVGED'
STARTTIME = Time('2014-01-03')
ENDTIME = Time('2014-01-04')
ZEROTIME_START = Time('2014-01-01')
ZEROTIME_END = Time('2014-01-02')

# Header defaults
TARG_RA = 345.0
TARG_DEC = -87.0
V2_REF = 200.0
V3_REF = -350.0
V3I_YANG = 42.0
VPARITY = -1

# Get the mock DB
db_path = os.path.join(os.path.dirname(__file__), 'data', 'engdb_mock.csv')
mock_db = Table.read(db_path)


# Some expected falues
Q_EXPECTED = np.asarray(
    [-0.36915286, 0.33763282, 0.05758533, 0.86395264]
)
J2FGS_MATRIX_EXPECTED = np.asarray(
    [-1.00444000e-03,   3.38145836e-03,   9.99993778e-01,
     9.99999496e-01,  -3.90000000e-14,   1.00444575e-03,
     3.39649146e-06,   9.99994283e-01,  -3.38145665e-03]
)
FSMCORR_EXPECTED = np.zeros((2,))
OBSTIME_EXPECTED = STARTTIME


def register_responses(mocker, response_db, starttime, endtime):
    request_url = ''.join([
        engdb_tools.ENGDB_BASE_URL,
        'Data/',
        '{mnemonic}',
        '?sTime={starttime}',
        '&eTime={endtime}'
    ])

    starttime_mil = int(starttime.unix * 1000)
    endtime_mil = int(endtime.unix * 1000)
    time_increment = (endtime_mil - starttime_mil) // len(response_db)

    response_generic = {
        'AllPoints': 1,
        'Count': 2,
        'ReqSTime': '/Date({:013d}+0000)/'.format(starttime_mil),
        'ReqETime': '/Date({:013d}+0000)/'.format(endtime_mil),
        'TlmMnemonic': None,
        'Data': [],
    }

    responses = {}
    for mnemonic in response_db.colnames:
        response = copy.deepcopy(response_generic)
        response['TlmMnemonic'] = mnemonic
        current_time = starttime_mil - time_increment
        for row in response_db:
            current_time += time_increment
            data = {}
            data['ObsTime'] = '/Date({:013d}+0000)/'.format(current_time)
            data['EUValue'] = row[mnemonic]
            response['Data'].append(data)
        mocker.get(
            request_url.format(
                mnemonic=mnemonic,
                starttime=starttime.iso,
                endtime=endtime.iso
            ),
            json=response
        )
        responses[mnemonic] = response

    return responses


@pytest.fixture
def eng_db():
    with requests_mock.Mocker() as rm:

        # Define response for aliveness
        url = ''.join([
            engdb_tools.ENGDB_BASE_URL,
            engdb_tools.ENGDB_METADATA
        ])
        rm.get(url, text='Success')

        # Define good responses
        good_responses = register_responses(
            rm,
            mock_db[1:2],
            STARTTIME,
            ENDTIME
        )

        # Define with zeros in the first row
        zero_responses = register_responses(
            rm,
            mock_db,
            ZEROTIME_START,
            ENDTIME
        )

        edb = engdb_tools.ENGDB_Service()

        # Test for good responses.
        for mnemonic in mock_db.colnames:
            r = edb.get_records(mnemonic, STARTTIME, ENDTIME)
            assert r == good_responses[mnemonic]

        # Test for zeros.
        for mnemonic in mock_db.colnames:
            r = edb.get_records(mnemonic, ZEROTIME_START, ENDTIME)
            assert r == zero_responses[mnemonic]

        yield edb


@pytest.fixture
def fits_file():
    hdu = fits.PrimaryHDU()
    header = hdu.header
    header['EXPSTART'] = STARTTIME.mjd
    header['EXPEND'] = ENDTIME.mjd
    header['TARG_RA'] = TARG_RA
    header['TARG_DEC'] = TARG_DEC
    header['V2_REF'] = V2_REF
    header['V3_REF'] = V3_REF
    header['V3I_YANG'] = V3I_YANG
    header['VPARITY'] = VPARITY
    hdul = fits.HDUList([hdu])

    with TemporaryDirectory() as path:
        file_path = os.path.join(path, 'fits.fits')
        hdul.writeto(file_path)
        yield file_path


def test_get_pointing_fail():
    with pytest.raises(Exception):
        q, j2fgs_matrix, fmscorr, obstime = stp.get_pointing(47892.0, 48256.0)


def test_get_pointing(eng_db):
        (q,
         j2fgs_matrix,
         fsmcorr,
         obstime) = stp.get_pointing(STARTTIME, ENDTIME)
        assert np.isclose(q, Q_EXPECTED).all()
        assert np.isclose(j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
        assert np.isclose(fsmcorr, FSMCORR_EXPECTED).all()
        assert obstime == STARTTIME


def test_get_pointing_list(eng_db):
        results = stp.get_pointing(STARTTIME, ENDTIME, result_type='all')
        assert isinstance(results, list)
        assert len(results) > 0
        assert np.isclose(results[0].q, Q_EXPECTED).all()
        assert np.isclose(results[0].j2fgs_matrix, J2FGS_MATRIX_EXPECTED).all()
        assert np.isclose(results[0].fsmcorr, FSMCORR_EXPECTED).all()
        assert results[0].obstime == STARTTIME


def test_get_pointing_with_zeros(eng_db):
    (q,
     j2fgs_matrix,
     fsmcorr,
     obstime) = stp.get_pointing(ZEROTIME_START, ENDTIME)
    assert j2fgs_matrix.any()
    (q_desired,
     j2fgs_matrix_desired,
     fsmcorr_desired,
     obstime) = stp.get_pointing(STARTTIME, ENDTIME)
    assert np.array_equal(q, q_desired)
    assert np.array_equal(j2fgs_matrix, j2fgs_matrix_desired)
    assert np.array_equal(fsmcorr, fsmcorr_desired)


def test_add_wcs_default(fits_file):
    try:
        stp.add_wcs(fits_file)
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            '\nException={}'.format(e)
        )

    hdul = fits.open(fits_file)
    header = hdul[0].header
    assert header['RA_V1'] == TARG_RA
    assert header['DEC_V1'] == TARG_DEC
    assert header['PA_V3'] == 0.
    assert header['CRVAL1'] == TARG_RA
    assert header['CRVAL2'] == TARG_DEC
    assert np.isclose(header['PC1_1'], -1.0)
    assert np.isclose(header['PC1_2'], 0.0)
    assert np.isclose(header['PC2_1'], 0.0)
    assert np.isclose(header['PC2_2'], 1.0)
    assert header['RA_REF'] == TARG_RA
    assert header['DEC_REF'] == TARG_DEC
    assert np.isclose(header['ROLL_REF'], 358.9045979379)
    assert header['WCSAXES'] == 0.


def test_add_wcs_fsmcorr_v1(fits_file):
    """Test with default value using FSM original correction"""
    try:
        stp.add_wcs(fits_file, fsmcorr_version='v1')
    except Exception as e:
        pytest.skip(
            'Live ENGDB service is not accessible.'
            '\nException={}'.format(e)
        )

    hdul = fits.open(fits_file)
    header = hdul[0].header
    assert header['RA_V1'] == TARG_RA
    assert header['DEC_V1'] == TARG_DEC
    assert header['PA_V3'] == 0.
    assert header['CRVAL1'] == TARG_RA
    assert header['CRVAL2'] == TARG_DEC
    assert np.isclose(header['PC1_1'], -1.0)
    assert np.isclose(header['PC1_2'], 0.0)
    assert np.isclose(header['PC2_1'], 0.0)
    assert np.isclose(header['PC2_2'], 1.0)
    assert header['RA_REF'] == TARG_RA
    assert header['DEC_REF'] == TARG_DEC
    assert np.isclose(header['ROLL_REF'], 358.9045979379)
    assert header['WCSAXES'] == 0.


def test_add_wcs_with_db(eng_db, fits_file):
    """Test using the database"""
    stp.add_wcs(fits_file)

    hdul = fits.open(fits_file)
    header = hdul[0].header
    assert np.isclose(header['RA_V1'], 348.9278669)
    assert np.isclose(header['DEC_V1'], -38.749239)
    assert np.isclose(header['PA_V3'], 50.1767077)
    assert np.isclose(header['CRVAL1'], 348.8776709)
    assert np.isclose(header['CRVAL2'], -38.854159)
    assert np.isclose(header['PC1_1'], 0.0385309)
    assert np.isclose(header['PC1_2'], 0.9992574)
    assert np.isclose(header['PC2_1'], 0.9992574)
    assert np.isclose(header['PC2_2'], -0.0385309)
    assert np.isclose(header['RA_REF'], 348.8776709)
    assert np.isclose(header['DEC_REF'], -38.854159)
    assert np.isclose(header['ROLL_REF'], 50.20832726650)
    assert header['WCSAXES'] == 0.


def test_add_wcs_with_db_fsmcorr_v1(eng_db, fits_file):
    """Test using the database with original FSM correction"""
    stp.add_wcs(fits_file, fsmcorr_version='v1')

    hdul = fits.open(fits_file)
    header = hdul[0].header
    assert np.isclose(header['RA_V1'], 348.9278669)
    assert np.isclose(header['DEC_V1'], -38.749239)
    assert np.isclose(header['PA_V3'], 50.1767077)
    assert np.isclose(header['CRVAL1'], 348.8776709)
    assert np.isclose(header['CRVAL2'], -38.854159)
    assert np.isclose(header['PC1_1'], 0.0385309)
    assert np.isclose(header['PC1_2'], 0.9992574)
    assert np.isclose(header['PC2_1'], 0.9992574)
    assert np.isclose(header['PC2_2'], -0.0385309)
    assert np.isclose(header['RA_REF'], 348.8776709)
    assert np.isclose(header['DEC_REF'], -38.854159)
    assert np.isclose(header['ROLL_REF'], 50.20832726650)
    assert header['WCSAXES'] == 0.
