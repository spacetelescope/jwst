"""set_telescope_pointing unit test using data from Mini-OTE-1 data

Mini-OTE-1 executed proposal jw00697, observation 13. As of 2021-11, this test
is considered the gold standard in pointing, where all known issues with the
pointing chain, from APT, through the OTB, through DMS processing, have been
resolved.

The test below uses the meta information from those exposures, and an offline
version of the engineering database telemetry, to compare current code
to the expected values from the Mini-OTE_1 test.

The full set of data covers 420 exposures. However, only a subset of exposures
are actually tested.
"""
import pytest

from pathlib import Path

import asdf
import numpy as np

import jwst.datamodels as dm
from jwst.lib import siafdb
from jwst.lib.file_utils import pushdir
import jwst.lib.set_telescope_pointing as stp
from jwst.lib.tests import engdb_mock

# Requires pysiaf
pytest.importorskip('pysiaf')

EXPS = ['jw00697013001_03101_00001_nrcblong_uncal',
        'jw00697013004_03103_00001_nrcb1_uncal',
        'jw00697013007_03101_00001_nrca4_uncal',
        'jw00697013009_03105_00001_nrcb1_uncal',
        'jw00697013002_03105_00001_nrcb3_uncal',
        'jw00697013008_03101_00001_nrcb2_uncal',
        'jw00697013008_03103_00001_nrcb3_uncal',
        'jw00697013002_03107_00001_nrcblong_uncal',
        'jw00697013007_03105_00001_nrcalong_uncal',
        'jw00697013003_03101_00001_nrcb1_uncal',
        'jw00697013006_03103_00001_nrca3_uncal']


@pytest.mark.bigdata
def test_aperture_ra(calc_wcs):
    """Compare aperture WCS information"""
    model, _,wcsinfo, _, _ = calc_wcs
    assert np.isclose(model.meta.wcsinfo.ra_ref, wcsinfo.ra)


@pytest.mark.bigdata
def test_aperture_dec(calc_wcs):
    """Compare aperture WCS information"""
    model, _,wcsinfo, _, _ = calc_wcs
    assert np.isclose(model.meta.wcsinfo.dec_ref, wcsinfo.dec)


@pytest.mark.bigdata
def test_aperture_pa(calc_wcs):
    """Compare aperture WCS information"""
    model, _,wcsinfo, _, _ = calc_wcs
    assert np.isclose(model.meta.aperture.position_angle, wcsinfo.pa)


@pytest.mark.bigdata
def test_v1_ra(calc_wcs):
    """Compare v1 WCS information"""
    model, _,_, vinfo, _ = calc_wcs
    assert np.isclose(model.meta.pointing.ra_v1, vinfo.ra)


@pytest.mark.bigdata
def test_v1_dec(calc_wcs):
    """Compare v1 WCS information"""
    model, _,_, vinfo, _ = calc_wcs
    assert np.isclose(model.meta.pointing.dec_v1, vinfo.dec)


@pytest.mark.bigdata
def test_v1_pa(calc_wcs):
    """Compare v1 WCS information"""
    model, _,_, vinfo, _ = calc_wcs
    assert np.isclose(model.meta.pointing.pa_v3, vinfo.pa)


@pytest.mark.bigdata
def test_roll_ref(calc_wcs):
    """Compare roll ref"""
    model, t_pars, wcsinfo, _, _ = calc_wcs
    roll_ref = stp.pa_to_roll_ref(wcsinfo.pa, t_pars.siaf)
    assert np.isclose(model.meta.wcsinfo.roll_ref, roll_ref)


@pytest.fixture(scope='module', params=EXPS)
def calc_wcs(databases, request):
    """Calculate the V1 and aperture """
    exposure = request.param
    siaf_db, metas = databases

    with dm.ImageModel((10, 10)) as model:
        model.update({'meta': metas[exposure]})

        # Calculate the pointing information
        t_pars = stp.t_pars_from_model(model, siaf_db=siaf_db, engdb_url='http://localhost')
        t_pars.update_pointing()
        wcsinfo, vinfo, transforms = stp.calc_wcs(t_pars)

        return model, t_pars, wcsinfo, vinfo, transforms


@pytest.fixture(scope='module')
def databases(rtdata_module):
    """Create the necessary databases needed to run pointing code

    The SIAF database needs to be open explicitly and requires `pysiaf`.

    The engineering database is handled as a context.

    Returns
    -------
    siaf_db, metas : `set_telescope_pointing.SiafDb`, dict
        Returns the tuple of the siaf database and all exposures meta information.
    """

    # Pin the PRD. Not testing changes in PRD
    siaf_db = siafdb.SiafDb(prd='PRDOPSSOC-055')

    # Get the exposures' meta information
    metas_path = rtdata_module.get_data('pointing/jw00697013_metas.asdf')
    with asdf.open(metas_path) as metas_asdf:
        metas = metas_asdf['metas']

    # Setup the engineering database
    engdb_path = Path('engdb')
    engdb_path.mkdir()
    with pushdir(engdb_path):
        paths = rtdata_module.data_glob('pointing/jw00697013_engdb')
        for path in paths:
            rtdata_module.get_data(path)

    # Pass on the database info.
    with engdb_mock.EngDB_Mocker(db_path=engdb_path):
        yield siaf_db, metas
