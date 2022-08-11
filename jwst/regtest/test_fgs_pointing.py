"""Test suite for ensuring correct FGS pointing"""
import logging
from pathlib import Path
import pytest

# Only run if `pysiaf` is installed.
pytest.importorskip('pysiaf')

from numpy import array                             # noqa: E402
from numpy import allclose, isclose                 # noqa: E402

from jwst.datamodels import Level1bModel            # noqa: E402
from jwst.lib import engdb_mast                     # noqa: E402
from jwst.lib import engdb_tools                    # noqa: E402
from jwst.lib import set_telescope_pointing as stp  # noqa: E402
from jwst.lib.file_utils import pushdir             # noqa: E402
from jwst.lib.siafdb import SiafDb                  # noqa: E402
from jwst.lib.tests.engdb_mock import EngDB_Mocker  # noqa: E402

# Set logging for the module to be tested.
logger = logging.getLogger(stp.__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)

# Get database paths.
DATA_PATH = Path(stp.__file__).parent / 'tests' / 'data'
SIAF_PATH = DATA_PATH / 'xml_data_siafxml' / 'SIAFXML'

# All the FGS GUIDER examples. Generated from proposal JW01029
FGS_ROOT_PATH = DATA_PATH / 'fgs_exposures'
FGS_PATHS = list(FGS_ROOT_PATH.glob('*.fits'))
WCS_ATTRIBUTES = ['crpix1', 'crpix2', 'crval1', 'crval2', 'pc_matrix']

FGS_TRUTHS = {
    'fgs_acq2-guider2': {
        'crpix1': 14.819122599592902, 'crpix2': 17.146828363188888,
        'crval1': 79.96718069128876, 'crval2': -69.46731727531973,
        'pc_matrix': [-0.9999870595413809, 0.005087312628765689, 0.005087312628765689, 0.9999870595413809]
    },
    'fgs_acq1-guider1': {
        'crpix1': 87.61603886804687, 'crpix2': 63.774429604665784,
        'crval1': 81.83366240613603, 'crval2': -69.9550484271187,
        'pc_matrix': [-0.9997653641994049, -0.02166140686177685, -0.02166140686177685, 0.9997653641994049]
    },
    'fgs_id-image-guider1': {
        'crpix1': 216.99326282453353, 'crpix2': 1078.9913232886374,
        'crval1': 81.83366240613603, 'crval2': -69.9550484271187,
        'pc_matrix': [-0.9997653641994049, -0.02166140686177685, -0.02166140686177685, 0.9997653641994049]
    },
    'fgs_track-guider1': {
        'crpix1': 15.69677730515923, 'crpix2': 15.056990925049945,
        'crval1': 81.83366240613603, 'crval2': -69.9550484271187,
        'pc_matrix': [-0.9997653641994049, -0.02166140686177685, -0.02166140686177685, 0.9997653641994049]
    },
    'fgs_id-stack-guider1': {
        'crpix1': 216.99326282453353, 'crpix2': 1078.9913232886374,
        'crval1': 81.83366240613603, 'crval2': -69.9550484271187,
        'pc_matrix': [-0.9997653641994049, -0.02166140686177685, -0.02166140686177685, 0.9997653641994049]
    },
    'fgs_fineguide-guider2': {
        'crpix1': 2.528106143013133, 'crpix2': 5.1871953296345055,
        'crval1': 79.96718069128876, 'crval2': -69.46731727531973,
        'pc_matrix': [-0.9999870595413809, 0.005087312628765689, 0.005087312628765689, 0.9999870595413809]
    },
    'fgs_id-stack-guider2': {
        'crpix1': 334.0005319984152, 'crpix2': 1046.0004786236082,
        'crval1': 79.96718069128876, 'crval2': -69.46731727531973,
        'pc_matrix': [-0.9999870595413809, 0.005087312628765689, 0.005087312628765689, 0.9999870595413809]
    },
    'fgs_fineguide-guider1': {
        'crpix1': 3.3286305599267507, 'crpix2': 3.32740499197007,
        'crval1': 81.83366240613603, 'crval2': -69.9550484271187,
        'pc_matrix': [-0.9997653641994049, -0.02166140686177685, -0.02166140686177685, 0.9997653641994049]
    },
    'fgs_track-guider2': {
        'crpix1': 15.285383252442102, 'crpix2': 17.29819212708435,
        'crval1': 79.96718069128876, 'crval2': -69.46731727531973,
        'pc_matrix': [-0.9999870595413809, 0.005087312628765689, 0.005087312628765689, 0.9999870595413809]
    },
    'fgs_id-image-guider2': {
        'crpix1': 334.0005319984152, 'crpix2': 1046.0004786236082,
        'crval1': 79.96718069128876, 'crval2': -69.46731727531973,
        'pc_matrix': [-0.9999870595413809, 0.005087312628765689, 0.005087312628765689, 0.9999870595413809]
    },
    'fgs_acq1-guider2': {
        'crpix1': 57.74390675801692, 'crpix2': 80.82385247510342,
        'crval1': 79.96718069128876, 'crval2': -69.46731727531973,
        'pc_matrix': [-0.9999870595413809, 0.005087312628765689, 0.005087312628765689, 0.9999870595413809]
    },
    'fgs_acq2-guider1': {
        'crpix1': 14.776268791950542, 'crpix2': 14.90005565778631,
        'crval1': 81.83366240613603, 'crval2': -69.9550484271187,
        'pc_matrix': [-0.9997653641994049, -0.02166140686177685, -0.02166140686177685, 0.9997653641994049]
    }
}

# Time frame
OBSSTART = '2022-05-22T18:07:34.000'
OBSEND = '2022-05-22T18:07:46.087'

# Define minimal model meta structure
META_ATTRIBUTES = ['pc1_1', 'pc1_2', 'pc2_1', 'pc2_2', 'crpix1', 'crpix2', 'crval1', 'crval2']
META_FGS1 = {
    'wcs': {
        'meta': {
            'aperture': {
                'name': 'FGS1_FULL',
            },
            'exposure': {
                'type': 'FGS_ACQ1',
            },
            'guidestar': {
                'gs_ra': 45.1234,
                'gs_dec': -45.1234,
            },
            'instrument': {
                'detector': 'GUIDER1',
            },
            'observation': {
                'date_beg': OBSSTART,
                'date_end': OBSEND,
                'date': '2017-01-01',
            },
        }
    },
    'expected': {
        'pc1_1': -0.9997653641994049,
        'pc1_2': -0.02166140686177685,
        'pc2_1': -0.02166140686177685,
        'pc2_2': 0.9997653641994049,
        'crpix1': 87.61603886804687,
        'crpix2': 63.774429604665784,
        'crval1': 45.1234,
        'crval2': -45.1234,
    }
}
META_FGS2 = {
    'wcs': {
        'meta': {
            'aperture': {
                'name': 'FGS2_FULL',
            },
            'exposure': {
                'type': 'FGS_ACQ1',
            },
            'guidestar': {
                'gs_ra': 45.1234,
                'gs_dec': -45.1234,
            },
            'instrument': {
                'detector': 'GUIDER2',
            },
            'observation': {
                'date_beg': OBSSTART,
                'date_end': OBSEND,
                'date': '2017-01-01',
            },
        }
    },
    'expected': {
        'pc1_1': -0.9999870595413809,
        'pc1_2': 0.005087312628765689,
        'pc2_1': 0.005087312628765689,
        'pc2_2': 0.9999870595413809,
        'crpix1': 86.80602146475417,
        'crpix2': 172.52832537903396,
        'crval1': 45.1234,
        'crval2': -45.1234,
    }
}


@pytest.mark.parametrize('attr', WCS_ATTRIBUTES)
def test_wcs_calc_guiding(get_guider_wcs, attr):
    """Test the WCS calculation"""
    exp_type, detector, wcs = get_guider_wcs

    assert allclose(wcs[attr], FGS_TRUTHS[f'{exp_type}-{detector}'][attr])


def get_guider_wcs_id(path):
    """Generate ids for get_guider_wcs"""
    id = path.stem
    return id


@pytest.fixture(params=FGS_PATHS, ids=get_guider_wcs_id, scope='module')
def get_guider_wcs(request, multi_engdb):
    siaf_db = SiafDb(SIAF_PATH)
    engdb = multi_engdb
    with Level1bModel(request.param) as model:
        exp_type = model.meta.exposure.type.lower()
        detector = model.meta.instrument.detector.lower()
        t_pars = stp.t_pars_from_model(model, siaf_db=siaf_db, engdb_url=engdb.base_url)
        wcs = stp.calc_wcs_guiding(model, t_pars)

    wcs = {key: value for key, value in zip(WCS_ATTRIBUTES, wcs)}
    return exp_type, detector, wcs


@pytest.mark.parametrize('attr', META_ATTRIBUTES)
def test_update_wcs(update_wcs, attr):
    """Test that meta info gets updated"""
    model, expected = update_wcs
    wcsinfo = model.meta.wcsinfo.instance

    assert isclose(wcsinfo[attr], expected[attr], atol=1e-15)


# ---------
# Utilities
# ---------
@pytest.fixture(scope='module')
def update_wcs(multi_engdb, make_level1b):
    """Update the model wcs info"""
    engdb = multi_engdb
    model, expected = make_level1b

    stp.update_wcs(model, engdb_url=engdb.base_url)

    return model, expected


@pytest.fixture(scope='module')
def engdb_jw01029(request, rtdata_module):
    """Setup the test engineering database"""
    rtdata = rtdata_module
    if not request.config.getoption('--bigdata'):
        pytest.skip('"--bigdata" not specified')

    db_path = Path('engdb_jw0109')
    db_path.mkdir()
    with pushdir(db_path):
        for path in rtdata.data_glob('fgs/pointing/engdb_jw01029'):
            rtdata.get_data(path)
    with EngDB_Mocker(db_path=db_path):
        engdb = engdb_tools.ENGDB_Service(base_url='http://localhost')
        yield engdb


@pytest.fixture(scope='module')
def mast():
    """Use the Mast database."""
    try:
        engdb = engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f'Live MAST Engineering Service not available: {exception}')
    yield engdb


@pytest.fixture(params=['engdb_jw01029', 'mast'], scope='module')
def multi_engdb(request):
    """Allow a test to use multiple database fixtures"""
    return request.getfixturevalue(request.param)


@pytest.fixture(scope='module', params=[META_FGS1, META_FGS2])
def make_level1b(request):
    model_meta = request.param
    data = array([1.])
    data.shape = (1, 1, 1, 1)
    model = Level1bModel(data)
    model.update(model_meta['wcs'])
    return model, model_meta['expected']
