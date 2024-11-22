import json
import logging
import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import polynomial

from jwst.extract_1d import extract as ex
from jwst.tests.helpers import LogWatcher


@pytest.fixture
def log_watcher(monkeypatch):
    # Set a log watcher to check for a log message at any level
    # in the extract_1d.extract module
    watcher = LogWatcher('')
    logger = logging.getLogger('jwst.extract_1d.extract')
    for level in ['debug', 'info', 'warning', 'error']:
        monkeypatch.setattr(logger, level, watcher)
    return watcher


@pytest.fixture()
def extract1d_ref_dict():
    apertures = [{'id': 'slit1'},
                 {'id': 'slit2', 'region_type': 'other'},
                 {'id': 'slit3', 'xstart': 10, 'xstop': 20, 'ystart': 10, 'ystop': 20},
                 {'id': 'slit4', 'bkg_coeff': [[10], [20]]},
                 {'id': 'slit5', 'bkg_coeff': None},
                 {'id': 'slit6', 'use_source_posn': True},
                 ]
    ref_dict = {'apertures': apertures}
    return ref_dict


@pytest.fixture()
def extract1d_ref_file(tmp_path, extract1d_ref_dict):
    filename = str(tmp_path / 'extract1d_ref.json')
    with open(filename, 'w') as fh:
        json.dump(extract1d_ref_dict, fh)
    return filename


@pytest.fixture()
def extract_defaults():
    default = {'bkg_coeff': None,
               'bkg_fit': None,
               'bkg_order': 0,
               'extract_width': None,
               'extraction_type': 'box',
               'independent_var': 'pixel',
               'match': 'exact match',
               'position_correction': 0,
               'smoothing_length': 0,
               'spectral_order': 1,
               'src_coeff': None,
               'subtract_background': False,
               'use_source_posn': False,
               'xstart': 0,
               'xstop': 49,
               'ystart': 0,
               'ystop': 49}
    return default


def test_read_extract1d_ref(extract1d_ref_dict, extract1d_ref_file):
    ref_dict = ex.read_extract1d_ref(extract1d_ref_file)
    assert ref_dict == extract1d_ref_dict


def test_read_extract1d_ref_bad_json(tmp_path):
    filename = str(tmp_path / 'bad_ref.json')
    with open(filename, 'w') as fh:
        fh.write('apertures: [bad,]\n')

    with pytest.raises(RuntimeError, match='Invalid JSON extract1d reference'):
        ex.read_extract1d_ref(filename)


def test_read_extract1d_ref_bad_type(tmp_path):
    filename = str(tmp_path / 'bad_ref.fits')
    with open(filename, 'w') as fh:
        fh.write('bad file\n')

    with pytest.raises(RuntimeError, match='must be JSON'):
        ex.read_extract1d_ref(filename)


def test_read_extract1d_ref_na():
    ref_dict = ex.read_extract1d_ref('N/A')
    assert ref_dict is None


def test_read_apcorr_ref():
    apcorr_model = ex.read_apcorr_ref(None, 'MIR_LRS-FIXEDSLIT')
    assert isinstance(apcorr_model, dm.MirLrsApcorrModel)


def test_get_extract_parameters_default(
        mock_nirspec_fs_one_slit, extract1d_ref_dict, extract_defaults):
    input_model = mock_nirspec_fs_one_slit

    # match a bare entry
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit1', 1, input_model.meta)

    # returned value has defaults except that use_source_posn
    # is switched on for NRS_FIXEDSLIT
    expected = extract_defaults
    expected['use_source_posn'] = True

    assert params == expected


def test_get_extract_parameters_na(mock_nirspec_fs_one_slit, extract_defaults):
    input_model = mock_nirspec_fs_one_slit

    # no reference input: defaults returned
    params = ex.get_extract_parameters(None, input_model, 'slit1', 1, input_model.meta)
    assert params == extract_defaults


@pytest.mark.parametrize('bgsub', [None, True])
@pytest.mark.parametrize('bgfit', ['poly', 'mean', None])
def test_get_extract_parameters_background(
        mock_nirspec_fs_one_slit, extract1d_ref_dict, bgsub, bgfit):
    input_model = mock_nirspec_fs_one_slit

    # match a slit with background defined
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit4', 1, input_model.meta,
        subtract_background=bgsub, bkg_fit=bgfit)

    # returned value has background switched on
    assert params['subtract_background'] is True
    assert params['bkg_coeff'] is not None
    assert params['bkg_fit'] == 'poly'
    assert params['bkg_order'] == 0


def test_get_extract_parameters_bg_ignored(mock_nirspec_fs_one_slit, extract1d_ref_dict):
    input_model = mock_nirspec_fs_one_slit

    # match a slit with background defined
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit4', 1, input_model.meta,
        subtract_background=False)

    # background parameters are ignored
    assert params['subtract_background'] is False
    assert params['bkg_fit'] is None


@pytest.mark.parametrize('slit', ['slit2', 'no_match'])
def test_get_extract_parameters_no_match(
        mock_nirspec_fs_one_slit, extract1d_ref_dict, slit):
    input_model = mock_nirspec_fs_one_slit

    # no slit with an appropriate region_type matched
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, slit, 1, input_model.meta)
    assert params == {'match': ex.NO_MATCH}


def test_get_extract_parameters_source_posn_exptype(
        mock_nirspec_bots, extract1d_ref_dict, extract_defaults):
    input_model = mock_nirspec_bots

    # match a bare entry
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit1', 1, input_model.meta,
        use_source_posn=None)

    # use_source_posn is switched off for NRS_BRIGHTOBJ
    assert params['use_source_posn'] is False


def test_get_extract_parameters_source_posn_from_ref(
        mock_nirspec_bots, extract1d_ref_dict, extract_defaults):
    input_model = mock_nirspec_bots

    # match an entry that explicity sets use_source_posn
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit6', 1, input_model.meta,
        use_source_posn=None)

    # returned value has use_source_posn switched off by default
    # for NRS_BRIGHTOBJ, but ref file overrides
    assert params['use_source_posn'] is True


def test_get_extract_parameters_smoothing(
        mock_nirspec_fs_one_slit, extract1d_ref_dict, extract_defaults):
    input_model = mock_nirspec_fs_one_slit

    # match an entry that explicity sets use_source_posn
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit1', 1, input_model.meta,
        smoothing_length=3)

    # returned value has input smoothing length
    assert params['smoothing_length'] == 3


def test_log_params(extract_defaults, log_watcher):
    log_watcher.message = 'Extraction parameters'

    # Defaults don't have dispaxis assigned yet - parameters are not logged
    ex.log_initial_parameters(extract_defaults)
    assert not log_watcher.seen

    # Add dispaxis: parameters are now logged
    extract_defaults['dispaxis'] = 1
    ex.log_initial_parameters(extract_defaults)
    log_watcher.assert_seen()


def test_create_poly():
    coeff = [1, 2, 3]
    poly = ex.create_poly(coeff)
    assert isinstance(poly, polynomial.Polynomial1D)
    assert poly.degree == 2
    assert poly(2) == 1 + 2 * 2 + 3 * 2**2


def test_create_poly_empty():
    coeff = []
    assert ex.create_poly(coeff) is None


def test_populate_time_keywords(mock_nirspec_bots, mock_10_spec):
    ex.populate_time_keywords(mock_nirspec_bots, mock_10_spec)

    # time keywords now added to output spectra
    for i, spec in enumerate(mock_10_spec.spec):
        assert spec.int_num == i + 1
        assert spec.start_time_mjd == mock_nirspec_bots.int_times['int_start_MJD_UTC'][i]
        assert spec.end_tdb == mock_nirspec_bots.int_times['int_end_BJD_TDB'][i]


def test_populate_time_keywords_no_table(mock_nirspec_fs_one_slit, mock_one_spec):
    ex.populate_time_keywords(mock_nirspec_fs_one_slit, mock_one_spec)

    # only int_num is added to spec
    assert mock_one_spec.spec[0].int_num == 1


def test_populate_time_keywords_multislit(mock_nirspec_mos, mock_10_spec):
    mock_nirspec_mos.meta.exposure.nints = 10
    ex.populate_time_keywords(mock_nirspec_mos, mock_10_spec)

    # no int_times - only int_num is added to spec
    # It is set to 1 for all spectra - no integrations in multislit data.
    assert mock_10_spec.spec[0].int_num == 1
    assert mock_10_spec.spec[9].int_num == 1


def test_populate_time_keywords_multislit_table(
        mock_nirspec_mos, mock_nirspec_bots, mock_10_spec, log_watcher):
    mock_nirspec_mos.meta.exposure.nints = 10
    mock_nirspec_mos.int_times = mock_nirspec_bots.int_times

    log_watcher.message = 'Not using INT_TIMES table'
    ex.populate_time_keywords(mock_nirspec_mos, mock_10_spec)
    log_watcher.assert_seen()

    # int_times present but not used - no update
    assert mock_10_spec.spec[0].int_num is None


def test_populate_time_keywords_averaged(
        mock_nirspec_fs_one_slit, mock_nirspec_bots, mock_10_spec, log_watcher):
    mock_nirspec_fs_one_slit.meta.exposure.nints = 10
    mock_nirspec_fs_one_slit.int_times = mock_nirspec_bots.int_times

    log_watcher.message = 'Not using INT_TIMES table'
    ex.populate_time_keywords(mock_nirspec_fs_one_slit, mock_10_spec)
    log_watcher.assert_seen()

    # int_times not used - no update
    assert mock_10_spec.spec[0].int_num is None


def test_populate_time_keywords_mismatched_table(mock_nirspec_bots, mock_10_spec, log_watcher):
    # mock 20 integrations - table has 10
    mock_nirspec_bots.data = np.vstack([mock_nirspec_bots.data, mock_nirspec_bots.data])
    log_watcher.message = 'Not using INT_TIMES table'
    ex.populate_time_keywords(mock_nirspec_bots, mock_10_spec)
    log_watcher.assert_seen()

    # int_times not used - no update
    assert mock_10_spec.spec[0].int_num is None


def test_populate_time_keywords_missing_ints(mock_nirspec_bots, mock_10_spec, log_watcher):
    mock_nirspec_bots.meta.exposure.integration_start = 20
    log_watcher.message = 'does not include rows'
    ex.populate_time_keywords(mock_nirspec_bots, mock_10_spec)
    log_watcher.assert_seen()

    # int_times not used - no update
    assert mock_10_spec.spec[0].int_num is None


def test_populate_time_keywords_ifu_table(
        mock_miri_ifu, mock_nirspec_bots, mock_10_spec, log_watcher):
    mock_miri_ifu.meta.exposure.nints = 10
    mock_miri_ifu.int_times = mock_nirspec_bots.int_times

    log_watcher.message = 'ignored for IFU'
    ex.populate_time_keywords(mock_miri_ifu, mock_10_spec)
    log_watcher.assert_seen()

    # int_times present but not used - no update
    assert mock_10_spec.spec[0].int_num is None


def test_populate_time_keywords_mismatched_spec(
        mock_nirspec_bots, mock_one_spec, log_watcher):
    log_watcher.message = "Don't understand n_output_spec"
    ex.populate_time_keywords(mock_nirspec_bots, mock_one_spec)
    log_watcher.assert_seen()
