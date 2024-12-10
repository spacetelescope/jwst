import json
import logging
import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import polynomial

from jwst.datamodels import ModelContainer
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
                 {'id': 'slit7', 'spectral_order': 20},
                 {'id': 'S200A1'}
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


@pytest.fixture()
def simple_profile():
    profile = np.zeros((50, 50), dtype=np.float32)
    profile[20:30, :] = 1.0
    return profile


@pytest.fixture()
def background_profile():
    profile = np.zeros((50, 50), dtype=np.float32)
    profile[:10, :] = 1.0
    profile[40:, :] = 1.0
    return profile


@pytest.fixture()
def create_extraction_inputs(mock_nirspec_fs_one_slit, extract1d_ref_dict):
    input_model = mock_nirspec_fs_one_slit
    slit = None
    output_model = dm.MultiSpecModel()
    ref_dict = extract1d_ref_dict
    slitname = 'S200A1'
    sp_order = 1
    exp_type = 'NRS_FIXEDSLIT'
    yield [input_model, slit, output_model, ref_dict,
           slitname, sp_order, exp_type]
    output_model.close()


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


@pytest.mark.parametrize('length', [3, 4, 2.8, 3.5])
def test_get_extract_parameters_smoothing(
        mock_nirspec_fs_one_slit, extract1d_ref_dict,
        extract_defaults, length):
    input_model = mock_nirspec_fs_one_slit

    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit1', 1, input_model.meta,
        smoothing_length=length)

    # returned value has input smoothing length, rounded to an
    # odd integer if necessary
    assert params['smoothing_length'] == 3


@pytest.mark.parametrize('length', [-1, 1, 2, 1.3])
def test_get_extract_parameters_smoothing_bad_value(
        mock_nirspec_fs_one_slit, extract1d_ref_dict,
        extract_defaults, length):
    input_model = mock_nirspec_fs_one_slit

    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, 'slit1', 1, input_model.meta,
        smoothing_length=length)

    # returned value has smoothing length 0
    assert params['smoothing_length'] == 0


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


def test_get_spectral_order(mock_nirspec_fs_one_slit):
    slit = mock_nirspec_fs_one_slit

    slit.meta.wcsinfo.spectral_order = 2
    assert ex.get_spectral_order(slit) == 2

    slit.meta.wcsinfo.spectral_order = None
    assert ex.get_spectral_order(slit) == 1

    slit.meta.wcsinfo = None
    assert ex.get_spectral_order(slit) == 1

    del slit.meta.wcsinfo
    assert ex.get_spectral_order(slit) == 1


def test_is_prism_false_conditions(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    assert ex.is_prism(model) is False

    model.meta.instrument.filter = None
    model.meta.instrument.grating = None
    assert ex.is_prism(model) is False

    model.meta.instrument.name = None
    assert ex.is_prism(model) is False


def test_is_prism_nirspec(mock_nirspec_fs_one_slit):
    mock_nirspec_fs_one_slit.meta.instrument.grating = 'PRISM'
    assert ex.is_prism(mock_nirspec_fs_one_slit) is True


def test_is_prism_miri(mock_miri_ifu):
    mock_miri_ifu.meta.instrument.filter = 'P750L'
    assert ex.is_prism(mock_miri_ifu) is True

def test_copy_keyword_info(mock_nirspec_fs_one_slit, mock_one_spec):
    expected = {'slitlet_id': 2,
                'source_id': 3,
                'source_name': '4',
                'source_alias': '5',
                'source_type': 'POINT',
                'stellarity': 0.5,
                'source_xpos': -0.5,
                'source_ypos': -0.5,
                'source_ra': 10.0,
                'source_dec': 10.0,
                'shutter_state': 'x'}
    for key, value in expected.items():
        setattr(mock_nirspec_fs_one_slit, key, value)
        assert not hasattr(mock_one_spec, key)

    ex.copy_keyword_info(mock_nirspec_fs_one_slit, 'slit_name', mock_one_spec)
    assert mock_one_spec.name == 'slit_name'
    for key, value in expected.items():
        assert getattr(mock_one_spec, key) == value


@pytest.mark.parametrize("partial", [True, False])
@pytest.mark.parametrize("lower,upper",
                         [(0, 19), (-1, 21),
                          (np.full(10, 0.0), np.full(10, 19.0)),
                          (np.linspace(-1, 0, 10), np.linspace(19, 20, 10)),])
def test_set_weights_from_limits_whole_array(lower, upper, partial):
    shape = (20, 10)
    profile = np.zeros(shape, dtype=float)
    yidx, _ = np.mgrid[:shape[0], :shape[1]]

    ex._set_weight_from_limits(profile, yidx, lower, upper, allow_partial=partial)
    assert np.all(profile == 1.0)


@pytest.mark.parametrize("lower,upper",
                         [(10, 12), (9.5, 12.5),
                          (np.linspace(9.5, 10, 10), np.linspace(12, 12.5, 10)),])
def test_set_weights_from_limits_whole_pixel(lower, upper):
    shape = (20, 10)
    profile = np.zeros(shape, dtype=np.float32)
    yidx, _ = np.mgrid[:shape[0], :shape[1]]

    ex._set_weight_from_limits(profile, yidx, lower, upper, allow_partial=False)
    assert np.all(profile[10:13] == 1.0)


@pytest.mark.parametrize("lower,upper",
                         [(9.5, 12.5),
                          (np.linspace(9.5, 10, 10), np.linspace(12, 12.5, 10)),])
def test_set_weights_from_limits_partial_pixel(lower, upper):
    shape = (20, 10)
    profile = np.zeros(shape, dtype=np.float32)
    yidx, _ = np.mgrid[:shape[0], :shape[1]]

    ex._set_weight_from_limits(profile, yidx, lower, upper, allow_partial=True)
    assert np.allclose(profile[10:13], 1.0)
    assert np.allclose(profile[9], 10 - lower)
    assert np.allclose(profile[13], upper - 12)


def test_set_weights_from_limits_overlap():
    shape = (20, 10)
    profile = np.zeros(shape, dtype=np.float32)
    yidx, _ = np.mgrid[:shape[0], :shape[1]]

    # Set an aperture with partial pixel edges
    ex._set_weight_from_limits(profile, yidx, 9.5, 10.5, allow_partial=True)
    assert np.allclose(profile[9], 0.5)
    assert np.allclose(profile[11], 0.5)
    assert np.allclose(profile[12], 0.0)

    # Set an overlapping region in the same profile
    ex._set_weight_from_limits(profile, yidx, 9.8, 11.5, allow_partial=True)

    # Higher weight from previous profile remains
    assert np.allclose(profile[9], 0.5)

    # Previous partial pixel is now fully included
    assert np.allclose(profile[11], 1.0)

    # New partial weight set on upper limit
    assert np.allclose(profile[12], 0.5)


def test_box_profile_horizontal(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = 1

    # Exclude 2 pixels at left and right
    params['xstart'] = 1.5
    params['xstop'] = 7.5

    # Exclude 2 pixels at top and bottom, set half pixel weights
    # for another pixel at top and bottom edges
    params['ystart'] = 2.5
    params['ystop'] = 6.5
    profile = ex.box_profile(shape, extract_defaults, wl_array)

    # ystart/stop sets partial weights, xstart/stop sets whole pixels only
    assert np.all(profile[2:3, 3:8] == 0.5)
    assert np.all(profile[7:8, 3:8] == 0.5)
    assert np.all(profile[3:7, 3:8] == 1.0)
    assert np.all(profile[:2] == 0.0)
    assert np.all(profile[8:] == 0.0)
    assert np.all(profile[:, :2] == 0.0)
    assert np.all(profile[:, 8:] == 0.0)


def test_box_profile_vertical(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = 2

    # Exclude 2 pixels at "left" and "right" - in transposed aperture
    params['ystart'] = 1.5
    params['ystop'] = 7.5

    # Exclude 2 pixels at "top" and "bottom", set half pixel weights
    # for another pixel at top and bottom edges
    params['xstart'] = 2.5
    params['xstop'] = 6.5
    profile = ex.box_profile(shape, extract_defaults, wl_array)

    # xstart/stop sets partial weights, ystart/stop sets whole pixels only
    assert np.all(profile[3:8, 2:3] == 0.5)
    assert np.all(profile[3:8, 7:8] == 0.5)
    assert np.all(profile[3:8, 3:7] == 1.0)

    assert np.all(profile[:2] == 0.0)
    assert np.all(profile[8:] == 0.0)
    assert np.all(profile[:, :2] == 0.0)
    assert np.all(profile[:, 8:] == 0.0)


@pytest.mark.parametrize('dispaxis', [1, 2])
def test_box_profile_bkg_coeff(extract_defaults, dispaxis):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = dispaxis

    # the definition for bkg_coeff is half a pixel off from start/stop definitions -
    # this will set the equivalent of start/stop 0-2, 7-9 -
    # 3 pixels at top and bottom of the array
    params['bkg_coeff'] = [[-0.5], [2.5], [6.5], [9.5]]

    profile, lower, upper = (
        ex.box_profile(shape, extract_defaults, wl_array,
                       coefficients='bkg_coeff', return_limits=True))
    if dispaxis == 2:
        profile = profile.T

    assert np.all(profile[:3] == 1.0)
    assert np.all(profile[7:] == 1.0)
    assert np.all(profile[3:7] == 0.0)
    assert lower == 0.0
    assert upper == 9.0


def test_box_profile_bkg_coeff_median(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = 1
    params['bkg_fit'] = 'median'

    # Attempt to set partial pixels at middle edges
    params['bkg_coeff'] = [[-0.5], [3.0], [6.0], [9.5]]

    profile, lower, upper = (
        ex.box_profile(shape, extract_defaults, wl_array,
                       coefficients='bkg_coeff', return_limits=True))

    # partial pixels are not allowed for fit type median - the profile is
    # set for whole pixels only
    assert np.all(profile[:3] == 1.0)
    assert np.all(profile[7:] == 1.0)
    assert np.all(profile[3:7] == 0.0)
    assert lower == 0.0
    assert upper == 9.0


@pytest.mark.parametrize('swap_order', [False, True])
def test_box_profile_bkg_coeff_poly(extract_defaults, swap_order):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = 1
    params['bkg_fit'] = 'poly'

    # Attempt to set partial pixels at middle edges
    if swap_order:
        # upper region first - this should make no difference.
        params['bkg_coeff'] = [[6.0], [9.5], [-0.5], [3.0]]
    else:
        params['bkg_coeff'] = [[-0.5], [3.0], [6.0], [9.5]]

    profile, lower, upper = (
        ex.box_profile(shape, extract_defaults, wl_array,
                       coefficients='bkg_coeff', return_limits=True))

    # partial pixels are allowed for fit type poly
    assert np.all(profile[:3] == 1.0)
    assert np.all(profile[7:] == 1.0)
    assert np.all(profile[3] == 0.5)
    assert np.all(profile[6] == 0.5)
    assert np.all(profile[4:6] == 0.0)
    assert lower == 0.0
    assert upper == 9.0


@pytest.mark.parametrize('independent_var', ['pixel', 'wavelength'])
def test_box_profile_src_coeff_constant(extract_defaults, independent_var):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = 1
    params['independent_var'] = independent_var

    # the definition for src_coeff is half a pixel off from start/stop definitions -
    # this will set the equivalent of start/stop 3-6, excluding
    # 3 pixels at top and bottom of the array
    params['src_coeff'] = [[2.5], [6.5]]

    profile, lower, upper = (
        ex.box_profile(shape, extract_defaults, wl_array,
                       coefficients='src_coeff', return_limits=True))
    assert np.all(profile[3:7] == 1.0)
    assert np.all(profile[:3] == 0.0)
    assert np.all(profile[7:] == 0.0)
    assert lower == 3.0
    assert upper == 6.0


@pytest.mark.parametrize('independent_var', ['pixel', 'wavelength'])
def test_box_profile_src_coeff_linear(extract_defaults, independent_var):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = 1
    params['independent_var'] = independent_var

    if independent_var == 'wavelength':
        slope = 1 / (wl_array[0, 1] - wl_array[0, 0])
        start = -0.5 - wl_array[0, 0] * slope
        stop = start + 4
    else:
        slope = 1.0
        start = -0.5
        stop = 3.5

    # Set linearly increasing upper and lower edges,
    # starting at the bottom of the array, with a width of 4 pixels
    params['src_coeff'] = [[start, slope], [stop, slope]]
    expected = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                [1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                [0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                [0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 1, 1, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 0],
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1]]

    profile, lower, upper = (
        ex.box_profile(shape, extract_defaults, wl_array,
                       coefficients='src_coeff', return_limits=True))
    assert np.allclose(profile, expected)

    # upper and lower limits are averages
    assert np.isclose(lower, 4.5)
    assert np.isclose(upper, 7.5)


def test_box_profile_mismatched_coeff(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = 1

    # set some mismatched coefficients limits
    params['src_coeff'] = [[2.5], [6.5], [9.5]]

    with pytest.raises(RuntimeError, match='must contain alternating lists'):
        ex.box_profile(shape, extract_defaults, wl_array)


@pytest.mark.parametrize('dispaxis', [1, 2])
def test_box_profile_from_width(extract_defaults, dispaxis):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params['dispaxis'] = dispaxis

    if dispaxis == 1:
        # Set ystart and ystop to center on pixel 4
        params['ystart'] = 2.0
        params['ystop'] = 6.0
    else:
        # Set xstart and xstop to center on pixel 4
        params['xstart'] = 2.0
        params['xstop'] = 6.0

    # Set a width to 6 pixels
    params['extract_width'] = 6.0

    profile = ex.box_profile(shape, extract_defaults, wl_array)
    if dispaxis == 2:
        profile = profile.T

    # Aperture is centered at pixel 4, to start at 1.5, end at 6.5
    assert np.all(profile[2:7] == 1.0)
    assert np.all(profile[1] == 0.5)
    assert np.all(profile[7] == 0.5)
    assert np.all(profile[0] == 0.0)
    assert np.all(profile[8:] == 0.0)


@pytest.mark.parametrize('middle', [None, 7])
@pytest.mark.parametrize('dispaxis', [1, 2])
def test_aperture_center(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    profile[1:4] = 1.0
    if dispaxis != 1:
        profile = profile.T
    slit_center, spec_center = ex.aperture_center(
        profile, dispaxis=dispaxis, middle_pix=middle)
    assert slit_center == 2.0
    if middle is None:
        assert spec_center == 4.5
    else:
        assert spec_center == middle


@pytest.mark.parametrize('middle', [None, 7])
@pytest.mark.parametrize('dispaxis', [1, 2])
def test_aperture_center_zero_weight(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    slit_center, spec_center = ex.aperture_center(
        profile, dispaxis=dispaxis, middle_pix=middle)
    assert slit_center == 4.5
    if middle is None:
        assert spec_center == 4.5
    else:
        assert spec_center == middle


@pytest.mark.parametrize('middle', [None, 7])
@pytest.mark.parametrize('dispaxis', [1, 2])
def test_aperture_center_variable_weight_by_slit(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    profile[1:4] = np.arange(10)
    if dispaxis != 1:
        profile = profile.T
    slit_center, spec_center = ex.aperture_center(
        profile, dispaxis=dispaxis, middle_pix=middle)
    assert slit_center == 2.0
    if middle is None:
        assert np.isclose(spec_center, 6.3333333)
    else:
        assert spec_center == middle


@pytest.mark.parametrize('middle', [None, 7])
@pytest.mark.parametrize('dispaxis', [1, 2])
def test_aperture_center_variable_weight_by_spec(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    profile[:, 1:4] = np.arange(10)[:, None]
    if dispaxis != 1:
        profile = profile.T
    slit_center, spec_center = ex.aperture_center(
        profile, dispaxis=dispaxis, middle_pix=middle)
    if middle is None:
        assert np.isclose(slit_center, 6.3333333)
        assert np.isclose(spec_center, 2.0)
    else:
        assert np.isclose(slit_center, 4.5)
        assert spec_center == middle


@pytest.mark.parametrize('resampled', [True, False])
@pytest.mark.parametrize('is_slit', [True, False])
@pytest.mark.parametrize('missing_bbox', [True, False])
def test_location_from_wcs_nirspec(
        monkeypatch, mock_nirspec_fs_one_slit, resampled, is_slit, missing_bbox):
    model = mock_nirspec_fs_one_slit

    # monkey patch in a transform for the wcs
    def slit2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 0.0, 1.0
        return return_one

    monkeypatch.setattr(model.meta.wcs, 'get_transform', slit2det)

    if not resampled:
        # also mock available frames, so it looks like unresampled cal data
        monkeypatch.setattr(model.meta.wcs, 'available_frames', ['gwa'])

    if missing_bbox:
        # also mock a missing bounding box - should have same results
        # for the test data
        monkeypatch.setattr(model.meta.wcs, 'bounding_box', None)

    if is_slit:
        middle, middle_wl, location = ex.location_from_wcs(model, model)
    else:
        middle, middle_wl, location = ex.location_from_wcs(model, None)

    # middle pixel is center of dispersion axis
    assert middle == int((model.data.shape[1] - 1) / 2)

    # middle wavelength is the wavelength at that point, from the mock wcs
    assert np.isclose(middle_wl, 7.74)

    # location is 1.0 - from the mocked transform function
    assert location == 1.0


@pytest.mark.parametrize('is_slit', [True, False])
def test_location_from_wcs_miri(monkeypatch, mock_miri_lrs_fs, is_slit):
    model = mock_miri_lrs_fs

    # monkey patch in a transform for the wcs
    def radec2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 1.0, 0.0
        return return_one

    monkeypatch.setattr(model.meta.wcs, 'backward_transform', radec2det())

    # Get the slit center from the WCS
    if is_slit:
        middle, middle_wl, location = ex.location_from_wcs(model, model)
    else:
        middle, middle_wl, location = ex.location_from_wcs(model, None)

    # middle pixel is center of dispersion axis
    assert middle == int((model.data.shape[0] - 1) / 2)

    # middle wavelength is the wavelength at that point, from the mock wcs
    assert np.isclose(middle_wl, 7.26)

    # location is 1.0 - from the mocked transform function
    assert location == 1.0


def test_location_from_wcs_missing_data(mock_miri_lrs_fs, log_watcher):
    # model is missing WCS information - None values are returned
    log_watcher.message = "Dithered pointing location not found"
    result = ex.location_from_wcs(mock_miri_lrs_fs, None)
    assert result == (None, None, None)
    log_watcher.assert_seen()


def test_location_from_wcs_wrong_exptype(mock_niriss_soss, log_watcher):
    # model is not a handled exposure type
    log_watcher.message = "Source position cannot be found for EXP_TYPE"
    result = ex.location_from_wcs(mock_niriss_soss, None)
    assert result == (None, None, None)
    log_watcher.assert_seen()


def test_location_from_wcs_bad_location(
        monkeypatch, mock_nirspec_fs_one_slit, log_watcher):
    model = mock_nirspec_fs_one_slit

    # monkey patch in a transform for the wcs
    def slit2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 0.0, np.nan
        return return_one

    monkeypatch.setattr(model.meta.wcs, 'get_transform', slit2det)

    # WCS transform returns NaN for the location
    log_watcher.message = "Source position could not be determined"
    result = ex.location_from_wcs(model, None)
    assert result == (None, None, None)
    log_watcher.assert_seen()


def test_location_from_wcs_location_out_of_range(
        monkeypatch, mock_nirspec_fs_one_slit, log_watcher):
    model = mock_nirspec_fs_one_slit

    # monkey patch in a transform for the wcs
    def slit2det(*args, **kwargs):
        def return_one(*args, **kwargs):
            return 0.0, 2000
        return return_one

    monkeypatch.setattr(model.meta.wcs, 'get_transform', slit2det)

    # WCS transform a value outside the bounding box
    log_watcher.message = "outside the bounding box"
    result = ex.location_from_wcs(model, None)
    assert result == (None, None, None)
    log_watcher.assert_seen()


def test_shift_by_source_location_horizontal(extract_defaults):
    location = 12.5
    nominal_location = 15.0
    offset = location - nominal_location

    extract_params = extract_defaults.copy()
    extract_params['dispaxis'] = 1

    ex.shift_by_source_location(location, nominal_location, extract_params)
    assert extract_params['xstart'] == extract_defaults['xstart']
    assert extract_params['xstop'] == extract_defaults['xstop']
    assert extract_params['ystart'] == extract_defaults['ystart'] + offset
    assert extract_params['ystop'] == extract_defaults['ystop'] + offset


def test_shift_by_source_location_vertical(extract_defaults):
    location = 12.5
    nominal_location = 15.0
    offset = location - nominal_location

    extract_params = extract_defaults.copy()
    extract_params['dispaxis'] = 2

    ex.shift_by_source_location(location, nominal_location, extract_params)
    assert extract_params['xstart'] == extract_defaults['xstart'] + offset
    assert extract_params['xstop'] == extract_defaults['xstop'] + offset
    assert extract_params['ystart'] == extract_defaults['ystart']
    assert extract_params['ystop'] == extract_defaults['ystop']


def test_shift_by_source_location_coeff(extract_defaults):
    location = 6.5
    nominal_location = 4.0
    offset = location - nominal_location

    extract_params = extract_defaults.copy()
    extract_params['dispaxis'] = 1
    extract_params['src_coeff'] = [[2.5, 1.0], [6.5, 1.0]]
    extract_params['bkg_coeff'] = [[-0.5], [3.0], [6.0], [9.5]]

    ex.shift_by_source_location(location, nominal_location, extract_params)
    assert extract_params['src_coeff'] == [[2.5 + offset, 1.0], [6.5 + offset, 1.0]]
    assert extract_params['bkg_coeff'] == [[-0.5 + offset], [3.0 + offset],
                                           [6.0 + offset], [9.5 + offset]]


@pytest.mark.parametrize('is_slit', [True, False])
def test_define_aperture_nirspec(mock_nirspec_fs_one_slit, extract_defaults, is_slit):
    model = mock_nirspec_fs_one_slit
    extract_defaults['dispaxis'] = 1
    if is_slit:
        slit = model
    else:
        slit = None
    exptype = 'NRS_FIXEDSLIT'
    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    ra, dec, wavelength, profile, bg_profile, limits = result
    assert np.isclose(ra, 45.05)
    assert np.isclose(dec, 45.1)
    assert wavelength.shape == (model.data.shape[1],)
    assert profile.shape == model.data.shape

    # Default profile is the full array
    assert np.all(profile == 1.0)
    assert limits == (0, model.data.shape[0] - 1, 0, model.data.shape[1] - 1)

    # Default bg profile is None
    assert bg_profile is None


@pytest.mark.parametrize('is_slit', [True, False])
def test_define_aperture_miri(mock_miri_lrs_fs, extract_defaults, is_slit):
    model = mock_miri_lrs_fs
    extract_defaults['dispaxis'] = 2
    if is_slit:
        slit = model
    else:
        slit = None
    exptype = 'MIR_LRS-FIXEDSLIT'
    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    ra, dec, wavelength, profile, bg_profile, limits = result
    assert np.isclose(ra, 45.05)
    assert np.isclose(dec, 45.1)
    assert wavelength.shape == (model.data.shape[1],)
    assert profile.shape == model.data.shape

    # Default profile is the full array
    assert np.all(profile == 1.0)
    assert limits == (0, model.data.shape[0] - 1, 0, model.data.shape[1] - 1)

    # Default bg profile is None
    assert bg_profile is None


def test_define_aperture_with_bg(mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults['dispaxis'] = 1
    slit = None
    exptype = 'NRS_FIXEDSLIT'

    extract_defaults['subtract_background'] = True
    extract_defaults['bkg_coeff'] = [[-0.5], [2.5]]

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    bg_profile = result[-2]

    # Bg profile has 1s in the first 3 rows
    assert bg_profile.shape == model.data.shape
    assert np.all(bg_profile[:3] == 1.0)
    assert np.all(bg_profile[3:] == 0.0)


def test_define_aperture_empty_aperture(mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults['dispaxis'] = 1
    slit = None
    exptype = 'NRS_FIXEDSLIT'

    # Set the extraction limits out of range
    extract_defaults['ystart'] = 2000
    extract_defaults['ystop'] = 3000

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    _, _, _, profile, _, limits = result

    assert np.all(profile == 0.0)
    assert limits == (2000, 3000, None, None)


def test_define_aperture_bad_wcs(monkeypatch, mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults['dispaxis'] = 1
    slit = None
    exptype = 'NRS_FIXEDSLIT'

    # Set a wavelength so wcs is not called to retrieve it
    model.wavelength = np.empty_like(model.data)
    model.wavelength[:] = np.linspace(3, 5, model.data.shape[1])

    # mock a bad wcs
    def return_nan(*args):
        return np.nan, np.nan, np.nan

    monkeypatch.setattr(model.meta, 'wcs', return_nan)

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    ra, dec = result[:2]

    # RA and Dec returned are none
    assert ra is None
    assert dec is None


def test_define_aperture_use_source(monkeypatch, mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults['dispaxis'] = 1
    slit = None
    exptype = 'NRS_FIXEDSLIT'

    # mock the source location function
    def mock_source_location(*args):
        return 24, 7.74, 9.5

    monkeypatch.setattr(ex, 'location_from_wcs', mock_source_location)

    # set parameters to extract a 6 pixel aperture, centered on source location
    extract_defaults['use_source_posn'] = True
    extract_defaults['extract_width'] = 6.0

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    _, _, _, profile, _, limits = result

    assert np.all(profile[:7] == 0.0)
    assert np.all(profile[7:13] == 1.0)
    assert np.all(profile[13:] == 0.0)


def test_extract_one_slit_horizontal(mock_nirspec_fs_one_slit, extract_defaults,
                                     simple_profile, background_profile):
    # update parameters to subtract background
    extract_defaults['dispaxis'] = 1
    extract_defaults['subtract_background'] = True
    extract_defaults['bkg_fit'] = 'poly'
    extract_defaults['bkg_order'] = 1

    # set a source in the profile region
    mock_nirspec_fs_one_slit.data[simple_profile != 0] += 1.0

    result = ex.extract_one_slit(mock_nirspec_fs_one_slit, -1, simple_profile,
                                 background_profile, extract_defaults)

    for data in result[:-1]:
        assert np.all(data > 0)
        assert data.shape == (mock_nirspec_fs_one_slit.data.shape[1],)

    # residuals from the 2D scene model should be zero - this simple case
    # is exactly modeled with a box profile
    scene_model = result[-1]
    assert scene_model.shape == mock_nirspec_fs_one_slit.data.shape
    assert np.allclose(np.abs(mock_nirspec_fs_one_slit.data - scene_model), 0)

    # flux should be 1.0 * npixels
    flux = result[0]
    npixels = result[-2]
    assert np.allclose(flux, npixels)

    # npixels is sum of profile
    assert np.all(npixels == np.sum(simple_profile, axis=0))


def test_extract_one_slit_vertical(mock_miri_lrs_fs, extract_defaults,
                                   simple_profile, background_profile):
    model = mock_miri_lrs_fs
    profile = simple_profile.T
    profile_bg = background_profile.T

    # update parameters to subtract background
    extract_defaults['dispaxis'] = 2
    extract_defaults['subtract_background'] = True
    extract_defaults['bkg_fit'] = 'poly'
    extract_defaults['bkg_order'] = 1

    # set a source in the profile region
    model.data[profile != 0] += 1.0

    result = ex.extract_one_slit(model, -1, profile, profile_bg, extract_defaults)

    for data in result[:-1]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[0],)

    # residuals from the 2D scene model should be zero - this simple case
    # is exactly modeled with a box profile
    scene_model = result[-1]
    assert scene_model.shape == model.data.shape
    assert np.allclose(np.abs(model.data - scene_model), 0)

    # flux should be 1.0 * npixels
    flux = result[0]
    npixels = result[-2]
    assert np.allclose(flux, npixels)

    # npixels is sum of profile
    assert np.all(npixels == np.sum(profile, axis=1))


def test_extract_one_slit_vertical_no_bg(mock_miri_lrs_fs, extract_defaults,
                                         simple_profile):
    model = mock_miri_lrs_fs
    profile = simple_profile.T
    extract_defaults['dispaxis'] = 2

    result = ex.extract_one_slit(model, -1, profile, None, extract_defaults)

    # flux and variances are nonzero
    for data in result[:4]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[0],)

    # background and variances are zero
    for data in result[4:8]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[0],)

    # npixels is the sum of the profile
    assert np.allclose(result[8], np.sum(simple_profile, axis=0))

    # scene model has 2D shape
    assert result[-1].shape == model.data.shape


def test_extract_one_slit_multi_int(mock_nirspec_bots, extract_defaults,
                                    simple_profile, log_watcher):
    model = mock_nirspec_bots
    extract_defaults['dispaxis'] = 1

    log_watcher.message = "Extracting integration 2"
    result = ex.extract_one_slit(model, 1, simple_profile, None, extract_defaults)
    log_watcher.assert_seen()

    # flux and variances are nonzero
    for data in result[:4]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[2],)

    # background and variances are zero
    for data in result[4:8]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[2],)

    # npixels is the sum of the profile
    assert np.allclose(result[8], np.sum(simple_profile, axis=0))

    # scene model has 2D shape
    assert result[-1].shape == model.data.shape[-2:]


def test_extract_one_slit_missing_var(mock_nirspec_fs_one_slit, extract_defaults,
                                      simple_profile):
    model = mock_nirspec_fs_one_slit
    extract_defaults['dispaxis'] = 1

    # Test that mismatched variances still extract okay.
    # This is probably only possible for var_flat, which is optional and
    # uninitialized if flat fielding is skipped, but the code has handling
    # for all 3 variance arrays.
    model.var_rnoise = np.zeros((10, 10))
    model.var_poisson = np.zeros((10, 10))
    model.var_flat = np.zeros((10, 10))

    result = ex.extract_one_slit(model, -1, simple_profile, None, extract_defaults)

    # flux is nonzero
    assert np.all(result[0] > 0)
    assert result[0].shape == (model.data.shape[1],)

    # variances are zero
    for data in result[1:4]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[1],)


def test_create_extraction_with_photom(create_extraction_inputs):
    model = create_extraction_inputs[0]
    model.meta.cal_step.photom = 'COMPLETE'

    ex.create_extraction(*create_extraction_inputs)

    output_model = create_extraction_inputs[2]
    assert output_model.spec[0].spec_table.columns['flux'].unit == 'Jy'


def test_create_extraction_without_photom(create_extraction_inputs):
    model = create_extraction_inputs[0]
    model.meta.cal_step.photom = 'SKIPPED'

    ex.create_extraction(*create_extraction_inputs)

    output_model = create_extraction_inputs[2]
    assert output_model.spec[0].spec_table.columns['flux'].unit == 'DN/s'


def test_create_extraction_missing_src_type(create_extraction_inputs):
    model = create_extraction_inputs[0]
    model.source_type = None
    model.meta.target.source_type = 'EXTENDED'

    ex.create_extraction(*create_extraction_inputs)

    output_model = create_extraction_inputs[2]
    assert output_model.spec[0].source_type == 'EXTENDED'


def test_create_extraction_no_match(create_extraction_inputs):
    create_extraction_inputs[4] = 'bad slitname'
    with pytest.raises(ValueError, match="Missing extraction parameters"):
        ex.create_extraction(*create_extraction_inputs)


def test_create_extraction_partial_match(create_extraction_inputs, log_watcher):
    # match a slit that has a mismatched spectral order specified
    create_extraction_inputs[4] = 'slit7'

    log_watcher.message = 'Spectral order 1 not found'
    with pytest.raises(ex.ContinueError):
        ex.create_extraction(*create_extraction_inputs)
    log_watcher.assert_seen()


def test_create_extraction_missing_dispaxis(create_extraction_inputs, log_watcher):
    create_extraction_inputs[0].meta.wcsinfo.dispersion_direction = None
    log_watcher.message = 'dispersion direction information is missing'
    with pytest.raises(ex.ContinueError):
        ex.create_extraction(*create_extraction_inputs)
    log_watcher.assert_seen()


def test_create_extraction_missing_wavelengths(create_extraction_inputs, log_watcher):
    model = create_extraction_inputs[0]
    model.wavelength = np.full_like(model.data, np.nan)
    log_watcher.message = 'Spectrum is empty; no valid data'
    with pytest.raises(ex.ContinueError):
        ex.create_extraction(*create_extraction_inputs)
    log_watcher.assert_seen()


def test_create_extraction_nrs_apcorr(create_extraction_inputs, nirspec_fs_apcorr,
                                      mock_nirspec_bots, log_watcher):
    model = mock_nirspec_bots
    model.source_type = 'POINT'
    model.meta.cal_step.photom = 'COMPLETE'
    create_extraction_inputs[0] = model

    log_watcher.message = 'Tabulating aperture correction'
    ex.create_extraction(*create_extraction_inputs, apcorr_ref_model=nirspec_fs_apcorr,
                         use_source_posn=False)
    log_watcher.assert_seen()


def test_create_extraction_one_int(create_extraction_inputs, mock_nirspec_bots, log_watcher):
    # Input model is a cube, but with only one integration
    model = mock_nirspec_bots
    model.data = model.data[0].reshape(1, *model.data.shape[-2:])
    create_extraction_inputs[0] = model

    log_watcher.message = '1 integration done'
    ex.create_extraction(*create_extraction_inputs, log_increment=1)
    output_model = create_extraction_inputs[2]
    assert len(output_model.spec) == 1
    log_watcher.assert_seen()


def test_create_extraction_log_increment(
        create_extraction_inputs, mock_nirspec_bots, log_watcher):
    create_extraction_inputs[0] = mock_nirspec_bots

    # all integrations are logged
    log_watcher.message = '... 9 integrations done'
    ex.create_extraction(*create_extraction_inputs, log_increment=1)
    log_watcher.assert_seen()


@pytest.mark.parametrize('use_source', [True, False, None])
@pytest.mark.parametrize('source_type', ['POINT', 'EXTENDED'])
def test_create_extraction_use_source(
        monkeypatch, create_extraction_inputs, mock_nirspec_fs_one_slit,
        use_source, source_type, log_watcher):
    model = mock_nirspec_fs_one_slit
    model.source_type = source_type
    create_extraction_inputs[0] = model

    # mock the source location function
    def mock_source_location(*args):
        return 24, 7.74, 9.5

    monkeypatch.setattr(ex, 'location_from_wcs', mock_source_location)

    if source_type != 'POINT' and use_source is None:
        # If not specified, source position should be used only if POINT
        log_watcher.message = 'Setting use_source_posn to False'
    elif use_source is True or (source_type == 'POINT' and use_source is None):
        # If explicitly set to True, or unspecified + source type is POINT,
        # source position is used
        log_watcher.message = 'Aperture start/stop: -15'
    else:
        # If False, source position is not used 
        log_watcher.message = 'Aperture start/stop: 0'
    ex.create_extraction(*create_extraction_inputs, use_source_posn=use_source)
    log_watcher.assert_seen()


def test_run_extract1d(mock_nirspec_mos):
    model = mock_nirspec_mos
    output_model, profile_model, scene_model = ex.run_extract1d(model)
    assert isinstance(output_model, dm.MultiSpecModel)
    assert profile_model is None
    assert scene_model is None
    output_model.close()


def test_run_extract1d_save_models(mock_niriss_wfss_l3):
    model = mock_niriss_wfss_l3
    output_model, profile_model, scene_model = ex.run_extract1d(
        model, save_profile=True, save_scene_model=True)
    assert isinstance(output_model, dm.MultiSpecModel)
    assert isinstance(profile_model, ModelContainer)
    assert isinstance(scene_model, ModelContainer)

    assert len(profile_model) == len(model)
    assert len(scene_model) == len(model)

    for pmodel in profile_model:
        assert isinstance(pmodel, dm.ImageModel)
    for smodel in scene_model:
        assert isinstance(smodel, dm.ImageModel)

    output_model.close()
    profile_model.close()
    scene_model.close()


def test_run_extract1d_save_cube_scene(mock_nirspec_bots):
    model = mock_nirspec_bots
    output_model, profile_model, scene_model = ex.run_extract1d(
        model, save_profile=True, save_scene_model=True)
    assert isinstance(output_model, dm.MultiSpecModel)
    assert isinstance(profile_model, dm.ImageModel)
    assert isinstance(scene_model, dm.CubeModel)

    assert profile_model.data.shape == model.data.shape[-2:]
    assert scene_model.data.shape == model.data.shape

    output_model.close()
    profile_model.close()
    scene_model.close()



def test_run_extract1d_tso(mock_nirspec_bots):
    model = mock_nirspec_bots
    output_model, _, _ = ex.run_extract1d(model)

    # time and integration keywords are populated
    for i, spec in enumerate(output_model.spec):
        assert spec.int_num == i + 1

    output_model.close()


@pytest.mark.parametrize('from_name_attr', [True, False])
def test_run_extract1d_slitmodel_name(mock_nirspec_fs_one_slit, from_name_attr):
    model = mock_nirspec_fs_one_slit
    slit_name = 'S200A1'
    if from_name_attr:
        model.name = slit_name
        model.meta.instrument.fixed_slit = None
    else:
        model.name = None
        model.meta.instrument.fixed_slit = 'S200A1'

    output_model, _, _ = ex.run_extract1d(model)
    assert output_model.spec[0].name == 'S200A1'

    output_model.close()


@pytest.mark.parametrize('from_name_attr', [True, False])
def test_run_extract1d_imagemodel_name(mock_miri_lrs_fs, from_name_attr):
    model = mock_miri_lrs_fs
    slit_name = 'test_slit_name'
    if from_name_attr:
        model.name = slit_name
    else:
        model.name = None

    output_model, _, _ = ex.run_extract1d(model)
    if from_name_attr:
        assert output_model.spec[0].name == 'test_slit_name'
    else:
        assert output_model.spec[0].name == 'MIR_LRS-FIXEDSLIT'
    output_model.close()


def test_run_extract1d_apcorr(mock_miri_lrs_fs, miri_lrs_apcorr_file, log_watcher):
    model = mock_miri_lrs_fs
    model.meta.target.source_type = 'POINT'

    log_watcher.message = 'Creating aperture correction'
    output_model, _, _ = ex.run_extract1d(model, apcorr_ref_name=miri_lrs_apcorr_file)
    log_watcher.assert_seen()

    output_model.close()


def test_run_extract1d_invalid():
    model = dm.MultiSpecModel()
    with pytest.raises(RuntimeError, match="Can't extract a spectrum"):
        ex.run_extract1d(model)


def test_run_extract1d_zeroth_order_slit(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    model.meta.wcsinfo.spectral_order = 0
    output_model, _, _ = ex.run_extract1d(model)

    # no spectra extracted for zeroth order
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_zeroth_order_image(mock_miri_lrs_fs):
    model = mock_miri_lrs_fs
    model.meta.wcsinfo.spectral_order = 0
    output_model, _, _ = ex.run_extract1d(model)

    # no spectra extracted for zeroth order
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_zeroth_order_multispec(mock_nirspec_mos):
    model = mock_nirspec_mos
    for slit in model.slits:
        slit.meta.wcsinfo.spectral_order = 0
    output_model, _, _ = ex.run_extract1d(model)

    # no spectra extracted for zeroth order
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_no_data(mock_niriss_wfss_l3):
    container = mock_niriss_wfss_l3
    for model in container:
        model.data = np.array([])
    output_model, _, _ = ex.run_extract1d(container)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_continue_error_slit(monkeypatch, mock_nirspec_fs_one_slit):
    def raise_continue_error(*args, **kwargs):
        raise ex.ContinueError('Test error')

    monkeypatch.setattr(ex, 'create_extraction', raise_continue_error)
    output_model, _, _ = ex.run_extract1d(mock_nirspec_fs_one_slit)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_continue_error_image(monkeypatch, mock_miri_lrs_fs):
    def raise_continue_error(*args, **kwargs):
        raise ex.ContinueError('Test error')

    monkeypatch.setattr(ex, 'create_extraction', raise_continue_error)
    output_model, _, _ = ex.run_extract1d(mock_miri_lrs_fs)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_continue_error_multislit(monkeypatch, mock_nirspec_mos):
    def raise_continue_error(*args, **kwargs):
        raise ex.ContinueError('Test error')

    monkeypatch.setattr(ex, 'create_extraction', raise_continue_error)
    output_model, _, _ = ex.run_extract1d(mock_nirspec_mos)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()
