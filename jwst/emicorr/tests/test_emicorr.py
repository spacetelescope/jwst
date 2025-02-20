"""Unit tests for EMI correction."""

import logging

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import RampModel, EmiModel

from jwst.emicorr import emicorr, emicorr_step
from jwst.tests.helpers import LogWatcher


subarray_example_case = {
    'MASK1550': {'frameclocks': 23968,
        'freqs': {
            'FAST': ["Hz390", "Hz10"],
            'SLOW': {
                'MIRIFULONG': ['Hz390', 'Hz10_slow_MIRIFULONG'],
                'MIRIFUSHORT': ['Hz390', 'Hz10_slow_MIRIFUSHORT'],
                'MIRIMAGE': ['Hz390', 'Hz10_slow_MIRIMAGE']}},
        'rowclocks': 82}
}

emimdl = {'frequencies': {
            "Hz390": {'frequency': 390.625,
                'phase_amplitudes': np.ones(20)+0.1},
            "Hz10": { 'frequency': 10.039216,
                'phase_amplitudes': np.ones(20)+0.5}},
        'subarray_cases': subarray_example_case
        }
emicorr_model = EmiModel(emimdl)


def mk_data_mdl(data, subarray, readpatt, detector):
    # create input_model
    input_model = RampModel(data=data)
    input_model.meta.instrument.name = 'MIRI'
    input_model.meta.instrument.detector = detector
    input_model.meta.exposure.type = 'MIR_4QPM'
    input_model.meta.subarray.name = subarray
    input_model.meta.exposure.readpatt = readpatt
    input_model.meta.subarray.xsize = 288
    input_model.meta.subarray.xstart = 1
    input_model.meta.exposure.nsamples = 1
    return input_model


@pytest.fixture
def log_watcher(monkeypatch):
    # Set a log watcher to check for a log message at any level
    # in the emicorr module
    watcher = LogWatcher('')
    logger = logging.getLogger('jwst.emicorr.emicorr')
    for level in ['debug', 'info', 'warning', 'error']:
        monkeypatch.setattr(logger, level, watcher)
    return watcher


def test_emicorrstep():
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, 'MASK1550', 'FAST', 'MIRIMAGE')
    expected_outmdl = input_model.copy()

    nirmdl = input_model.copy()
    nirmdl.meta.instrument.name = 'NIRISS'

    stp = emicorr_step.EmiCorrStep()
    nir_result = stp.call(nirmdl)

    # expect no change because instrument is not MIRI
    assert expected_outmdl.data.all() == nir_result.data.all()


@pytest.mark.parametrize('data_case', ['flat', 'linear'])
@pytest.mark.parametrize('algorithm', ['sequential', 'joint'])
def test_apply_emicorr(data_case, algorithm):
    data = np.ones((1, 5, 20, 20))
    if data_case == 'linear':
        linear_ramp = np.arange(1, 6, dtype=float)
        data[0, :, ...] = linear_ramp[:, None, None]
    input_model = mk_data_mdl(data, 'MASK1550', 'FAST', 'MIRIMAGE')
    pars = {
        'save_onthefly_reffile': None,
        'nints_to_phase': None,
        'nbins': None,
        'scale_reference': True,
        'onthefly_corr_freq': None,
        'use_n_cycles': None,
        'algorithm': algorithm
    }
    outmdl = emicorr.apply_emicorr(input_model.copy(), emicorr_model, **pars)
    assert isinstance(outmdl, RampModel)

    # flat or linear ramp data shows no correction
    assert np.allclose(outmdl.data, input_model.data)


@pytest.mark.parametrize('algorithm', ['sequential', 'joint'])
@pytest.mark.parametrize('emicorr_model', ['bad_model', None])
@pytest.mark.parametrize('output_ext', ['fits', 'asdf'])
def test_apply_emicorr_with_freq(tmp_path, algorithm, log_watcher, emicorr_model, output_ext):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, 'MASK1550', 'FAST', 'MIRIMAGE')

    # emicorr_model is ignored if user frequencies are provided
    onthefly_corr_freq = [218.3]

    output_name = tmp_path / f'otf_reffile.{output_ext}'
    expected_output_name = tmp_path / 'otf_reffile.asdf'

    log_watcher.message = "'joint' algorithm cannot be used"
    outmdl = emicorr.apply_emicorr(
        input_model.copy(), emicorr_model, onthefly_corr_freq=onthefly_corr_freq,
        save_onthefly_reffile=str(output_name), nints_to_phase=None,
        nbins=None, scale_reference=True, algorithm=algorithm)

    # flat data has no correction
    assert np.allclose(outmdl.data, input_model.data)

    # joint model should warn and fall back to sequential for user frequencies
    if algorithm == 'joint':
        log_watcher.assert_seen()
    else:
        log_watcher.assert_not_seen()

    # reference file is saved to an asdf file
    assert expected_output_name.exists()


def test_get_subarcase():
    subarray, readpatt, detector = 'MASK1550', 'FAST', 'MIRIMAGE'
    subarray_info_r = emicorr.get_subarcase(emicorr_model, subarray, readpatt, detector)
    subname_r, rowclocks_r, frameclocks_r, freqs2correct_r = subarray_info_r

    # test if we get the right configuration
    compare_real = ['MASK1550', 82, 23968, ["Hz390", "Hz10"]]
    subname_real, rowclocks_real, frameclocks_real, freqs2correct_real = compare_real

    assert subname_real == subname_r
    assert rowclocks_real == rowclocks_r
    assert frameclocks_real == frameclocks_r
    assert freqs2correct_real[0] == freqs2correct_r[0]


def test_sloper():
    data = np.ones((5, 5, 5))
    outarray, intercept = emicorr.sloper(data)

    compare_arr, compare_intercept = np.zeros((5, 5)), np.ones((5, 5))
    assert compare_arr.all() == outarray.all()
    assert compare_intercept.all() == intercept.all()


def test_minmed():
    data = np.ones((5, 5, 5))
    compare_arr = data.copy()
    medimg = emicorr.minmed(data)
    assert compare_arr.all() == medimg.all()


def test_rebin():
    data = np.ones(10)
    data[1] = 0.55
    data[5] = 1.55
    data[9] = 2.0
    rebinned_data = emicorr.rebin(data, [7])
    compare_arr = np.array([1., 0.55, 1., 1., 1.55, 1., 1.])
    assert compare_arr.all() == rebinned_data.all()
