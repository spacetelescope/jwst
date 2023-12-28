"""

Unit tests for EMI correction

"""


import numpy as np
from jwst.emicorr import emicorr, emicorr_step
from stdatamodels.jwst.datamodels import Level1bModel, EmiModel

default_subarray_cases = emicorr.default_subarray_cases

def mk_data_mdl(data, subarray, readpatt, detector):
    # create input_model
    input_model = Level1bModel(data=data)
    input_model.meta.instrument.name = 'MIRI'
    input_model.meta.instrument.detector = detector
    input_model.meta.exposure.type = 'MIR_4QPM'
    input_model.meta.subarray.name = subarray
    input_model.meta.exposure.readpatt = readpatt
    input_model.meta.subarray.xsize = 288
    input_model.meta.subarray.xstart = 1
    input_model.meta.exposure.nsamples = 1
    return input_model


def test_EmiCorrStep():
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, 'MASK1550', 'FAST', 'MIRIMAGE')
    expected_outmdl = input_model.copy()

    nirmdl = input_model.copy()
    nirmdl.meta.instrument.name = 'NIRISS'

    stp = emicorr_step.EmiCorrStep()
    nir_result = stp.call(nirmdl)

    # expect no change because instrument is not MIRI
    assert expected_outmdl.data.all() == nir_result.data.all()


def test_do_correction():
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, 'MASK1550', 'FAST', 'MIRIMAGE')
    pars = {
        'save_intermediate_results': False,
        'user_supplied_reffile': None,
        'nints_to_phase': None,
        'nbins': None,
        'scale_reference': True
    }
    emimdl = {'frequencies': {
                "Hz390": {'frequency': 390.625,
                    'phase_amplitudes': np.ones(20)+0.1},
                "Hz10": { 'frequency': 10.039216,
                    'phase_amplitudes': np.ones(20)+0.5}},
            'subarray_cases': default_subarray_cases
    }
    save_onthefly_reffile = None
    emicorr_model = EmiModel(emimdl)
    outmdl = emicorr.do_correction(input_model, emicorr_model, save_onthefly_reffile, **pars)

    assert outmdl is not None


def test_apply_emicorr():
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, 'MASK1550', 'FAST', 'MIRIMAGE')
    emicorr_model, save_onthefly_reffile = None, None
    outmdl = emicorr.apply_emicorr(input_model, emicorr_model, save_onthefly_reffile,
                save_intermediate_results=False, user_supplied_reffile=None,
                nints_to_phase=None, nbins_all=None, scale_reference=True)

    assert outmdl is not None


def test_get_subarcase():
    subarray, readpatt, detector = 'MASK1550', 'FAST', 'MIRIMAGE'
    subarray_info_r = emicorr.get_subarcase(default_subarray_cases, subarray, readpatt, detector)
    subname_r, rowclocks_r, frameclocks_r, freqs2correct_r = subarray_info_r

    # set up a fake configuration
    subarray_info_f = emicorr.get_subarcase(default_subarray_cases, 'FULL', 'kkkkk', detector)
    subname_f, rowclocks_f, frameclocks_f, freqs2correct_f = subarray_info_f

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


def test_iter_stat_sig_clip():
    data = np.ones((5, 5, 5, 5))
    data[0, 0, 0, 0] = 0.55
    data[1, 1, 1, 1] = 1.11
    data[2, 2, 2, 2] = 2.22
    data[3, 3, 3, 3] = 3.33
    data[4, 4, 4, 4] = 4.44
    u = np.where(data > 1)
    dmean, dmedian, dsigma, dmask = emicorr.iter_stat_sig_clip(data[u])

    compare_mean, compare_median = 2.7750000000000004, 0.0
    compare_sigma, compare_mask = 1.326703756361177, np.ones(np.size(data), dtype='b')+1

    assert compare_mean == dmean
    assert compare_median == dmedian
    assert compare_sigma == dsigma
    assert compare_mask.all() == dmask.all()


def test_rebin():
    data = np.ones(10)
    data[1] = 0.55
    data[5] = 1.55
    data[9] = 2.0
    rebinned_data = emicorr.rebin(data, [7])

    compare_arr = np.array([1., 0.55, 1. , 1., 1.55, 1., 1.])

    assert compare_arr.all() == rebinned_data.all()
