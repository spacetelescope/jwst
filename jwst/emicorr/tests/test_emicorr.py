"""Unit tests for EMI correction."""

import warnings

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import RampModel, EmiModel

from jwst.pipeline import Detector1Pipeline
from jwst.emicorr import emicorr, emicorr_step


@pytest.fixture()
def subarray_example_case():
    example = {
        "MASK1550": {
            "frameclocks": 23968,
            "freqs": {
                "FAST": ["Hz390", "Hz10"],
                "SLOW": {
                    "MIRIFULONG": ["Hz390", "Hz10_slow_MIRIFULONG"],
                    "MIRIFUSHORT": ["Hz390", "Hz10_slow_MIRIFUSHORT"],
                    "MIRIMAGE": ["Hz390", "Hz10_slow_MIRIMAGE"],
                },
            },
            "rowclocks": 82,
        },
        "FULL_FAST": {
            "frameclocks": 277504,
            "freqs": {
                "FAST": ["Hz10"],
                "SLOW": {
                    "MIRIFULONG": ["Hz10_slow_MIRIFULONG"],
                    "MIRIFUSHORT": ["Hz10_slow_MIRIFUSHORT"],
                    "MIRIMAGE": ["Hz10_slow_MIRIMAGE"],
                },
            },
            "rowclocks": 271,
        },
        "FULL_SLOW": {
            "frameclocks": 2388992,
            "freqs": {
                "SLOW": {
                    "MIRIFULONG": ["Hz10_slow_MIRIFULONG"],
                    "MIRIFUSHORT": ["Hz10_slow_MIRIFUSHORT"],
                    "MIRIMAGE": ["Hz10_slow_MIRIMAGE"],
                },
            },
            "rowclocks": 2333,
        },
    }
    return example


@pytest.fixture()
def emicorr_model(subarray_example_case):
    frequencies = {
        "Hz390": {"frequency": 390.625, "phase_amplitudes": np.ones(20) + 0.1},
        "Hz10": {"frequency": 10.039216, "phase_amplitudes": np.ones(20) + 0.5},
    }

    emi_model = EmiModel()
    emi_model.frequencies = frequencies
    emi_model.subarray_cases = subarray_example_case

    emi_model.meta.reftype = "emicorr"
    emi_model.meta.author = "JWST calibration pipeline"
    emi_model.meta.description = "EMI correction file"
    emi_model.meta.pedigree = "GROUND"
    emi_model.meta.useafter = "2024-01-01T00:00:00"

    return emi_model


def mk_data_mdl(data, subarray, readpatt, detector):
    # create input_model
    input_model = RampModel(data=data)
    input_model.meta.instrument.name = "MIRI"
    input_model.meta.instrument.detector = detector
    input_model.meta.exposure.type = "MIR_4QPM"
    input_model.meta.subarray.name = subarray
    input_model.meta.exposure.readpatt = readpatt
    input_model.meta.subarray.xsize = 288
    input_model.meta.subarray.xstart = 1
    input_model.meta.exposure.nsamples = 1
    input_model.meta.observation.date = "2024-01-01"
    input_model.meta.observation.time = "00:00:00"
    return input_model


def add_emi(data):
    amp = 0.1
    phase = 1 / 3

    freq = 10.0
    rowclocks = 2000
    period = (1.0 / freq) / 10.0e-6

    ni, ng, ny, nx = data.shape
    extra_rowclocks = (1024.0 - ny) * (4 + 3.0)
    readtimes = np.zeros((ng, ny, nx))
    for i in range(ng):
        for j in range(ny):
            for k in range(nx // 4):
                readtimes[i, j, k * 4 : k * 4 + 4] = i * extra_rowclocks + j * rowclocks + k

    emi = amp * np.sin(2 * np.pi * readtimes / period + phase)
    data = data + emi[None, :, :, :]
    return data


@pytest.fixture()
def data_without_emi():
    nint = 1
    ngroup = 10
    ny = 400
    nx = 200
    data = np.ones((nint, ngroup, ny, nx))
    linear_ramp = np.arange(1, ngroup + 1, dtype=float)
    data[:, :, ...] = linear_ramp[None, :, None, None]
    return data


@pytest.fixture()
def data_without_emi_3int():
    nint = 3
    ngroup = 10
    ny = 400
    nx = 200
    data = np.ones((nint, ngroup, ny, nx))
    linear_ramp = np.arange(1, ngroup + 1, dtype=float)
    data[:, :, ...] = linear_ramp[None, :, None, None]
    return data


@pytest.fixture()
def data_with_emi(data_without_emi):
    return add_emi(data_without_emi)


@pytest.fixture()
def data_with_emi_3int(data_without_emi_3int):
    return add_emi(data_without_emi_3int)


@pytest.fixture()
def model_with_emi(emicorr_model):
    freq = 10.0
    pa = np.sin(2 * np.pi * np.arange(500) / 500)

    subarray_cases = {
        "FULL_FAST": {
            "frameclocks": 200000,
            "freqs": {
                "FAST": ["test"],
            },
            "rowclocks": 2000,
        }
    }
    frequencies = {
        "test": {"frequency": freq, "phase_amplitudes": pa},
    }

    emicorr_model.frequencies = frequencies
    emicorr_model.subarray_cases = subarray_cases
    return emicorr_model


@pytest.fixture()
def module_log_watcher(monkeypatch):
    # Set a log watcher to check for a log message at any level
    # in the emicorr module
    watcher = LogWatcher("")
    logger = logging.getLogger("jwst.emicorr.emicorr")
    for level in ["debug", "info", "warning", "error"]:
        monkeypatch.setattr(logger, level, watcher)
    return watcher


def test_emicorrstep_skip_default():
    step = emicorr_step.EmiCorrStep()
    # Default is that step is skipped
    assert step.skip == True


@pytest.mark.parametrize("skip", [True, False])
def test_run_in_pipeline(skip):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "FAST", "MIRIMAGE")

    pipeline_instance = Detector1Pipeline()

    assert pipeline_instance.steps["emicorr"]["skip"] == True

    # Run the pipeline, omitting steps that are incompatible with our datamodel
    cleaned = pipeline_instance.call(
        input_model,
        steps={
            "emicorr": {"skip": skip},
            "dq_init": {"skip": True},
            "saturation": {"skip": True},
            "ipc": {"skip": True},
            "reset": {"skip": True},
            "linearity": {"skip": True},
            "rscd": {"skip": True},
            "dark_current": {"skip": True},
            "jump": {"skip": True},
            "ramp_fit": {"skip": True},
        },
    )

    if skip:
        assert cleaned.meta.cal_step.emicorr == "SKIPPED"
    else:
        assert cleaned.meta.cal_step.emicorr == "COMPLETE"

    input_model.close()
    cleaned.close()


def test_emicorrstep_skip_instrument(log_watcher):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "FAST", "MIRIMAGE")

    nirmdl = input_model.copy()
    nirmdl.meta.instrument.name = "NIRISS"

    watcher = log_watcher("stpipe.EmiCorrStep", message="not implemented for instrument")
    step = emicorr_step.EmiCorrStep()
    nir_result = step.call(nirmdl, skip=False)
    watcher.assert_seen()

    # expect no change because instrument is not MIRI
    assert np.all(input_model.data == nir_result.data)
    assert nir_result.meta.cal_step.emicorr == "SKIPPED"


def test_emicorrstep_skip_readpatt(log_watcher):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "ANY", "MIRIMAGE")

    watcher = log_watcher("stpipe.EmiCorrStep", message="not implemented for read pattern")
    step = emicorr_step.EmiCorrStep()
    result = step.call(input_model, skip=False)
    watcher.assert_seen()

    # expect no change because read pattern is not supported
    assert np.all(input_model.data == result.data)
    assert result.meta.cal_step.emicorr == "SKIPPED"


def test_emicorrstep_skip_no_reffile(monkeypatch, log_watcher):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "FAST", "MIRIMAGE")

    # mock no reference file found
    monkeypatch.setattr(emicorr_step.EmiCorrStep, "get_reference_file", lambda *args: "N/A")
    step = emicorr_step.EmiCorrStep()

    watcher = log_watcher("stpipe.EmiCorrStep", message="No reference file")
    result = step.call(input_model, skip=False)
    watcher.assert_seen()

    # expect no change because step is skipped
    assert np.all(input_model.data == result.data)
    assert result.meta.cal_step.emicorr == "SKIPPED"


def test_emicorrstep_skip_for_failure(monkeypatch, log_watcher):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "FAST", "MIRIMAGE")

    # mock a failure internal to the correction algorithm
    monkeypatch.setattr(emicorr, "apply_emicorr", lambda *args, **kwargs: None)
    step = emicorr_step.EmiCorrStep()

    watcher = log_watcher("stpipe.EmiCorrStep", message="Step skipped")
    result = step.call(input_model, skip=False)
    watcher.assert_seen()

    # expect no change because step is skipped
    assert np.all(input_model.data == result.data)
    assert result.meta.cal_step.emicorr == "SKIPPED"


def test_emicorrstep_skip_for_small_groups(log_watcher):
    data = np.ones((1, 2, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "FAST", "MIRIMAGE")

    watcher = log_watcher("jwst.emicorr.emicorr", message="cannot be performed for ngroups=2")
    step = emicorr_step.EmiCorrStep()
    result = step.call(input_model, skip=False)
    watcher.assert_seen()

    # step is skipped
    assert np.all(input_model.data == result.data)
    assert result.meta.cal_step.emicorr == "SKIPPED"


@pytest.mark.parametrize("algorithm", ["sequential", "joint"])
@pytest.mark.parametrize("subarray", ["MASK1550", "FULL"])
def test_emicorrstep_succeeds(algorithm, subarray):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, subarray, "FAST", "MIRIMAGE")

    step = emicorr_step.EmiCorrStep()
    result = step.call(input_model, skip=False, algorithm=algorithm)

    # step completes but we expect no change for flat data
    assert np.all(input_model.data == result.data)
    assert result.meta.cal_step.emicorr == "COMPLETE"


@pytest.mark.parametrize("readpatt", ["FAST", "SLOW"])
def test_emicorrstep_user_freq(tmp_path, readpatt):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "FULL", readpatt, "MIRIMAGE")
    input_model.meta.filename = "test.fits"

    corr_freq = [218.3]
    step = emicorr_step.EmiCorrStep()
    result = step.call(
        input_model,
        skip=False,
        onthefly_corr_freq=corr_freq,
        save_intermediate_results=True,
        output_dir=str(tmp_path),
    )

    # step completes but we expect no change for flat data
    assert np.all(input_model.data == result.data)
    assert result.meta.cal_step.emicorr == "COMPLETE"

    # on-the-fly reference file is saved
    assert (tmp_path / "test_emi_ref_waves.asdf").exists()


def test_emicorrstep_user_reffile(tmp_path, emicorr_model):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "FAST", "MIRIMAGE")

    model_name = str(tmp_path / "emicorr.asdf")
    emicorr_model.save(model_name)

    step = emicorr_step.EmiCorrStep()
    with warnings.catch_warnings():
        # Warning emitted for numpy 1.26 but not numpy 2.2
        warnings.filterwarnings("ignore", message="Polyfit may be poorly conditioned")
        result = step.call(input_model, skip=False, user_supplied_reffile=model_name)

    # step completes but we expect no change for flat data
    assert np.all(input_model.data == result.data)
    assert result.meta.cal_step.emicorr == "COMPLETE"


@pytest.mark.parametrize("data_case", ["flat", "linear"])
@pytest.mark.parametrize("algorithm", ["sequential", "joint"])
@pytest.mark.parametrize("readpatt", ["FAST", "FASTR1"])
def test_apply_emicorr_noiseless(data_case, algorithm, emicorr_model, readpatt):
    data = np.ones((1, 5, 20, 20))
    if data_case == "linear":
        linear_ramp = np.arange(1, 6, dtype=float)
        data[0, :, ...] = linear_ramp[:, None, None]
    input_model = mk_data_mdl(data, "MASK1550", readpatt, "MIRIMAGE")
    pars = {
        "save_onthefly_reffile": None,
        "nints_to_phase": None,
        "nbins": None,
        "scale_reference": True,
        "onthefly_corr_freq": None,
        "use_n_cycles": None,
        "algorithm": algorithm,
    }
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Polyfit may be poorly conditioned")
        outmdl = emicorr.apply_emicorr(input_model.copy(), emicorr_model, **pars)
    assert isinstance(outmdl, RampModel)

    # flat or linear ramp data shows no correction
    assert np.allclose(outmdl.data, input_model.data)


@pytest.mark.parametrize("algorithm,accuracy", [("sequential", 0.1), ("joint", 0.0001)])
def test_apply_emicorr(data_without_emi, data_with_emi, model_with_emi, algorithm, accuracy):
    input_model = mk_data_mdl(data_with_emi, "FULL", "FAST", "MIRIMAGE")
    expected_model = mk_data_mdl(data_without_emi, "FULL", "FAST", "MIRIMAGE")

    outmdl = emicorr.apply_emicorr(input_model.copy(), model_with_emi, algorithm=algorithm)

    # Corrected ramp should be close to ramp without noise
    assert np.allclose(outmdl.data, expected_model.data, rtol=accuracy)


def test_apply_emicorr_separate_ints(data_without_emi_3int, data_with_emi_3int, model_with_emi):
    input_model = mk_data_mdl(data_with_emi_3int, "FULL", "FAST", "MIRIMAGE")
    expected_model = mk_data_mdl(data_without_emi_3int, "FULL", "FAST", "MIRIMAGE")

    outmdl = emicorr.apply_emicorr(
        input_model.copy(), model_with_emi, algorithm="joint", fit_ints_separately=True
    )

    # Corrected ramp should be close to ramp without noise
    assert np.allclose(outmdl.data, expected_model.data, rtol=1e-2)


@pytest.mark.parametrize("algorithm", ["sequential", "joint"])
@pytest.mark.parametrize("model_placeholder", ["bad_model", None])
@pytest.mark.parametrize("output_ext", ["fits", "asdf"])
def test_apply_emicorr_with_freq(tmp_path, algorithm, log_watcher, model_placeholder, output_ext):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "MASK1550", "FAST", "MIRIMAGE")

    onthefly_corr_freq = [218.3]

    output_name = tmp_path / f"otf_reffile.{output_ext}"
    expected_output_name = tmp_path / "otf_reffile.asdf"

    watcher = log_watcher("jwst.emicorr.emicorr", message="'joint' algorithm cannot be used")

    # emicorr model is ignored if user frequencies are provided - any placeholder will work
    outmdl = emicorr.apply_emicorr(
        input_model.copy(),
        model_placeholder,
        onthefly_corr_freq=onthefly_corr_freq,
        save_onthefly_reffile=str(output_name),
        nints_to_phase=None,
        nbins=None,
        scale_reference=True,
        algorithm=algorithm,
    )

    # flat data has no correction
    assert np.allclose(outmdl.data, input_model.data)

    # joint model should warn and fall back to sequential for user frequencies
    if algorithm == "joint":
        watcher.assert_seen()
    else:
        watcher.assert_not_seen()

    # reference file is saved to an asdf file
    assert expected_output_name.exists()


def test_apply_emicorr_no_config_found(log_watcher, emicorr_model):
    data = np.ones((1, 5, 20, 20))
    input_model = mk_data_mdl(data, "SUB64", "FAST", "MIRIMAGE")

    watcher = log_watcher("jwst.emicorr.emicorr", message="No correction match")
    result = emicorr.apply_emicorr(input_model, emicorr_model)
    assert result is None
    watcher.assert_seen()


@pytest.mark.parametrize("detector", ["MIRIMAGE", "MIRIFULONG"])
@pytest.mark.parametrize("subarray", ["FULL", "MASK1550"])
@pytest.mark.parametrize("readpatt", ["FAST", "SLOW"])
def test_get_subarcase(emicorr_model, subarray_example_case, detector, subarray, readpatt):
    subarray_info_r = emicorr.get_subarcase(emicorr_model, subarray, readpatt, detector)
    subname_r, rowclocks_r, frameclocks_r, freqs2correct_r = subarray_info_r

    # test if we get the right configuration
    if subarray == "FULL":
        if readpatt == "FAST":
            subname_real = "FULL_FAST"
        else:
            subname_real = "FULL_SLOW"
    else:
        subname_real = subarray
    rowclocks_real = subarray_example_case[subname_real]["rowclocks"]
    frameclocks_real = subarray_example_case[subname_real]["frameclocks"]
    freqs2correct_real = subarray_example_case[subname_real]["freqs"][readpatt]
    if readpatt == "SLOW":
        freqs2correct_real = freqs2correct_real[detector]

    assert subname_real == subname_r
    assert rowclocks_real == rowclocks_r
    assert frameclocks_real == frameclocks_r
    assert freqs2correct_real == freqs2correct_r


def test_get_subarcase_bad_readpatt(emicorr_model, log_watcher):
    subarray = "FULL"
    detector = "MIRIMAGE"

    # readpattern not recognized
    readpatt = "ANY"
    watcher = log_watcher("jwst.emicorr.emicorr", message="does not include expected string")
    subarray_info_r = emicorr.get_subarcase(emicorr_model, subarray, readpatt, detector)
    watcher.assert_seen()

    assert subarray_info_r == (None, None, None, None)


def test_get_subarcase_missing_subarray(emicorr_model, log_watcher):
    detector = "MIRIMAGE"
    readpatt = "FAST"

    # subarray not recognized
    subarray = "ANY"
    watcher = log_watcher("jwst.emicorr.emicorr", message="Subarray ANY not found")
    subarray_info_r = emicorr.get_subarcase(emicorr_model, subarray, readpatt, detector)
    watcher.assert_seen()

    assert subarray_info_r == (subarray, None, None, None)


def test_get_frequency_info(emicorr_model):
    freqname = "Hz10"
    expected = (10.039216, np.full(20, 1.5))
    freq, pa = emicorr.get_frequency_info(emicorr_model, freqname)
    assert freq == expected[0]
    assert np.all(pa == expected[1])

    freqname = "Hz390"
    expected = (390.625, np.full(20, 1.1))
    freq, pa = emicorr.get_frequency_info(emicorr_model, freqname)
    assert freq == expected[0]
    assert np.all(pa == expected[1])

    freqname = "Hz218"
    with pytest.raises(AttributeError, match=f"No attribute '{freqname}'"):
        emicorr.get_frequency_info(emicorr_model, freqname)


def test_sloper():
    data = np.ones((5, 5, 5))
    outarray, intercept = emicorr.sloper(data)

    compare_arr, compare_intercept = np.zeros((5, 5)), np.ones((5, 5))
    assert np.all(compare_arr == outarray)
    assert np.all(compare_intercept == intercept)


def test_minmed_ones():
    data = np.ones((5, 5, 5))
    compare_arr = data.copy()
    medimg = emicorr.minmed(data)
    assert np.all(compare_arr == medimg)


def test_minmed_2_arrays():
    shape = (5, 5)
    data = np.array([np.full(shape, 0.5), np.full(shape, 0.2)])

    # takes the minimum for two arrays
    medimg = emicorr.minmed(data)
    assert np.all(medimg == 0.2)


def test_minmed_3_arrays():
    shape = (5, 5)
    data = np.array([np.full(shape, 0.5), np.full(shape, 0.3), np.full(shape, 0.2)])

    # takes the median for >2 arrays
    medimg = emicorr.minmed(data)
    assert np.all(medimg == 0.3)


def test_minmed_ignore_nans():
    shape = (5, 5)
    data = np.array([np.full(shape, 0.5), np.full(shape, 0.3), np.full(shape, 0.2)])
    data[0, 0, 0] = np.nan
    data[1, 0, 1] = np.nan
    data[2, 0, 2] = np.nan

    # takes the median for >2 arrays, ignores NaNs
    medimg = emicorr.minmed(data)

    assert medimg[0, 0] == 0.25
    assert medimg[0, 1] == 0.35
    assert medimg[0, 2] == 0.4
    assert np.all(medimg.flat[3:] == 0.3)


def test_rebin_shrink():
    data = np.ones(10)
    data[1] = 0.55
    data[5] = 1.55
    data[9] = 2.0
    rebinned_data = emicorr.rebin(data, (7,))
    compare_arr = np.array([1.0, 0.55, 1.0, 1.0, 1.55, 1.0, 1.0])
    assert np.all(compare_arr == rebinned_data)


def test_rebin_grow():
    data = np.ones(10)
    data[1] = 0.55
    data[5] = 1.55
    data[9] = 2.0
    rebinned_data = emicorr.rebin(data, (15,))
    compare_arr = np.array(
        [1.0, 1.0, 0.55, 1.0, 1.0, 1.0, 1.0, 1.0, 1.55, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0]
    )
    assert np.all(compare_arr == rebinned_data)


def test_rebin_same():
    data = np.ones(10)
    data[1] = 0.55
    data[5] = 1.55
    data[9] = 2.0
    rebinned_data = emicorr.rebin(data, (10,))
    assert np.all(data == rebinned_data)


def test_rebin_error():
    data = np.ones(10)
    with pytest.raises(ValueError, match="dimensions must match"):
        emicorr.rebin(data, (10, 10))
