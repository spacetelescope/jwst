import os

import pytest
import shutil
import json
import numpy as np
from pathlib import Path

from stdatamodels.jwst.datamodels.dqflags import pixel
from stdatamodels.jwst import datamodels
from jwst.stpipe import Step
from jwst.assign_wcs import AssignWcsStep
from jwst.background import BackgroundStep
from jwst.background.background_sub_wfss import (
    subtract_wfss_bkg,
    _mask_from_source_cat,
    _sufficient_background_pixels,
    _ScalingFactorComputer,
)

BKG_SCALING = 0.123
DETECTOR_SHAPE = (2048, 2048)
INITIAL_NAN_FRACTION = 1e-4
INITIAL_OUTLIER_FRACTION = 1e-3


@pytest.fixture(scope="module")
def known_bkg():
    """Make a simplified version of the reference background model data."""

    ny, nx = DETECTOR_SHAPE
    y, x = np.mgrid[:ny, :nx]
    gradient = x * y / (nx * ny)
    gradient = gradient - np.mean(gradient)
    return gradient + 1


@pytest.fixture(scope="module")
def mock_data(known_bkg):
    """Synthetic data with NaNs, noise, and the known background structure
    but rescaled. Later tests will ensure we can retrieve the proper scaling."""

    err_scaling = 0.05
    nan_fraction = INITIAL_NAN_FRACTION

    # make random data and error arrays
    rng = np.random.default_rng(seed=42)
    data = rng.normal(0, 1, DETECTOR_SHAPE)
    # ensure all errors are positive and not too close to zero
    err = err_scaling * (1 + rng.normal(0, 1, DETECTOR_SHAPE) ** 2)

    # add NaNs
    num_nans = int(data.size * nan_fraction)
    nan_indices = np.unravel_index(rng.choice(data.size, num_nans), data.shape)
    data[nan_indices] = np.nan
    err[nan_indices] = np.nan
    original_data_mean = np.nanmean(data)

    # add some outliers
    num_outliers = int(data.size * INITIAL_OUTLIER_FRACTION)
    outlier_indices = np.unravel_index(rng.choice(data.size, num_outliers), data.shape)
    data[outlier_indices] = rng.normal(100, 1, num_outliers)

    data[nan_indices] = np.nan
    err[nan_indices] = np.nan

    # also add a small background to the data with same structure
    # as the known reference background to see if it will get removed
    data += known_bkg * BKG_SCALING

    return data, err, original_data_mean


@pytest.fixture(scope="module")
def make_wfss_datamodel(data_path, mock_data):
    """Generate WFSS Observation"""
    wcsinfo = {
        "dec_ref": -27.79156387419731,
        "ra_ref": 53.16247756038121,
        "roll_ref": 0.04254766236781744,
        "v2_ref": -290.1,
        "v3_ref": -697.5,
        "v3yangle": 0.56987,
        "vparity": -1,
    }

    observation = {"date": "2023-01-05", "time": "8:59:37"}

    exposure = {
        "duration": 11.805952,
        "end_time": 58119.85416,
        "exposure_time": 11.776,
        "frame_time": 0.11776,
        "group_time": 0.11776,
        "groupgap": 0,
        "integration_time": 11.776,
        "nframes": 1,
        "ngroups": 8,
        "nints": 1,
        "nresets_between_ints": 0,
        "nsamples": 1,
        "sample_time": 10.0,
        "start_time": 58668.72509857639,
        "zero_frame": False,
    }

    subarray = {"xsize": DETECTOR_SHAPE[0], "ysize": DETECTOR_SHAPE[1], "xstart": 1, "ystart": 1}

    instrument = {"filter_position": 1, "pupil_position": 1}

    image = datamodels.ImageModel(DETECTOR_SHAPE)

    image.meta.wcsinfo._instance.update(wcsinfo)
    image.meta.instrument._instance.update(instrument)
    image.meta.observation._instance.update(observation)
    image.meta.subarray._instance.update(subarray)
    image.meta.exposure._instance.update(exposure)

    image.data = mock_data[0]
    image.err = mock_data[1]
    image.original_data_mean = mock_data[2]  # just add this here for convenience
    image.dq = np.isnan(image.data)

    image.meta.source_catalog = str(data_path / "test_cat.ecsv")

    return image


@pytest.fixture
def make_nrc_wfss_datamodel(make_wfss_datamodel):
    """Make a NIRCAM WFSS datamodel and call AssignWCS to populate its WCS"""
    data = make_wfss_datamodel.copy()
    data.meta.instrument.filter = "F250M"
    data.meta.instrument.pupil = "GRISMC"
    data.meta.instrument.detector = "NRCALONG"
    data.meta.instrument.channel = "LONG"
    data.meta.instrument.name = "NIRCAM"
    data.meta.exposure.type = "NRC_WFSS"
    data.meta.instrument.module = "A"
    result = AssignWcsStep.call(data)

    return result


@pytest.fixture
def make_nis_wfss_datamodel(make_wfss_datamodel):
    """Make a NIRISS WFSS datamodel and call AssignWCS to populate its WCS"""
    data = make_wfss_datamodel.copy()
    data.meta.instrument.filter = "GR150C"
    data.meta.instrument.pupil = "F090W"
    data.meta.instrument.detector = "NIS"
    data.meta.instrument.name = "NIRISS"
    data.meta.exposure.type = "NIS_WFSS"
    result = AssignWcsStep.call(data)

    return result


@pytest.fixture
def make_nis_wfss_sub64(make_wfss_datamodel):
    """Make a NIRISS WFSS datamodel with subarray 64x2048"""
    model = make_wfss_datamodel.copy()
    model.meta.instrument.filter = "GR150C"
    model.meta.instrument.pupil = "F090W"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.name = "NIRISS"
    model.meta.exposure.type = "NIS_WFSS"

    model.meta.subarray.xsize = 2048
    model.meta.subarray.ysize = 64
    model.meta.subarray.xstart = 1
    model.meta.subarray.ystart = 1985
    model.meta.subarray.name = "WFSS64C"

    model.data = model.data[-64:, :2048]  # simulate subarray by slicing data
    model.err = model.err[-64:, :2048]
    model.dq = model.dq[-64:, :2048]

    result = AssignWcsStep.call(model)

    return result


@pytest.fixture()
def bkg_file(tmp_cwd, make_wfss_datamodel, known_bkg):
    """Mock background reference file"""

    bkg_fname = "ref_bkg.fits"
    bkg_image = make_wfss_datamodel.copy()
    bkg_image.data = known_bkg
    bkg_image.save(tmp_cwd / Path(bkg_fname))

    return bkg_fname


def shared_tests(sci, mask, original_data_mean):
    """
    Tests that are common to all WFSS modes.

    Note that NaN fraction test in test_nrc_wfss_background and test_nis_wfss_background
    cannot be applied to the full run tests because the background reference files contain
    NaNs in some cases (specifically for NIRISS)
    """

    # re-mask data so "real" sources are ignored here
    sci[~mask] = np.nan

    # test that the background has been subtracted from the data to within some fraction of
    # the noise in the data. There's probably a first-principles way to determine the tolerance,
    # but this is ok for the purposes of this test.
    # ignore the outliers here too
    sci[sci > 50] = np.nan
    tol = 0.01 * np.nanstd(sci)
    assert np.isclose(np.nanmean(sci), original_data_mean, atol=tol)


def test_nrc_wfss_background(make_nrc_wfss_datamodel, bkg_file):
    """Test background subtraction for NIRCAM WFSS modes."""
    data = make_nrc_wfss_datamodel.copy()

    # Get References
    wavelenrange = Step().get_reference_file(data, "wavelengthrange")

    # do the subtraction
    result = subtract_wfss_bkg(data, bkg_file, wavelenrange)
    sci = result.data.copy()

    # ensure NaN fraction did not increase. Rejecting outliers during determination
    # of factor should not have carried over into result.
    nan_frac = np.sum(np.isnan(sci)) / sci.size
    assert np.isclose(nan_frac, INITIAL_NAN_FRACTION, rtol=1e-2)

    # re-compute mask to ignore "real" sources for tests
    mask = _mask_from_source_cat(result, wavelenrange)

    shared_tests(sci, mask, data.original_data_mean)


@pytest.mark.parametrize("subarray", [None, "WFSS64C"])
def test_nis_wfss_background(subarray, make_nis_wfss_datamodel, make_nis_wfss_sub64, bkg_file):
    """Test background subtraction for NIRISS WFSS modes."""
    if subarray == "WFSS64C":
        data = make_nis_wfss_sub64.copy()
    else:
        data = make_nis_wfss_datamodel.copy()

    # Get References
    wavelenrange = Step().get_reference_file(data, "wavelengthrange")

    # do the subtraction
    result = subtract_wfss_bkg(data, bkg_file, wavelenrange)
    sci = result.data.copy()

    # ensure NaN fraction did not increase. Rejecting outliers during determination
    # of factor should not have carried over into result.
    nan_frac = np.sum(np.isnan(sci)) / sci.size
    assert np.isclose(nan_frac, INITIAL_NAN_FRACTION, rtol=1e-2)

    mask = _mask_from_source_cat(result, wavelenrange)
    shared_tests(sci, mask, data.original_data_mean)


# test both filters because they have opposite dispersion directions
@pytest.mark.parametrize("pupil", ["GRISMC", "GRISMR"])
def test_nrc_wfss_full_run(pupil, make_nrc_wfss_datamodel):
    """
    Test full run of NIRCAM WFSS background subtraction.

    The residual structure in the background will not look as nice as in
    test_nis_wfss_background because here it's taken from a reference file,
    so the bkg has real detector imperfections
    while the data is synthetic and just has a mock gradient
    """
    data = make_nrc_wfss_datamodel.copy()
    data.meta.instrument.pupil = pupil

    # do the subtraction. set all options to ensure they are at least recognized
    result = BackgroundStep.call(
        data,
        None,
        wfss_maxiter=3,
        wfss_outlier_percent=0.5,
        wfss_rms_stop=0,
    )

    sci = result.data.copy()
    # re-derive mask to ignore "real" sources for tests
    wavelenrange = Step().get_reference_file(data, "wavelengthrange")
    mask = _mask_from_source_cat(result, wavelenrange)
    shared_tests(sci, mask, data.original_data_mean)
    assert isinstance(result.meta.background.scaling_factor, float)


@pytest.mark.parametrize("filt", ["GR150C", "GR150R"])
def test_nis_wfss_full_run(filt, make_nis_wfss_datamodel):
    """
    Test full run of NIRISS WFSS background subtraction.

    The residual structure in the background will not look as nice as in
    test_nis_wfss_background because here it's taken from a reference file,
    so the bkg has real detector imperfections
    while the data is synthetic and just has a mock gradient
    """
    data = make_nis_wfss_datamodel.copy()
    data.meta.instrument.filter = filt

    # do the subtraction. set all options to ensure they are at least recognized
    result = BackgroundStep.call(
        data,
        None,
        wfss_maxiter=3,
        wfss_outlier_percent=0.5,
        wfss_rms_stop=0,
    )

    sci = result.data.copy()
    # re-derive mask to ignore "real" sources for tests
    wavelenrange = Step().get_reference_file(data, "wavelengthrange")
    mask = _mask_from_source_cat(result, wavelenrange)
    shared_tests(sci, mask, data.original_data_mean)
    assert isinstance(result.meta.background.scaling_factor, float)


def test_sufficient_background_pixels():
    model = datamodels.ImageModel(data=np.zeros((2048, 2048)), dq=np.zeros((2048, 2048)))
    refpix_flags = pixel["DO_NOT_USE"] | pixel["REFERENCE_PIXEL"]
    model.dq[:4, :] = refpix_flags
    model.dq[-4:, :] = refpix_flags
    model.dq[:, :4] = refpix_flags
    model.dq[:, -4:] = refpix_flags

    bkg = np.ones((2048, 2048), dtype=float)

    bkg_mask = np.ones((2048, 2048), dtype=bool)
    # With full array minus refpix available for bkg, should be sufficient
    assert _sufficient_background_pixels(model.dq, bkg_mask, bkg)

    bkg_mask[4:-4, :] = 0
    bkg_mask[:, 4:-4] = 0
    # Now mask out entire array, mocking full source coverage of detector -
    # no pixels should be available for bkg
    assert not _sufficient_background_pixels(model.dq, bkg_mask, bkg)


def test_sufficient_background_pixels_nonoverlapping():
    """Valid pixels in bkg, and valid pixels in bkg_mask, but these do not overlap."""
    dq = np.zeros((2048, 2048), dtype=int)
    bkg_mask = np.zeros((2048, 2048), dtype=bool)
    bkg = np.zeros((2048, 2048), dtype=float)
    bkg_mask[:1024, :1024] = 1
    bkg[1024:, 1024:] = 1.0
    assert not _sufficient_background_pixels(dq, bkg_mask, bkg)


def test_weighted_mean(make_wfss_datamodel, bkg_file):
    sci = make_wfss_datamodel.data
    var = make_wfss_datamodel.err**2
    with datamodels.open(bkg_file) as bkg_model:
        bkg = bkg_model.data

    # put 0.1% zero values in variance to ensure coverage of previous bug where zero-valued
    # variances in real data caused factor = 1/np.inf = 0
    rng = np.random.default_rng(seed=42)
    n_bad = int(var.size / 1000)
    bad_i = rng.choice(var.size - 1, n_bad)
    var[np.unravel_index(bad_i, var.shape)] = 0.0

    # instantiate scaling factor computer
    rescaler = _ScalingFactorComputer()

    # just get the weighted mean without iteration
    # to check it's as expected, mask outliers
    sci[sci > 50] = np.nan
    factor = rescaler.err_weighted_mean(sci, bkg, var)
    original_data_mean = make_wfss_datamodel.original_data_mean
    expected_factor = BKG_SCALING + original_data_mean
    assert np.isclose(factor, expected_factor, atol=1e-3)

    # ensure it still works after iteration
    for niter in [1, 2, 5]:
        for p in [2, 0.5, 0.1]:
            rescaler = _ScalingFactorComputer(p=p, maxiter=niter)
            assert (
                rescaler.delta_rms_thresh == 0
            )  # check rms_thresh=None input sets thresh properly

            factor, mask_out = rescaler(sci, bkg, var)
            mask_fraction = np.sum(mask_out) / mask_out.size
            max_mask_fraction = p * niter * 2 + INITIAL_NAN_FRACTION

            assert np.isclose(factor, expected_factor, atol=1e-3)
            assert mask_fraction <= max_mask_fraction
            assert mask_fraction > INITIAL_NAN_FRACTION

    # test that variance stopping criterion works
    # tune the RMS thresh to take roughly half the iterations
    # need lots of significant digits here because iterating makes little difference
    # for this test case
    maxiter = 10
    delta_rms_thresh = 1e-4
    p = 100 * INITIAL_OUTLIER_FRACTION / 2
    rescaler = _ScalingFactorComputer(
        p=p, dispersion_axis=1, delta_rms_thresh=delta_rms_thresh, maxiter=maxiter
    )
    factor, mask_out = rescaler(sci, bkg, var)
    assert rescaler._iters_run_last_call < maxiter

    # test putting mask=None works ok, and that maxiter=0 just gives you err weighted mean
    rescaler = _ScalingFactorComputer(maxiter=0)
    factor, mask_out = rescaler(sci, bkg, var)
    assert np.all(mask_out == 0)
    assert factor == rescaler.err_weighted_mean(sci, bkg, var)

    # test invalid inputs
    with pytest.raises(ValueError):
        rescaler = _ScalingFactorComputer(dispersion_axis=5, delta_rms_thresh=1)

    with pytest.raises(ValueError):
        rescaler = _ScalingFactorComputer(dispersion_axis=None, delta_rms_thresh=1)


@pytest.fixture()
def mock_asn_and_data(tmp_path_factory, data_path, make_nis_wfss_datamodel):
    # Create temp dir and copy the catalog in there
    tmp_path = tmp_path_factory.mktemp("asn_input")
    shutil.copy(str(data_path / "test_cat.ecsv"), str(tmp_path / "test_cat.ecsv"))
    # Save the datmodel into a rate file but remove the catalog to make sure it is
    # added back in by the asn_intake module
    make_nis_wfss_datamodel.meta.source_catalog = None
    ratefile = tmp_path / "jw01000001001_test_00001_nis_rate.fits"
    make_nis_wfss_datamodel.save(str(ratefile))
    # Pretend this is also the direct image, save with a diff name
    i2dfile = tmp_path / "jw01000-o001_t001_niriss_i2d.fits"
    make_nis_wfss_datamodel.save(str(i2dfile))
    # Pretend this is also the segment, save with a diff name
    segmfile = tmp_path / "jw01000-o001_t001_niriss_segm.fits"
    make_nis_wfss_datamodel.save(str(segmfile))

    data = {
        "asn_type": "spec2",
        "asn_rule": "Asn_Lv2WFSS",
        "program": "01000",
        "asn_pool": "jw010000_pool.csv",
        "products": [
            {
                "name": "jw01000001001_test_00001_nis",
                "members": [
                    {
                        "expname": "jw01000001001_test_00001_nis_rate.fits",
                        "exptype": "science",
                        "exposerr": "null",
                    },
                    {"expname": "jw01000-o001_t001_niriss_i2d.fits", "exptype": "direct_image"},
                    {"expname": "test_cat.ecsv", "exptype": "sourcecat"},
                    {"expname": "jw01000-o001_t001_niriss_segm.fits", "exptype": "segmap"},
                ],
            }
        ],
    }

    asn_name = str(tmp_path / "jw010000-wfss_test_spec2_00001_asn.json")
    with open(asn_name, "w") as asn:
        json.dump(data, asn)

    return [asn_name, ratefile, i2dfile, segmfile]


def test_wfss_asn_input(mock_asn_and_data):
    # get the file name of asn and other file objects
    asn_name, ratefile = mock_asn_and_data[0], mock_asn_and_data[1]
    i2dfile, segmfile = mock_asn_and_data[2], mock_asn_and_data[3]
    # change the working directory into the temp so it can find all files
    cwd = os.getcwd()
    os.chdir(ratefile.parents[0])
    result = BackgroundStep.call(asn_name)
    # return to previous working dir
    os.chdir(cwd)

    assert result.meta.source_catalog == "test_cat.ecsv"
    assert result.meta.direct_image == i2dfile.name
    assert result.meta.segmentation_map == segmfile.name
    assert result.meta.cal_step.bkg_subtract == "COMPLETE"
