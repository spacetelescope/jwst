import pytest
import numpy as np
from pathlib import Path

from stdatamodels.jwst.datamodels.dqflags import pixel
from stdatamodels.jwst import datamodels
from jwst.stpipe import Step
from jwst.assign_wcs import AssignWcsStep
from jwst.background.background_sub_wfss import (subtract_wfss_bkg, 
                                            _mask_from_source_cat,
                                            _sufficient_background_pixels,
                                            _ScalingFactorComputer)

BKG_SCALING = 0.123
DETECTOR_SHAPE = (2048, 2048)
INITIAL_NAN_FRACTION = 1e-4

@pytest.fixture(scope="module")
def known_bkg():
    """Make a simplified version of the reference background model data."""

    ny, nx = DETECTOR_SHAPE
    y, x = np.mgrid[:ny, :nx]
    gradient = x * y / (nx*ny)
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
    err = err_scaling*(1 + rng.normal(0, 1, DETECTOR_SHAPE)**2)

    # add NaNs
    num_nans = int(data.size * nan_fraction)
    nan_indices = np.unravel_index(rng.choice(data.size, num_nans), data.shape)
    data[nan_indices] = np.nan
    err[nan_indices] = np.nan
    original_data_mean = np.nanmean(data)

    # also add a small background to the data with same structure
    # as the known reference background to see if it will get removed
    data += known_bkg*BKG_SCALING

    return data, err, original_data_mean


@pytest.fixture(scope='module')
def make_wfss_datamodel(data_path, mock_data, known_bkg):

    """Generate WFSS Observation"""
    wcsinfo = {
        'dec_ref': -27.79156387419731,
        'ra_ref': 53.16247756038121,
        'roll_ref': 0.04254766236781744,
        'v2_ref': -290.1,
        'v3_ref': -697.5,
        'v3yangle': 0.56987,
        'vparity': -1}

    observation = {
        'date': '2023-01-05',
        'time': '8:59:37'}

    exposure = {
        'duration': 11.805952,
        'end_time': 58119.85416,
        'exposure_time': 11.776,
        'frame_time': 0.11776,
        'group_time': 0.11776,
        'groupgap': 0,
        'integration_time': 11.776,
        'nframes': 1,
        'ngroups': 8,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'sample_time': 10.0,
        'start_time': 58668.72509857639,
        'zero_frame': False}

    subarray = {'xsize': DETECTOR_SHAPE[0],
                'ysize': DETECTOR_SHAPE[1],
                'xstart': 1,
                'ystart': 1}

    instrument = {
        'filter_position': 1,
        'pupil_position': 1}

    image = datamodels.ImageModel(DETECTOR_SHAPE)

    image.meta.wcsinfo._instance.update(wcsinfo)
    image.meta.instrument._instance.update(instrument)
    image.meta.observation._instance.update(observation)
    image.meta.subarray._instance.update(subarray)
    image.meta.exposure._instance.update(exposure)

    image.data = mock_data[0]
    image.err = mock_data[1]
    image.original_data_mean = mock_data[2] #just add this here for convenience
    image.dq = np.isnan(image.data)

    image.meta.source_catalog = str(data_path / "test_cat.ecsv")

    return image


@pytest.fixture()
def bkg_file(tmp_cwd, make_wfss_datamodel, known_bkg):
    """Mock background reference file"""

    bkg_fname = "ref_bkg.fits"
    bkg_image = make_wfss_datamodel.copy()
    bkg_image.data = known_bkg
    bkg_image.save(tmp_cwd / Path(bkg_fname))
    
    return bkg_fname


def test_make_wfss_ref_bkg(make_wfss_ref_bkg):
    bkg_model = datamodels.open(make_wfss_ref_bkg)
    bkg = bkg_model.data
    import matplotlib.pyplot as plt
    plt.imshow(bkg, origin='lower')
    plt.show()


filter_list = ['F250M', 'F277W', 'F335M', 'F356W', 'F460M',
               'F356W', 'F410M', 'F430M', 'F444W']  # + ['F480M', 'F322W2', 'F300M']


@pytest.mark.parametrize("pupils", ['GRISMC', 'GRISMR'])
@pytest.mark.parametrize("filters", filter_list)
@pytest.mark.parametrize("detectors", ['NRCALONG', 'NRCBLONG'])
def test_nrc_wfss_background(tmp_cwd, filters, pupils, detectors, make_wfss_datamodel, bkg_file):
    """Test background subtraction for NIRCAM WFSS modes."""
    data = make_wfss_datamodel

    data.meta.instrument.filter = filters
    data.meta.instrument.pupil = pupils
    data.meta.instrument.detector = detectors
    data.meta.instrument.channel = 'LONG'
    data.meta.instrument.name = 'NIRCAM'
    data.meta.exposure.type = 'NRC_WFSS'

    if data.meta.instrument.detector == 'NRCALONG':
        data.meta.instrument.module = 'A'
    elif data.meta.instrument.detector == 'NRCBLONG':
        data.meta.instrument.module = 'B'

    # Get References
    wcs_corrected = AssignWcsStep.call(data)
    wavelenrange = Step().get_reference_file(wcs_corrected, "wavelengthrange")

    # do the subtraction
    result = subtract_wfss_bkg(wcs_corrected, bkg_file, wavelenrange)

    # ensure NaN fraction did not increase. Rejecting outliers during determination
    # of factor should not have carried over into result.
    nan_frac = np.sum(np.isnan(result.data))/result.data.size
    rtol = 1/(result.data.size*INITIAL_NAN_FRACTION)
    assert np.isclose(nan_frac, INITIAL_NAN_FRACTION, rtol=rtol)

    # re-mask data so "real" sources are ignored here
    mask = _mask_from_source_cat(result, wavelenrange)
    result.data[~mask] = np.nan

    # test that the background has been subtracted from the data to within some fraction of
    # the noise in the data. There's probably a first-principles way to determine the tolerance,
    # but this is ok for the purposes of this test.
    tol = 0.01*np.nanstd(result.data)
    assert np.isclose(np.nanmean(result.data), result.original_data_mean, atol=tol)


@pytest.mark.parametrize("filters", ['GR150C', 'GR150R'])
@pytest.mark.parametrize("pupils", ['F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W'])
def test_nis_wfss_background(filters, pupils, make_wfss_datamodel, bkg_file):
    """Test background subtraction for NIRISS WFSS modes."""
    data = make_wfss_datamodel

    data.meta.instrument.filter = filters
    data.meta.instrument.pupil = pupils
    data.meta.instrument.detector = 'NIS'
    data.meta.instrument.name = 'NIRISS'
    data.meta.exposure.type = 'NIS_WFSS'

    # Get References
    wcs_corrected = AssignWcsStep.call(data)
    wavelenrange = Step().get_reference_file(wcs_corrected, "wavelengthrange")

    # do the subtraction
    result = subtract_wfss_bkg(wcs_corrected, bkg_file, wavelenrange)

    # ensure NaN fraction did not increase. Rejecting outliers during determination
    # of factor should not have carried over into result.
    nan_frac = np.sum(np.isnan(result.data))/result.data.size
    rtol = 1/(result.data.size*INITIAL_NAN_FRACTION)
    assert np.isclose(nan_frac, INITIAL_NAN_FRACTION, rtol=rtol)

    # re-mask data so "real" sources are ignored here
    mask = _mask_from_source_cat(result, wavelenrange)
    data = result.data[mask]

    tol = 0.01*np.nanstd(result.data)
    assert np.isclose(np.nanmean(result.data), result.original_data_mean, atol=tol)


def test_sufficient_background_pixels():
    model = datamodels.ImageModel(data=np.zeros((2048, 2048)),
                                  dq=np.zeros((2048, 2048)))
    refpix_flags = pixel['DO_NOT_USE'] | pixel['REFERENCE_PIXEL']
    model.dq[:4, :] = refpix_flags
    model.dq[-4:, :] = refpix_flags
    model.dq[:, :4] = refpix_flags
    model.dq[:, -4:] = refpix_flags

    bkg_mask = np.ones((2048, 2048), dtype=bool)
    # With full array minux refpix available for bkg, should be sufficient
    assert _sufficient_background_pixels(model.dq, bkg_mask)

    bkg_mask[4: -4, :] = 0
    bkg_mask[:, 4: -4] = 0
    # Now mask out entire array, mocking full source coverage of detector -
    # no pixels should be available for bkg
    assert not _sufficient_background_pixels(model.dq, bkg_mask)


def test_weighted_mean(make_wfss_datamodel, bkg_file):
    
    sci = make_wfss_datamodel.data
    var = make_wfss_datamodel.err**2
    with datamodels.open(bkg_file) as bkg_model:
        bkg = bkg_model.data

    # instantiate scaling factor computer
    rescaler = _ScalingFactorComputer()

    # just get the weighted mean without iteration
    factor = rescaler.err_weighted_mean(sci, bkg, var)
    original_data_mean = make_wfss_datamodel.original_data_mean
    expected_factor = BKG_SCALING+original_data_mean
    assert np.isclose(factor, expected_factor, atol=1e-3)

    # ensure it still works after iteration
    for niter in [1,2,5]:
        for p in [2, 0.5, 0.1]:
            rescaler = rescaler = _ScalingFactorComputer(p=p, maxiter=niter)
            factor, mask_out = rescaler(sci, bkg, var)
            mask_fraction = np.sum(mask_out)/mask_out.size
            max_mask_fraction = p*niter*2 + INITIAL_NAN_FRACTION

            assert np.isclose(factor, expected_factor, atol=1e-3)
            assert mask_fraction <= max_mask_fraction
            assert mask_fraction > INITIAL_NAN_FRACTION
    
    # TODO: test that this works with variance stopping criterion
    # TODO: test invalid inputs
    # TODO: test passing in all the different options

