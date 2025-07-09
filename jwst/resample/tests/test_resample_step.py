import pytest
import warnings
from copy import deepcopy

import asdf
import numpy as np
from astropy.io import fits

from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from stdatamodels.jwst.datamodels import ImageModel, dqflags
from stcal.resample.utils import compute_mean_pixel_area, build_driz_weight

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.assign_wcs import AssignWcsStep
from jwst.assign_wcs.util import compute_scale
from jwst.exp_to_source import multislit_to_container
from jwst.extract_2d import Extract2dStep
from jwst.resample import ResampleSpecStep, ResampleStep
from jwst.resample.resample import input_jwst_model_to_dict
from jwst.resample.resample_step import GOOD_BITS
from jwst.resample.resample_spec import ResampleSpec, compute_spectral_pixel_scale
from jwst.resample.resample_utils import load_custom_wcs


_FLT32_EPS = np.finfo(np.float32).eps


def _set_photom_kwd(im):
    xmin = im.meta.subarray.xstart - 1
    xmax = xmin + im.meta.subarray.xsize
    ymin = im.meta.subarray.ystart - 1
    ymax = ymin + im.meta.subarray.ysize

    im.meta.wcs.array_shape = im.data.shape

    if im.meta.wcs.bounding_box is None:
        bb = ((xmin - 0.5, xmax - 0.5), (ymin - 0.5, ymax - 0.5))
        im.meta.wcs.bounding_box = bb

    mean_pixel_area = compute_mean_pixel_area(
        im.meta.wcs,
        shape=im.data.shape,
    )

    if mean_pixel_area:
        im.meta.photometry.pixelarea_steradians = mean_pixel_area
        im.meta.photometry.pixelarea_arcsecsq = mean_pixel_area * np.rad2deg(3600) ** 2


def miri_rate_model():
    xsize = 72
    ysize = 416
    sregion = (
        "POLYGON ICRS  10.355323877 -22.353560934 10.355437846 -22.353464295 "
        + "10.354477543 -22.352498313 10.354363599 -22.352595345"
    )

    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.data += 5
    im.var_rnoise += 1
    im.meta.wcsinfo = {
        "dec_ref": 40,
        "ra_ref": 100,
        "roll_ref": 0.0,
        "v2_ref": -453.5134,
        "v3_ref": -373.4826,
        "v3yangle": 0.0,
        "vparity": -1,
        "s_region": sregion,
    }
    im.meta.instrument = {"detector": "MIRIMAGE", "filter": "P750L", "name": "MIRI"}
    im.meta.observation = {"date": "2019-01-01", "time": "17:00:00"}
    im.meta.subarray = {
        "fastaxis": 1,
        "name": "SLITLESSPRISM",
        "slowaxis": 2,
        "xsize": xsize,
        "xstart": 1,
        "ysize": ysize,
        "ystart": 529,
    }
    im.meta.exposure = {
        "duration": 11.805952,
        "end_time": 58119.85416,
        "exposure_time": 11.776,
        "measurement_time": 11.65824,
        "frame_time": 0.11776,
        "group_time": 0.11776,
        "groupgap": 0,
        "integration_time": 11.776,
        "nframes": 1,
        "ngroups": 100,
        "nints": 1,
        "nresets_between_ints": 0,
        "nsamples": 1,
        "readpatt": "FAST",
        "sample_time": 10.0,
        "start_time": 58119.8333,
        "type": "MIR_LRS-SLITLESS",
        "zero_frame": False,
    }
    return im


@pytest.fixture
def miri_rate():
    im = miri_rate_model()
    yield im
    im.close()


@pytest.fixture
def miri_cal(miri_rate):
    im = AssignWcsStep.call(miri_rate)
    _set_photom_kwd(im)

    # Add non-zero values to check flux conservation
    im.data += 1.0

    yield im
    im.close()


@pytest.fixture
def miri_rate_zero_crossing():
    xsize = 1032
    ysize = 1024
    sregion = (
        "POLYGON ICRS  10.355323877 -22.353560934 10.355437846 -22.353464295 "
        + "10.354477543 -22.352498313 10.354363599 -22.352595345"
    )
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise = np.random.random(shape)
    im.meta.wcsinfo = {
        "dec_ref": 2.16444343946559e-05,
        "ra_ref": -0.00026031780056776,
        "roll_ref": 0.0,
        "v2_ref": -415.0690466121227,
        "v3_ref": -400.575920398547,
        "v3yangle": 0.0,
        "vparity": -1,
        "s_region": sregion,
    }
    im.meta.instrument = {"detector": "MIRIMAGE", "filter": "P750L", "name": "MIRI"}
    im.meta.observation = {"date": "2019-01-01", "time": "17:00:00"}
    im.meta.subarray = {
        "fastaxis": 1,
        "name": "FULL",
        "slowaxis": 2,
        "xsize": xsize,
        "xstart": 1,
        "ysize": ysize,
        "ystart": 1,
    }
    im.meta.exposure = {
        "duration": 11.805952,
        "end_time": 58119.85416,
        "exposure_time": 11.776,
        "frame_time": 0.11776,
        "group_time": 0.11776,
        "groupgap": 0,
        "integration_time": 11.776,
        "nframes": 1,
        "ngroups": 100,
        "nints": 1,
        "nresets_between_ints": 0,
        "nsamples": 1,
        "readpatt": "FAST",
        "sample_time": 10.0,
        "start_time": 58119.8333,
        "type": "MIR_LRS-FIXEDSLIT",
        "zero_frame": False,
    }

    yield im
    im.close()


@pytest.fixture
def miri_rate_pair(miri_rate_zero_crossing):
    im1 = miri_rate_zero_crossing
    # Create a nodded version
    im2 = im1.copy()
    im2.meta.wcsinfo.ra_ref = 0.00026308279776455
    im2.meta.wcsinfo.dec_ref = -2.1860888891293e-05
    im1 = AssignWcsStep.call(im1)
    im2 = AssignWcsStep.call(im2)

    yield im1, im2
    im1.close()
    im2.close()


@pytest.fixture
def nircam_rate():
    xsize = 204
    ysize = 204
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise += 0
    im.meta.wcsinfo = {
        "ctype1": "RA---TAN",
        "ctype2": "DEC--TAN",
        "dec_ref": 11.99875540218638,
        "ra_ref": 22.02351763251896,
        "roll_ref": 0.005076934167039675,
        "v2_ref": 86.039011,
        "v3_ref": -493.385704,
        "v3yangle": -0.07385127,
        "vparity": -1,
        "wcsaxes": 2,
    }
    im.meta.instrument = {
        "channel": "LONG",
        "detector": "NRCALONG",
        "filter": "F444W",
        "lamp_mode": "NONE",
        "module": "A",
        "name": "NIRCAM",
        "pupil": "CLEAR",
    }
    im.meta.subarray = {
        "fastaxis": -1,
        "name": "FULL",
        "slowaxis": 2,
        "xsize": xsize,
        "xstart": 1,
        "ysize": ysize,
        "ystart": 1,
    }
    im.meta.observation = {
        "activity_id": "01",
        "date": "2021-10-25",
        "exposure_number": "00001",
        "obs_id": "V42424001001P0000000001101",
        "observation_label": "nircam_ptsrc_only",
        "observation_number": "001",
        "program_number": "42424",
        "sequence_id": "1",
        "time": "16:58:27.258",
        "visit_group": "01",
        "visit_id": "42424001001",
        "visit_number": "001",
    }
    im.meta.exposure = {
        "duration": 161.05155,
        "end_time": 59512.70899968495,
        "exposure_time": 150.31478,
        "measurement_time": 139.57801,
        "frame_time": 10.73677,
        "group_time": 21.47354,
        "groupgap": 1,
        "integration_time": 150.31478,
        "mid_time": 59512.70812980775,
        "nframes": 1,
        "ngroups": 7,
        "nints": 1,
        "nresets_at_start": 1,
        "nresets_between_ints": 1,
        "readpatt": "BRIGHT1",
        "sample_time": 10,
        "start_time": 59512.70725993055,
        "type": "NRC_IMAGE",
    }
    im.meta.photometry = {
        "pixelarea_steradians": 1e-13,
        "pixelarea_arcsecsq": 4e-3,
    }
    yield im
    im.close()


@pytest.fixture
def nirspec_rate():
    ysize = 2048
    xsize = 2048
    shape = (ysize, xsize)
    im = ImageModel(shape)
    im.var_rnoise += 1
    im.meta.target = {"ra": 100.1237, "dec": 39.86}
    im.meta.wcsinfo = {
        "dec_ref": 40,
        "ra_ref": 100,
        "roll_ref": 0,
        "v2_ref": -453.5134,
        "v3_ref": -373.4826,
        "v3yangle": 0.0,
        "vparity": -1,
    }
    im.meta.instrument = {
        "detector": "NRS1",
        "filter": "CLEAR",
        "grating": "PRISM",
        "name": "NIRSPEC",
        "gwa_tilt": 37.0610,
        "gwa_xtilt": 0.0001,
        "gwa_ytilt": 0.0001,
        "fixed_slit": "S200A1",
    }
    im.meta.subarray = {
        "fastaxis": 1,
        "name": "SUBS200A1",
        "slowaxis": 2,
        "xsize": 72,
        "xstart": 1,
        "ysize": 416,
        "ystart": 529,
    }
    im.meta.observation = {"program_number": "1234", "date": "2016-09-05", "time": "8:59:37"}
    im.meta.exposure = {
        "duration": 11.805952,
        "end_time": 58119.85416,
        "exposure_time": 11.776,
        "measurement_time": 11.65824,
        "frame_time": 0.11776,
        "group_time": 0.11776,
        "groupgap": 0,
        "integration_time": 11.776,
        "nframes": 1,
        "ngroups": 100,
        "nints": 1,
        "nresets_between_ints": 0,
        "nsamples": 1,
        "readpatt": "NRSRAPID",
        "sample_time": 10.0,
        "start_time": 58119.8333,
        "type": "NRS_FIXEDSLIT",
        "zero_frame": False,
    }

    yield im
    im.close()


@pytest.fixture
def nirspec_cal(nirspec_rate):
    im = AssignWcsStep.call(nirspec_rate)

    # Since the ra_targ, and dec_targ are flux-weighted, we need non-zero
    # flux values.
    im.data += 1.0

    im = Extract2dStep.call(im)
    for slit in im.slits:
        _set_photom_kwd(slit)
    yield im
    im.close()


@pytest.fixture
def nirspec_cal_pair(nirspec_rate):
    # copy the rate model to make files with different filters
    rate1 = nirspec_rate
    rate1.meta.instrument.grating = "G140H"
    rate1.meta.instrument.filter = "F070LP"
    rate2 = nirspec_rate.copy()
    rate2.meta.instrument.grating = "G140H"
    rate2.meta.instrument.filter = "F100LP"

    im1 = AssignWcsStep.call(nirspec_rate)
    im2 = AssignWcsStep.call(rate2)

    # Since the ra_targ, and dec_targ are flux-weighted, we need non-zero
    # flux values.
    im1.data += 1.0
    im2.data += 1.0

    im1 = Extract2dStep.call(im1)
    im2 = Extract2dStep.call(im2)
    for slit in im1.slits:
        _set_photom_kwd(slit)
    for slit in im2.slits:
        _set_photom_kwd(slit)
    yield im1, im2
    im1.close()
    im2.close()


@pytest.fixture
def nirspec_lamp(nirspec_rate):
    nirspec_rate.meta.exposure.type = "NRS_LAMP"
    nirspec_rate.meta.instrument.lamp_mode = "FIXEDSLIT"
    nirspec_rate.meta.instrument.lamp_state = "FLAT"
    im = AssignWcsStep.call(nirspec_rate)
    im.data += 1.0

    im = Extract2dStep.call(im)
    for slit in im.slits:
        _set_photom_kwd(slit)

    yield im
    im.close()


def test_nirspec_wcs_roundtrip(nirspec_cal):
    im = ResampleSpecStep.call(nirspec_cal)

    for slit in im.slits:
        x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
        ra, dec, lam = slit.meta.wcs(x, y)
        xp, yp = slit.meta.wcs.invert(ra, dec, lam)

        assert_allclose(x, xp, rtol=0, atol=1e-8)
        assert_allclose(y, yp, rtol=0, atol=3e-4)
    im.close()


def test_nirspec_lamp_wcs_roundtrip(nirspec_lamp):
    # RuntimeWarning: invalid value encountered in sqrt
    with np.errstate(invalid="ignore"):
        im = ResampleSpecStep.call(nirspec_lamp)

    for slit in im.slits:
        x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
        ra, dec, lam = slit.meta.wcs(x, y)
        xp, yp = slit.meta.wcs.invert(ra, dec, lam)

        assert_allclose(x, xp, rtol=0, atol=1e-8)
        assert_allclose(y, yp, rtol=0, atol=3e-4)
    im.close()


def test_miri_wcs_roundtrip(miri_cal):
    im = ResampleSpecStep.call(miri_cal)
    x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
    ra, dec, lam = im.meta.wcs(x, y)
    xp, yp = im.meta.wcs.invert(ra, dec, lam)

    assert_allclose(x, xp, atol=1e-8)
    assert_allclose(y, yp, atol=1e-8)
    im.close()


def test_single_image_file_input(nircam_rate, tmp_cwd):
    """Ensure step can be run on a single image file."""
    # Create a temporary file with the input data
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    im.meta.filename = "test_input.fits"
    im.save("test_input.fits")

    # Run the step on the file
    result_from_memory = ResampleStep.call(im)
    result_from_file = ResampleStep.call("test_input.fits")

    # Check that the output is as expected
    assert_allclose(result_from_file.data, result_from_memory.data, equal_nan=True)

    result_from_file.close()
    result_from_memory.close()
    im.close()


@pytest.mark.parametrize("ratio", [0.5, 0.7, 1.0])
def test_pixel_scale_ratio_imaging(nircam_rate, ratio):
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    _set_photom_kwd(im)
    im.data += 5
    result1 = ResampleStep.call(im)
    result2 = ResampleStep.call(im, pixel_scale_ratio=ratio)

    assert_allclose(
        np.array(result1.data.shape), np.array(result2.data.shape) * ratio, rtol=1, atol=1
    )

    # Make sure the photometry keywords describing the solid angle of a pixel
    # are updated
    area1 = result1.meta.photometry.pixelarea_steradians
    area2 = result2.meta.photometry.pixelarea_steradians
    assert_allclose(area1 * ratio**2, area2, rtol=1e-6)

    assert result1.meta.resample.pixel_scale_ratio == 1.0
    assert result2.meta.resample.pixel_scale_ratio == ratio

    im.close()
    result1.close()
    result2.close()


@pytest.mark.parametrize("units", ["MJy", "MJy/sr"])
@pytest.mark.parametrize("ratio", [0.7, 1.0, 1.3])
def test_pixel_scale_ratio_spec_miri(miri_cal, ratio, units):
    miri_cal.meta.bunit_data = units

    # Make an input pixel scale equivalent to the specified ratio
    input_scale = compute_spectral_pixel_scale(miri_cal.meta.wcs, disp_axis=2)
    pscale = 3600.0 * input_scale / ratio

    result1 = ResampleSpecStep.call(miri_cal)
    result2 = ResampleSpecStep.call(miri_cal, pixel_scale_ratio=ratio)
    result3 = ResampleSpecStep.call(miri_cal, pixel_scale=pscale)

    # pixel_scale and pixel_scale_ratio should be equivalent
    nn = np.isnan(result2.data) | np.isnan(result3.data)
    assert_allclose(result2.data[~nn], result3.data[~nn], rtol=1e-6)

    # Check result2 for expected results

    # wavelength size does not change
    assert result1.data.shape[0] == result2.data.shape[0]

    # spatial dimension is scaled
    assert_allclose(result1.data.shape[1], result2.data.shape[1] / ratio, atol=1)

    # data is non-trivial
    assert np.nansum(result1.data) > 0.0
    assert np.nansum(result2.data) > 0.0

    # flux is conserved
    if "sr" not in units:
        # flux density conservation: sum over pixels in each row
        # needs to be about the same, other than the edges
        # Check the maximum sums, to avoid edges.
        assert_allclose(
            np.max(np.nansum(result1.data, axis=1)),
            np.max(np.nansum(result2.data, axis=1)),
            rtol=0.05,
        )
    else:
        # surface brightness conservation: weighted total values are the same
        assert np.allclose(
            np.nansum(result1.data * result1.wht),
            np.nansum(result2.data * result2.wht),
            rtol=5.0e-3,
            equal_nan=True,
        )
        assert np.allclose(
            np.nansum(result1.data * result1.wht, axis=1),
            np.nansum(result2.data * result2.wht, axis=1),
            rtol=5.0e-3,
            equal_nan=True,
        )

    # output area is updated either way
    area1 = result1.meta.photometry.pixelarea_steradians
    area2 = result2.meta.photometry.pixelarea_steradians
    area3 = result2.meta.photometry.pixelarea_steradians
    assert_allclose(area1 / area2, ratio)
    assert_allclose(area1 / area3, ratio)

    assert result1.meta.resample.pixel_scale_ratio == 1.0
    assert result2.meta.resample.pixel_scale_ratio == ratio
    assert_allclose(result3.meta.resample.pixel_scale_ratio, ratio)

    result1.close()
    result2.close()
    result3.close()


@pytest.mark.parametrize("units", ["MJy", "MJy/sr"])
@pytest.mark.parametrize("ratio", [0.7, 1.0, 1.3])
def test_pixel_scale_ratio_1spec_miri_pair(miri_rate_pair, ratio, units):
    im1, im2 = miri_rate_pair
    _set_photom_kwd(im1)
    _set_photom_kwd(im2)
    im1.meta.bunit_data = units
    im2.meta.bunit_data = units
    im1.meta.filename = "file1.fits"
    im2.meta.filename = "file2.fits"
    im1.data += 1.0
    im2.data += 1.0

    # Make an input pixel scale equivalent to the specified ratio
    input_scale = compute_spectral_pixel_scale(im1.meta.wcs, disp_axis=2)
    pscale = 3600.0 * input_scale / ratio

    result1 = ResampleSpecStep.call([im1, im2])
    result2 = ResampleSpecStep.call([im1, im2], pixel_scale_ratio=ratio)
    result3 = ResampleSpecStep.call([im1, im2], pixel_scale=pscale)

    # pixel_scale and pixel_scale_ratio should be equivalent
    nn = ~(np.isnan(result2.data) | np.isnan(result3.data))
    assert_allclose(result2.data[nn], result3.data[nn], rtol=2.000001 * _FLT32_EPS)

    # Check result2 for expected results

    # wavelength size does not change
    assert result1.data.shape[0] == result2.data.shape[0]

    # spatial dimension is scaled
    assert_allclose(result1.data.shape[1], result2.data.shape[1] / ratio, rtol=0, atol=1)

    # data is non-trivial
    assert np.nansum(result1.data) > 0.0
    assert np.nansum(result2.data) > 0.0

    # flux is conserved
    if "sr" not in units:
        # flux density conservation: sum over pixels in each row
        # needs to be about the same, other than the edges
        # Check the maximum sums, to avoid edges.
        assert np.allclose(
            np.max(np.nansum(result1.data, axis=1)),
            np.max(np.nansum(result2.data, axis=1)),
            rtol=0.05,
        )
    else:
        # surface brightness conservation: weighted total values are the same
        assert np.allclose(
            np.sum(result1.data * result1.wht),
            np.sum(result2.data * result2.wht),
            rtol=1.0e-5,
            equal_nan=True,
        )

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=RuntimeWarning, message="All-NaN slice encountered"
            )
            assert np.allclose(
                np.nanmedian(result1.data, axis=1),
                np.nanmedian(result2.data, axis=1),
                rtol=1.0e-5,
                equal_nan=True,
            )

    # output area is updated either way
    area1 = result1.meta.photometry.pixelarea_steradians
    area2 = result2.meta.photometry.pixelarea_steradians
    area3 = result2.meta.photometry.pixelarea_steradians
    assert_allclose(area1 / area2, ratio)
    assert_allclose(area1 / area3, ratio)

    assert result1.meta.resample.pixel_scale_ratio == 1.0
    assert result2.meta.resample.pixel_scale_ratio == ratio
    assert_allclose(result3.meta.resample.pixel_scale_ratio, ratio)

    result1.close()
    result2.close()
    result3.close()


@pytest.mark.parametrize("units", ["MJy", "MJy/sr"])
@pytest.mark.parametrize("ratio", [0.7, 1.0, 1.3])
def test_pixel_scale_ratio_spec_nirspec(nirspec_cal, ratio, units):
    for slit in nirspec_cal.slits:
        slit.meta.bunit_data = units

    # Make an input pixel scale equivalent to the specified ratio
    input_scale = compute_spectral_pixel_scale(nirspec_cal.slits[0].meta.wcs, disp_axis=1)
    pscale = 3600.0 * input_scale / ratio

    result1 = ResampleSpecStep.call(nirspec_cal)
    result2 = ResampleSpecStep.call(nirspec_cal, pixel_scale_ratio=ratio)
    result3 = ResampleSpecStep.call(nirspec_cal, pixel_scale=pscale)

    for slit1, slit2, slit3 in zip(result1.slits, result2.slits, result3.slits):
        # pixel_scale and pixel_scale_ratio should be equivalent
        nn = np.isnan(slit2.data) | np.isnan(slit3.data)
        assert_allclose(slit2.data[~nn], slit3.data[~nn], rtol=1e-6)

        # Check result2 for expected results

        # wavelength size does not change
        assert slit1.data.shape[1] == slit2.data.shape[1]

        # spatial dimension is scaled
        assert_allclose(slit1.data.shape[0], slit2.data.shape[0] / ratio, atol=1)

        # data is non-trivial
        assert np.nansum(slit1.data) > 0.0
        assert np.nansum(slit2.data) > 0.0

        # flux is conserved
        if "sr" not in units:
            # flux density conservation: sum over pixels in each column
            # needs to be about the same, other than edge effects.
            # Check the maximum sums, to avoid edges.
            assert_allclose(
                np.max(np.nansum(slit1.data, axis=0)),
                np.max(np.nansum(slit2.data, axis=0)),
                rtol=0.05,
            )
        else:
            # surface brightness conservation: weighted total values are the same
            assert np.allclose(
                np.nansum(slit1.data * slit1.wht, axis=0),
                np.nansum(slit2.data * slit2.wht, axis=0),
                rtol=1.0e-5,
                equal_nan=True,
            )

        # output area is updated either way
        area1 = slit1.meta.photometry.pixelarea_steradians
        area2 = slit2.meta.photometry.pixelarea_steradians
        area3 = slit3.meta.photometry.pixelarea_steradians
        assert_allclose(area1 / area2, ratio)
        assert_allclose(area1 / area3, ratio)

    assert result1.meta.resample.pixel_scale_ratio == 1.0
    assert result2.meta.resample.pixel_scale_ratio == ratio
    assert_allclose(result3.meta.resample.pixel_scale_ratio, ratio)

    result1.close()
    result2.close()
    result3.close()


def test_weight_type(nircam_rate, tmp_cwd):
    """Check that weight_type of exptime and ivm work"""
    im1 = AssignWcsStep.call(nircam_rate, sip_approx=False)
    _set_photom_kwd(im1)
    im1.var_rnoise[:] = 0
    im2 = im1.copy()
    im3 = im1.copy()
    im1.data += 10
    im2.data += 5
    im3.data += 5
    im1.var_rnoise += 1 / 10
    im2.var_rnoise += 1 / 5
    im3.var_rnoise += 1 / 5
    im2.meta.observation.sequence_id = "2"
    im3.meta.observation.sequence_id = "3"

    c = ModelLibrary([im1, im2, im3])
    assert len(c.group_names) == 3

    result1 = ResampleStep.call(c, weight_type="ivm", blendheaders=False, save_results=True)

    # assert_allclose(result1.data, result2.data)
    # assert_allclose(result1.wht, result2.wht)
    assert_allclose(result1.data[100:105, 100:105], 7.5, rtol=1e-2)
    assert_allclose(result1.wht[100:105, 100:105], 19.5, rtol=1e-2)

    result2 = ResampleStep.call(c, weight_type="exptime", blendheaders=False)

    assert_allclose(result2.data[100:105, 100:105], 6.667, rtol=1e-2)
    expectation_value = 407.0
    assert_allclose(result2.wht[100:105, 100:105], expectation_value, rtol=1e-2)

    # remove measurement time to force use of exposure time
    # this also implicitly shows that measurement time was indeed used above
    expected_ratio = im1.meta.exposure.exposure_time / im1.meta.exposure.measurement_time
    with c:
        for j, im in enumerate(c):
            del im.meta.exposure.measurement_time
            c.shelve(im, j)

    result3 = ResampleStep.call(c, weight_type="exptime", blendheaders=False)
    assert_allclose(result3.data[100:105, 100:105], 6.667, rtol=1e-2)
    assert_allclose(result3.wht[100:105, 100:105], expectation_value * expected_ratio, rtol=1e-2)

    im1.close()
    im2.close()
    im3.close()
    result1.close()
    result2.close()
    result3.close()


def test_sip_coeffs_do_not_propagate(nircam_rate):
    im = AssignWcsStep.call(nircam_rate, sip_degree=2)
    _set_photom_kwd(im)

    # Check some SIP keywords produced above
    assert im.meta.wcsinfo.cd1_1 is not None
    assert im.meta.wcsinfo.ctype1 == "RA---TAN-SIP"

    # Make sure no PC matrix stuff is there
    assert im.meta.wcsinfo.pc1_1 is None

    result = ResampleStep.call(im)

    # Verify that SIP-related keywords do not propagate to resampled output
    assert result.meta.wcsinfo.cd1_1 is None
    assert result.meta.wcsinfo.ctype1 == "RA---TAN"

    # Make sure we have a PC matrix
    assert result.meta.wcsinfo.pc1_1 is not None

    im.close()
    result.close()


def test_build_interpolated_output_wcs(miri_rate_pair):
    im1, im2 = miri_rate_pair

    driz = ResampleSpec(ModelContainer([im1, im2]))
    output_wcs = driz.build_interpolated_output_wcs([im1, im2])

    # Make sure that all RA, Dec values in the input image have a location in
    # the output frame
    grid = grid_from_bounding_box(im2.meta.wcs.bounding_box)
    ra, dec, lam = im2.meta.wcs(*grid)
    x, y = output_wcs.invert(ra, dec, lam)
    nn = ~(np.isnan(x) | np.isnan(y))
    assert np.sum(nn) > 0
    assert np.all(x[nn] > -1)

    # Make sure the output slit size is larger than the input slit size
    # for this nodded data
    assert output_wcs.array_shape[1] > ra.shape[1]


def test_build_nirspec_output_wcs(nirspec_cal_pair):
    im1, im2 = nirspec_cal_pair
    containers = multislit_to_container([im1, im2])
    driz = ResampleSpec(containers["1"])
    output_wcs = driz.build_nirspec_output_wcs(containers["1"])

    # Make sure that all slit values in the input images have a
    # location in the output frame, in both RA/Dec and slit units
    output_s2d = output_wcs.get_transform("slit_frame", "detector")
    for im in [im1, im2]:
        grid = grid_from_bounding_box(im.slits[0].meta.wcs.bounding_box)

        # check slit values
        input_d2s = im.slits[0].meta.wcs.get_transform("detector", "slit_frame")
        sx, sy, lam = input_d2s(*grid)
        x, y = output_s2d(np.full_like(sy, 0), sy, lam * 1e6)
        nn = ~(np.isnan(x) | np.isnan(y))
        assert np.sum(nn) > 0
        assert np.all(y[nn] > -1)

        # check RA, Dec, lam
        ra, dec, lam = im.slits[0].meta.wcs(*grid)
        x, y = output_wcs.invert(ra, dec, lam)
        nn = ~(np.isnan(x) | np.isnan(y))
        assert np.sum(nn) > 0
        assert np.all(y[nn] > -1)

    # Make a WCS for each input individually
    containers = multislit_to_container([im1])
    driz = ResampleSpec(containers["1"])
    compare_wcs_1 = driz.build_nirspec_output_wcs(containers["1"])

    containers = multislit_to_container([im2])
    driz = ResampleSpec(containers["1"])
    compare_wcs_2 = driz.build_nirspec_output_wcs(containers["1"])

    # The output shape should be the larger of the two
    assert output_wcs.array_shape[0] == max(
        compare_wcs_1.array_shape[0], compare_wcs_2.array_shape[0]
    )
    assert output_wcs.array_shape[1] == max(
        compare_wcs_1.array_shape[1], compare_wcs_2.array_shape[1]
    )


def test_wcs_keywords(nircam_rate):
    """Make sure certain wcs keywords are removed after resample"""
    im = AssignWcsStep.call(nircam_rate)
    result = ResampleStep.call(im)

    assert result.meta.wcsinfo.v2_ref is None
    assert result.meta.wcsinfo.v3_ref is None
    assert result.meta.wcsinfo.ra_ref is None
    assert result.meta.wcsinfo.dec_ref is None
    assert result.meta.wcsinfo.roll_ref is None
    assert result.meta.wcsinfo.v3yangle is None
    assert result.meta.wcsinfo.vparity is None

    im.close()
    result.close()


@pytest.mark.parametrize(
    "n_images,weight_type",
    [
        (1, "ivm"),
        (2, "ivm"),
        (3, "ivm"),
        (9, "ivm"),
        (1, "exptime"),
        (2, "exptime"),
        (3, "exptime"),
        (9, "exptime"),
    ],
)
def test_resample_variance(nircam_rate, n_images, weight_type):
    """Test that resampled variance and error arrays are computed properly"""
    err = 0.02429
    var_rnoise = 0.00034
    var_poisson = 0.00025
    im = AssignWcsStep.call(nircam_rate)
    _set_photom_kwd(im)
    im.var_rnoise += var_rnoise
    im.var_poisson += var_poisson
    im.err += err
    im.meta.filename = "foo.fits"

    c = ModelLibrary([im.copy() for _ in range(n_images)])

    result = ResampleStep.call(c, blendheaders=False, weight_type=weight_type)

    # Verify that the combined uncertainty goes as 1 / sqrt(N)
    assert_allclose(result.err[5:-5, 5:-5].mean(), err / np.sqrt(n_images), atol=1e-5)
    assert_allclose(result.var_rnoise[5:-5, 5:-5].mean(), var_rnoise / n_images, atol=1e-7)
    assert_allclose(result.var_poisson[5:-5, 5:-5].mean(), var_poisson / n_images, atol=1e-7)

    im.close()
    result.close()


@pytest.mark.parametrize("enable_ctx", [True, False])
@pytest.mark.parametrize("report_var", [True, False])
@pytest.mark.parametrize("enable_err", [True, False])
def test_resample_variance_context_disable(
    nircam_rate, enable_ctx, report_var, enable_err, tmp_cwd
):
    """Test that con, var, and err arrays respect their control flags."""
    n_images = 3
    err = 0.02429
    var_rnoise = 0.00034
    var_poisson = 0.00025
    im = AssignWcsStep.call(nircam_rate)
    _set_photom_kwd(im)
    im.var_rnoise += var_rnoise
    im.var_poisson += var_poisson
    im.err += err
    im.meta.filename = "foo.fits"
    # Adding bunit_err covers a bug where the ERR extension was being created just to hold that
    im.meta.bunit_err = "MJy/sr"

    c = ModelLibrary([im.copy() for _ in range(n_images)])

    result = ResampleStep.call(
        c,
        enable_ctx=enable_ctx,
        report_var=report_var,
        enable_err=enable_err,
        blendheaders=False,
    )

    if report_var and enable_err:
        # Verify that the combined uncertainty goes as 1 / sqrt(N)
        assert_allclose(result.err[5:-5, 5:-5].mean(), err / np.sqrt(n_images), atol=1e-5)
        assert_allclose(result.var_rnoise[5:-5, 5:-5].mean(), var_rnoise / n_images, atol=1e-7)
        assert_allclose(result.var_poisson[5:-5, 5:-5].mean(), var_poisson / n_images, atol=1e-7)
    elif enable_err and not report_var:
        assert result.hasattr("err")
        assert_allclose(result.err[5:-5, 5:-5].mean(), err / np.sqrt(n_images), atol=1e-5)
        assert not result.hasattr("var_flat")
        assert not result.hasattr("var_rnoise")
        assert not result.hasattr("var_poisson")
    else:
        assert not result.hasattr("err")
        assert not result.hasattr("var_flat")
        assert not result.hasattr("var_rnoise")
        assert not result.hasattr("var_poisson")
    if enable_ctx:
        assert result.con.shape == (1, 210, 210)
        bitmap = [1, 1, 1]  # 3 images, all contributing.
        expected_con = int("".join(map(str, bitmap)), 2)
        isin = result.con == expected_con
        # The majority of the image should have all 3 images contributing
        assert np.sum(isin) / isin.size > 0.95
    else:
        assert not result.hasattr("con")

    im.close()

    result.save("tmp.fits")
    result.close()

    # Cover a bug where the ERR extension was being created even when enable_err=False
    # just to hold the bunit_err keyword.
    with fits.open("tmp.fits") as hdul:
        if not enable_err:
            assert "ERR" not in hdul


@pytest.mark.parametrize("shape", [(0,), (10, 1)])
def test_resample_undefined_variance(nircam_rate, shape):
    """Test that resampled variance and error arrays are computed properly"""
    im = AssignWcsStep.call(nircam_rate)
    im.var_rnoise = np.ones(shape, dtype=im.var_rnoise.dtype.type)
    im.var_poisson = np.ones(shape, dtype=im.var_poisson.dtype.type)
    im.var_flat = np.ones(shape, dtype=im.var_flat.dtype.type)
    im.meta.filename = "foo.fits"
    c = ModelLibrary([im])

    with pytest.warns(RuntimeWarning, match="'var_rnoise' array not available"):
        result = ResampleStep.call(c, blendheaders=False)

    # no valid variance - output error and variance are all NaN
    assert_allclose(result.err, np.nan)
    assert_allclose(result.var_rnoise, np.nan)
    assert_allclose(result.var_poisson, np.nan)
    assert_allclose(result.var_flat, np.nan)

    im.close()
    result.close()


@pytest.mark.parametrize("ratio", [0.7, 1.2])
@pytest.mark.parametrize("rotation", [0, 15, 135])
@pytest.mark.parametrize("crpix", [(256, 488), (700, 124)])
@pytest.mark.parametrize("crval", [(22.04019, 11.98262), (22.0404, 11.983)])
@pytest.mark.parametrize("shape", [(1020, 1010)])
def test_custom_wcs_resample_imaging(nircam_rate, ratio, rotation, crpix, crval, shape):
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    im.data += 5

    result = ResampleStep.call(
        im, output_shape=shape, crpix=crpix, crval=crval, rotation=rotation, pixel_scale_ratio=ratio
    )

    t = result.meta.wcs.forward_transform

    # make sure results are nontrivial
    assert not np.all(np.isnan(result.data))

    # test rotation
    pc = t["pc_rotation_matrix"].matrix.value
    orientation = np.rad2deg(np.arctan2(pc[0, 1], pc[1, 1]))
    assert_allclose(rotation, orientation)

    # test CRPIX
    assert_allclose((-t["crpix1"].offset.value, -t["crpix2"].offset.value), crpix)

    # test CRVAL
    assert_allclose(t(*crpix), crval)

    # test output image shape
    assert result.data.shape == shape[::-1]

    im.close()
    result.close()


@pytest.mark.parametrize(
    "output_shape2, match", [((1205, 1100), True), ((1222, 1111), False), (None, True)]
)
def test_custom_refwcs_resample_imaging(nircam_rate, output_shape2, match, tmp_path):
    # make some data with a WCS and some random values
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    rng = np.random.default_rng(seed=77)
    im.data[:, :] = rng.random(im.data.shape)

    if output_shape2 is None:
        crpix = None
    else:
        crpix = (output_shape2[-1] // 2, output_shape2[-2] // 2)
    crval = tuple(np.mean(im.meta.wcs.footprint(), axis=0))
    rotation = 15
    ratio = 0.7

    result = ResampleStep.call(
        im,
        output_shape=output_shape2,
        crpix=crpix,
        crval=crval,
        rotation=rotation,
        pixel_scale_ratio=ratio,
    )

    # make sure results are nontrivial
    data1 = result.data
    weight1 = result.wht

    assert not np.all(np.isnan(data1))

    if crpix is not None:
        assert_allclose(result.meta.wcs(*crpix), crval, rtol=1e-12, atol=0)

    refwcs = str(tmp_path / "resample_refwcs.asdf")
    asdf.AsdfFile({"wcs": result.meta.wcs, "array_shape": data1.shape}).write_to(refwcs)

    result = ResampleStep.call(im, output_shape=output_shape2, output_wcs=refwcs)

    data2 = result.data
    weight2 = result.wht
    assert not np.all(np.isnan(data2))

    if output_shape2 is not None:
        assert data2.shape == output_shape2[::-1]

    if match:
        # test output image shape
        assert data1.shape == data2.shape
        assert_allclose(data1, data2, equal_nan=True, rtol=1.0e-7, atol=1e-7)

    # make sure pixel values are similar, accounting for scale factor
    # (assuming inputs are in surface brightness units)
    input_mean = np.nanmean(im.data)
    total_weight = np.sum(weight1)
    output_mean_1 = np.nansum(data1 * weight1) / total_weight
    output_mean_2 = np.nansum(data2 * weight2) / total_weight
    # rtol and atol values are from np.isclose default settings.
    assert_allclose(input_mean, output_mean_1, rtol=1e-5, atol=1e-8)
    assert_allclose(input_mean, output_mean_2, rtol=1e-5, atol=1e-8)

    im.close()
    result.close()


@pytest.mark.parametrize("weight_type", ["ivm", "exptime"])
def test_custom_refwcs_pixel_shape_imaging(nircam_rate, tmp_path, weight_type):
    # make some data with a WCS and some random values
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    rng = np.random.default_rng(seed=77)
    im.data[:, :] = 0.0
    im.data[5:-5, 5:-5] = rng.random(tuple(i - 10 for i in im.data.shape))
    im.dq[:, :] = 1
    im.dq[5:-5, 5:-5] = 0

    data0 = deepcopy(im.data)

    crpix = (600, 550)
    crval = (22.04019, 11.98262)
    rotation = 15
    ratio = 0.7

    im.meta.group_id = "1"
    im_dict = input_jwst_model_to_dict(
        im, weight_type=weight_type, enable_var=True, compute_err=True
    )

    in_weight = build_driz_weight(
        im_dict, weight_type=weight_type, good_bits=GOOD_BITS, flag_name_map=dqflags.pixel
    )

    # first pass - create a reference output WCS:
    assert np.all(np.isfinite(im.data))
    result = ResampleStep.call(
        im,
        output_shape=(1205, 1100),
        crpix=crpix,
        crval=crval,
        rotation=rotation,
        pixel_scale_ratio=ratio,
        weight_type=weight_type,
    )
    # TODO: next assert would fail. Why does resample step modify input data???
    # assert np.all(np.isfinite(im.data))

    data1 = result.data
    wht1 = result.wht

    # remove the bounding box so shape is set from pixel_shape
    # and also set a top-level pixel area
    pixel_area = 1e-13
    refwcs = str(tmp_path / "resample_refwcs.asdf")
    result.meta.wcs.bounding_box = None
    asdf.AsdfFile({"wcs": result.meta.wcs, "pixel_area": pixel_area}).write_to(refwcs)
    result = ResampleStep.call(im, output_wcs=refwcs, weight_type=weight_type)

    data2 = result.data
    wht2 = result.wht
    assert not np.all(np.isnan(data2))

    # test output image shape
    assert data1.shape == data2.shape
    assert_allclose(data1, data2, equal_nan=True)

    # make sure pixel values are similar, accounting for scale factor
    # (assuming inputs are in surface brightness units)
    assert np.isclose(
        np.sum(data0 * in_weight) / np.sum(in_weight),
        np.nansum(data1 * wht1) / np.sum(wht1),
        atol=1e-4,
    )
    assert np.isclose(np.sum(data0 * in_weight), np.nansum(data2 * wht2), atol=1e-4)

    # check that output pixel area is set from input
    # rtol and atol values are from np.isclose default settings.
    assert_allclose(result.meta.photometry.pixelarea_steradians, pixel_area, rtol=1e-5, atol=1e-8)

    im.close()
    result.close()


@pytest.mark.parametrize("ratio", [0.7, 1.0, 1.3])
def test_custom_refwcs_resample_miri(miri_cal, tmp_path, ratio):
    im = miri_cal
    miri_cal.meta.bunit_data = "MJy"

    # mock a spectrum by giving the first slit some random
    # values at the center
    rng = np.random.default_rng(seed=77)
    new_values = rng.random(im.data.shape)

    center = im.data.shape[1] // 2
    im.data[:] = 0.0
    im.data[:, center - 2 : center + 2] = new_values[:, center - 2 : center + 2]

    # first pass: create a reference output WCS with a custom pixel scale
    result = ResampleSpecStep.call(im, pixel_scale_ratio=ratio)

    # make sure results are nontrivial
    data1 = result.data
    assert not np.all(np.isnan(data1))

    # save the wcs from the output
    refwcs = str(tmp_path / "resample_refwcs.asdf")
    asdf.AsdfFile({"wcs": result.meta.wcs}).write_to(refwcs)

    # run again, this time using the created WCS as input
    result = ResampleSpecStep.call(im, output_wcs=refwcs)
    data2 = result.data
    assert not np.all(np.isnan(data2))

    # check output data against first pass
    assert data1.shape == data2.shape
    assert_allclose(data1, data2, equal_nan=True, rtol=1e-4)

    # make sure flux is conserved: sum over spatial dimension
    # should be same in input and output
    # (assuming inputs are in flux density units)
    input_sum = np.nanmean(np.nansum(im.data, axis=1))
    output_sum_1 = np.nanmean(np.nansum(data1, axis=1))
    output_sum_2 = np.nanmean(np.nansum(data2, axis=1))
    assert_allclose(input_sum, output_sum_1, rtol=0.005)
    assert_allclose(input_sum, output_sum_2, rtol=0.005)

    im.close()
    result.close()


@pytest.mark.parametrize("ratio", [0.7, 1.0, 1.3])
def test_custom_refwcs_resample_nirspec(nirspec_cal, tmp_path, ratio):
    im = nirspec_cal
    for slit in im.slits:
        slit.meta.bunit_data = "MJy"

    # mock a spectrum by giving the first slit some random
    # values at the center
    rng = np.random.default_rng(seed=77)
    new_values = rng.random(im.slits[0].data.shape)

    center = im.slits[0].data.shape[0] // 2
    im.slits[0].data[:] = 0.0
    im.slits[0].data[center - 2 : center + 2, :] = new_values[center - 2 : center + 2, :]

    # first pass: create a reference output WCS with a custom pixel scale
    result = ResampleSpecStep.call(im, pixel_scale_ratio=ratio)

    # make sure results are nontrivial
    data1 = result.slits[0].data
    assert not np.all(np.isnan(data1))

    # save the wcs from the output
    refwcs = str(tmp_path / "resample_refwcs.asdf")
    asdf.AsdfFile({"wcs": result.slits[0].meta.wcs}).write_to(refwcs)

    # run again, this time using the created WCS as input
    result = ResampleSpecStep.call(im, output_wcs=refwcs)

    data2 = result.slits[0].data
    assert not np.all(np.isnan(data2))

    # check output data against first pass
    assert data1.shape == data2.shape
    assert_allclose(data1, data2, equal_nan=True, rtol=1e-4)

    # make sure flux is conserved: sum over spatial dimension
    # should be same in input and output
    # (assuming inputs are in flux density units)
    input_sum = np.nanmean(np.nansum(im.slits[0].data, axis=0))
    output_sum_1 = np.nanmean(np.nansum(data1, axis=0))
    output_sum_2 = np.nanmean(np.nansum(data2, axis=0))
    assert_allclose(input_sum, output_sum_1, rtol=0.005)
    assert_allclose(input_sum, output_sum_2, rtol=0.005)

    im.close()
    result.close()


def test_custom_refwcs_pixel_shape_nirspec(nirspec_cal, tmp_path):
    im = nirspec_cal
    for slit in im.slits:
        slit.meta.bunit_data = "MJy/sr"

    # mock a spectrum by giving the first slit some random
    # values at the center
    rng = np.random.default_rng(seed=77)
    new_values = rng.random(im.slits[0].data.shape)

    center = im.slits[0].data.shape[0] // 2
    im.slits[0].data[:] = 0.0
    im.slits[0].data[center - 2 : center + 2, :] = new_values[center - 2 : center + 2, :]

    # first pass: create a reference output WCS with a custom pixel scale
    ratio = 0.7
    result = ResampleSpecStep.call(im, pixel_scale_ratio=ratio)

    # make sure results are nontrivial
    data1 = result.slits[0].data
    assert not np.all(np.isnan(data1))

    # remove the bounding box from the WCS so shape is set from pixel_shape
    # and also set a top-level pixel area
    pixel_area = 1e-13
    refwcs = str(tmp_path / "resample_refwcs.asdf")
    asdf.AsdfFile({"wcs": result.slits[0].meta.wcs, "pixel_area": pixel_area}).write_to(refwcs)

    # run again, this time using the created WCS as input
    result = ResampleSpecStep.call(im, output_wcs=refwcs)

    data2 = result.slits[0].data
    assert not np.all(np.isnan(data2))

    # check output data against first pass
    assert data1.shape == data2.shape
    assert_allclose(data1, data2, equal_nan=True, rtol=1e-4)

    # check that output pixel area is set from output_wcs
    assert_allclose(result.slits[0].meta.photometry.pixelarea_steradians, pixel_area)

    im.close()
    result.close()


@pytest.mark.parametrize("ratio", [1.3, 1])
def test_custom_wcs_pscale_resample_imaging(nircam_rate, ratio):
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    im.data += 5

    crval = (22.04019, 11.98262)
    input_scale = compute_scale(wcs=im.meta.wcs, fiducial=crval)
    result = ResampleStep.call(im, pixel_scale_ratio=ratio, pixel_scale=3600 * input_scale * 0.75)
    output_scale = compute_scale(wcs=result.meta.wcs, fiducial=crval)

    # test scales are close
    assert_allclose(output_scale, input_scale * 0.75)

    im.close()
    result.close()


@pytest.mark.parametrize("ratio", [1.3, 1])
def test_custom_wcs_pscale_resample_miri(miri_cal, ratio):
    im = miri_cal

    # pass both ratio and direct scale: ratio is ignored in favor of scale
    input_scale = compute_spectral_pixel_scale(im.meta.wcs, disp_axis=2)
    result = ResampleSpecStep.call(
        im, pixel_scale_ratio=ratio, pixel_scale=3600 * input_scale * 0.75
    )
    output_scale = compute_spectral_pixel_scale(result.meta.wcs, disp_axis=2)

    # test scales are close to scale specified, regardless of ratio
    # rtol and atol values are from np.allclose default settings.
    assert_allclose(output_scale, input_scale * 0.75, rtol=1e-05, atol=1e-08)

    result.close()


@pytest.mark.parametrize("ratio", [1.3, 1])
def test_custom_wcs_pscale_resample_nirspec(nirspec_cal, ratio):
    im = nirspec_cal.slits[0]

    # pass both ratio and direct scale: ratio is ignored in favor of scale
    input_scale = compute_spectral_pixel_scale(im.meta.wcs, disp_axis=1)
    result = ResampleSpecStep.call(
        nirspec_cal, pixel_scale_ratio=ratio, pixel_scale=3600 * input_scale * 0.75
    )
    output_scale = compute_spectral_pixel_scale(result.slits[0].meta.wcs, disp_axis=1)

    # test scales are close to scale specified, regardless of ratio
    # rtol and atol values are from np.allclose default settings.
    assert_allclose(output_scale, input_scale * 0.75, rtol=1e-05, atol=1e-08)

    result.close()


@pytest.mark.parametrize("wcs_attr", ["pixel_shape", "array_shape", "bounding_box"])
def test_custom_wcs_input(tmp_path, nircam_rate, wcs_attr):
    # make a valid WCS
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    wcs = im.meta.wcs

    # store values in a dictionary
    wcs_dict = {
        "array_shape": im.data.shape,
        "pixel_shape": im.data.shape[::-1],
        "bounding_box": wcs.bounding_box,
    }

    # Set all attributes to None
    for attr in ["pixel_shape", "array_shape", "bounding_box"]:
        setattr(wcs, attr, None)

    # Set the attribute to the correct value
    setattr(wcs, wcs_attr, wcs_dict[wcs_attr])

    # write the WCS to an asdf file
    refwcs = str(tmp_path / "test_wcs.asdf")
    asdf.AsdfFile({"wcs": wcs, "pixel_area": 7919}).write_to(refwcs)

    # load the WCS from the asdf file
    loaded_wcs = load_custom_wcs(refwcs)
    assert loaded_wcs["pixel_area"] == 7919

    # check that the loaded WCS has the correct values
    for attr in ["pixel_shape", "array_shape"]:
        assert_allclose(getattr(loaded_wcs["wcs"], attr), wcs_dict[attr])


@pytest.mark.parametrize(
    "override,value", [("pixel_shape", (300, 400)), ("array_shape", (400, 300))]
)
def test_custom_wcs_input_overrides(tmp_path, nircam_rate, override, value):
    # make a valid WCS
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    wcs = im.meta.wcs

    # remove existing shape keys if testing shape overrides
    if override != "pixel_area":
        wcs.pixel_shape = None
        wcs.bounding_box = None

    expected_array_shape = im.data.shape
    expected_pixel_shape = im.data.shape[::-1]

    # write the WCS to an asdf file with a top-level override
    refwcs = str(tmp_path / "test_wcs.asdf")
    asdf.AsdfFile({"wcs": wcs, override: value}).write_to(refwcs)

    # check for expected values when read back in
    keys = ["pixel_shape", "array_shape"]
    loaded_wcs = load_custom_wcs(refwcs)["wcs"]

    for key in keys:
        if key == override:
            assert_allclose(getattr(loaded_wcs, key), value)
        elif key == "pixel_shape":
            if override == "array_shape":
                assert_allclose(getattr(loaded_wcs, key), value[::-1])
            else:
                assert_allclose(getattr(loaded_wcs, key), expected_pixel_shape)
        elif key == "array_shape":
            if override == "pixel_shape":
                assert_allclose(getattr(loaded_wcs, key), value[::-1])
            else:
                assert_allclose(getattr(loaded_wcs, key), expected_array_shape)


def test_custom_wcs_input_error(tmp_path, nircam_rate):
    # make a valid WCS
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    wcs = im.meta.wcs

    # remove shape settings
    wcs.pixel_shape = None
    wcs.array_shape = None
    wcs.bounding_box = None

    # write the WCS to an asdf file
    refwcs = str(tmp_path / "test_wcs.asdf")
    asdf.AsdfFile({"wcs": wcs}).write_to(refwcs)

    # loading the file without shape info should produce an error
    with pytest.raises(ValueError, match="'output_shape' is required"):
        loaded_wcs = load_custom_wcs(refwcs)["wcs"]

    # providing an output shape should succeed
    output_shape = (300, 400)
    loaded_wcs = load_custom_wcs(refwcs, output_shape=output_shape)["wcs"]

    # array shape is opposite of input values (numpy convention)
    assert np.all(loaded_wcs.array_shape == output_shape[::-1])

    # pixel shape matches
    assert np.all(loaded_wcs.pixel_shape == output_shape)

    # bounding box is not set
    assert loaded_wcs.bounding_box is None


def test_pixscale(nircam_rate):
    # check that if both 'pixel_scale_ratio' and 'pixel_scale' are passed in,
    # that 'pixel_scale' overrides correctly
    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    _set_photom_kwd(im)
    pixarea = im.meta.photometry.pixelarea_arcsecsq

    # check when both pixel_scale and pixel_scale_ratio are passed in
    res = ResampleStep.call(im, pixel_scale=0.04, pixel_scale_ratio=0.7)
    assert_allclose(res.meta.resample.pixel_scale_ratio, 0.04 / np.sqrt(pixarea))

    # just pixel_scale
    res = ResampleStep.call(im, pixel_scale=0.04)
    assert_allclose(res.meta.resample.pixel_scale_ratio, 0.04 / np.sqrt(pixarea))

    # just pixel_scale_ratio
    res = ResampleStep.call(im, pixel_scale_ratio=0.7)
    assert res.meta.resample.pixel_scale_ratio == 0.7

    im.close()
    res.close()


def test_phot_keywords(nircam_rate):
    # test that resample keywords agree with photometry keywords after step is run

    im = AssignWcsStep.call(nircam_rate, sip_approx=False)
    _set_photom_kwd(im)

    orig_pix_area_sr = im.meta.photometry.pixelarea_steradians
    orig_pix_area_arcsec = im.meta.photometry.pixelarea_arcsecsq

    # first run by setting `pixel_scale`
    res = ResampleStep.call(im, pixel_scale=0.04)
    new_psr = res.meta.resample.pixel_scale_ratio

    assert_allclose(
        res.meta.resample.pixel_scale_ratio,
        0.04 / np.sqrt(orig_pix_area_arcsec),
        atol=0,
        rtol=1e-12,
    )
    assert_allclose(
        res.meta.photometry.pixelarea_steradians, orig_pix_area_sr * new_psr**2, atol=0, rtol=1e-12
    )
    assert_allclose(
        res.meta.photometry.pixelarea_arcsecsq,
        orig_pix_area_arcsec * new_psr**2,
        atol=0,
        rtol=1e-12,
    )

    im.close()
    res.close()


def test_missing_nominal_area(miri_cal, tmp_path):
    # remove nominal pixel area
    miri_cal.meta.photometry.pixelarea_steradians = None
    miri_cal.meta.photometry.pixelarea_arcsecsq = None

    # result should still process okay
    result = ResampleSpecStep.call(miri_cal)

    # no area keywords in output
    assert result.meta.photometry.pixelarea_steradians is None
    assert result.meta.photometry.pixelarea_arcsecsq is None

    # direct pixel scale setting is not supported
    result2 = ResampleSpecStep.call(miri_cal, pixel_scale=0.5)
    assert_allclose(result2.data, result.data, equal_nan=True)
    assert result2.meta.resample.pixel_scale_ratio == 1.0

    # setting pixel_scale_ratio is still allowed,
    # output area is still None
    result3 = ResampleSpecStep.call(miri_cal, pixel_scale_ratio=0.5)
    assert result3.data.shape[0] == result.data.shape[0]
    assert result3.data.shape[1] < result.data.shape[1]
    assert result3.meta.resample.pixel_scale_ratio == 0.5
    assert result3.meta.photometry.pixelarea_steradians is None
    assert result3.meta.photometry.pixelarea_arcsecsq is None

    # specifying a custom WCS without nominal area sets output pixel
    # scale ratio to 1, since direct scale cannot be computed
    # save the wcs from the output
    refwcs = str(tmp_path / "resample_refwcs.asdf")
    asdf.AsdfFile({"wcs": result3.meta.wcs}).write_to(refwcs)
    result4 = ResampleSpecStep.call(miri_cal, output_wcs=refwcs)
    assert result4.data.shape == result3.data.shape
    assert result4.meta.resample.pixel_scale_ratio == 1.0
    assert result4.meta.photometry.pixelarea_steradians is None
    assert result4.meta.photometry.pixelarea_arcsecsq is None

    result.close()
    result2.close()
    result3.close()
    result4.close()


def test_nirspec_lamp_pixscale(nirspec_lamp, tmp_path):
    # RuntimeWarning: invalid value encountered in sqrt
    with np.errstate(invalid="ignore"):
        result = ResampleSpecStep.call(nirspec_lamp)

    # output data should have the same wavelength size,
    # spatial size is close
    assert_allclose(result.slits[0].data.shape[0], nirspec_lamp.slits[0].data.shape[0], atol=5)
    assert result.slits[0].data.shape[1] == nirspec_lamp.slits[0].data.shape[1]

    # test pixel scale setting: will not work without sky-based WCS
    # RuntimeWarning: invalid value encountered in sqrt
    with np.errstate(invalid="ignore"):
        result2 = ResampleSpecStep.call(nirspec_lamp, pixel_scale=0.5)
    assert_allclose(result2.slits[0].data, result.slits[0].data, equal_nan=True)
    assert result2.meta.resample.pixel_scale_ratio == 1.0

    # setting pixel_scale_ratio is still allowed
    # RuntimeWarning: invalid value encountered in sqrt
    with np.errstate(invalid="ignore"):
        result3 = ResampleSpecStep.call(nirspec_lamp, pixel_scale_ratio=0.5)
    assert result3.slits[0].data.shape[0] < result.slits[0].data.shape[0]
    assert result3.slits[0].data.shape[1] == result.slits[0].data.shape[1]
    assert result3.meta.resample.pixel_scale_ratio == 0.5

    # specifying a custom WCS should work, but output ratio is 1.0,
    # since output scale cannot be determined
    refwcs = str(tmp_path / "resample_refwcs.asdf")
    asdf.AsdfFile({"wcs": result3.slits[0].meta.wcs}).write_to(refwcs)
    # RuntimeWarning: invalid value encountered in sqrt
    with np.errstate(invalid="ignore"):
        result4 = ResampleSpecStep.call(nirspec_lamp, output_wcs=refwcs)
    assert result4.slits[0].data.shape == result3.slits[0].data.shape
    assert result4.meta.resample.pixel_scale_ratio == 1.0

    result.close()
    result2.close()
    result3.close()
    result4.close()
