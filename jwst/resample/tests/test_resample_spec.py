from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose
import numpy as np
import pytest

from jwst.datamodels import ImageModel
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst.resample import ResampleSpecStep, ResampleStep

@pytest.fixture
def nirspec_rate():
    wcsinfo = {
        'dec_ref': 40,
        'ra_ref': 100,
        'roll_ref': 0,
        'v2_ref': -453.5134,
        'v3_ref': -373.4826,
        'v3yangle': 0.0,
        'vparity': -1}

    instrument = {
        'detector': 'NRS1',
        'filter': 'CLEAR',
        'grating': 'PRISM',
        'name': 'NIRSPEC',
        'gwa_tilt': 37.0610,
        'gwa_xtilt': 0.0001,
        'gwa_ytilt': 0.0001}

    subarray = {
        'fastaxis': 1,
        'name': 'SUBS200A1',
        'slowaxis': 2,
        'xsize': 72,
        'xstart': 1,
        'ysize': 416,
        'ystart': 529}

    observation = {
        'date': '2016-09-05',
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
        'ngroups': 100,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'readpatt': 'NRSRAPID',
        'sample_time': 10.0,
        'start_time': 58119.8333,
        'type': 'NRS_FIXEDSLIT',
        'zero_frame': False}

    im = ImageModel((2048, 2048))

    im.meta.wcsinfo = wcsinfo
    im.meta.instrument = instrument
    im.meta.observation = observation
    im.meta.exposure = exposure
    im.meta.subarray = subarray

    return im


@pytest.fixture
def miri_rate():
    wcsinfo = {
        'dec_ref': 40,
        'ra_ref': 100,
        'roll_ref': 0.0,
        'v2_ref': -453.5134,
        'v3_ref': -373.4826,
        'v3yangle': 0.0,
        'vparity': -1}

    instrument = {
        'detector': 'MIRIMAGE',
        'filter': 'P750L',
        'name': 'MIRI'}

    observation = {
        'date': '2019-01-01',
        'time': '17:00:00'}

    subarray = {
        'fastaxis': 1,
        'name': 'SLITLESSPRISM',
        'slowaxis': 2,
        'xsize': 72,
        'xstart': 1,
        'ysize': 416,
        'ystart': 529}

    exposure = {
        'duration': 11.805952,
        'end_time': 58119.85416,
        'exposure_time': 11.776,
        'frame_time': 0.11776,
        'group_time': 0.11776,
        'groupgap': 0,
        'integration_time': 11.776,
        'nframes': 1,
        'ngroups': 100,
        'nints': 1,
        'nresets_between_ints': 0,
        'nsamples': 1,
        'readpatt': 'FAST',
        'sample_time': 10.0,
        'start_time': 58119.8333,
        'type': 'MIR_LRS-SLITLESS',
        'zero_frame': False}

    im = ImageModel((416, 72))
    im.data += 5.

    im.meta.wcsinfo = wcsinfo
    im.meta.instrument = instrument
    im.meta.observation = observation
    im.meta.exposure = exposure
    im.meta.subarray = subarray

    return im


@pytest.fixture
def nircam_rate():
    wcsinfo = {
        'ctype1': 'RA---TAN',
        'ctype2': 'DEC--TAN',
        'dec_ref': 11.99875540218638,
        'ra_ref': 22.02351763251896,
        'roll_ref': 0.005076934167039675,
        'v2_ref': 86.039011,
        'v3_ref': -493.385704,
        'v3yangle': -0.07385127,
        'vparity': -1,
        'wcsaxes': 2}
    instrument = {
        'channel': 'LONG',
        'detector': 'NRCALONG',
        'filter': 'F444W',
        'lamp_mode': 'NONE',
        'module': 'A',
        'name': 'NIRCAM',
        'pupil': 'CLEAR'}
    subarray = {
        'fastaxis': -1,
        'name': 'FULL',
        'slowaxis': 2,
        'xsize': 2048,
        'xstart': 1,
        'ysize': 2048,
        'ystart': 1}
    observation = {
        'activity_id': '01',
        'date': '2021-10-25',
        'exposure_number': '00001',
        'obs_id': 'V42424001001P0000000001101',
        'observation_label': 'nircam_ptsrc_only',
        'observation_number': '001',
        'program_number': '42424',
        'sequence_id': '1',
        'time': '16:58:27.258',
        'visit_group': '01',
        'visit_id': '42424001001',
        'visit_number': '001'}
    exposure = {
        'duration': 161.05155,
        'end_time': 59512.70899968495,
        'exposure_time': 150.31478,
        'frame_time': 10.73677,
        'group_time': 21.47354,
        'groupgap': 1,
        'integration_time': 150.31478,
        'mid_time': 59512.70812980775,
        'nframes': 1,
        'ngroups': 7,
        'nints': 1,
        'nresets_at_start': 1,
        'nresets_between_ints': 1,
        'readpatt': 'BRIGHT1',
        'sample_time': 10,
        'start_time': 59512.70725993055,
        'type': 'NRC_IMAGE'}

    im = ImageModel((2048, 2048))
    im.data += 5.

    im.meta.wcsinfo = wcsinfo
    im.meta.instrument = instrument
    im.meta.observation = observation
    im.meta.exposure = exposure
    im.meta.subarray = subarray

    return im


def test_nirspec_wcs_roundtrip(nirspec_rate):
    im = AssignWcsStep.call(nirspec_rate)
    im = Extract2dStep.call(im)
    im = ResampleSpecStep.call(im)

    for slit in im.slits:
        x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
        ra, dec, lam = slit.meta.wcs(x, y)
        xp, yp = slit.meta.wcs.invert(ra, dec, lam)

        assert_allclose(x, xp, atol=1e-8)
        assert_allclose(y, yp, atol=1e-8)


def test_miri_wcs_roundtrip(miri_rate):
    im = AssignWcsStep.call(miri_rate)
    im = ResampleSpecStep.call(im)

    x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
    ra, dec, lam = im.meta.wcs(x, y)
    xp, yp = im.meta.wcs.invert(ra, dec, lam)
    assert_allclose(x, xp, atol=1e-8)
    assert_allclose(y, yp, atol=1e-8)


@pytest.mark.parametrize("ratio", [0.5, 0.7, 1.0])
def test_pixel_scale_ratio_spec(miri_rate, ratio):
    im = AssignWcsStep.call(miri_rate)
    result1 = ResampleSpecStep.call(im)
    result2 = ResampleSpecStep.call(im, pixel_scale_ratio=ratio)
    assert_allclose(np.array(result1.data.shape), np.array(result2.data.shape) * ratio, rtol=1, atol=1)


@pytest.mark.parametrize("ratio", [0.5, 0.7, 1.0])
def test_pixel_scale_ratio_imaging(nircam_rate, ratio):
    im = AssignWcsStep.call(nircam_rate)
    result1 = ResampleStep.call(im)
    result2 = ResampleStep.call(im, pixel_scale_ratio=ratio)
    assert_allclose(np.array(result1.data.shape), np.array(result2.data.shape) * ratio, rtol=1, atol=1)
    e1 = 50
    e2 = int(50 / ratio)
    assert np.mean(result1.data[e1:-e1,e1:-e1]) == np.mean(result2.data[e2:-e2,e2:-e2])
