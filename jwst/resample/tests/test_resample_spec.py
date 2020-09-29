from gwcs.wcstools import grid_from_bounding_box
from numpy.testing import assert_allclose

from jwst.datamodels import ImageModel
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst.resample import ResampleSpecStep


def test_spatial_transform_nirspec():
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

    im = AssignWcsStep.call(im)
    im = Extract2dStep.call(im)
    im = ResampleSpecStep.call(im)

    for slit in im.slits:
        x, y = grid_from_bounding_box(slit.meta.wcs.bounding_box)
        ra, dec, lam = slit.meta.wcs(x, y)
        xp, yp = slit.meta.wcs.invert(ra, dec, lam)

        assert_allclose(x, xp, atol=1e-8)
        assert_allclose(y, yp, atol=1e-8)


def test_spatial_transform_miri():
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

    im.meta.wcsinfo = wcsinfo
    im.meta.instrument = instrument
    im.meta.observation = observation
    im.meta.exposure = exposure
    im.meta.subarray = subarray

    im = AssignWcsStep.call(im)
    im = ResampleSpecStep.call(im)

    x, y = grid_from_bounding_box(im.meta.wcs.bounding_box)
    ra, dec, lam = im.meta.wcs(x, y)
    xp, yp = im.meta.wcs.invert(ra, dec, lam)
    assert_allclose(x, xp, atol=1e-8)
    assert_allclose(y, yp, atol=1e-8)
