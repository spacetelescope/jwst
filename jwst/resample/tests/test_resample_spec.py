import numpy as np
from numpy.testing import assert_allclose

from ...datamodels import ImageModel
from jwst.assign_wcs import AssignWcsStep
from jwst.extract_2d import Extract2dStep
from jwst.resample import ResampleSpecStep

from gwcs.wcstools import grid_from_bounding_box


def test_spatial_transform_nirspec():
    wcsinfo = {
        'dec_ref': -0.00601415671349804,
        'ra_ref': -0.02073605215697509,
        'roll_ref': -0.0,
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

    im = ImageModel()
    im.data = np.random.rand(2048, 2048)
    im.error = np.random.rand(2048, 2048)
    im.dq = np.random.rand(2048, 2048)

    im.meta.wcsinfo._instance.update(wcsinfo)
    im.meta.instrument._instance.update(instrument)
    im.meta.observation._instance.update(observation)
    im.meta.exposure._instance.update(exposure)
    im.meta.subarray._instance.update(subarray)
    im = AssignWcsStep.call(im)
    im = Extract2dStep.call(im)
    im = ResampleSpecStep.call(im)

    for slit in im.slits:
        x, y =grid_from_bounding_box(slit.meta.wcs.bounding_box)
        ra, dec, lam = slit.meta.wcs(x, y)

        ra1 = np.where(ra < 0, 360 + ra, ra)
        assert_allclose(slit.meta.wcs.invert(ra, dec, lam), slit.meta.wcs.invert(ra1, dec, lam))


def test_spatial_transform_miri():
    wcsinfo = {
        'dec_ref': -0.00601415671349804,
        'ra_ref': -0.02073605215697509,
        'roll_ref': -0.0,
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

    im = ImageModel()
    im.data = np.random.rand(416, 72)
    im.error = np.random.rand(416, 72)
    im.dq = np.random.rand(416, 72)

    im.meta.wcsinfo._instance.update(wcsinfo)
    im.meta.instrument._instance.update(instrument)
    im.meta.observation._instance.update(observation)
    im.meta.exposure._instance.update(exposure)
    im.meta.subarray._instance.update(subarray)

    out = AssignWcsStep.call(im)
    out = ResampleSpecStep.call(out)
    x, y =grid_from_bounding_box(out.meta.wcs.bounding_box)
    ra, dec, lam = out.meta.wcs(x, y)
    ra1 = np.where(ra < 0, 360 + ra, ra)
    assert_allclose(out.meta.wcs.invert(ra, dec, lam), out.meta.wcs.invert(ra1, dec, lam))
