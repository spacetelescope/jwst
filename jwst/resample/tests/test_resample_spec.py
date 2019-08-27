import numpy as np
from numpy.testing import assert_allclose

from gwcs.wcstools import grid_from_bounding_box

from ...datamodels import ImageModel
from jwst.assign_wcs import AssignWcsStep
from jwst.resample import ResampleSpecStep


wcsinfo = {
    'dec_ref': -0.00601415671349804,
    'ra_ref': -0.02073605215697509,
    'roll_ref': -0.0,
    'v2_ref': -453.5134,
    'v3_ref': -373.4826,
    'v3yangle': 0.0,
    'vparity': -1
}


instrument = {
    'detector': 'MIRIMAGE',
    'filter': 'P750L',
    'name': 'MIRI'
}


observation = {
    'date': '2019-01-01',
    'time': '17:00:00'}


subarray = {
    'fastaxis': 1,
    'name': 'SUBPRISM',
    'slowaxis': 2,
    'xsize': 72,
    'xstart': 1,
    'ysize': 416,
    'ystart': 529
}


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


def test_spatial_transform():
    """
    Calling the backwards WCS transform gives the same results
    for ``negative RA`` and ``negative RA + 360``.
    """
    im = ImageModel()
    im.meta.wcsinfo._instance.update(wcsinfo)
    im.meta.instrument._instance.update(instrument)
    im.meta.exposure._instance.update(exposure)
    im.meta.observation._instance.update(observation)
    im.meta.subarray._instance.update(subarray)

    im = AssignWcsStep.call(im)
    im.data = np.random.rand(416, 72)
    im.error = np.random.rand(416, 72)
    im.dq = np.random.rand(416, 72)

    im = ResampleSpecStep.call(im)
    x, y =grid_from_bounding_box(im.meta.wcs.bounding_box)
    ra, dec, lam = im.meta.wcs(x, y)
    ra1 = np.where(ra < 0, 360 + ra, ra)
    assert_allclose(im.meta.wcs.invert(ra, dec, lam), im.meta.wcs.invert(ra1, dec, lam))
