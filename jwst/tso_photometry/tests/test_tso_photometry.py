import math

import numpy as np
import pytest

from jwst import datamodels
from jwst.tso_photometry.tso_photometry import tso_aperture_photometry

shape = (7, 100, 150)
xcenter = 75.
ycenter = 50.


def mk_data_array(shape, value, background, xcenter, ycenter, radius):
    """Create a data array"""

    data = np.zeros(shape, dtype=np.float32) + background

    xlow = int(math.floor(xcenter - radius))
    xhigh = int(math.ceil(xcenter + radius)) + 1
    ylow = int(math.floor(ycenter - radius))
    yhigh = int(math.ceil(ycenter + radius)) + 1
    xlow = max(xlow, 0)
    xhigh = min(xhigh, shape[-1])
    ylow = max(ylow, 0)
    yhigh = min(yhigh, shape[-2])

    radius2 = radius**2
    for j in range(ylow, yhigh):
        for i in range(xlow, xhigh):
            dist2 = (float(i) - xcenter)**2 + (float(j) - ycenter)**2
            if dist2 < radius2:
                data[:, j, i] = (value + background)

    return data


def dummy_wcs(x, y):
    """Placeholder WCS"""

    global xcenter, ycenter

    crpix1 = xcenter
    crpix2 = ycenter
    cdelt1 = 0.031
    cdelt2 = 0.031
    crval1 = 15. * (((34.4147 / 60. + 48.) / 60.) + 15.)    # 15h 48m 34.4147s
    crval2 = ((24.295 / 60. + 9.) / 60.) + 28.              # 28d 09' 24.295"

    dec = (y + 1. - crpix2) * cdelt2 + crval2
    cosdec = math.cos(dec * math.pi / 180.)
    ra = (x + 1. - crpix1) * cdelt1 / cosdec + crval1

    return ra, dec


def set_meta(datamodel, sub64p=False):
    """Assign some metadata"""

    datamodel.meta.exposure.nints = datamodel.data.shape[0]
    datamodel.meta.exposure.integration_time = 10.7
    datamodel.meta.exposure.start_time = 58704.62847

    datamodel.meta.exposure.integration_start = 5
    datamodel.meta.exposure.integration_end = (
        datamodel.meta.exposure.integration_start + datamodel.meta.exposure.nints - 1)

    datamodel.meta.instrument.name = 'NIRCAM'
    datamodel.meta.instrument.detector = 'NRCA2'
    datamodel.meta.instrument.channel = 'SHORT'
    datamodel.meta.instrument.filter = 'F210M'
    datamodel.meta.target.catalog_name = 'R Coronae Borealis'
    if sub64p:
        datamodel.meta.instrument.pupil = 'WLP8'
        datamodel.meta.subarray.name = 'SUB64P'
    else:
        datamodel.meta.instrument.pupil = 'CLEAR'
        datamodel.meta.subarray.name = 'FULL'

    datamodel.meta.bunit_data = 'MJy/sr'
    datamodel.meta.bunit_err = 'MJy/sr'
    # NOTE: this is a dummy value that leaves the mock test data values
    # unchanged during the unit conversion in tso_photometry.
    datamodel.meta.photometry.pixelarea_steradians = 1.0e-6

    datamodel.meta.wcs = dummy_wcs


def include_int_times(datamodel):
    """Create an int_times table and copy it into the data model."""

    int_start = datamodel.meta.exposure.integration_start

    nrows = int_start + shape[0] + 2            # create a few extra rows
    # integration_number and time_arr are one_indexed
    integration_number = np.arange(1, nrows + 1, dtype=np.float32)
    time_arr = np.arange(1, nrows + 1, dtype=np.float64) + 58700.
    dummy = np.arange(nrows, dtype=np.float64)

    it_dtype = [('integration_number', '<i4'),
                ('int_start_MJD_UTC', '<f8'),
                ('int_mid_MJD_UTC', '<f8'),
                ('int_end_MJD_UTC', '<f8'),
                ('int_start_BJD_TDB', '<f8'),
                ('int_mid_BJD_TDB', '<f8'),
                ('int_end_BJD_TDB', '<f8')]

    otab = np.array(list(zip(integration_number,
                             dummy, time_arr, dummy,
                             dummy, dummy, dummy)),
                    dtype=it_dtype)
    datamodel.int_times = otab.copy()

    return time_arr


def test_tso_phot_1():

    global shape, xcenter, ycenter

    # pupil = 'CLEAR' and subarray = 'FULL'

    value = 17.
    background = 0.8
    radius = 5.
    radius_inner = 8.
    radius_outer = 11.
    data = mk_data_array(shape, value, background, xcenter, ycenter, radius)
    datamodel = datamodels.CubeModel(data)
    set_meta(datamodel, sub64p=False)

    # Use a larger radius than was used for creating the data.
    catalog = tso_aperture_photometry(datamodel, xcenter, ycenter,
                                      radius + 1., radius_inner,
                                      radius_outer)

    assert catalog.meta['instrument'] == datamodel.meta.instrument.name
    assert catalog.meta['filter'] == datamodel.meta.instrument.filter
    assert catalog.meta['subarray'] == datamodel.meta.subarray.name
    assert catalog.meta['detector'] == datamodel.meta.instrument.detector
    assert catalog.meta['channel'] == datamodel.meta.instrument.channel
    assert catalog.meta['target_name'] == datamodel.meta.target.catalog_name
    assert math.isclose(catalog.meta['xcenter'], xcenter, abs_tol=0.01)
    assert math.isclose(catalog.meta['ycenter'], ycenter, abs_tol=0.01)

    assert np.allclose(catalog['aperture_sum'].value, 1263.4778, rtol=1.e-7)
    assert np.allclose(catalog['aperture_sum_err'].value, 0., atol=1.e-7)
    assert np.allclose(catalog['net_aperture_sum'].value, 1173., rtol=1.e-7)
    assert np.allclose(catalog['annulus_sum'].value, 143.256627, rtol=1.e-7)
    assert np.allclose(catalog['annulus_sum_err'].value, 0., atol=1.e-7)
    assert np.allclose(catalog['annulus_mean'].value, background, rtol=1.e-7)

    assert np.allclose(catalog['annulus_mean'].value, 0.8, rtol=1.e-6)
    assert np.allclose(catalog['annulus_mean_err'].value, 0., rtol=1.e-7)
    assert np.allclose(catalog['net_aperture_sum'].value, 1173., rtol=1.e-7)
    assert np.allclose(catalog['net_aperture_sum_err'].value, 0., atol=1.e-7)


def test_tso_phot_2():

    # pupil = 'WLP8' and subarray = 'SUB64P'

    shape = (7, 64, 64)
    xcenter = 31.
    ycenter = 31.

    value = 17.
    background = 0.8
    radius = 50.
    radius_inner = None
    radius_outer = None
    data = mk_data_array(shape, value, background, xcenter, ycenter, radius)
    datamodel = datamodels.CubeModel(data)
    set_meta(datamodel, sub64p=True)

    catalog = tso_aperture_photometry(datamodel, xcenter, ycenter,
                                      radius, radius_inner,
                                      radius_outer)

    assert catalog.meta['instrument'] == datamodel.meta.instrument.name
    assert catalog.meta['filter'] == datamodel.meta.instrument.filter
    assert catalog.meta['subarray'] == datamodel.meta.subarray.name
    assert catalog.meta['pupil'] == datamodel.meta.instrument.pupil
    assert math.isclose(catalog.meta['xcenter'], xcenter, abs_tol=0.01)
    assert math.isclose(catalog.meta['ycenter'], ycenter, abs_tol=0.01)

    assert np.allclose(catalog['aperture_sum'].value,
                       (value + background) * shape[-1] * shape[-2],
                       rtol=1.e-6)
    assert np.allclose(catalog['aperture_sum_err'].value, 0., atol=1.e-7)


def test_tso_phot_3():

    global shape, xcenter, ycenter

    value = 17.
    background = 0.8
    radius = 5.
    radius_inner = 8.
    radius_outer = 11.
    data = mk_data_array(shape, value, background, xcenter, ycenter, radius)
    datamodel = datamodels.CubeModel(data)
    set_meta(datamodel, sub64p=False)

    int_times = include_int_times(datamodel)

    catalog = tso_aperture_photometry(datamodel, xcenter, ycenter,
                                      radius + 1., radius_inner,
                                      radius_outer)

    offset = datamodel.meta.exposure.integration_start - 1
    slc = slice(offset, offset + shape[0])
    assert np.allclose(catalog['MJD'], int_times[slc], atol=1.e-8)


def test_tso_phot_4():

    global shape, xcenter, ycenter

    value = 17.
    background = 0.8
    radius = 5.
    radius_inner = 8.
    radius_outer = 11.
    data = mk_data_array(shape, value, background, xcenter, ycenter, radius)
    datamodel = datamodels.CubeModel(data)
    set_meta(datamodel, sub64p=False)

    int_times = include_int_times(datamodel)
    # Modify the column of integration numbers so that they extend outside
    # the range of integration numbers in the science data.  This shouldn't
    # happen, i.e. this is to test an edge case.
    int_start = datamodel.meta.exposure.integration_start
    datamodel.int_times['integration_number'] += (2 * int_start)

    catalog = tso_aperture_photometry(datamodel, xcenter, ycenter,
                                      radius + 1., radius_inner,
                                      radius_outer)

    int_times = np.array([58704.62853, 58704.628655, 58704.62878,
                          58704.62890, 58704.6290, 58704.6291511,
                          58704.629275])
    assert np.allclose(catalog['MJD'], int_times, rtol=1.e-8)


def test_tso_phot_5():

    global shape, xcenter, ycenter

    value = 17.
    background = 0.8
    radius = 5.
    radius_inner = 8.
    radius_outer = 11.
    data = mk_data_array(shape, value, background, xcenter, ycenter, radius)
    datamodel = datamodels.CubeModel(data)
    set_meta(datamodel, sub64p=False)
    datamodel.meta.bunit_data = 'MJy'  # unexpected data unit

    with pytest.raises(ValueError):
        tso_aperture_photometry(datamodel, xcenter, ycenter, radius + 1., radius_inner,
                                radius_outer)
