"""Unit tests for the non-linearity correction.

Authors:
    M. Cracraft
"""
from stdatamodels.jwst.datamodels import dqflags, LinearityModel, RampModel

from jwst.linearity import LinearityStep
from jwst.linearity.linearity import do_correction as lincorr
import numpy as np


def test_coeff_dq():
    """Test linearity algorithm with random data ramp (does algorithm match expected algorithm)
    also test a variety of dq flags and expected output """

    # size of integration
    nints = 1
    ngroups = 160
    xsize = 103
    ysize = 102

    # create the data and groupdq arrays
    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)
    im.meta.instrument.detector = 'MIRIMAGE'

    # Create reference file
    numcoeffs = 5

    ref_model = LinearityModel((numcoeffs, ysize, xsize))
    ref_model.dq = np.zeros((ysize, xsize), dtype=int)

    ref_model.meta.instrument.name = 'MIRI'
    ref_model.meta.instrument.detector = 'MIRIMAGE'
    ref_model.meta.subarray.xstart = 1
    ref_model.meta.subarray.ystart = 1
    ref_model.meta.subarray.xsize = xsize
    ref_model.meta.subarray.ysize = ysize

    # Set coefficient values in reference file to check the algorithm
    # Equation is DNcorr = L0 + L1*DN(i) + L2*DN(i)^2 + L3*DN(i)^3 + L4*DN(i)^4
    # DN(i) = signal in pixel, Ln = coefficient from ref file
    # L0 = 0 for all pixels for CDP6

    coeffs = np.asarray([0.0e+00, 0.85, 4.62e-06, -6.16e-11, 7.23e-16], dtype=float)

    # pixels to test using default coeffs
    ref_model.coeffs[:, 30, 50] = coeffs
    ref_model.coeffs[:, 35, 35] = coeffs
    ref_model.coeffs[:, 35, 36] = coeffs
    L0 = 0
    L1 = 0.85
    L2 = 4.62E-6
    L3 = -6.16E-11
    L4 = 7.23E-16

    # check behavior with NaN coefficients: should not alter pixel values
    coeffs2 = np.asarray([L0, np.nan, L2, L3, L4], dtype=float)
    ref_model.coeffs[:, 20, 50] = coeffs2
    im.data[0, 50, 20, 50] = 500.0

    # test case where all coefficients are zero, the linearity reference file may
    # not mark these pixels as NO_LIN_CORR. The code will mark these pixels as
    # NO_LIN_CORR
    ref_model.coeffs[:, 25, 25] = 0.0
    im.data[0, 50, 25, 25] = 600.0

    tgroup = 2.775

    # set pixel values (DN) for specific pixels up the ramp
    im.data[0, :, 30, 50] = np.arange(ngroups) * 100 * tgroup

    scival = 40000.0
    im.data[0, 45, 30, 50] = scival  # to check linearity multiplication is done correctly
    im.data[0, 30, 35, 36] = 35  # pixel to check that dq=2 meant no correction was applied

    # check if dq flags in pixeldq are correctly populated in output
    im.pixeldq[50, 40] = dqflags.pixel['DO_NOT_USE']
    im.pixeldq[50, 41] = dqflags.pixel['SATURATED']
    im.pixeldq[50, 42] = dqflags.pixel['DEAD']
    im.pixeldq[50, 43] = dqflags.pixel['HOT']

    # set dq flags in DQ of reference file
    ref_model.dq[35, 35] = dqflags.pixel['DO_NOT_USE']
    ref_model.dq[35, 36] = dqflags.pixel['NO_LIN_CORR']
    ref_model.dq[30, 50] = dqflags.pixel['GOOD']
    ref_model.dq[25, 25] = dqflags.pixel['GOOD']  # Testing the linerity sets this to NO_LIN_CORR

    # run through Linearity pipeline
    outfile = lincorr(im, ref_model)

    # check that multiplication of polynomial was done correctly for specified pixel
    outval = L0 + (L1 * scival) + (L2 * scival**2) + (L3 * scival**3) + (L4 * scival**4)
    assert np.isclose(outfile.data[0, 45, 30, 50], outval, rtol=0.00001)

    # check that dq value was handled correctly
    assert outfile.pixeldq[35, 35] == dqflags.pixel['DO_NOT_USE']
    assert outfile.pixeldq[35, 36] == dqflags.pixel['NO_LIN_CORR']
    assert outfile.pixeldq[25, 25] == dqflags.pixel['NO_LIN_CORR']
    # NO_LIN_CORR, sci value should not change
    assert outfile.data[0, 30, 35, 36] == 35
    # NaN coefficient should not change data value
    assert outfile.data[0, 50, 20, 50] == 500.0
    # coefficients all zero should not change data value
    assert outfile.data[0, 50, 25, 25] == 600.0


def test_saturation():
    """Check that correction is not applied for groups flagged as SATURATED in GROUPDQ"""

    # size of integration
    nints = 1
    ngroups = 16
    xsize = 1032
    ysize = 1024

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)
    im.meta.instrument.detector = 'MIRIMAGE'

    # set groupdq pixels to saturated
    im.groupdq[0, 10:, 200, 150] = dqflags.pixel['SATURATED']  # saturated dq flag
    im.data[0, 15, 200, 150] = 1000.0  # value shouldn't change

    # run through Linearity pipeline
    outfile = LinearityStep.call(im)

    # pixel flagged as saturated, shouldn't change
    assert outfile.data[0, 15, 200, 150] == 1000.0


def test_nolincorr():
    """Check that correction is not applied for pixels flagged NO_LIN_CORR in the DQ array
    of the reference file."""

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)
    im.meta.instrument.detector = 'MIRIMAGE'

    # set data value
    im.data[0, 5, 500, 500] = 35

    # Create reference file
    dq = np.zeros((ysize, xsize), dtype=int)
    numcoeffs = 3

    # set reference file DQ to 'NO_LIN_CORR'
    dq[500, 500] = dqflags.pixel['NO_LIN_CORR']

    ref_model = LinearityModel((numcoeffs, 1024, 1032))
    ref_model.dq = dq

    ref_model.meta.instrument.name = 'MIRI'
    ref_model.meta.instrument.detector = 'MIRIMAGE'
    ref_model.meta.subarray.xstart = 1
    ref_model.meta.subarray.xsize = xsize
    ref_model.meta.subarray.ystart = 1
    ref_model.meta.subarray.ysize = ysize

    # run through pipeline (saturation and linearity steps)
    outfile = lincorr(im, ref_model)

    assert outfile.pixeldq[500, 500] == dqflags.pixel['NO_LIN_CORR']
    assert outfile.data[0, 5, 500, 500] == 35  # NO_LIN_CORR, sci value should not change


def test_pixeldqprop():
    """Check that linearity reference file DQ values are populated into the PIXELDQ array of the data
    This does not fully test mapping, as we have not provided a dq_def to define the
    mapping of dq flags in the reference file. it's a straightforward mapping (1-1)"""

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)
    im.meta.instrument.detector = 'MIRIMAGE'

    # Create reference file
    dq = np.zeros((ysize, xsize), dtype=int)
    numcoeffs = 5

    # set PIXELDQ to 'NO_LIN_CORR'
    dq[500, 500] = dqflags.pixel['NO_LIN_CORR']
    dq[550, 550] = dqflags.pixel['DO_NOT_USE']
    dq[560, 550] = dqflags.pixel['HOT']
    dq[550, 560] = dqflags.pixel['DEAD']
    dq[500, 300] = np.bitwise_or(dqflags.pixel['HOT'], dqflags.pixel['DO_NOT_USE'])

    ref_model = LinearityModel((numcoeffs, ysize, xsize))
    ref_model.dq = dq

    coeffs = np.asarray([0.0e+00, 0.85, 4.62e-06, -6.16e-11, 7.23e-16], dtype=float)

    # pixels to test using default coeffs.
    ref_model.coeffs[:, 550, 550] = coeffs
    ref_model.coeffs[:, 560, 550] = coeffs
    ref_model.coeffs[:, 550, 560] = coeffs
    ref_model.coeffs[:, 500, 300] = coeffs

    ref_model.meta.instrument.name = 'MIRI'
    ref_model.meta.instrument.detector = 'MIRIMAGE'
    ref_model.meta.subarray.xstart = 1
    ref_model.meta.subarray.xsize = xsize
    ref_model.meta.subarray.ystart = 1
    ref_model.meta.subarray.ysize = ysize

    # run through linearity correction
    outfile = lincorr(im, ref_model)

    assert outfile.pixeldq[500, 500] == dqflags.pixel['NO_LIN_CORR']
    assert outfile.pixeldq[550, 550] == dqflags.pixel['DO_NOT_USE']
    assert outfile.pixeldq[560, 550] == dqflags.pixel['HOT']
    assert outfile.pixeldq[550, 560] == dqflags.pixel['DEAD']
    assert outfile.pixeldq[500, 300] == np.bitwise_or(
        dqflags.pixel['HOT'], dqflags.pixel['DO_NOT_USE'])


def test_lin_subarray():
    """Test that the pipeline properly extracts the subarray from the reference file.
    put dq flags in specific pixels and make sure they match in the output subarray file"""

    # create input data
    # create model of data with 0 value array
    ngroups = 50
    ysize = 224
    xsize = 288

    # create a JWST datamodel for MIRI data
    im = RampModel((1, ngroups, ysize, xsize))
    im.data += 1

    im.meta.instrument.name = 'MIRI'
    im.meta.instrument.detector = 'MIRIMAGE'
    im.meta.instrument.filter = 'F1500W'
    im.meta.instrument.band = 'N/A'
    im.meta.observation.date = '2016-06-01'
    im.meta.observation.time = '00:00:00'
    im.meta.exposure.type = 'MIR_IMAGE'
    im.meta.subarray.name = 'MASK1550'
    im.meta.subarray.xstart = 1
    im.meta.subarray.xsize = xsize
    im.meta.subarray.ystart = 467
    im.meta.subarray.ysize = ysize

    # Read in reference file

    dq = np.zeros((1024, 1032), dtype=int)
    numcoeffs = 3

    # place dq flags in dq array that would be in subarray
    # MASK1550 file has colstart=1, rowstart=467
    dq[542, 100:105] = 1

    ref_model = LinearityModel((numcoeffs, 1024, 1032))
    # set all the linear terms =1, so it does not trip the check if
    # the linear terms = 0, which results in DQ of NO_LIN_CORR
    ref_model.coeffs[1,:,:] = 1
    ref_model.dq = dq

    ref_model.meta.instrument.name = 'MIRI'
    ref_model.meta.instrument.detector = 'MIRIMAGE'
    ref_model.meta.description = "MIRI LINEARITY Correction"
    ref_model.meta.reftype = "LINEARITY"
    ref_model.meta.author = "Monty Pytest"
    ref_model.meta.pedigree = "GROUND"
    ref_model.meta.useafter = '2015-08-01T00:00:00'
    ref_model.meta.subarray.xstart = 1
    ref_model.meta.subarray.xsize = 1032
    ref_model.meta.subarray.ystart = 1
    ref_model.meta.subarray.ysize = 1024

    # run through pipeline
    outfile = lincorr(im, ref_model)

    # read dq array
    outpixdq = outfile.pixeldq

    # check for dq flag in pixeldq of subarray image
    assert (outpixdq[76, 100] == 1)
    assert (outpixdq[76, 104] == 1)


def test_err_array():
    """Test that the error array is not changed by the linearity step"""

    # size of integration
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create a JWST datamodel for MIRI data
    im = RampModel((1, ngroups, ysize, xsize))
    im.data += 1
    im.err += 2
    # set file header values
    im.meta.instrument.detector = 'MIRIMAGE'
    im.meta.instrument.name = 'MIRI'
    im.meta.observation.date = '2018-01-01'
    im.meta.observation.time = '00:00:00'
    im.meta.subarray.xstart = 1
    im.meta.subarray.xsize = xsize
    im.meta.subarray.ystart = 1
    im.meta.subarray.ysize = ysize

    # run pipeline
    outfile = LinearityStep.call(im)

    # check output of error array
    # test that the science data are not changed
    np.testing.assert_allclose(im.err, outfile.err)


def make_rampmodel(nints, ngroups, ysize, xsize):
    """Function to provide ramp model to tests"""

    dm_ramp = RampModel((nints, ngroups, ysize, xsize))
    dm_ramp.data += 1

    dm_ramp.meta.instrument.name = 'MIRI'
    dm_ramp.meta.observation.date = '2018-01-01'
    dm_ramp.meta.observation.time = '00:00:00'
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp
