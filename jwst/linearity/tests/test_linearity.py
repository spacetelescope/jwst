"""Unit tests for the non-linearity correction.

Authors:
    M. Cracraft
"""
from jwst.linearity import LinearityStep
from jwst.linearity.linearity import do_correction as lincorr
from jwst.datamodels import dqflags, LinearityModel, MIRIRampModel
from astropy.io import fits
import numpy as np


def test_coeff_dq():
    """test linearity algorithm with random data ramp (does algorithm match expected algorithm)
    also test a variety of dq flags and expected output """

    # size of integration
    nints = 1
    ngroups = 160
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)
    im.meta.instrument.detector = 'MIRIMAGE'

    # Create reference file
    numcoeffs = 5

    ref_model = LinearityModel()
    ref_model.data = np.zeros(shape=(numcoeffs, ysize, xsize), dtype=np.float32)
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
    L0 = 0
    L1 = 0.85
    L2 = 4.62E-6
    L3 = -6.16E-11
    L4 = 7.23E-16

    coeffs = np.asfarray([0.0e+00,  0.85, 4.62e-06,  -6.16e-11, 7.23e-16])

    ref_model.data[:, 300, 500] = coeffs
    print(ref_model.data[:, 300, 500])

    # check behavior with NaN coefficients: should not alter pixel values
    coeffs2 = np.asfarray([L0, np.nan, L2, L3, L4])
    ref_model.data[:, 200, 500] = coeffs2
    im.data[0, 50, 200, 500] = 500.0

    tgroup = 2.775

    # set pixel values (DN) for specific pixels up the ramp
    im.data[0, :, 300, 500] = np.arange(ngroups) * 100 * tgroup

    scival = 40000.0
    im.data[0, 45, 300, 500] = scival  # to check linearity multiplication is done correctly
    im.data[0, 30, 350, 360] = 35  # pixel to check that dq=2 meant no correction was applied

    # check if dq flags in pixeldq are correctly populated in output
    im.pixeldq[50, 40] = 1
    im.pixeldq[50, 41] = 2
    im.pixeldq[50, 42] = 1024
    im.pixeldq[50, 43] = 2048

    # set dq flags in DQ of reference file
    ref_model.dq[350, 350] = 1  # DO_NOT_USE
    ref_model.dq[350, 360] = dqflags.pixel['NO_LIN_CORR']  # NO_LIN_CORR
    ref_model.dq[300, 500] = 0  # good pixel

    # run through Linearity pipeline
    outfile = lincorr(im, ref_model)

    # check that multiplication of polynomial was done correctly for specified pixel
    outval = L0+(L1*scival)+(L2*scival**2)+(L3*scival**3)+(L4*scival**4)

    print(outfile.data[0, 45, 300, 500])
    assert(np.isclose(outfile.data[0, 45, 300, 500], outval, rtol=0.1))

    # check that dq value was handled correctly

    assert(outfile.pixeldq[350, 350] == 1)
    assert(outfile.pixeldq[350, 360] == dqflags.pixel['NO_LIN_CORR'])  # NO_LIN_CORR flag
    assert(outfile.data[0, 30, 350, 360] == 35)  # NO_LIN_CORR, sci value should not change
    assert(outfile.data[0, 50, 200, 500] == 500.0)  # NaN coefficient should not change data value


def test_saturation():
    """Check that correction is not applied for groups flagged as SATURATED in GROUPDQ"""

    # size of integration
    nints = 1
    ngroups = 160
    xsize = 1032
    ysize = 1024

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)
    im.meta.instrument.detector = 'MIRIMAGE'

    # set groupdq pixels to saturated
    im.groupdq[0, 100:, 200, 150] = dqflags.pixel['SATURATED']  # saturated dq flag
    im.data[0, 150, 200, 150] = 1000.0  # value shouldn't change

    # run through Linearity pipeline
    outfile = LinearityStep.call(im)

    assert(outfile.data[0, 150, 200, 150] == 1000.0)  # pixel flagged as saturated, shouldn't change


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
    dq = np.zeros((1024, 1032), dtype=int)
    numcoeffs = 3

    # set reference file DQ to 'NO_LIN_CORR'
    dq[500, 500] = dqflags.pixel['NO_LIN_CORR']

    ref_model = LinearityModel()
    ref_model.data = np.zeros(shape=(numcoeffs, 1024, 1032), dtype=np.float32)
    ref_model.dq = dq

    ref_model.meta.instrument.name = 'MIRI'
    ref_model.meta.instrument.detector = 'MIRIMAGE'
    ref_model.meta.subarray.xstart = 1
    ref_model.meta.subarray.xsize = 1032
    ref_model.meta.subarray.ystart = 1
    ref_model.meta.subarray.ysize = 1024

    # run through pipeline (saturation and linearity steps)
    outfile = lincorr(im, ref_model)

    assert(outfile.pixeldq[500, 500] == dqflags.pixel['NO_LIN_CORR'])
    assert(outfile.data[0, 5, 500, 500] == 35)  # NO_LIN_CORR, sci value should not change


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
    dq = np.zeros((1024, 1032), dtype=int)
    numcoeffs = 3

    # set PIXELDQ to 'NO_LIN_CORR'
    dq[500, 500] = dqflags.pixel['NO_LIN_CORR']
    dq[550, 550] = dqflags.pixel['DO_NOT_USE']
    dq[560, 550] = dqflags.pixel['HOT']
    dq[550, 560] = dqflags.pixel['DEAD']
    dq[500, 300] = 2049  # Hot and DO_NOT_USE

    ref_model = LinearityModel()
    ref_model.data = np.zeros(shape=(numcoeffs, 1024, 1032), dtype=np.float32)
    ref_model.dq = dq

    ref_model.meta.instrument.name = 'MIRI'
    ref_model.meta.instrument.detector = 'MIRIMAGE'
    ref_model.meta.subarray.xstart = 1
    ref_model.meta.subarray.xsize = 1032
    ref_model.meta.subarray.ystart = 1
    ref_model.meta.subarray.ysize = 1024

    # run through linearity correction
    outfile = lincorr(im, ref_model)

    assert(outfile.pixeldq[500, 500] == dqflags.pixel['NO_LIN_CORR'])
    assert(outfile.pixeldq[550, 550] == dqflags.pixel['DO_NOT_USE'])
    assert(outfile.pixeldq[560, 550] == dqflags.pixel['HOT'])
    assert(outfile.pixeldq[550, 560] == dqflags.pixel['DEAD'])
    assert(outfile.pixeldq[500, 300] == 2049)


def test_lin_subarray():
    """Test that the pipeline properly extracts the subarray from the reference file.
    put dq flags in specific pixels and make sure they match in the output subarray file"""

    # create input data
    # create model of data with 0 value array
    ngroups = 50
    ysize = 224
    xsize = 288

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    im = MIRIRampModel(data=data, pixeldq=pixeldq, groupdq=groupdq)

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

    ref_model = LinearityModel()
    ref_model.data = np.zeros(shape=(numcoeffs, 1024, 1032), dtype=np.float32)
    ref_model.dq = dq

    ref_model.meta.instrument.name = 'MIRI'
    ref_model.meta.instrument.detector = 'MIRIMAGE'
    ref_model.meta.subarray.xstart = 1
    ref_model.meta.subarray.xsize = 1032
    ref_model.meta.subarray.ystart = 1
    ref_model.meta.subarray.ysize = 1024

    # run through pipeline
    outfile = lincorr(im, ref_model)

    # read dq array
    outpixdq = outfile.pixeldq

    # check for dq flag in pixeldq of subarray image
    assert(outpixdq[76, 100] == 1)
    assert(outpixdq[76, 104] == 1)


def test_wave_f2100w():
    """make sure correct filter (F2100W) reference file is called"""

    # create input data
    # create model of data with 0 value array
    nints = 1
    ngroups = 50
    ysize = 1024
    xsize = 1032

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)

    im.meta.instrument.detector = 'MIRIMAGE'
    im.meta.instrument.filter = 'F2100W'
    im.meta.instrument.band = 'N/A'
    im.meta.exposure.type = 'MIR_IMAGE'

    # set filter
    filter = 'F2100W'

    # run pipeline
    outfile = LinearityStep.call(im)

    # get linearity reference file from header
    linfile = outfile.meta.ref_file.linearity.name
    print('Linfile', linfile)

    # parse linfile name to discard crds://
    file = linfile.split("/")[2]
    print('split file', file)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/'+file
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    filt = refhead['FILTER']
    print(filt)
    print(filter)
   
    # test if filters match
    assert(filt == filter)


def test_wave_f2300c():
    """make sure correct filter (F2300C) reference file is called"""

    # create input data
    # create model of data with 0 value array
    nints = 1
    ngroups = 50
    ysize = 1024
    xsize = 1032

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)

    im.meta.instrument.detector = 'MIRIMAGE'
    im.meta.instrument.filter = 'F2300C'
    im.meta.instrument.band = 'N/A'
    im.meta.exposure.type = 'MIR_IMAGE'

    # set filter
    filter = 'F2300C'

    # run pipeline
    outfile = LinearityStep.call(im)

    # get linearity reference file from header
    linfile = outfile.meta.ref_file.linearity.name
    print('Linfile', linfile)

    # parse linfile name to discard crds://
    file = linfile.split("/")[2]
    print('split file', file)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + file
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    filt = refhead['FILTER']
    print(filt)
    print(filter)

    # test if filters match
    assert (filt == filter)


def test_wave_f2550w():
    """make sure correct filter (F2550W) reference file is called"""

    # create input data
    # create model of data with 0 value array
    nints = 1
    ngroups = 50
    ysize = 1024
    xsize = 1032

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)

    im.meta.instrument.detector = 'MIRIMAGE'
    im.meta.instrument.filter = 'F2550W'
    im.meta.instrument.band = 'N/A'
    im.meta.exposure.type = 'MIR_IMAGE'

    # set filter
    filter = 'F2550W'

    # run pipeline
    outfile = LinearityStep.call(im)

    # get linearity reference file from header
    linfile = outfile.meta.ref_file.linearity.name
    print('Linfile', linfile)

    # parse linfile name to discard crds://
    file = linfile.split("/")[2]
    print('split file', file)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + file
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    filt = refhead['FILTER']
    print(filt)
    print(filter)

    # test if filters match
    assert (filt == filter)


def test_err_array():
    """test that the error array is not changed by the linearity step"""

    # size of integration
    ngroups = 20
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)
    err = np.full(csize, 2.0)

    # create a JWST datamodel for MIRI data
    im = MIRIRampModel(data=data, pixeldq=pixeldq, groupdq=groupdq, err=err)
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

    np.testing.assert_array_equal(np.full(csize, 2.0, dtype=float),
                                  outfile.err, err_msg='no changes should be seen in array ')


def test_ifulong_short():
    """make sure correct (MIRIFULONG-SHORT) reference file is called"""

    # create input data
    # create model of data with 0 value array
    nints = 1
    ngroups = 50
    ysize = 1024
    xsize = 1032

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)

    im.meta.instrument.detector = 'MIRIFULONG'
    im.meta.instrument.filter = 'F2100W'
    im.meta.instrument.band = 'SHORT'
    im.meta.exposure.type = 'MIR_MRS'

    # set band, det
    set_band = 'SHORT'
    set_det = 'MIRIFULONG'

    # run pipeline
    outfile = LinearityStep.call(im)

    # get linearity reference file from header
    linfile = outfile.meta.ref_file.linearity.name
    print('Linfile', linfile)

    # parse linfile name to discard crds://
    file = linfile.split("/")[2]
    print('split file', file)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + file
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    file_det = refhead['DETECTOR']
    file_band = refhead['BAND']

    # test if values match
    assert set_band == file_band
    assert set_det == file_det


def test_ifulong_med():
    """make sure correct (MIRIFULONG-MEDIUM) reference file is called"""

    # create input data
    # create model of data with 0 value array
    nints = 1
    ngroups = 50
    ysize = 1024
    xsize = 1032

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)

    im.meta.instrument.detector = 'MIRIFULONG'
    im.meta.instrument.filter = 'F2100W'
    im.meta.instrument.band = 'MEDIUM'
    im.meta.exposure.type = 'MIR_MRS'

    # set band, det
    set_band = 'MEDIUM'
    set_det = 'MIRIFULONG'

    # run pipeline
    outfile = LinearityStep.call(im)

    # get linearity reference file from header
    linfile = outfile.meta.ref_file.linearity.name
    print('Linfile', linfile)

    # parse linfile name to discard crds://
    file = linfile.split("/")[2]
    print('split file', file)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + file
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    file_det = refhead['DETECTOR']
    file_band = refhead['BAND']

    # test if values match
    assert set_band == file_band
    assert set_det == file_det


def test_ifulong_long():
    """make sure correct (MIRIFULONG-LONG) reference file is called"""

    # create input data
    # create model of data with 0 value array
    nints = 1
    ngroups = 50
    ysize = 1024
    xsize = 1032

    # create a JWST datamodel for MIRI data
    im = make_rampmodel(nints, ngroups, ysize, xsize)

    im.meta.instrument.detector = 'MIRIFULONG'
    im.meta.instrument.filter = 'F2100W'
    im.meta.instrument.band = 'LONG'
    im.meta.exposure.type = 'MIR_MRS'

    # set band, det
    set_band = 'LONG'
    set_det = 'MIRIFULONG'

    # run pipeline
    outfile = LinearityStep.call(im)

    # get linearity reference file from header
    linfile = outfile.meta.ref_file.linearity.name
    print('Linfile', linfile)

    # parse linfile name to discard crds://
    file = linfile.split("/")[2]
    print('split file', file)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + file
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    file_det = refhead['DETECTOR']
    file_band = refhead['BAND']

    # test if values match
    assert set_band == file_band
    assert set_det == file_det


def make_rampmodel(nints, ngroups, ysize, xsize):
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data, pixeldq=pixeldq, groupdq=groupdq)

    dm_ramp.meta.instrument.name = 'MIRI'
    dm_ramp.meta.observation.date = '2018-01-01'
    dm_ramp.meta.observation.time = '00:00:00'
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp
