import pytest
from astropy.io import fits
from jwst.dq_init import DQInitStep
from jwst.dq_init.dq_initialization import do_dqinit
from jwst.datamodels import MIRIRampModel, MaskModel
from jwst.datamodels import ImageModel, dqflags
import numpy as np
import numpy.testing as nptest
import setup_models as mod

# jwst_miri_mask_0022.fits is MIRIFUSHORT
# jwst_miri_mask_0023.fits is MIRIMAGE
# jwst_miri_mask_0021.fits is MIRIFULONG


def test_dq_im():
    # Check that PIXELDQ is initialized with the information from the reference file.
    # test that a flagged value in the reference file flags the PIXELDQ array
    # read in reference mask file
    #hduref = fits.open('/grp/crds/jwst/references/jwst/jwst_miri_mask_0023.fits')
    #filedir = '/ifs/jwst/wit/miri/pipelinetests/cracraft/build7_x/'
    
    with MaskModel('/grp/crds/jwst/references/jwst/jwst_miri_mask_0023.fits') as ref_data:
        pass
    # size of integration
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = mod.make_rawramp(nints, ngroups, ysize, xsize)

    #print(dqflags.pixel['DEAD'])
    # create a MaskModel for the dq input mask
    #ref_data = mod.make_maskmodel(ysize, xsize)

    # edit reference file with known bad pixel values
    #ref_data = hduref['DQ'].data
    ref_data.dq[100, 100] = 2   # Dead pixel
    ref_data.dq[200, 100] = 4   # Hot pixel
    ref_data.dq[300, 100] = 8   # Unreliable_slope
    ref_data.dq[400, 100] = 16  # RC
    ref_data.dq[500, 100] = 1   # Do_not_use
    ref_data.dq[100, 200] = 3   # Dead pixel + do not use
    ref_data.dq[200, 200] = 5   # Hot pixel + do not use
    ref_data.dq[300, 200] = 9   # Unreliable slope + do not use
    ref_data.dq[400, 200] = 17  # RC + do not use

    # write out reference file with modifications
    #hduref.writeto(filedir+'dq_image_testref.fits', overwrite=True)

    # run pipeline step. Pipeline step can take in science data model, but needs a string
    # to read in the modified reference file to override the mask, so file must be written 
    # first and then read in as an override.
    #outfile = DQInitStep.call(dm_ramp, override_mask=filedir+'dq_image_testref.fits')
    outfile = do_dqinit(dm_ramp, ref_data)
    print(outfile.err.shape)
    dqdata = outfile.pixeldq

    # assert that the pixels read back in match the mapping from ref data to science data
    assert(dqdata[100, 100] == dqflags.pixel['DEAD'])
    assert(dqdata[200, 100] == dqflags.pixel['HOT'])
    assert(dqdata[300, 100] == dqflags.pixel['UNRELIABLE_SLOPE'])
    assert(dqdata[400, 100] == dqflags.pixel['RC'])
    assert(dqdata[500, 100] == dqflags.pixel['DO_NOT_USE'])
    assert(dqdata[100, 200] == 1025)
    assert(dqdata[200, 200] == 2049)
    assert(dqdata[300, 200] == 16777217)
    assert(dqdata[400, 200] == 16385)


def test_groupdq():

    # Check that GROUPDQ extension is added to the data and all values are initialized to zero.
    # size of integration
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = mod.make_rawramp(nints, ngroups, ysize, xsize)

    # create a MaskModel for the dq input mask
    ref_data = mod.make_maskmodel(ysize, xsize)

    # run the correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # check that GROUPDQ was created and initialized to zero
    groupdq = outfile.groupdq

    np.testing.assert_array_equal(np.full((1, ngroups, ysize, xsize),
                                          0,
                                          dtype=int),
                                  groupdq,
                                  err_msg='groupdq not initialized to zero')


def test_err():
    # Check that a 4-D ERR array is initialized and all values are zero.
    # size of integration
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = mod.make_rawramp(nints, ngroups, ysize, xsize)

    # create a MaskModel for the dq input mask
    ref_data = mod.make_maskmodel(ysize, xsize)

    # run correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # check that ERR array was created and initialized to zero
    errarr = outfile.err  # should these dimensions match the data dimensions?
                          # errarr = array([], shape=(0, 0, 0, 0), dtype=float32)

    #print(len(errarr))
    assert(errarr.ndim == 4)  # check that output err array is 4-D
    assert(np.all(errarr == 0))  # check that values are 0


def test_dq_ifushort():
    # test that a flagged value in the reference file flags the PIXELDQ array
    hduref = fits.open('/grp/crds/jwst/references/jwst/jwst_miri_mask_0022.fits')
    filedir = '/ifs/jwst/wit/miri/pipelinetests/cracraft/build7_x/'
    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    
    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data)
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    # write reference file with known bad pixel values
    ref_data = hduref['DQ'].data
    ref_data[100, 100] = 2   # Dead pixel
    ref_data[200, 100] = 4   # Hot pixel
    ref_data[300, 100] = 8   # Unreliable_slope
    ref_data[400, 100] = 16  # RC
    ref_data[500, 100] = 1   # Do_not_use
    ref_data[100, 200] = 3   # Dead pixel + do not use
    ref_data[200, 200] = 5   # Hot pixel + do not use
    ref_data[300, 200] = 9   # Unreliable slope + do not use
    ref_data[400, 200] = 17  # RC + do not use

    # write out modified reference file    
    hduref.writeto(filedir+'dq_ifushort_testref.fits', overwrite=True)

    # run pipeline step
    outfile = DQInitStep.call(dm_ramp, override_mask=filedir+'dq_ifushort_testref.fits')
    dqdata = outfile.pixeldq

    assert(dqdata[100, 100] == 1024)
    assert(dqdata[200, 100] == 2048)
    assert(dqdata[300, 100] == 16777216)
    assert(dqdata[400, 100] == 16384)
    assert(dqdata[500, 100] == 1)
    assert(dqdata[100, 200] == 1025)
    assert(dqdata[200, 200] == 2049)
    assert(dqdata[300, 200] == 16777217)
    assert(dqdata[400, 200] == 16385)

    # check that GROUPDQ was created and initialized to zero
    groupdq = outfile.groupdq

    np.testing.assert_array_equal(np.full((csize),
                                          0,
                                          dtype=int),
                                  groupdq,
                                  err_msg='groupdq not initialized to zero')


def test_dq_ifulong():
    # test that a flagged value in the reference file flags the PIXELDQ array
    hduref = fits.open('/grp/crds/jwst/references/jwst/jwst_miri_mask_0021.fits')
    filedir = '/ifs/jwst/wit/miri/pipelinetests/cracraft/build7_x/'
    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data)
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    # write reference file with known bad pixel values
    ref_data = hduref['DQ'].data
    ref_data[100, 100] = 2   # Dead pixel
    ref_data[200, 100] = 4   # Hot pixel
    ref_data[300, 100] = 8   # Unreliable_slope
    ref_data[400, 100] = 16  # RC
    ref_data[500, 100] = 1   # Do_not_use
    ref_data[100, 200] = 3   # Dead pixel + do not use
    ref_data[200, 200] = 5   # Hot pixel + do not use
    ref_data[300, 200] = 9   # Unreliable slope + do not use
    ref_data[400, 200] = 17  # RC + do not use

    # write out modified reference file
    hduref.writeto(filedir+'dq_ifulong_testref.fits', overwrite=True)
    
    # run pipeline step
    outfile = DQInitStep.call(dm_ramp, override_mask=filedir+'dq_ifulong_testref.fits')
    dqdata = outfile.pixeldq

    assert(dqdata[100, 100] == 1024)
    assert(dqdata[200, 100] == 2048)
    assert(dqdata[300, 100] == 16777216)
    assert(dqdata[400, 100] == 16384)
    assert(dqdata[500, 100] == 1)
    assert(dqdata[100, 200] == 1025)
    assert(dqdata[200, 200] == 2049)
    assert(dqdata[300, 200] == 16777217)
    assert(dqdata[400, 200] == 16385)

    # check that GROUPDQ was created and initialized to zero
    groupdq = outfile.groupdq

    np.testing.assert_array_equal(np.full((csize),
                                          0,
                                          dtype=int),
                                  groupdq,
                                  err_msg='groupdq not initialized to zero')


def test_dq_subarray():

    # Test that the pipeline properly extracts the subarray from the reference file.
    # put dq flags in specific pixels and make sure they match in the output subarray file

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
    refhdu = fits.open('/grp/crds/jwst/references/jwst/jwst_miri_mask_0023.fits')
    dq = refhdu['DQ'].data

    # place dq flags in dq array that would be in subarray
    # MASK1550 file has colstart=1, rowstart=467
    dq[542, 100] = 2
    dq[550, 100] = 1
    dq[580, 80] = 4

    # write ref file
    refhdu.writeto('testmaskref.fits', overwrite=True)

    # run through pipeline saturation step
    outfile = DQInitStep.call(im, override_mask='testmaskref.fits')

    # read dq array
    outpixdq = outfile.pixeldq

    # check for dq flag in pixeldq of subarray image
    assert(outpixdq[76, 100] == 1024)
    assert(outpixdq[84, 100] == 1)
    assert(outpixdq[114, 80] == 2048)  # check that pixel was flagged 'NO_SAT_CHECK'


@pytest.mark.skip("Getting ValueError: Array has wrong number of dimensions")
def test_dq_wrong_modeltype():

    # size of integration

    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (ysize, xsize)
    data = np.full(csize, 1.0)
    #groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    image = ImageModel(data=data)

    nptest.assert_raises(ValueError, DQInitStep.call(image)), 'Not TypeError raised'


def test_dq_add1_groupdq():
    """
    Test if the dq_init code set the groupdq flag on the first
    group to 'do_not_use' by adding 1 to the flag, not overwriting to 1
    Also test whether two flags on the same pixel are added together.
    """
    hduref = fits.open('/grp/crds/jwst/references/jwst/jwst_miri_mask_0021.fits')

    # size of integration
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data, pixeldq=pixeldq, groupdq=groupdq)

    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    #print(dm_ramp.pixeldq.shape())

    # write reference file with known bad pixel values
    ref_data = hduref['DQ'].data
    ref_data[505, 505] = 1   # Do_not_use
    ref_data[400, 500] = 3  # do_not_use and dead pixel

    # set a flag in the pixel dq - would an uncal file already contain a pixel dq array?
    dm_ramp.pixeldq[505, 505] = 4

    # write ref file
    hduref.writeto('testmaskref.fits', overwrite=True)

    # run the dq init correction step
    dm_ramp_dq = DQInitStep.call(dm_ramp, override_mask='testmaskref.fits')

    # test if pixels in pixeldq were incremented in value by 1
    assert(dm_ramp_dq.pixeldq[505, 505] == 5)  # check that previous dq flag is added to mask value
    assert(dm_ramp_dq.pixeldq[400, 500] == 1025)  # check two flags propagate correctly
