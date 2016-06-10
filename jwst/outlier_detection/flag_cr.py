from __future__ import absolute_import, unicode_literals, division, print_function
import importlib

import numpy as np
from scipy import ndimage
from stsci.tools import bitmask

from .. import datamodels
from .. import assign_wcs
# from ..stpipe import Step

from . import quickDeriv

CRBIT = np.uint32(datamodels.dqflags.pixel.get('JUMP_DET', 4))

def do_detection(input_models, blot_models, ref_filename, **pars):
    """
    Flags DQ array for cosmic rays in input images.

    The science frame in each ImageModel in input_models is compared to
    the corresponding blotted median image in blot_models.  The result is
    an updated DQ array in each ImageModel in input_models.

    Parameters
    ----------
    input_models: JWST ModelContainer object
        data model container holding science ImageModels, modified in place

    blot_models : JWST ModelContainer object
        data model container holding ImageModels of the median output frame
        blotted back to the wcs and frame of the ImageModels in input_models

    ref_filename : dict
        Contains JWST ModelContainers for 'gain' and 'readnoise' reference files

    Returns
    -------
    None
        The dq array in each input model is modified in place
    """

    #gain_models = build_reffile_container(input_models, 'gain')
    #rn_models = build_reffile_container(input_models, 'readnoise')

    gain_models = ref_filename['gain']
    rn_models = ref_filename['readnoise']

    for image, blot, gain, rn in zip(input_models, blot_models, gain_models,
        rn_models):
        flag_cr(image, blot, gain, rn, **pars)


#def build_reffile_container_func(input_models, reftype):
#    """
#    Return a ModelContainer of reference file datamodels.
#
#    Parameters
#    ----------
#
#    input_models: ModelContainer
#        the science data, ImageModels in a ModelContainer
#
#    reftype: string
#        type of reference file
#
#    Returns
#    -------
#
#    a ModelContainer with corresponding reference files for each input model
#    """
#
#    """
#    model_mapping = {'gain':'GainModel', 'readnoise':'ReadnoiseModel'}
#    module = importlib.import_module('..datamodels')
#    model_name = model_mapping[reftype]
#    ref_model = getattr(module, model_name)
#    print("ref_model defined as: {} of type {}".format(ref_model, type(ref_model)))
#    """
#    q = Step()
#    reffiles = [q.get_reference_file(im, reftype) for im in input_models]
#
#    # Check if all the ref files are the same.  If so build it by reading
#    # the reference file just once.
#    if len(set(reffiles)) <= 1:
#        length = len(input_models)
#        ref_list = [datamodels.open(reffiles[0])] * length
#        #ref_list = [ref_model(reffiles[0])] * length
#    else:
#        ref_list = [datamodels.open(ref) for ref in reffiles]
#        #ref_list = [ref_model(ref) for ref in reffiles]
#    return datamodels.ModelContainer(ref_list)


def buildMask(dqarr, bitvalue):
    """ Builds a bit-mask from an input DQ array and a bitvalue flag"""

    bitvalue = bitmask.interpret_bits_value(bitvalue)
    if bitvalue is None:
        return (np.ones(dqarr.shape, dtype=np.uint32))
    return np.logical_not(np.bitwise_and(dqarr, ~bitvalue)).astype(np.uint32)


def flag_cr(sci_image, blot_image, gain_image, readnoise_image, **pars):
    """
    Masks outliers in science image

    Mask blemishes in dithered data by comparing a science image
    with a model image and the derivative of the model image.

    Parameters
    ----------
    sci_image : ImageModel
        the science data

    blot_image : ImageModel
        the blotted median image of the dithered science frames

    gain_image : GainModel
        the 2-D gain array

    readnoise_image : ReadnoiseModel
        the 2-D read noise array

    pars : dict
        the user parameters for Outlier Detection

    Default parameters:

    grow     = 1               # Radius to mask [default=1 for 3x3]
    ctegrow  = 0               # Length of CTE correction to be applied
    snr      = "4.0 3.0"       # Signal-to-noise ratio
    scale    = "0.5 0.4"       # scaling factor applied to the derivative
    backg    = 0               # Background value
    """

    grow = pars.get('grow', 1)
    ctegrow = pars.get('ctegrow', 0)
    snr = pars.get('snr', '4.0 3.0').split()
    scale = pars.get('scale', '0.5 0.4').split()
    backg = pars.get('backg', 0)
    # Get necessary parameters from the meta tree
    try:
        subtracted_background = sci_image.meta.skybg
    except AttributeError:
        subtracted_background = backg
    try:
        exptime = float(sci_image.meta.exposure.exposure_time)
    except AttributeError:
        exptime = 100.

    input_image = sci_image.data * exptime
    blot_data = blot_image.data * exptime
    blot_deriv = quickDeriv.qderiv(blot_data)

    # # This mask can take into account any crbits values
    # # specified by the user to be ignored.
    # dq_mask = buildMask(sci_image.dq, CRBIT)

    #parse out the SNR and scaling information
    snr1=float(snr[0])
    snr2=float(snr[1])
    mult1 = float(scale[0])
    mult2 = float(scale[1])

    gain = gain_image.data
    read_noise = readnoise_image.data

    # Define output cosmic ray mask to populate
    cr_mask = np.zeros(sci_image.shape, dtype=np.uint8)

    # Set scaling factor to 1 since scaling has already been accounted for
    # in blotted image
    exp_mult = 1.

    ##################   COMPUTATION PART I    ###################
    # Create a temporary array mask
    __t1 = np.absolute(input_image - blot_data)
    __ta = np.sqrt(gain * np.absolute(blot_data * exp_mult +
        subtracted_background * exp_mult) + read_noise ** 2)
    __tb = ( mult1 * blot_deriv + snr1 * __ta / gain )
    del __ta
    __t2 = __tb / exp_mult
    del __tb
    __tmp1 = np.logical_not(np.greater(__t1, __t2))
    del __t1
    del __t2

    # Create a convolution kernel that is 3 x 3 of 1's
    kernel = np.ones((3,3), dtype=np.uint8)
    # Create an output tmp file the same size as the input temp mask array
    __tmp2 = np.zeros(__tmp1.shape, dtype=np.int16)
    # Convolve the mask with the kernel
    ndimage.convolve(__tmp1, kernel, output=__tmp2, mode='nearest', cval=0)
    del kernel
    del __tmp1

    ##################   COMPUTATION PART II    ###################
    # Create the CR Mask
    __xt1 = np.absolute(input_image - blot_data)
    __xta = np.sqrt(gain * np.absolute(blot_data * exp_mult +
        subtracted_background * exp_mult) + read_noise * read_noise)
    __xtb = ( mult2 *blot_deriv + snr2 * __xta / gain )
    del __xta
    __xt2 = __xtb / exp_mult
    del __xtb
    # Must use a bitwise 'and' to create the mask with numarray objects.
    np.logical_not(np.greater(__xt1, __xt2) & np.less(__tmp2, 9), cr_mask)

    del __xt1
    del __xt2
    del __tmp2

    ##################   COMPUTATION PART III    ###################
    # Flag additional cte 'radial' and 'tail' pixels surrounding CR
    # pixels as CRs

    # In both the 'radial' and 'length' kernels below, 0=good and
    # 1=bad, so that upon convolving the kernels with cr_mask, the
    # convolution output will have low->bad and high->good from which
    # 2 new arrays are created having 0->bad and 1->good. These 2 new
    # arrays are then AND'ed to create a new cr_mask.

    # recast cr_mask to int for manipulations below; will recast to
    # Bool at end
    cr_mask_orig_bool = cr_mask.copy()
    cr_mask = cr_mask_orig_bool.astype(np.int8)

    # make radial convolution kernel and convolve it with original cr_mask
    cr_grow_kernel = np.ones((grow, grow))
    cr_grow_kernel_conv = cr_mask.copy()
    ndimage.convolve(cr_mask, cr_grow_kernel, output=cr_grow_kernel_conv)

    # make tail convolution kernel and (shortly) convolve it with original cr_mask
    cr_ctegrow_kernel = np.zeros((2*ctegrow+1,2*ctegrow+1))
    cr_ctegrow_kernel_conv = cr_mask.copy()

    # which pixels are masked by tail kernel depends on readout direction
    # We could put useful info in here for CTE masking if needed.  Code
    # remains below.  For now, we set to zero, which turns off CTE masking.
    ctedir = 0
    if (ctedir == 1):  # HRC: amp C or D ; WFC: chip = sci,1 ; WFPC2
        cr_ctegrow_kernel[0:ctegrow, ctegrow] = 1    #  'positive' direction
    if (ctedir == -1): # HRC: amp A or B ; WFC: chip = sci,2
        cr_ctegrow_kernel[ctegrow+1:2*ctegrow+1, ctegrow] = 1 #'negative' direction
    if (ctedir == 0):  # NICMOS: no cte tail correction
        pass

    # finally do the tail convolution
    ndimage.convolve(cr_mask, cr_ctegrow_kernel, output=cr_ctegrow_kernel_conv)

    # select high pixels from both convolution outputs; then 'and' them to
    # create new cr_mask
    where_cr_grow_kernel_conv = np.where(cr_grow_kernel_conv < grow*grow, 0, 1)
    where_cr_ctegrow_kernel_conv = np.where(cr_ctegrow_kernel_conv < ctegrow, 0, 1)

    # combine masks and cast back to Bool
    np.logical_and(where_cr_ctegrow_kernel_conv, where_cr_grow_kernel_conv, cr_mask)
    cr_mask = cr_mask.astype(bool)

    # Update the DQ array in the input image
    np.bitwise_or(sci_image.dq, np.invert(cr_mask) * CRBIT, sci_image.dq)

    del cr_mask_orig_bool
    del cr_grow_kernel
    del cr_grow_kernel_conv
    del cr_ctegrow_kernel
    del cr_ctegrow_kernel_conv
    del where_cr_grow_kernel_conv
    del where_cr_ctegrow_kernel_conv

    # write out the updated file to disk to preserve the changes
    sci_image.save(sci_image.meta.filename)

    # # write out the dq array as a separate file
    # outfilename = sci_image.meta.filename.split('.')[0] + '_dq.fits'
    # out_dq = datamodels.ImageModel()
    # out_dq.data = result_dq
    # out_dq.to_fits(outfilename, clobber=True)

    # ######## Save the cosmic ray mask file to disk
    # _cr_file = np.zeros(input_image.shape, np.uint32)
    # _cr_file = np.where(cr_mask, 1, 0).astype(np.uint32)

    # _pf = util.createFile(_cr_file, outfile=outfile, header = None)
