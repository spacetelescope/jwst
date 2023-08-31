#
#  Module for applying MIRI EMI correction
#

import numpy as np
import logging
from astropy.io import fits

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, emicorr_model, **pars):
    """
    EMI-correct a JWST data model using an emicorr model

    Parameters
    ----------
    input_model : JWST data model
        input science data model to be emi-corrected

    emicorr_model : JWST data model
        data model containing emi correction

    pars : dict
        Optional user-specified parameters to modify how outlier_detection
        will operate.  Valid parameters include:
            save_intermediate_results - saves the output into a file and the
                                        reference file (if created on-the-fly)
            user_supplied_reffile - reference file supplied by the user

    Returns
    -------
    output_model : JWST data model
        emi-corrected science data model

    """
    save_intermediate_results = pars['save_intermediate_results']
    user_supplied_reffile = pars['user_supplied_reffile']
    output_model = apply_emicorr(input_model, emicorr_model,
                        save_intermediate_results=save_intermediate_results,
                        user_supplied_reffile=user_supplied_reffile)
    output_model.meta.cal_step.emicorr = 'COMPLETE'

    return output_model


def apply_emicorr(input_model, emicorr_model, save_intermediate_results=False,
        user_supplied_reffile=None):
    """
    EMI-corrects data and error arrays by the following procedure:
         [repeat recursively for each discrete emi frequency desired]
    1) Read image data.
    2) Make very crude slope image and fixed pattern "super"bias for each
       integration, ignoring everything (nonlin, saturation, badpix, etc).
    3) Subtract scaled slope image and bias from each frame of each integration.
    4) Calculate phase of every pixel in the image at the desired emi frequency
       (e.g. 390 Hz) relative to the first pixel in the image.
    5) Make a binned, phased amplitude (pa) wave from the cleaned data (plot
       cleaned vs phase, then bin by phase).
    6)* Measure the phase shift between the binned pa wave and the input
       reference wave
    7)* Use look-up table to get the aligned reference wave value for each pixel
       (the noise array corresponding to the input image).

     6a-7a) Alternately, use the binned pa wave instead of the reference wave to
       "self-correct"

    8) Subtract the noise array from the input image and return the cleaned result.

         [repeat for next frequency using cleaned output as input]
      [removing highest amplitude emi first will probably give best results]

    Parameters
    ----------
    input_model : JWST data model
        input science data model to be emi-corrected

    emicorr_model : JWST data model
         ImageModel of emi

    save_intermediate_results : bool
        Saves the output into a file and the reference file (if created on-the-fly)

    user_supplied_reffile : str
        Reference file supplied by the user

    Returns
    -------
    output_model : JWST data model
        input science data model which has been emi-corrected
    """

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

    # Read image data
    data = input_model.data.copy()

    # Make very crude slope image and fixed pattern "super"bias for each
    # integration, ignoring everything (nonlin, saturation, badpix, etc)
    img_slope()
    fix_pattern_superbias()

    # Subtract scaled slope image and bias from each frame of each integration
    subtract_scaled_slope()

    # Calculate phase of every pixel in the image at the desired emi frequency
    # (e.g. 390 Hz) relative to the first pixel in the image
    calc_phase()

    # Make a binned, phased amplitude (pa) wave from the cleaned data (plot
    # cleaned vs phase, then bin by phase)
    bin_pa_wave()

    # Measure the phase shift between the binned pa wave and the input
    # reference wave
    measure_phase_shift()

    # Use look-up table to get the aligned reference wave value for each pixel
    # (the noise array corresponding to the input image)
    get_ref_wave()

    # Subtract the noise array from the input image and return the cleaned result
    subtract_img_noise()

    # DQ values remain the same as the input
    log.info('DQ values in the reference file NOT used to update the output DQ.')

    return output_model


def mk_reffile_waveform(input_model, emicorr_ref_filename):
    """ Create the reference file from the input science data.

    Parameters
    ----------
    input_model: JWST data model
        Input science data model to be emi-corrected

    emicorr_ref_filename : str
        Name of the reference file

    Returns
    -------
        Nothing
    """
    # initialize the reference fits file
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())

    # image
    #e = fits.ImageHDU(arr, name='extension_name')

    # table
    #col_names = ['integration_number', 'int_start_MJD_UTC', 'int_mid_MJD_UTC',
    #             'int_end_MJD_UTC', 'int_start_BJD_TDB', 'int_mid_BJD_TDB',
    #             'int_end_BJD_TDB']
    #c1 = fits.Column(name=col_names[0], array=np.array([]), format='D')
    # e = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7], name='INT_TIMES')

    # add extension
    #hdulist.append(e)

    # use the reference file model to format and write file
    model = datamodels.EmiModel(hdulist)
    model.save(emicorr_ref_filename)
    log.debug('On-the-fly reference file written as: %s', emicorr_ref_filename)


