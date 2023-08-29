#
#  Module for applying MIRI EMI correction
#

import numpy as np
import logging

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


def reffile_waveform(input_model, emicorr_ref_file, save_intermediary_resutls,
        ret_phased_wav=True):
    """ Read or create the reference file from the input science data.

    Parameters
    ----------
    input_model: JWST data model
        input science data model to be emi-corrected

    emicorr_ref_file : str
        Path and name of the reference file

    save_intermediate_results : bool
        Saves the output into a file and the reference file (if created on-the-fly)

    ret_phased_wav : bool
        If True, the phased wave from this data set will be returned

    Returns
    -------
    emicorr_ref_file: string
        Path and name of the reference file just created

    phased_wav : float
        The phased wave from this data set (only if ret_phased_wav=True)
    """

    if emicorr_ref_file != 'N/A':
        # read the reference file
    else:
        log.info('No emicorr CRDS reference file found, one will be created from the data.')

