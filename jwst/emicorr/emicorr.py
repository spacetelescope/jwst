#
#  Module for applying MIRI EMI correction
#

import numpy as np
import logging
from astropy.io import fits

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


subarray_cases = {

    # 390Hz out-of-phase - these all need 390hz correction

    "SLITLESSPRISM": {
        "rowclocks": 28,
        "frameclocks": 15904
    },

    "MASKLYOT": {
        "rowclocks": 90,
        "frameclocks": 32400
    },

    "SUB64": {
        "rowclocks": 28,
        "frameclocks": 8512
    },

    "SUB128": {
        "rowclocks": 44,
        "frameclocks": 11904
    },

    "MASK1140": {
        "rowclocks": 82,
        "frameclocks": 23968
    },

    "MASK1550": {
        "rowclocks": 82,
        "frameclocks": 23968
    },

    # 390Hz already in-phase for these, but may need corr for other
    # frequencies (e.g. 10Hz heater noise)

    "FULL_FASTR1": {
        "rowclocks": 271,
        "frameclocks": 277504
    },

    "FULL_SLOWR1": {
        "rowclocks": 2333,
        "frameclocks": 2388992
    },

    "BRIGHTSKY": {
        "rowclocks": 162,
        "frameclocks": 86528
    },

    "SUB256": {
        "rowclocks": 96,
        "frameclocks": 29952
    }

}


def do_correction(input_model, emicorr_model, **pars):
    """
    EMI-correct a JWST data model using an emicorr model

    Parameters
    ----------
    input_model : JWST data model
        Input science data model to be emi-corrected

    emicorr_model : JWST data model
        Data model containing emi correction

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
    # get the subarray case
    subarray = input_model.meta.subarray.name
    readpatt = input_model.meta.exposure.readpatt
    subname, rowclocks, frameclocks = get_subarcase(subarray, readpatt)
    if rowclocks is None:
        return subname

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

    # Read image data and set up some variables
    data = output_model.data
    nints, ngroups, ny, nx = np.shape(data)
    dd_all = np.ones((nints, ngroups, ny, nx/4))
    # s[0] = 4 in idl
    # s[1] = nx
    # s[2] = ny
    # s[3] = ngroups
    # s[4] = nints
    times = np.ones((nints, ngroups, ny, nx/4), dtype='uint64')

    # Find the phase
    counter = 0
    for inti in range(nints):
        log.info('Working on integration: {}'.format(inti))

        # Remove source signal and fixed bias from each integration ramp
        # (linear is good enough for phase finding)

        # do linear fit for source + sky
        s0 = sloper(data[inti, ngroups-2, :, :], intercept=mm0)

        # subtract source+sky from each frame of this ramp
        for groupi in ngroups:
            data[inti, groupi, ...] = data[inti, groupi, ...] - (s0 * groupi)

        # make a self-superbias
        m0 = minmed(data[inti, ngroups-2, :, :])

        # subtract self-superbias from each frame of this ramp
        for groupi in ngroups:
            data[inti, groupi, ...] = data[inti, groupi, ...] - m0

            # de-interleave each frame into the 4 separate output channels and
            # average (or median) them together for S/N
            tmp0 = np.ones((inti, groupi, :, (nx/4)*4+0))
            d0 = data[tmp0]
            tmp1 = np.ones((inti, groupi, :, (nx/4)*4+1))
            d1 = data[tmp1]
            tmp2 = np.ones((inti, groupi, :, (nx/4)*4+2))
            d2 = data[tmp2]
            tmp3 = np.ones((inti, groupi, :, (nx/4)*4+3))
            d3 = data[tmp3]
            dd = (d0 + d1 + d2 + d3)/4.

            # fix a bad ref col

    # DQ values remain the same as the input
    log.info('DQ values remain the same as input data.')

    return output_model


def sloper():
    pass


def minmed():
    pass


def get_subarcase(subarray, readpatt):
    """ Get the rowclocks and frameclocks values for the given configuration.

    Parameters
    ----------
    subarray : str
        Keyword value

    readpatt : str
        Keyword value

    Returns
    -------
    subname : str
        Modified subarray name

    rowclocks : int

    frameclocks : int

    """
    subname = subarray
    if subname == 'FULL':
        subname = subname + '_' + readpatt
    if subname in subarray_cases:
        rowclocks = subarray_cases[subname]["rowclocks"]
        frameclocks = subarray_cases[subname]["frameclocks"]
        return subname, rowclocks, frameclocks
    else:
        return subname, None, None


def mk_reffile_waveform(input_model, emicorr_ref_filename, save_mdl=False):
    """ Create the reference file from the input science data.

    Parameters
    ----------
    input_model : JWST data model
        Input science data model to be emi-corrected

    emicorr_ref_filename : str or None
        Name of the reference file

    save_mdl : bool
        Save the on-the-fly reference file

    Returns
    -------
    return_mdl: JWST data model
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
    if save_mdl:
        model.save(emicorr_ref_filename)
        log.debug('On-the-fly reference file written as: %s', emicorr_ref_filename)

    return model


