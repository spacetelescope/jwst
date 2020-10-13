# Module for calculating pathloss correction for science data sets

import logging
import math
import numpy as np

from gwcs import wcstools

from jwst.assign_wcs import nirspec, util
import jwst.datamodels as datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# There are 30 slices in the NIRSpec IFU, numbered from 0 to 29
NIRSPEC_IFU_SLICES = np.arange(30)


def get_center(exp_type, input):
    """Get the center of the target in the aperture.
    (0.0, 0.0) is the aperture center.  Coordinates go
    from -0.5 to 0.5.
    """
    if exp_type == "NRS_IFU":

        # Currently assume IFU sources are centered
        return (0.0, 0.0)

    elif exp_type in ["NRS_MSASPEC", "NRS_FIXEDSLIT", "NRS_BRIGHTOBJ"]:

        # MSA centering is specified in the MiltiSlit model
        # "input" treated as a slit object
        try:
            xcenter = input.source_xpos
            ycenter = input.source_ypos
        except AttributeError:
            log.warning("Unable to get source center from model")
            log.warning("Using 0.0, 0.0")
            xcenter = 0.0
            ycenter = 0.0
        return (xcenter, ycenter)

    else:
        log.warning(f'No method to get centering for exp_type {exp_type}')
        log.warning("Using (0.0, 0.0)")
        return (0.0, 0.0)


def get_aperture_from_model(input_model, match):
    """Figure out the correct aperture based on the value of the 'match'
    parameter.  For MSA, match is the number of shutters, for fixed slit,
    match is the name of the slit
    """
    if input_model.meta.exposure.type == 'NRS_MSASPEC':
        for aperture in input_model.apertures:
            if aperture.shutters == match:
                return aperture
    elif input_model.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ',
                                            'NIS_SOSS']:
        for aperture in input_model.apertures:
            log.debug(aperture.name)
            if aperture.name == match:
                return aperture
    else:
        log.warning(f'Unable to get aperture from exp_type {input_model.meta.exposure.type}')

    # If nothing matches, return None
    return None


def calculate_pathloss_vector(pathloss_refdata,
                              pathloss_wcs,
                              xcenter,
                              ycenter):
    """
    Calculate the pathloss vectors from the pathloss model using the
    coordinates of the center of the target to interpolate the
    pathloss value as a function of wavelength at that location

    Parameters:
    -----------
    pathloss_refdata : numpy ndarray
        The input pathloss data array

    pathloss_wcs : wcs attribute from model

    xcenter : float
        The x-center of the target (-0.5 to 0.5)

    ycenter : float
        The y-center of the target (-0.5 to 0.5)

    Returns:
    --------
    wavelength : numpy ndarray
        The 1-d wavelength array

    pathloss : numpy ndarray
        The corresponding 1-d pathloss array

    is_inside_slitlet : bool
        Returns True if the object source position is inside the slitlet,
        otherwise returns False

    """
    is_inside_slitlet = True
    if len(pathloss_refdata.shape) == 3:
        wavesize, nrows, ncols = pathloss_refdata.shape
    else:
        wavesize = pathloss_refdata.shape[0]
    wavelength = np.zeros(wavesize, dtype=np.float32)
    pathloss_vector = np.zeros(wavesize, dtype=np.float32)

    # uniformsource.data is 1-d. We just return it, along with
    # a vector of wavelengths calculated using the WCS.
    # Uniform source is always inside the slitlet
    if len(pathloss_refdata.shape) == 1:
        crpix1 = pathloss_wcs.crpix1
        crval1 = pathloss_wcs.crval1
        cdelt1 = pathloss_wcs.cdelt1
        for i in np.arange(wavesize):
            wavelength[i] = crval1 + (float(i+1) - crpix1)*cdelt1
        return wavelength, pathloss_refdata, True

    # pointsource.data is 3-d, so we have to extract a wavelength vector
    # at the specified location.  We do this using bilinear interpolation
    else:
        crpix3 = pathloss_wcs.crpix3
        crval3 = pathloss_wcs.crval3
        cdelt3 = pathloss_wcs.cdelt3
        for i in np.arange(wavesize):
            wavelength[i] = crval3 +(float(i+1) - crpix3) * cdelt3

        # Calculate python index of object center
        crpix1 = pathloss_wcs.crpix1
        crval1 = pathloss_wcs.crval1
        cdelt1 = pathloss_wcs.cdelt1
        crpix2 = pathloss_wcs.crpix2
        crval2 = pathloss_wcs.crval2
        cdelt2 = pathloss_wcs.cdelt2
        object_colindex = crpix1 + (xcenter - crval1) / cdelt1 - 1
        object_rowindex = crpix2 + (ycenter - crval2) / cdelt2 - 1
        if (object_colindex < 0 or object_colindex >= (ncols - 1) or
            object_rowindex < 0 or object_rowindex >= (nrows - 1)):
            is_inside_slitlet = False
        else:

            # Do bilinear interpolation to get the array of
            # path loss vs wavelength
            dx1 = object_colindex - int(object_colindex)
            dx2 = 1.0 - dx1
            dy1 = object_rowindex - int(object_rowindex)
            dy2 = 1.0 - dy1
            a11 = dx1*dy1
            a12 = dx1*dy2
            a21 = dx2*dy1
            a22 = dx2*dy2
            j, i = int(object_colindex), int(object_rowindex)
            pathloss_vector = (a22*pathloss_refdata[:, i, j]
                               + a12*pathloss_refdata[:, i+1, j]
                               + a21*pathloss_refdata[:, i, j+1]
                               + a11*pathloss_refdata[:, i+1, j+1])

        return wavelength, pathloss_vector, is_inside_slitlet


def do_correction(input_model, pathloss_model=None, inverse=False, source_type=None, correction_pars=None):
    """
    Short Summary
    -------------
    Execute all tasks for Path Loss Correction

    Parameters
    ----------
    input_model : data model object
        science data to be corrected

    pathloss_model : pathloss model object or None
        pathloss correction data

    inverse : boolean
        Invert the math operations used to apply the flat field.

    source_type : str or None
        Force processing using the specified source type.

    correction_pars : dict or None
        Correction parameters to use instead of recalculation.

    Returns
    -------
    output_model, corrections : jwst.datamodel.DataModel, jwst.datamodel.datamodel
        2-tuple of the corrected science data with pathloss extensions added, and a
        model of the correction arrays.

    """
    if not pathloss_model and not correction_pars:
        raise RuntimeError(
            'Neither a PathLossModel nor PathLossStep correction parameters given.'
            ' One needs to be specified.'
        )
    exp_type = input_model.meta.exposure.type
    log.info(f'Input exposure type is {exp_type}')
    output_model = input_model.copy()

    if exp_type == 'NRS_MSASPEC':
        corrections = do_correction_mos(output_model, pathloss_model,
                                        inverse, source_type, correction_pars)
    elif exp_type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
        corrections = do_correction_fixedslit(output_model, pathloss_model,
                                              inverse, source_type, correction_pars)
    elif exp_type == 'NRS_IFU':
        corrections = do_correction_ifu(output_model, pathloss_model,
                                        inverse, source_type, correction_pars)
    elif exp_type == 'NIS_SOSS':
        if correction_pars:
            log.warning('Use of correction_pars with NIS_SOSS is not implemented. Skipping')
            output_model.meta.cal_step.pathloss = 'SKIPPED'
            corrections = None
        elif inverse:
            log.warning('Use of inversion with NIS_SOSS is not implemented. Skipping')
            output_model.meta.cal_step.pathloss = 'SKIPPED'
            corrections = None
        elif source_type is not None:
            log.warning('Forcing of source type with NIS_SOSS is not implemented. Skipping')
            output_model.meta.cal_step.pathloss = 'SKIPPED'
            corrections = None
        elif inverse:
            log.warning('Use of inversion with NIS_SOSS is not implemented. Skipping')
            output_model.meta.cal_step.pathloss = 'SKIPPED'
            corrections = None
        elif source_type is not None:
            log.warning('Forcing of source type with NIS_SOSS is not implemented. Skipping')
            output_model.meta.cal_step.pathloss = 'SKIPPED'
            corrections = None
        else:
            corrections = do_correction_soss(output_model, pathloss_model)

    return output_model, corrections


def interpolate_onto_grid(wavelength_grid, wavelength_vector, pathloss_vector):
    """
    Get the value of pathloss by interpolating each non-NaN element of
    wavelength_grid into pathloss_vector using the index lookup of
    wavelength_vector.  Pixels with wavelengths outside the range of the
    reference file should have a correction of NaN.

    Parameters:
    -----------

    wavelength_grid: numpy ndarray (2-d)

    The grid of wavelengths for each science data pixel

    wavelength_vector: numpy ndarray (1-d)

    Vector of wavelengths

    pathloss_vector:  numpy ndarray (1-d)

    Corresponding vector of pathloss values

    Returns:
    --------

    grid of pathloss corrections for each non-Nan pixel

    """

    # Need to set the pathloss correction of pixels whose wavelength is outside
    # the wavelength range of the reference file to NaN.  This trick will acomplish
    # that while still allowing the use of array linear interpolation
    #
    # Pad out the wavelength and pathloss vectors by adding another element
    # at the beginning and end, and put NaN in the new first and last elements
    # of the extended pathloss vector
    extended_pathloss_vector = np.zeros(len(pathloss_vector) + 2)
    extended_pathloss_vector[1:-1] = pathloss_vector
    extended_pathloss_vector[0] = np.nan
    extended_pathloss_vector[-1] = np.nan
    extended_wavelength_vector = np.zeros(len(wavelength_vector) + 2)
    extended_wavelength_vector[1:-1] = wavelength_vector
    extended_wavelength_vector[0] = wavelength_vector[0] - 0.1
    extended_wavelength_vector[-1] = wavelength_vector[-1] + 0.1

    # Find the indices in the original wavelength array that correspond
    # to the lower of 2 adjacent wavelength values spanning the wavelength
    # values in the wavelength grid.  NaNs and values > max wavelength will
    # return an index to an element 1 past the array, values below min wavelength
    # will return 0
    upper_indices = np.searchsorted(wavelength_vector,
                                    wavelength_grid)

    # Move these indices so they correspond to the extended arrays
    lower_indices = upper_indices
    upper_indices = upper_indices + 1

    # Now we can just proceed without worrying about values outside the wavelength
    # array
    numerator = wavelength_grid - extended_wavelength_vector[lower_indices]
    denominator = (extended_wavelength_vector[upper_indices]
                   - extended_wavelength_vector[lower_indices])
    fraction = numerator / denominator

    pathloss_grid = wavelength_grid * 0.0

    pathloss_grid = (extended_pathloss_vector[lower_indices]
                     + fraction*(extended_pathloss_vector[upper_indices]
                                 - extended_pathloss_vector[lower_indices]))

    return pathloss_grid


def is_pointsource(srctype):
    """Returns True if srctype is POINT"""
    if srctype is None:
        return False
    elif srctype.upper() == 'POINT':
        return True
    else:
        return False


def do_correction_mos(data, pathloss, inverse=False, source_type=None, correction_pars=None):
    """Path loss correction for NIRSpec MOS

    Data is modified in-place.

    Parameters
    ----------
    data : jwst.datamodel.DataModel
        The NIRSpec MOS data to be corrected.

    pathloss : jwst.datamodel.DataModel or None
        The pathloss reference data.

    inverse : boolean
        Invert the math operations used to apply the flat field.

    source_type : str or None
        Force processing using the specified source type.

    correction_pars : jwst.datamodels.MultiSlitModel or None
        The precomputed pathloss to apply instead of recalculation.

    Returns
    -------
    corrections : jwst.datamodel.MultiSlitModel
        The pathloss corrections applied.
    """
    exp_type = data.meta.exposure.type

    # Loop over all MOS slitlets
    corrections = datamodels.MultiSlitModel()
    for slit_number, slit in enumerate(data.slits):
        log.info(f'Working on slit {slit_number}')

        if correction_pars:
            correction = correction_pars.slits[slit_number]
        else:
            correction = _corrections_for_mos(slit, pathloss, exp_type, source_type)
        corrections.slits.append(correction)

        # Apply the correction
        if not correction:
            log.warning(f'No correction provided for slit {slit_number}. Skipping')
            continue

        if not inverse:
            slit.data /= correction.data
        else:
            slit.data *= correction.data
        slit.err /= correction.data
        slit.var_poisson /= correction.data**2
        slit.var_rnoise /= correction.data**2
        if slit.var_flat is not None and np.size(slit.var_flat) > 0:
            slit.var_flat /= correction.data**2
        slit.pathloss_point = correction.pathloss_point
        slit.pathloss_uniform = correction.pathloss_uniform

    # Set step status to complete
    data.meta.cal_step.pathloss = 'COMPLETE'

    return corrections


def do_correction_fixedslit(data, pathloss, inverse=False, source_type=None, correction_pars=None):
    """Path loss correction for NIRSpec fixed-slit modes

    Data is modified in-place.

    Parameters
    ----------
    data : jwst.datamodel.DataModel
        The NIRSpec fixed-slit data to be corrected.

    pathloss : jwst.datamodel.DataModel
        The pathloss reference data.

    inverse : boolean
        Invert the math operations used to apply the flat field.

    source_type : str or None
        Force processing using the specified source type.

    correction_pars : jwst.datamodels.MultiSlitModel or None
        The precomputed pathloss to apply instead of recalculation.

    Returns
    -------
    corrections : jwst.datamodel.MultiSlitModel
        The pathloss corrections applied.
    """
    exp_type = data.meta.exposure.type

    # Loop over all slits contained in the input
    corrections = datamodels.MultiSlitModel()
    for slit_number, slit in enumerate(data.slits):
        log.info(f'Working on slit {slit.name}')


        if correction_pars:
            correction = correction_pars.slits[slit_number]
        else:
            correction = _corrections_for_fixedslit(slit, pathloss, exp_type, source_type)
        corrections.slits.append(correction)

        # Apply the correction
        if not correction:
            log.warning(f'No correction provided for slit {slit_number}. Skipping')
            continue

        if not inverse:
            slit.data /= correction.data
        else:
            slit.data *= correction.data
        slit.err /= correction.data
        slit.var_poisson /= correction.data**2
        slit.var_rnoise /= correction.data**2
        if slit.var_flat is not None and np.size(slit.var_flat) > 0:
            slit.var_flat /= correction.data**2
        slit.pathloss_point = correction.pathloss_point
        slit.pathloss_uniform = correction.pathloss_uniform

        slit.data /= correction.data
        slit.err /= correction.data
        slit.var_poisson /= correction.data**2
        slit.var_rnoise /= correction.data**2
        if slit.var_flat is not None and np.size(slit.var_flat) > 0:
            slit.var_flat /= correction.data**2
        slit.pathloss_point = correction.pathloss_point
        slit.pathloss_uniform = correction.pathloss_uniform

    # Set step status to complete
    data.meta.cal_step.pathloss = 'COMPLETE'

    return corrections


def do_correction_ifu(data, pathloss, inverse=False, source_type=None, correction_pars=None):
    """Path loss correction for NIRSpec IFU

    Data is modified in-place.

    Parameters
    ----------
    data : jwst.datamodel.DataModel
        The NIRSpec fixed-slit data to be corrected.

    pathloss : jwst.datamodel.DataModel
        The pathloss reference data.

    inverse : boolean
        Invert the math operations used to apply the flat field.

    source_type : str or None
        Force processing using the specified source type.

    correction_pars : jwst.datamodels.MultiSlitModel or None
        The precomputed pathloss to apply instead of recalculation.

    Returns
    -------
    corrections : jwst.datamodel.MultiSlitModel
        The pathloss corrections applied.
    """
    if correction_pars:
        correction = correction_pars
    else:
        correction = _corrections_for_ifu(data, pathloss, source_type)

    if not inverse:
        data.data /= correction.data
    else:
        data.data *= correction.data
    data.err /= correction.data
    data.var_poisson /= correction.data**2
    data.var_rnoise /= correction.data**2
    if data.var_flat is not None and np.size(data.var_flat) > 0:
        data.var_flat /= correction.data**2
    data.pathloss_point = correction.pathloss_point
    data.pathloss_uniform = correction.pathloss_uniform

    # This might be useful to other steps
    data.wavelength = correction.wavelength

    # Set the step status to complete
    data.meta.cal_step.pathloss = 'COMPLETE'

    return correction


def do_correction_soss(data, pathloss):
    """Path loss correction for NIRSpec SOSS

    NIRISS SOSS pathloss correction is basically a correction for the
    flux from the 2nd and 3rd order dispersion that falls outside the
    subarray aperture.  The correction depends on the pupil wheel position
    and column number (or wavelength).  The simple option is to do the
    correction by column number, then the only interpolation needed is a
    1-d interpolation into the pupil wheel position dimension.

    Data is modified in-place.

    Parameters
    ----------
    data : jwst.datamodel.DataModel
        The NIRSpec SOSS data to be corrected.

    pathloss : jwst.datamodel.DataModel
        The pathloss reference data.
    """
    # Omit correction if this is a TSO observation
    if data.meta.visit.tsovisit:
        log.warning("NIRISS SOSS TSO observations skip the pathloss step")
        data.meta.cal_step.pathloss = 'SKIPPED'
        return

    # Get the pupil wheel position
    pupil_wheel_position = data.meta.instrument.pupil_position
    if pupil_wheel_position is None:
        log.warning('Unable to get pupil wheel position from PWCPOS keyword '
                    f'for {data.meta.filename}')
        log.warning("Pathloss correction skipped")
        data.meta.cal_step.pathloss = 'SKIPPED'
        return

    # Get the aperture from the reference file that matches the subarray
    subarray = data.meta.subarray.name
    aperture = get_aperture_from_model(pathloss, subarray)
    if aperture is None:
        log.warning('Unable to get Aperture from reference file '
                    f'for subarray {subarray}')
        log.warning("Pathloss correction skipped")
        data.meta.cal_step.pathloss = 'SKIPPED'
        return

    else:
        log.info(f'Aperture {aperture.name} selected from reference file')

    # Set up pathloss correction array
    pathloss_array = aperture.pointsource_data[0]
    nrows, ncols = pathloss_array.shape
    _, data_ncols = data.data.shape
    correction = np.ones(data_ncols, dtype=np.float32)
    crpix1 = aperture.pointsource_wcs.crpix1
    crval1 = aperture.pointsource_wcs.crval1
    cdelt1 = aperture.pointsource_wcs.cdelt1
    pupil_wheel_index = crpix1 + (pupil_wheel_position - crval1) / cdelt1 - 1

    if pupil_wheel_index < 0 or pupil_wheel_index > (ncols - 2):
        log.warning("Pupil Wheel position outside reference file coverage")
        log.warning("Setting pathloss correction to 1.0")
    else:
        ix = int(pupil_wheel_index)
        dx = pupil_wheel_index - ix
        crpix2 = aperture.pointsource_wcs.crpix2
        crval2 = aperture.pointsource_wcs.crval2
        cdelt2 = aperture.pointsource_wcs.cdelt2
        for row in range(data_ncols):
            row_1indexed = row + 1
            refrow_index = math.floor(crpix2 + (row_1indexed - crval2) / cdelt2 - 0.5)
            if refrow_index < 0 or refrow_index > (nrows - 1):
                correction[row] = 1.0
            else:
                correction[row] = (1.0 - dx) * pathloss_array[refrow_index, ix] + \
                                  dx * pathloss_array[refrow_index, ix + 1]

    # Create and apply the 2D correction
    pathloss_2d = np.broadcast_to(correction, data.data.shape)
    data.data /= pathloss_2d
    data.err /= pathloss_2d
    data.var_poisson /= pathloss_2d**2
    data.var_rnoise /= pathloss_2d**2
    if data.var_flat is not None and np.size(data.var_flat) > 0:
        data.var_flat /= pathloss_2d**2
    data.pathloss_point = pathloss_2d

    # Set step status to complete
    data.meta.cal_step.pathloss = 'COMPLETE'


def _corrections_for_mos(slit, pathloss, exp_type, source_type=None):
    """Calculate the correction arrasy for MOS slit

    Parameters
    ----------
    slit : jwst.datamodels.SlitModel
        The slit being operated on.

    pathloss : jwst.datamodels.DataModel
        The pathloss reference data

    exp_type : str
        Exposure type

    source_type : str or None
        Force processing using the specified source type.

    Returns
    -------
    correction : jwst.datamodels.SlitModel
        The correction arrays
    """
    correction = None
    size = slit.data.size

    # Only work on slits with data.size > 0
    if size > 0:

        # Get centering
        xcenter, ycenter = get_center(exp_type, slit)
        # Calculate the 1-d wavelength and pathloss vectors
        # for the source position
        # Get the aperture from the reference file that matches the slit
        nshutters = util.get_num_msa_open_shutters(slit.shutter_state)
        aperture = get_aperture_from_model(pathloss, nshutters)
        if aperture is not None:
            (wavelength_pointsource,
             pathloss_pointsource_vector,
             is_inside_slitlet) = calculate_pathloss_vector(aperture.pointsource_data,
                                                            aperture.pointsource_wcs,
                                                            xcenter, ycenter)
            (wavelength_uniformsource,
             pathloss_uniform_vector,
             dummy) = calculate_pathloss_vector(aperture.uniform_data,
                                                aperture.uniform_wcs,
                                                xcenter, ycenter)
            if is_inside_slitlet:

                # Wavelengths in the reference file are in meters,
                # need them to be in microns
                wavelength_pointsource *= 1.0e6
                wavelength_uniformsource *= 1.0e6

                wavelength_array = slit.wavelength

                # Compute the point source pathloss 2D correction
                pathloss_2d_ps = interpolate_onto_grid(
                    wavelength_array,
                    wavelength_pointsource,
                    pathloss_pointsource_vector)

                # Compute the uniform source pathloss 2D correction
                pathloss_2d_un = interpolate_onto_grid(
                    wavelength_array,
                    wavelength_uniformsource,
                    pathloss_uniform_vector)

                # Use the appropriate correction for this slit
                if is_pointsource(source_type or slit.source_type):
                    pathloss_2d = pathloss_2d_ps
                else:
                    pathloss_2d = pathloss_2d_un

                # Save the corrections. The `data` portion is the correction used.
                # The individual ones will be saved in the respective attributes.
                correction = datamodels.SlitModel(data=pathloss_2d)
                correction.pathloss_point = pathloss_2d_ps
                correction.pathloss_uniform = pathloss_2d_un
            else:
                log.warning("Source is outside slit.")
        else:
            log.warning("Cannot find matching pathloss model for slit with"
                        f"{nshutters} shutters")
    else:
        log.warning(f"Slit has data size = {size}")

    return correction


def _corrections_for_fixedslit(slit, pathloss, exp_type, source_type):
    """Calculate the correction arrasy for Fixed-slit slit

    Parameters
    ----------
    slit : jwst.datamodels.SlitModel
        The slit being operated on.

    pathloss : jwst.datamodels.DataModel
        The pathloss reference data

    exp_type : str
        Exposure type

    source_type : str or None
        Force processing using the specified source type.

    Returns
    -------
    correction : jwst.datamodels.SlitModel
        The correction arrays
    """
    correction = None

    # Get centering
    xcenter, ycenter = get_center(exp_type, slit)
    # Calculate the 1-d wavelength and pathloss vectors for the source position
    # Get the aperture from the reference file that matches the slit
    aperture = get_aperture_from_model(pathloss, slit.name)

    if aperture is not None:
        log.info(f'Using aperture {aperture.name}')
        (wavelength_pointsource,
         pathloss_pointsource_vector,
         is_inside_slit) = calculate_pathloss_vector(aperture.pointsource_data,
                                                     aperture.pointsource_wcs,
                                                     xcenter, ycenter)
        (wavelength_uniformsource,
         pathloss_uniform_vector,
         dummy) = calculate_pathloss_vector(aperture.uniform_data,
                                            aperture.uniform_wcs,
                                            xcenter, ycenter)
        if is_inside_slit:

            # Wavelengths in the reference file are in meters,
            # need them to be in microns
            wavelength_pointsource *= 1.0e6
            wavelength_uniformsource *= 1.0e6

            wavelength_array = slit.wavelength

            # Compute the point source pathloss 2D correction
            pathloss_2d_ps = interpolate_onto_grid(
                wavelength_array,
                wavelength_pointsource,
                pathloss_pointsource_vector)

            # Compute the uniform source pathloss 2D correction
            pathloss_2d_un = interpolate_onto_grid(
                wavelength_array,
                wavelength_uniformsource,
                pathloss_uniform_vector)

            # Use the appropriate correction for this slit
            if is_pointsource(source_type or slit.source_type):
                pathloss_2d = pathloss_2d_ps
            else:
                pathloss_2d = pathloss_2d_un

            # Save the corrections. The `data` portion is the correction used.
            # The individual ones will be saved in the respective attributes.
            correction = datamodels.SlitModel(data=pathloss_2d)
            correction.pathloss_point = pathloss_2d_ps
            correction.pathloss_uniform = pathloss_2d_un

        else:
            log.warning('Source is outside slit. Skipping '
                        f'pathloss correction for slit {slit.name}')
    else:
        log.warning(f'Cannot find matching pathloss model for {slit.name}')
        log.warning('Skipping pathloss correction for this slit')

    return correction


def _corrections_for_ifu(data, pathloss, source_type):
    """Calculate the correction arrasy for MOS slit

    Parameters
    ----------
    data : jwst.datamodels.SlitModel
        The data being operated on.

    pathloss : jwst.datamodels.DataModel
        The pathloss reference data

    source_type : str or None
        Force processing using the specified source type.

    Returns
    -------
    correction : jwst.datamodels.SlitModel
        The correction arrays
    """

    # IFU targets are always inside slit
    # Get centering
    xcenter, ycenter = get_center(data.meta.exposure.type, None)

    # Calculate the 1-d wavelength and pathloss vectors for the source position
    aperture = pathloss.apertures[0]
    (wavelength_pointsource,
     pathloss_pointsource_vector,
     dummy) = calculate_pathloss_vector(aperture.pointsource_data,
                                        aperture.pointsource_wcs,
                                        xcenter, ycenter)
    (wavelength_uniformsource,
     pathloss_uniform_vector,
     dummy) = calculate_pathloss_vector(aperture.uniform_data,
                                        aperture.uniform_wcs,
                                        xcenter, ycenter)
    # Wavelengths in the reference file are in meters;
    # need them to be in microns
    wavelength_pointsource *= 1.0e6
    wavelength_uniformsource *= 1.0e6

    # Create the 2-d wavelength arrays, initialize with NaNs
    wavelength_array = np.zeros(data.shape, dtype=np.float32)
    wavelength_array.fill(np.nan)
    for slice in NIRSPEC_IFU_SLICES:
        slice_wcs = nirspec.nrs_wcs_set_input(data, slice)
        x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
        ra, dec, wavelength = slice_wcs(x, y)
        valid = ~np.isnan(wavelength)
        x = x[valid]
        y = y[valid]
        wavelength_array[y.astype(int), x.astype(int)] = wavelength[valid]

    # Compute the point source pathloss 2D correction
    pathloss_2d_ps = interpolate_onto_grid(
        wavelength_array,
        wavelength_pointsource,
        pathloss_pointsource_vector)

    # Compute the uniform source pathloss 2D correction
    pathloss_2d_un = interpolate_onto_grid(
        wavelength_array,
        wavelength_uniformsource,
        pathloss_uniform_vector)

    # Use the appropriate correction for the source type
    if is_pointsource(source_type or data.meta.target.source_type):
        pathloss_2d = pathloss_2d_ps
    else:
        pathloss_2d = pathloss_2d_un

    # Save the corrections. The `data` portion is the correction used.
    # The individual ones will be saved in the respective attributes.
    correction = type(data)(data=pathloss_2d)
    correction.pathloss_point = pathloss_2d_ps
    correction.pathloss_uniform = pathloss_2d_un
    correction.wavelength = wavelength_array

    return correction
