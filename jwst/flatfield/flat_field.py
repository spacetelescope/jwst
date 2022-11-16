#
#  Module for applying flat fielding
#

import logging
import math

import numpy as np
from gwcs.wcstools import grid_from_bounding_box

from .. import datamodels
from ..datamodels import dqflags
from ..lib import reffile_utils
from ..assign_wcs import nirspec

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

MICRONS_100 = 1.e-4  # 100 microns, in meters

# This is for NIRSpec.
FIXED_SLIT_TYPES = ["NRS_LAMP", "NRS_AUTOWAVE", "NRS_BRIGHTOBJ", "NRS_FIXEDSLIT"]
NIRSPEC_SPECTRAL_EXPOSURES = ['NRS_AUTOWAVE', 'NRS_BRIGHTOBJ', 'NRS_FIXEDSLIT', 'NRS_IFU', 'NRS_LAMP', 'NRS_MSASPEC']

# Dispersion direction, predominantly horizontal or vertical.  These values
# are to be compared with keyword DISPAXIS from the input header.
HORIZONTAL = 1
VERTICAL = 2


def do_correction(input_model,
                  flat=None, fflat=None, sflat=None, dflat=None, user_supplied_flat=None,
                  inverse=False):
    """Flat-field a JWST data model using a flat-field model

    Parameters
    ----------
    input_model : JWST data model
        Input science data model to be flat-fielded.

    flat : JWST data model, or None
        Data model containing flat-field for all instruments other than
        NIRSpec spectrographic data.

    fflat : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.  Used only for NIRSpec data.

    sflat : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.  Used only for NIRSpec data.

    dflat : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.  Used only for NIRSpec data.

    user_supplied_flat : ~jwst.datamodels.DataModel
        If supplied, all other reference flats and flat creation are
        ignored in favor of the specified flat.

    inverse : boolean
        Invert the math operations used to apply the flat field.

    Returns
    -------
    output_model : data model
        The data model for the flat-fielded science data.

    flat_applied : ~jwst.datamodels.MultiSlitModel, ~jwst.datamodels.ImageModel
        Data model containing the interpolated flat fields (NIRSpec data only), or
        just the input flat.
    """

    # Initialize the output model as a copy of the input
    output_model = input_model.copy()

    # NIRSpec spectrographic data are processed differently from other
    # types of data (including NIRSpec imaging).  The test on flat is
    # needed because NIRSpec imaging data are processed by do_flat_field().
    if ((input_model.meta.exposure.type in NIRSPEC_SPECTRAL_EXPOSURES) and
            (input_model.meta.instrument.lamp_mode != 'IMAGE')):
        flat_applied = do_nirspec_flat_field(output_model, fflat, sflat, dflat,
                                             user_supplied_flat=user_supplied_flat,
                                             inverse=inverse)
    else:
        if user_supplied_flat is not None:
            flat = user_supplied_flat
        if flat is None:
            log.warning("No flat found or supplied; step will be skipped.")
            output_model.meta.cal_step.flat_field = 'SKIPPED'
            flat_applied = None
        else:
            do_flat_field(output_model, flat, inverse=inverse)
            flat_applied = flat

    return output_model, flat_applied


#
# These functions are for non-NIRSpec flat fielding, or for NIRSpec imaging.
#

def do_flat_field(output_model, flat_model, inverse=False):
    """Apply flat-fielding for non-NIRSpec modes, updating the output model.

    Parameters
    ----------
    output_model : JWST data model
        flat-fielded input science data model, modified in-place

    flat_model : JWST data model
        data model containing flat-field

    inverse : boolean
        Invert the math operations used to apply the flat field.
    """
    if output_model.meta.instrument.name == "NIRSPEC":
        log.debug("Flat field correction for NIRSpec imaging data.")
    else:
        log.debug("Flat field correction for non-NIRSpec modes.")

    any_updated = False  # will set True if any flats applied

    # Check to see if flat data array is smaller than science data
    if ((output_model.data.shape[-1] > flat_model.data.shape[-1]) or
            (output_model.data.shape[-2] > flat_model.data.shape[-2])):
        log.warning('Reference data array is smaller than science data')
        log.warning('Step will be skipped')

    # Apply flat to MultiSlits
    elif isinstance(output_model, datamodels.MultiSlitModel):

        # Apply flat to each slit contained in the input
        for slit in output_model.slits:
            log.debug('Applying flat to slit %s' % (slit.name))
            apply_flat_field(slit, flat_model, inverse=inverse)
            any_updated = True

    # Apply flat to all other models
    else:
        apply_flat_field(output_model, flat_model, inverse=inverse)
        any_updated = True

    if any_updated:
        output_model.meta.cal_step.flat_field = 'COMPLETE'
    else:
        output_model.meta.cal_step.flat_field = 'SKIPPED'


def apply_flat_field(science, flat, inverse=False):
    """Flat field the data and error arrays.

    Extended summary
    ----------------
    The science data and error arrays will be divided by the flat field.
    The data quality array will be updated based on bad pixels in flat
    field arrays. Applies portion of flat field corresponding to science
    image subarray.

    Parameters
    ----------
    science : JWST data model
        input science data model

    flat : JWST data model
        flat field data model

    inverse : boolean
        Invert the math operations used to apply the flat field.
    """

    # Extract subarray from reference data, if necessary
    if reffile_utils.ref_matches_sci(science, flat):
        flat_data = flat.data
        flat_dq = flat.dq
        flat_err = flat.err
    else:
        log.info("Extracting matching subarray from flat")
        sub_flat = reffile_utils.get_subarray_model(science, flat)
        flat_data = sub_flat.data.copy()
        flat_dq = sub_flat.dq.copy()
        flat_err = sub_flat.err.copy()
        sub_flat.close()

    # Find pixels in the flat that have a value of NaN and set
    # DQ = DO_NOT_USE + NO_FLAT_FIELD
    bad_flag = dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['NO_FLAT_FIELD']

    flat_nan = np.isnan(flat_data)
    flat_dq[flat_nan] = np.bitwise_or(flat_dq[flat_nan], bad_flag)

    # Find pixels in the flat that have a value of zero, and set
    # DQ = DO_NOT_USE + NO_FLAT_FIELD
    flat_zero = np.where(flat_data == 0.)
    flat_dq[flat_zero] = np.bitwise_or(flat_dq[flat_zero], bad_flag)

    # Find pixels in the flat that have the DQ bit for NO_FLAT_FIELD set
    # Set the DO_NOT_USE flag for such pixels
    flat_noflat = np.where(np.bitwise_and(flat_dq, dqflags.pixel['NO_FLAT_FIELD']))
    flat_dq[flat_noflat] = np.bitwise_or(flat_dq[flat_noflat], dqflags.pixel['DO_NOT_USE'])

    # Find all pixels in the flat that have a DQ value of DO_NOT_USE
    flat_bad = np.bitwise_and(flat_dq, dqflags.pixel['DO_NOT_USE'])

    # Reset the flat value of all bad pixels to 1.0, so that no
    # correction is made
    flat_data[np.where(flat_bad)] = 1.0

    # Now let's apply the correction to science data and error arrays.  Rely
    # on array broadcasting to handle the cubes
    if not inverse:
        science.data /= flat_data
    else:
        science.data *= flat_data

    # Update the variances using BASELINE algorithm.  For guider data, it has
    # not gone through ramp fitting so there is no Poisson noise or readnoise
    if not isinstance(science, datamodels.GuiderCalModel):
        flat_data_squared = flat_data ** 2
        science.var_poisson /= flat_data_squared
        science.var_rnoise /= flat_data_squared
        science.var_flat = science.data ** 2 / flat_data_squared * flat_err ** 2
        science.err = np.sqrt(science.var_poisson + science.var_rnoise + science.var_flat)
    else:
        flat_data_squared = flat_data ** 2
        science.var_flat = science.data ** 2 / flat_data_squared * flat_err ** 2

        # Set the output ERR to be the combined input ERR plus flatfield ERR, summed in quadrature
        science.err = np.sqrt(science.err**2 + science.var_flat)

    # Combine the science and flat DQ arrays
    science.dq = np.bitwise_or(science.dq, flat_dq)

    # Find all pixels in the flat that have a DQ value of NON_SCIENCE
    # add the DO_NOT_USE flag to these pixels. We don't want to use these pixels
    # in further steps
    flag_nonsci = np.bitwise_and(science.dq, dqflags.pixel['NON_SCIENCE']).astype(bool)
    science.dq[flag_nonsci] = np.bitwise_or(science.dq[flag_nonsci], dqflags.pixel['DO_NOT_USE'])


#
# The following functions are for NIRSpec spectrographic data.
#


def do_nirspec_flat_field(output_model, f_flat_model, s_flat_model, d_flat_model,
                          user_supplied_flat=None, inverse=False):
    """Apply flat-fielding for NIRSpec data, updating in-place.

    Calls one of 3 functions depending on whether the data is 1) NIRSpec IFU,
    2) NIRSpec BRIGHTOBJ or 3) NIRSpec MSA or Fixed-slit.

    Parameters
    ----------
    output_model : JWST data model
        Science data model, modified (flat fielded) in-place.

    f_flat_model : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    user_supplied_flat : ~jwst.datamodels.DataModel or None
        If provided, override all other calculated or reference-file-retrieved
        flat information and use this data.

    inverse : boolean
        Invert the math operations used to apply the flat field.

    Returns
    -------
    ~jwst.datamodels.MultiSlitModel or ~jwst.datamodels.ImageModel
        The interpolated flat field(s).
    """

    log.debug("Flat field correction for NIRSpec spectrographic data.")

    exposure_type = output_model.meta.exposure.type
    try:
        dispaxis = output_model.meta.wcsinfo.dispersion_direction
    except AttributeError:
        if len(output_model.slits) > 0:
            dispaxis = output_model.slits[0].meta.wcsinfo.dispersion_direction
        else:
            dispaxis = None
    if dispaxis is None:
        log.warning("Can't determine dispaxis, assuming horizontal.")
        dispaxis = HORIZONTAL

    if exposure_type == "NRS_BRIGHTOBJ":
        if not isinstance(output_model, datamodels.SlitModel):
            log.error("NIRSpec BRIGHTOBJ data is not a SlitModel; "
                      "don't know how to process it.")
            raise RuntimeError("Input is {}; expected SlitModel"
                               .format(type(output_model)))
        return nirspec_brightobj(output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis,
                                 user_supplied_flat=user_supplied_flat, inverse=inverse)

    # We expect NIRSpec IFU data to be an IFUImageModel, but it's conceivable
    # that the slices have been copied out into a MultiSlitModel, so
    # check for that case.
    if not hasattr(output_model, "slits"):
        if exposure_type == "NRS_IFU" or (
                exposure_type in ["NRS_AUTOWAVE", "NRS_LAMP"] and output_model.meta.instrument.lamp_mode == 'IFU'):
            if not isinstance(output_model, datamodels.IFUImageModel):
                log.error("NIRSpec IFU data is not an IFUImageModel; "
                          "don't know how to process it.")
                raise RuntimeError("Input is {}; expected IFUImageModel"
                                   .format(type(output_model)))
            return nirspec_ifu(output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis,
                               user_supplied_flat=user_supplied_flat, inverse=inverse)
        else:
            raise RuntimeError(f'No flat field algorithm exists for handling data {output_model}')

    # For datamodels with slits, MSA and Fixed slit modes:
    else:
        return nirspec_fs_msa(output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis,
                              user_supplied_flat=user_supplied_flat, inverse=inverse)


def nirspec_fs_msa(output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis,
                   user_supplied_flat=None, inverse=False):
    """Apply flat-fielding for NIRSpec fixed slit and MSA data, in-place

    Parameters
    ----------
    output_model : `~jwst.datamodels.MultiSlitModel`
        MultiSlitModel, modified (flat fielded) slit-by-slit, in-place.

    f_flat_model : `~jwst.datamodels.NirspecFlatModel` or None
        Flat field for the fore optics.

    s_flat_model : `~jwst.datamodels.NirspecFlatModel` or None
        Flat field for the spectrograph.

    d_flat_model : `~jwst.datamodels.NirspecFlatModel` or None
        Flat field for the detector.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    user_supplied_flat : ~jwst.datamodels.DataModel or None
        If provided, override all other calculated or reference-file-retrieved
        flat information and use this data.

    inverse : boolean
        Invert the math operations used to apply the flat field.

    Returns
    -------
    interpolated_flats: `~jwst.datamodels.MultiSlitModel`
        The interpolated flat field, one for each slit.
    """

    exposure_type = output_model.meta.exposure.type
    primary_slit = output_model.meta.instrument.fixed_slit

    # Create a list to hold the list of slits.  This will eventually be used
    # to extend the MultiSlitModel.slits attribute.  We do it this way to
    # postpone validation until the end, which is faster.
    flat_slits = []

    # A flag to make sure at least one slit was flatfielded, so we can set
    # "COMPLETE", otherwise we set "SKIP"
    any_updated = False

    for slit_idx, slit in enumerate(output_model.slits):
        log.info("Working on slit %s", slit.name)
        if exposure_type == "NRS_MSASPEC":
            slit_nt = slit  # includes quadrant info
        else:
            slit_nt = None

        if user_supplied_flat is not None:
            slit_flat = user_supplied_flat.slits[slit_idx]
        else:
            if exposure_type == "NRS_FIXEDSLIT" and \
                    slit.name == primary_slit and slit.source_type.upper() == "POINT":

                # For fixed-slit exposures, if this is the primary slit
                # and it contains a point source, compute the flat-field
                # corrections for both uniform (without wavecorr) and point
                # source (with wavecorr) modes, applying only the point
                # source version to the data.

                # First compute a flat appropriate for a uniform source,
                # which means NOT using corrected wavelengths
                slit_flat = flat_for_nirspec_slit(
                    slit, f_flat_model, s_flat_model, d_flat_model,
                    dispaxis, exposure_type, slit_nt, output_model.meta.subarray,
                    use_wavecorr=False
                )

                # Store the result for uniform source
                slit.flatfield_uniform = slit_flat.data

                # Now compute a flat appropriate for a point source,
                # which means using corrected wavelengths
                slit_flat = flat_for_nirspec_slit(
                    slit, f_flat_model, s_flat_model, d_flat_model,
                    dispaxis, exposure_type, slit_nt, output_model.meta.subarray,
                    use_wavecorr=True
                )

                # Store the result for point source; this will be
                # the version actually applied to the data below.
                slit.flatfield_point = slit_flat.data

            else:
                # Build the flat for this slit the normal way, without any
                # specification for whether we want to use corrected wavelengths
                slit_flat = flat_for_nirspec_slit(
                    slit, f_flat_model, s_flat_model, d_flat_model,
                    dispaxis, exposure_type, slit_nt, output_model.meta.subarray,
                    use_wavecorr=None
                )
                if slit_flat is None:
                    log.debug(f'Slit {slit} flat field could not be determined.')
                    continue

            # Append the SlitDataModel to the list of slits
            flat_slits.append(slit_flat)

        # Now let's apply the correction to science data and error arrays.  Rely
        # on array broadcasting to handle the cubes
        if not inverse:
            slit.data /= slit_flat.data
        else:
            slit.data *= slit_flat.data

        # Update the variances using BASELINE algorithm
        flat_data_squared = slit_flat.data ** 2
        slit.var_poisson /= flat_data_squared
        slit.var_rnoise /= flat_data_squared
        slit.var_flat = slit.data ** 2 / flat_data_squared * slit_flat.err ** 2
        slit.err = np.sqrt(slit.var_poisson + slit.var_rnoise + slit.var_flat)

        # Combine the science and flat DQ arrays
        slit.dq |= slit_flat.dq

        any_updated = True

    if any_updated:
        output_model.meta.cal_step.flat_field = 'COMPLETE'
    else:
        output_model.meta.cal_step.flat_field = 'SKIPPED'

    # Create an output model for the interpolated flat fields.
    if user_supplied_flat:
        interpolated_flat = user_supplied_flat
    else:
        interpolated_flat = datamodels.MultiSlitModel()
        interpolated_flat.update(output_model, only="PRIMARY")
        interpolated_flat.slits.extend(flat_slits)

    return interpolated_flat


def nirspec_brightobj(output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis,
                      user_supplied_flat=None, inverse=False):
    """Apply flat-fielding for NIRSpec BRIGHTOBJ data, in-place

    Parameters
    ----------
    output_model : JWST data model
        CubeModel, modified (flat fielded) plane by plane, in-place.

    f_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the fore optics.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    user_supplied_flat : ~jwst.datamodels.ImageModel or None
        A pre-computed flat to use directly. If supplied,
        all other inputs are ignored

    inverse : boolean
        Invert the math operations used to apply the flat field.

    Returns
    -------
    ~jwst.datamodels.ImageModel
        The interpolated flat field.
    """

    if user_supplied_flat is not None:
        log.info(f'Pre-computed flat {user_supplied_flat} provided. Using the flat directly')
        interpolated_flat = user_supplied_flat
    else:
        interpolated_flat = flat_for_nirspec_brightobj(
            output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis
        )

    if not inverse:
        output_model.data /= interpolated_flat.data
    else:
        output_model.data *= interpolated_flat.data
    output_model.dq |= interpolated_flat.dq

    # Update the variances and uncertainty array using BASELINE algorithm
    flat_data_squared = interpolated_flat.data ** 2
    output_model.var_poisson /= flat_data_squared
    output_model.var_rnoise /= flat_data_squared
    output_model.var_flat = output_model.data ** 2 / flat_data_squared * interpolated_flat.err ** 2
    output_model.err = np.sqrt(
        output_model.var_poisson + output_model.var_rnoise + output_model.var_flat
    )

    output_model.meta.cal_step.flat_field = 'COMPLETE'

    return interpolated_flat


def nirspec_ifu(output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis,
                user_supplied_flat=None, inverse=False):
    """Apply flat-fielding for NIRSpec IFU data, in-place

    Parameters
    ----------
    output_model : JWST data model
        Science data model, modified (flat fielded) in-place.

    f_flat_model : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    user_supplied_flat : ~jwst.datamodels.ImageModel or None
        A pre-computed flat to use directly. If supplied,
        all other inputs are ignored

    inverse : boolean
        Invert the math operations used to apply the flat field.

    Returns
    -------
    ~jwst.datamodels.ImageModel
        The interpolated flat field.
    """

    if user_supplied_flat is not None:
        log.info(f'Pre-computed flat {user_supplied_flat} provided. Using the flat directly')
        flat = user_supplied_flat.data
        flat_dq = user_supplied_flat.dq
        flat_err = user_supplied_flat.err
        any_updated = True
    else:
        flat, flat_dq, flat_err, any_updated = flat_for_nirspec_ifu(
            output_model, f_flat_model, s_flat_model, d_flat_model, dispaxis
        )

    if any_updated:
        if not inverse:
            output_model.data /= flat
        else:
            output_model.data *= flat
        output_model.dq |= flat_dq

        # Update the variances and uncertainty array using BASELINE algorithm
        flat_data_squared = flat ** 2
        output_model.var_poisson /= flat_data_squared
        output_model.var_rnoise /= flat_data_squared
        output_model.var_flat = output_model.data ** 2 / flat_data_squared * flat_err ** 2
        output_model.err = np.sqrt(
            output_model.var_poisson + output_model.var_rnoise + output_model.var_flat
        )

        output_model.meta.cal_step.flat_field = 'COMPLETE'

        # Create an output model for the interpolated flat fields.
        interpolated_flats = datamodels.ImageModel(data=flat, dq=flat_dq, err=flat_err)
        interpolated_flats.update(output_model, only="PRIMARY")
    else:
        output_model.meta.cal_step.flat_field = 'SKIPPED'
        interpolated_flats = None

    return interpolated_flats


def create_flat_field(wl, f_flat_model, s_flat_model, d_flat_model,
                      xstart, xstop, ystart, ystop,
                      exposure_type, dispaxis, slit_name, slit_nt=None):
    """Extract and combine flat field components for NIRSpec

    Parameters
    ----------
    wl : 2-D ndarray
        Wavelength at each pixel of the 2-D slit array.  This array has
        shape (ystop - ystart, xstop - xstart).

    f_flat_model : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    xstart, ystart : int
        Starting pixel numbers (zero indexed) for the slice containing
        the data for the current slit.

    xstop, ystop : int
        End of the slice containing the data for the current slit.  The
        start and stop values are Python slice notation, i.e. the region
        to be extracted is [ystart:ystop, xstart:xstop].

    exposure_type : str
        The exposure type refers to fixed_slit, IFU, or using the
        micro-shutter array.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    slit_name : str
        The name of the slit currently being processed.

    slit_nt : namedtuple or None
        For MSA data only, info about the current slit.

    Returns
    -------
    flat_2d : ndarray, 2-D, float
        The flat field, interpolated over wavelength, same shape as `wl`.
        Divide the 2-D extracted spectrum by this array to correct for
        flat-field variations.

    flat_dq : ndarray, 2-D, uint32
        The data quality array corresponding to flat_2d.

    flat_err : ndarray, 2-D, float
        The error array corresponding to flat_2d.
    """

    f_flat, f_flat_dq, f_flat_err = fore_optics_flat(
        wl, f_flat_model, exposure_type, dispaxis,
        slit_name, slit_nt)

    s_flat, s_flat_dq, s_flat_err = spectrograph_flat(
        wl, s_flat_model, xstart, xstop, ystart, ystop,
        exposure_type, dispaxis, slit_name)
    d_flat, d_flat_dq, d_flat_err = detector_flat(
        wl, d_flat_model, xstart, xstop, ystart, ystop,
        exposure_type, dispaxis, slit_name)

    flat_2d = f_flat * s_flat * d_flat
    flat_dq = combine_dq(f_flat_dq, s_flat_dq, d_flat_dq,
                         default_shape=flat_2d.shape)

    # Combine the uncertainty arrays, excluding the ones that are None
    sum_var = np.zeros_like(flat_2d)
    if f_flat_err is not None:
        np.place(f_flat, f_flat == 0, 1.)
        sum_var += f_flat_err ** 2 / f_flat ** 2
    if s_flat_err is not None:
        np.place(s_flat, s_flat == 0, 1.)
        sum_var += s_flat_err ** 2 / s_flat ** 2
    if d_flat_err is not None:
        np.place(d_flat, d_flat == 0, 1.)
        sum_var += d_flat_err ** 2 / d_flat ** 2
    flat_err = flat_2d * np.sqrt(sum_var)

    mask = np.bitwise_and(flat_dq, dqflags.pixel['DO_NOT_USE'])
    flat_2d[np.where(mask)] = 1.

    return flat_2d, flat_dq, flat_err


def fore_optics_flat(wl, f_flat_model, exposure_type, dispaxis,
                     slit_name, slit_nt):
    """Extract the flat for the fore optics part.

    Parameters
    ----------
    wl : 2-D ndarray
        Wavelength at each pixel of the 2-D slit array.

    f_flat_model : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.

    exposure_type : str
        The exposure type refers to fixed_slit, IFU, or MSA.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    slit_name : str
        The name of the slit currently being processed.

    slit_nt : namedtuple or None
        For MOS data (only), this is used to get the quadrant number and
        the indices of the current shutter in the Y and X directions.

    Returns
    -------
    f_flat : ndarray, 2d, float32
        The computed flat field for the fore optics.

    f_flat_dq : ndarray, 2d, uint32, or None
        The associated data quality array.

    f_flat_err : ndarray, 2d, float32, or None
        The associated error array.
    """

    if f_flat_model is None:
        f_flat = np.ones(wl.shape, dtype=np.float32)
        f_flat_dq = None
        f_flat_err = None
        return f_flat, f_flat_dq, f_flat_err

    if slit_nt is None:
        quadrant = None
    else:
        quadrant = slit_nt.quadrant - 1  # convert to zero indexed

    (tab_wl, tab_flat) = read_flat_table(f_flat_model, exposure_type,
                                         slit_name, quadrant)
    if tab_wl.max() < MICRONS_100:
        log.warning("Wavelengths in f_flat table appear to be in meters")

    # While there actually is a slowly varying flat field for the MSA mode,
    # it's a 1-D array, not 2-D.  This array will be applied by incorporating
    # it into tab_flat.  So even for the MSA, there will not be any 2-D
    # image flat, so set the variable to 1.
    flat_2d = 1.

    f_flat_dq = None
    f_flat_err = None

    if exposure_type == "NRS_MSASPEC":
        # The MOS "image" is in MSA coordinates (shutter index in x and y),
        # not detector pixel coordinates.
        # This is an example to show what xcen and ycen mean:
        #       shutter_id = xcen + (ycen - 1) * 365
        msa_y, msa_x = slit_nt.ycen, slit_nt.xcen
        msa_y -= 1  # convert to zero indexed
        msa_x -= 1
        full_array_flat = f_flat_model.quadrants[quadrant].data
        full_array_err = f_flat_model.quadrants[quadrant].err

        # Get the wavelength corresponding to each plane in the "image".
        image_wl = read_image_wl(f_flat_model, quadrant)
        if image_wl.max() < MICRONS_100:
            log.warning("Wavelengths in f_flat image appear to be in meters.")
        one_d_flat = full_array_flat[:, msa_y, msa_x]
        one_d_err = full_array_err[:, msa_y, msa_x]

        # The wavelengths and flat-field values read from the reference
        # table are tab_wl and tab_flat respectively.  We need to combine
        # the 1-D MSA flat from the reference "image" with the table
        # values, so interpolate the 1-D MSA flat at each wavelength in
        # tab_wl, and multiply into tab_flat.
        #     numpy.interp(x, xp, fp, left, right)
        #       x:  values at which to interpolate
        #       xp:  array of independent-variable values (must be increasing)
        #       fp:  array of data values
        #       left, right:  values to return for out-of-bounds x
        tab_flat *= np.interp(tab_wl, image_wl, one_d_flat, 1., 1.)
        f_flat_err = np.interp(wl, image_wl, one_d_err, 0., 0.)

    # The shape of the output array is obtained from `wl`.

    f_flat, f_flat_dq = combine_fast_slow(wl, flat_2d, f_flat_dq,
                                          tab_wl, tab_flat, dispaxis)

    # Find pixels in the flat that have a value of NaN and add to
    # DQ mask, DO_NOT_USE + NO_FLAT_FIELD
    bad_flag = dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['NO_FLAT_FIELD']
    flat_nan = np.isnan(f_flat)
    f_flat_dq[flat_nan] = np.bitwise_or(f_flat_dq[flat_nan], bad_flag)

    # Find pixels in the flat have have a value of zero, and add to
    # DQ mask,  DO_NOT_USE + NO_FLAT_FIELD
    flat_zero = np.where(f_flat == 0.)
    f_flat_dq[flat_zero] = np.bitwise_or(f_flat_dq[flat_zero], bad_flag)

    # Find all pixels in the flat that have a DQ value of DO_NOT_USE
    flat_bad = np.bitwise_and(f_flat_dq, dqflags.pixel['DO_NOT_USE'])

    # Reset the flat value of all bad pixels to 1.0, so that no
    # correction is made
    f_flat[np.where(flat_bad)] = 1.0
    return f_flat, f_flat_dq, f_flat_err


def spectrograph_flat(wl, s_flat_model,
                      xstart, xstop, ystart, ystop,
                      exposure_type, dispaxis, slit_name):
    """Extract the flat for the spectrograph part.

    Parameters
    ----------
    wl : 2-D ndarray
        Wavelength at each pixel of the 2-D slit array.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    xstart, ystart : int
        Starting pixel numbers (zero indexed) for the slice containing
        the data for the current slit.

    xstop, ystop : int
        End of the slice containing the data for the current slit.  The
        start and stop values are Python slice notation, i.e. the region
        to be extracted is [ystart:ystop, xstart:xstop].

    exposure_type : str
        The exposure type refers to fixed_slit, IFU, or using the
        micro-shutter array.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    slit_name : str
        The name of the slit currently being processed.

    Returns
    -------
    s_flat : ndarray, 2d, float32
        The computed flat field for the spectrograph.

    s_flat_dq : ndarray, 2d, uint32, or None
        The associated data quality array.

    s_flat_err : ndarray, 2d, float32, or None
        The associated error array.
    """

    if s_flat_model is None:
        s_flat = np.ones(wl.shape, dtype=np.float32)
        s_flat_dq = None
        s_flat_err = None
        return s_flat, s_flat_dq, s_flat_err

    quadrant = None

    if xstart >= xstop or ystart >= ystop:
        return 1., None

    (tab_wl, tab_flat) = read_flat_table(s_flat_model, exposure_type,
                                         slit_name, quadrant)
    if tab_wl.max() < MICRONS_100:
        log.warning("Wavelengths in s_flat table appear to be in meters")

    full_array_flat = s_flat_model.data

    full_array_dq = s_flat_model.dq
    full_array_err = s_flat_model.err

    if len(full_array_flat.shape) == 3:
        # MSA data
        image_flat = full_array_flat[:, ystart:ystop, xstart:xstop]
        image_dq = full_array_dq[:, ystart:ystop, xstart:xstop]
        image_err = full_array_err[:, ystart:ystop, xstart:xstop]
        # Get the wavelength corresponding to each plane in the image.
        image_wl = read_image_wl(s_flat_model, quadrant)

        if image_wl.max() < MICRONS_100:
            log.warning("Wavelengths in s_flat image appear to be in meters")
        flat_2d, s_flat_dq, s_flat_err = interpolate_flat(image_flat, image_dq,
                                                          image_err, image_wl, wl)
    else:
        flat_2d = full_array_flat[ystart:ystop, xstart:xstop]
        s_flat_dq = full_array_dq[ystart:ystop, xstart:xstop]
        s_flat_err = full_array_err[ystart:ystop, xstart:xstop]

    # Find pixels in the flat that have a value of NaN and add to
    # DQ mask, DO_NOT_USE + NO_FLAT_FIELD
    bad_flag = dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['NO_FLAT_FIELD']

    flat_nan = np.isnan(flat_2d)
    s_flat_dq[flat_nan] = np.bitwise_or(s_flat_dq[flat_nan], bad_flag)

    # Find pixels in the flat have have a value of zero, and add to
    # DQ mask, DO_NOT_USE + NO_FLAT_FIELD
    flat_zero = np.where(flat_2d == 0.)
    s_flat_dq[flat_zero] = np.bitwise_or(s_flat_dq[flat_zero], bad_flag)

    # Find all pixels in the flat that have a DQ value of DO_NOT_USE
    flat_bad = np.bitwise_and(s_flat_dq, dqflags.pixel['DO_NOT_USE'])

    # Reset the flat value of all bad pixels to 1.0, so that no
    # correction is made
    flat_2d[np.where(flat_bad)] = 1.0

    s_flat, s_flat_dq = combine_fast_slow(wl, flat_2d, s_flat_dq,
                                          tab_wl, tab_flat, dispaxis)

    return s_flat, s_flat_dq, s_flat_err


def detector_flat(wl, d_flat_model,
                  xstart, xstop, ystart, ystop,
                  exposure_type, dispaxis, slit_name):
    """Extract the flat for the detector part.

    Parameters
    ----------
    wl : 2-D ndarray
        Wavelength at each pixel of the 2-D slit array.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    xstart, ystart : int
        Starting pixel numbers (zero indexed) for the slice containing
        the data for the current slit.

    xstop, ystop : int
        End of the slice containing the data for the current slit.  The
        start and stop values are Python slice notation, i.e. the region
        to be extracted is [ystart:ystop, xstart:xstop].

    exposure_type : str
        The exposure type refers to fixed_slit, IFU, or using the
        micro-shutter array.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    slit_name : str
        The name of the slit currently being processed.

    Returns
    -------
    d_flat : ndarray, 2d, float32
        The computed flat field for the detector.

    d_flat_dq : ndarray, 2d, uint32, or None
        The associated data quality array.

    d_flat_err : ndarray, 2d, float32, or None
        The associated error array.
    """

    if d_flat_model is None:
        d_flat = np.ones(wl.shape, dtype=np.float32)
        d_flat_dq = None
        d_flat_err = None
        return d_flat, d_flat_dq, d_flat_err

    quadrant = None

    if xstart >= xstop or ystart >= ystop:
        return 1., None

    (tab_wl, tab_flat) = read_flat_table(d_flat_model, exposure_type,
                                         slit_name, quadrant)
    if tab_wl.max() < MICRONS_100:
        log.warning("Wavelengths in d_flat table appear to be in meters.")

    full_array_flat = d_flat_model.data
    full_array_dq = d_flat_model.dq
    full_array_err = d_flat_model.err
    image_flat = full_array_flat[:, ystart:ystop, xstart:xstop]
    image_dq = full_array_dq[..., ystart:ystop, xstart:xstop]
    image_err = full_array_err[..., ystart:ystop, xstart:xstop]
    # Get the wavelength corresponding to each plane in the image.
    image_wl = read_image_wl(d_flat_model, quadrant)
    if image_wl.max() < MICRONS_100:
        log.warning("Wavelengths in d_flat image appear to be in meters.")

    flat_2d, d_flat_dq, d_flat_err = interpolate_flat(image_flat, image_dq,
                                                      image_err, image_wl, wl)

    # Find pixels in the flat that have a value of NaN and add to
    # DQ mask, DO_NOT_USE + NO_FLAT_FIELD
    bad_flag = dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['NO_FLAT_FIELD']

    flat_nan = np.isnan(flat_2d)
    d_flat_dq[flat_nan] = np.bitwise_or(d_flat_dq[flat_nan], bad_flag)

    # Find pixels in the flat have have a value of zero, and add to
    # DQ mask, DO_NOT_USE + NO_FLAT_FIELD
    flat_zero = np.where(flat_2d == 0.)
    d_flat_dq[flat_zero] = np.bitwise_or(d_flat_dq[flat_zero], bad_flag)

    # Find all pixels in the flat that have a DQ value of DO_NOT_USE
    flat_bad = np.bitwise_and(d_flat_dq, dqflags.pixel['DO_NOT_USE'])

    # Reset the flat value of all bad pixels to 1.0, so that no
    # correction is made
    flat_2d[np.where(flat_bad)] = 1.0

    d_flat, d_flat_dq = combine_fast_slow(wl, flat_2d, d_flat_dq,
                                          tab_wl, tab_flat, dispaxis)
    return d_flat, d_flat_dq, d_flat_err


def combine_dq(f_flat_dq, s_flat_dq, d_flat_dq, default_shape):
    """Combine non-None DQ arrays via bitwise_or.

    Parameters
    ----------
    f_flat_dq : ndarray
        The DQ array for the fore optics component.

    s_flat_dq : ndarray
        The DQ array for the spectrograph component.

    d_flat_dq : ndarray
        The DQ array for the detector component.

    default_shape : tuple
        If all three of the DQ arrays (see above) are None, use this shape
        to create a DQ array filled with zero.

    Returns
    -------
    flat_dq : ndarray, 2-D, uint32
        The DQ array resulting from combining the input DQ arrays via
        bitwise OR.
    """
    BADFLAT = dqflags.pixel['NO_FLAT_FIELD'] | dqflags.pixel['DO_NOT_USE']
    BADFLAT = BADFLAT | dqflags.pixel['UNRELIABLE_FLAT']

    dq_list = []
    if f_flat_dq is not None:
        dq_list.append(f_flat_dq)
    if s_flat_dq is not None:
        dq_list.append(s_flat_dq)
    if d_flat_dq is not None:
        dq_list.append(d_flat_dq)
    n_dq = len(dq_list)

    # Combine the component flat dq arrays.  If there are none, make a
    # dq array with all BADFLAT bits set
    flat_dq = np.zeros(default_shape, dtype=np.uint32)
    if n_dq == 0:
        flat_dq = np.bitwise_or(flat_dq, BADFLAT)
    else:
        for dq_component in dq_list:
            flat_dq = np.bitwise_or(flat_dq, dq_component)

    # Flag DO_NOT_USE where some or all of the flats had NO_FLAT_FIELD set
    do_not_use_loc = np.where(np.bitwise_and(flat_dq, dqflags.pixel['NO_FLAT_FIELD']))
    flat_dq[do_not_use_loc] = np.bitwise_or(flat_dq[do_not_use_loc], dqflags.pixel['DO_NOT_USE'])

    # Flag DO_NOT_USE, NO_FLAT_FIELD and UNRELIABLE_FLAT where some or all the
    # flats had DO_NOT_USE set
    iloc = np.where(np.bitwise_and(flat_dq, dqflags.pixel['DO_NOT_USE']))
    flat_dq[iloc] = np.bitwise_or(flat_dq[iloc], BADFLAT)

    return flat_dq


def read_image_wl(flat_model, quadrant=None):
    """Read wavelengths for the image planes.

    Extended summary
    ----------------
    This function should only be called if the SCI array in `flat_model`
    is 3-D.  The purpose is to get the wavelength for each of the planes
    in the SCI array.

    Parameters
    ----------
    flat_model : ~jwst.datamodels.NirspecFlatModel or ~jwst.datamodels.NirspecQuadFlatModel
        Flat field for the current component.

    quadrant : int (0, 1, 2, or 3)
        The quadrant of the micro-shutter array.  This is only needed for
        fore-optics for MSA (MOS) data.

    Returns
    -------
    wavelength : ndarray, 1-D
        An array of wavelengths, one for each plane of the SCI array.
    """

    if quadrant is not None:  # NRS_MSASPEC
        wavelength = flat_model.quadrants[quadrant].wavelength["wavelength"]
    else:
        wavelength = flat_model.wavelength["wavelength"]

    if len(wavelength.shape) > 1:
        n = wavelength.shape[-1]
        try:
            wl = wavelength.reshape((n,))
        except ValueError:
            log.error("Image wavelength array has shape %s; "
                      "don't know how to interpret that.",
                      str(wavelength.shape))
            raise ValueError("Expected either a scalar column or just one row.")
        wavelength = wl
    # The assumption here is that any NaN or non-positive wavelengths will
    # only be at the end of the array.  If there are embedded NaNs or zero
    # or negative wavelengths, the following will result in the wavelengths
    # becoming out of synch with the flat-field data.  In this case, it
    # will be necessary to also filter (along the first image axis) and
    # return the flat-field data.
    filter1 = np.logical_not(np.isnan(wavelength))  # skip NaNs
    wavelength = wavelength[filter1]
    filter2 = (wavelength > 0.)
    wavelength = wavelength[filter2]

    return wavelength


def read_flat_table(flat_model, exposure_type, slit_name=None, quadrant=None):
    """Read the table (the "fast" variation).

    Parameters
    ----------
    flat_model : NIRSpec flat-field object
        This contains the flat field table from which we will read the
        "fast" variation flat-field data.

    exposure_type : str
        The exposure type refers to fixed_slit, IFU, or using the
        micro-shutter array.  In this function we just need to check for
        fixed-slit types.

    slit_name : str
        The name of the slit.  This is only needed for fixed-slit data, in
        which case it is used for selecting the relevant row of the table.

    quadrant : int (0, 1, 2, or 3)
        The quadrant of the micro-shutter array.  This is only needed for
        fore-optics for MSA (MOS) data.

    Returns
    -------
    tab_wl : ndarray, 1-D, float
        The column of wavelengths read from the fast-variation table.

    tab_flat : ndarray, 1-D, float
        The column of flat_field values read from the fast-variation table.
        `tab_wl` and `tab_flat` should be the same length.
    """

    if quadrant is not None:  # NRS_MSASPEC
        data = flat_model.quadrants[quadrant].flat_table
    else:
        data = flat_model.flat_table

    try:
        slit_col = data["slit_name"]
    except KeyError:
        slit_col = None
    try:
        nelem_col = data["nelem"]
    except KeyError:
        nelem_col = None
    wl_col = data["wavelength"]
    flat_col = data["data"]

    nrows = len(wl_col)
    row = None
    # Note that it's only for fixed-slit data that we need to select the
    # row based on the slit name.
    if exposure_type in FIXED_SLIT_TYPES and slit_col is not None and slit_name is not None:
        slit_name_lc = slit_name.lower()
        for i in range(nrows):
            # Note:  The .strip() is a workaround.  As of the time of
            # writing, the value of a text string may have trailing blanks
            # if there is only one row in the table.
            column_value = slit_col[i].lower().strip()
            if column_value == "any" or column_value == slit_name_lc:
                row = i
                break
        if row is None:
            log.error("Slit name %s not found in flat field table", slit_name)
            raise ValueError("{} not found in SLIT_NAME column (nor was 'ANY')"
                             .format(slit_name))

    nelem = None  # initial value
    if row is not None:
        # Table contains arrays; use the row that was found above.
        tab_wl = wl_col[row].copy()
        tab_flat = flat_col[row].copy()
        if nelem_col is not None:
            nelem = nelem_col[row]
    else:
        # There was no SLIT_NAME column, or the data are not fixed-slit.
        if len(wl_col.shape) > 1:
            # Table contains arrays, but there should be only one row.
            tab_wl = wl_col[0].copy()
            tab_flat = flat_col[0].copy()
            if nelem_col is not None:
                nelem = nelem_col[0]
        else:
            # Table contains scalar columns.
            tab_wl = wl_col.copy()
            tab_flat = flat_col.copy()
            if nelem_col is not None:
                nelem = nelem_col[0]  # arbitrary choice of row
    if nelem is not None:
        if len(tab_wl) < nelem:
            log.error("The fast_variation array size %d in "
                      "the data model is too small, and table data were "
                      "truncated.", len(tab_wl))
            nelem = len(tab_wl)  # truncated!
        else:
            tab_wl = tab_wl[:nelem]
            tab_flat = tab_flat[:nelem]
    else:
        nelem = len(tab_wl)

    # Trailing dummy rows should have been taken care of via nelem above,
    # but only if an nelem column was present and was set correctly.
    filter1 = np.logical_not(np.isnan(tab_wl))  # skip NaNs
    filter2 = np.logical_not(np.isnan(tab_flat))
    filter = np.logical_and(filter1, filter2)
    n1 = filter.sum(dtype=np.intp)
    if n1 != nelem:
        log.debug("The table wavelength or flat-field data array contained "
                  "%d NaNs; these have been skipped.", nelem - n1)
        tab_wl = tab_wl[filter]
        tab_flat = tab_flat[filter]
    del filter1, filter2, filter
    # Skip zero or negative wavelengths, and skip zero flat-field values.
    filter1 = (tab_wl > 0.)
    filter2 = (tab_flat != 0.)
    filter = np.logical_and(filter1, filter2)
    n2 = filter.sum(dtype=np.intp)
    if n2 != n1:
        log.debug("The table wavelength or flat-field data array contained "
                  "%d zero or negative values; these have been skipped.",
                  n1 - n2)
        tab_wl = tab_wl[filter]
        tab_flat = tab_flat[filter]
    del filter1, filter2, filter

    # Check that the wavelengths are increasing.  This is a requirement
    # for using np.searchsorted (see wl_interpolate).
    if len(tab_wl) > 1:
        diff = tab_wl[1:] - tab_wl[0:-1]
        if np.any(diff <= 0.):
            log.warning("Wavelengths in the fast-variation table "
                        "must be strictly increasing.")

    return tab_wl, tab_flat


def combine_fast_slow(wl, flat_2d, flat_dq, tab_wl, tab_flat, dispaxis):
    """Multiply the image by the tabular values.

    Parameters
    ----------
    wl : 2-D ndarray
        Wavelength at each pixel of the 2-D slit array.

    flat_2d : 2-D ndarray or a scalar (float)
        The flat field derived from the image part of the reference file,
        or a scalar (e.g. 1) if there is no image part in the current
        reference file.

    flat_dq : ndarray or None
        If not None, the data quality array corresponding to `flat_2d`.
        A copy of this will be updated with flags which may be set based
        on the fast-variation component, and the updated array will be
        returned.

    tab_wl : ndarray, 1-D
        Wavelengths corresponding to `tab_flat`.

    tab_flat : ndarray, 1-D
        The flat field from the table part of the reference file.  This
        is the "fast" variation of the flat, i.e. fast with respect to
        wavelength.

    dispaxis : int
        1 is horizontal, 2 is vertical.

    Returns
    -------
    flat_2d * values : ndarray, 2-D, float32
        The product of `flat_2d` and the values in `tab_flat` interpolated
        to the wavelengths of the science image, i.e. `wl`.

    combined_dq : ndarray, 2-D, uint32
        The updated data quality array corresponding to `flat_2d`.  If a
        pixel wavelength is less than or equal to zero, or if it's not
        within the range of `tab_wl`, NO_FLAT_FIELD will be used to flag
        this condition, and the fast-variation flat field value at that
        pixel will be set to 1.
    """

    wl_c = clean_wl(wl, dispaxis)
    dwl = np.zeros_like(wl_c)

    if flat_dq is None:
        combined_dq = np.zeros(wl.shape, dtype=np.uint32)
    else:
        combined_dq = flat_dq.copy()

    if dispaxis == HORIZONTAL:
        dwl[:, 0:-1] = wl_c[:, 1:] - wl_c[:, 0:-1]
        dwl[:, -1] = dwl[:, -2]
    elif dispaxis == VERTICAL:
        dwl[0:-1, :] = wl_c[1:, :] - wl_c[0:-1, :]
        dwl[-1, :] = dwl[-2, :]

    # Values averaged within tab_flat.
    values = np.zeros_like(wl_c)
    (ny, nx) = wl_c.shape
    # Abscissas and weights for 3-point Gaussian integration, but taking
    # the width of the interval to be 1, so the result will be the average
    # over the interval.
    d = math.sqrt(0.6) / 2.
    dx = np.array([-d, 0., d])
    wgt = np.array([5., 8., 5.]) / 18.
    for j in range(ny):
        for i in range(nx):
            if wl[j, i] <= 0.:  # note:  wl, not wl_c
                values[j, i] = 1.
            else:
                # average the tabular data over the range of wavelengths
                temp = g_average(wl_c[j, i], dwl[j, i],
                                 tab_wl, tab_flat, dx, wgt)
                if temp is None:
                    values[j, i] = 1.
                    combined_dq[j, i] |= dqflags.pixel['NO_FLAT_FIELD']
                    combined_dq[j, i] |= dqflags.pixel['DO_NOT_USE']
                else:
                    values[j, i] = temp

    return flat_2d * values, combined_dq


def clean_wl(wl, dispaxis):
    """Replace zeros and/or NaNs in the wavelength array.

    Parameters
    ----------
    wl : 2-D ndarray
        Wavelength at each pixel of the 2-D slit array.

    dispaxis : int
        1 is horizontal, 2 is vertical.

    Returns
    -------
    wl_c : 2-D array
        A copy of `wl`, but with zero and negative values replaced with
        an average wavelength.  For each column (row) in the dispersion
        direction, the average to find a replacement value is taken along
        the cross-dispersion direction.
    """

    wl_c = wl.copy()  # so we can replace zeros
    wl_c[wl_c <= 0.] = np.nan
    shape = wl_c.shape

    if dispaxis == HORIZONTAL:
        i0 = None
        i1 = None
        for i in range(shape[1]):
            temp = wl_c[:, i]
            if np.any(np.isfinite(temp)):
                if i0 is None:
                    i0 = i  # first i with some non-zero wl
                i1 = i  # last i (so far) with some non-zero wl
                replacement_value = np.nanmean(temp)
                temp[np.isnan(temp)] = replacement_value
        if i0 is not None and i0 > 0:
            temp = wl_c[:, i0].reshape(shape[0], 1)
            wl_c[:, 0:i0] = temp.copy()
        if i1 is not None and i1 < shape[1] - 1:
            temp = wl_c[:, i1].reshape(shape[0], 1)
            wl_c[:, i1:] = temp.copy()
    elif dispaxis == VERTICAL:
        j0 = None
        j1 = None
        for j in range(shape[0]):
            temp = wl_c[j, :]
            if np.any(np.isfinite(temp)):
                if j0 is None:
                    j0 = j  # first j with some non-zero wl
                j1 = j  # last j (so far) with some non-zero wl
                replacement_value = np.nanmean(temp)
                temp[np.isnan(temp)] = replacement_value
        if j0 is not None and j0 > 0:
            temp = wl_c[j0, :].reshape(1, shape[1])
            wl_c[0:j0, :] = temp.copy()
        if j1 is not None and j1 < shape[0] - 1:
            temp = wl_c[j1, :].reshape(1, shape[1])
            wl_c[j1:, :] = temp.copy()
    else:
        wl_c = wl.copy()

    return wl_c


def g_average(wl0, dwl0, tab_wl, tab_flat, dx, wgt):
    """Gaussian integration.

    Parameters
    ----------
    wl0 : float
        Wavelength at the center of the current pixel.

    dwl0 : float
        Width (in wavelength units) of the current pixel.

    tab_wl : ndarray, 1-D
        Array of wavelengths corresponding to `tab_flat` flat-field values.

    tab_flat : ndarray, 1-D
        Array of flat-field values.

    dx : ndarray, 1-D
        Array of offsets within a pixel, e.g. -0.3873, 0.0, +0.3873

    wgt : ndarray, 1-D
        Array of weights, e.g. 5/18, 8/18, 5/18

    Returns
    -------
    float or None
        The average value of `tab_flat` over the current pixel.  None
        will be returned if any of the wavelengths used for computing
        the average is outside the range of wavelengths in `tab_wl`.
    """

    npts = len(dx)
    wavelengths = wl0 + dwl0 * dx
    sum = 0.
    for k in range(npts):
        value = wl_interpolate(wavelengths[k], tab_wl, tab_flat)
        if value is None:  # wavelengths[k] was out of bounds
            return None
        sum += (value * wgt[k])

    return sum


def wl_interpolate(wavelength, tab_wl, tab_flat):
    """Interpolate the flat field at the specified wavelength.

    Extended summary
    ----------------
    Linear interpolation is used.

    Parameters
    ----------
    wavelength : float
        The wavelength (microns) at which to find the flat-field value.

    tab_wl : ndarray, 1-D
        Array of wavelengths corresponding to `tab_flat` flat-field values.
        These are assumed to be strictly increasing.

    tab_flat : ndarray, 1-D
        Array of flat-field values.

    Returns
    -------
    float or None
        The flat-field value (from `tab_flat`) at `wavelength`.
        None will be returned if `wavelength` is not positive or is
        outside the range of `tab_wl`.
    """

    if wavelength <= 0. or wavelength < tab_wl[0] or wavelength > tab_wl[-1]:
        return None

    n0 = np.searchsorted(tab_wl, wavelength) - 1
    p = (wavelength - tab_wl[n0]) / (tab_wl[n0 + 1] - tab_wl[n0])
    q = 1. - p

    return q * tab_flat[n0] + p * tab_flat[n0 + 1]


def interpolate_flat(image_flat, image_dq, image_err, image_wl, wl):
    """Interpolate within the 3-D flat field image to get a 2-D flat.

    Parameters
    ----------
    image_flat : ndarray, 3-D
        This is slice [:, ystart:ystop, xstart:xstop] of the flat field
        reference image.  This slice covers the spatial extent of the
        extracted 2-D spectrum and includes all of the wavelength axis
        (the first axis) of the reference image.

    image_dq : ndarray, 2-D or 3-D
        This is slice [..., ystart:ystop, xstart:xstop] of the data quality
        array for the flat field reference image.

    image_err : ndarray, 2-D or 3-D
        This is slice [..., ystart:ystop, xstart:xstop] of the data quality
        array for the flat field reference image.

    image_wl : ndarray, 1-D
        The wavelength for each plane of the flat field reference image.

    wl : ndarray, 2-D
        The wavelength at each pixel of the 2-D extracted science spectrum.

    Returns
    -------
    flat_2d : ndarray, 2-D, float32
        The flat field, interpolated over wavelength, same shape as `wl`.
        Divide the 2-D extracted spectrum by this array to correct for
        flat-field variations.

    flat_dq : ndarray, 2-D, uint32
        The data quality array corresponding to `flat_2d`.

    flat_err : ndarray, 2-D, float32
        The error array corresponding to `flat_2d`.
    """

    if len(image_flat.shape) < 3:
        return image_flat, image_dq, image_err

    (nz, ysize, xsize) = image_flat.shape
    if nz == 1:
        if len(image_dq.shape) == 2:
            # dq and err arrays are the same size, so treat them the same
            return image_flat.reshape((ysize, xsize)), image_dq, image_err
        else:
            return (image_flat.reshape((ysize, xsize)),
                    image_dq.reshape((ysize, xsize)),
                    image_err.reshape((ysize, xsize)))

    grid = np.indices((ysize, xsize), dtype=np.intp)
    ixpixel = grid[1]
    iypixel = grid[0]

    # The initial value of -1 is a flag to indicate that elements have not
    # been assigned valid values yet.
    k = np.zeros(wl.shape, dtype=np.intp) - 1

    # Truncate the index for wavelengths that are outside the range of
    # image_wl.  The indices need to be assigned harmless values to avoid
    # indexing out of bounds.
    #   Why do we set the upper limit of k to nz - 2?
    #   Because we interpolate using elements k and k + 1.

    # for wavelengths < lower limit (image_wl[0]) set k to 0
    k[:, :] = np.where(wl <= image_wl[0], 0, k)

    # for wavelengths > upper limit (image_wl[nz-1]) set k to nz-2
    k[:, :] = np.where(wl >= image_wl[nz - 1], nz - 2, k)

    # Look for the correct interval for linear interpolation.
    for k_test in range(nz - 1):
        test1 = np.logical_and(wl >= image_wl[k_test],
                               wl < image_wl[k_test + 1])

        # If an element of k is not -1, it has already been assigned, and
        # I don't want to clobber it.
        test2 = np.logical_and(k == -1, test1)
        k[:, :] = np.where(test2, k_test, k)

        if np.all(k >= 0):
            break
    # Use linear interpolation within the 3-D flat field to get a 2-D
    # flat field.
    denom = image_wl[k + 1] - image_wl[k]
    zero_denom = (denom == 0.)
    denom = np.where(zero_denom, 1., denom)

    # linear interpolation equation
    # p is the linear interpolation wavelength scaling factor
    # flat = flat[k] + (flat[k+1] - flat[k])*p
    # flat = flat[k] - flat[k]*p + flat[k+1]*p
    # flat = (1-p)*flat[k] + p*flat[k+1]

    p = np.where(zero_denom, 0., (wl - image_wl[k]) / denom)
    q = 1. - p
    flat_2d = q * image_flat[k, iypixel, ixpixel] + \
        p * image_flat[k + 1, iypixel, ixpixel]
    if len(image_err.shape) == 2:
        flat_err = image_err.copy()
    else:
        flat_err = q * image_err[k, iypixel, ixpixel] + \
            p * image_err[k + 1, iypixel, ixpixel]

    if len(image_dq.shape) == 2:
        flat_dq = image_dq.copy()
    else:
        flat_dq = np.where(p == 0.,
                           image_dq[k, iypixel, ixpixel],
                           np.bitwise_or(image_dq[k, iypixel, ixpixel],
                                         image_dq[k + 1, iypixel, ixpixel]))

        flat_bad = np.bitwise_and(flat_dq, dqflags.pixel['DO_NOT_USE'])
        # Reset the flat value of all bad pixels to 1.0, so that no
        # correction is made
        flat_2d[np.where(flat_bad)] = 1.0

    # If the wavelength at a pixel is outside the range of wavelengths
    # for the reference image, flag the pixel as bad.  Note that this will
    # also result in the computed flat field being set to 1.
    mask = (wl < image_wl[0])
    flat_dq[mask] = np.bitwise_or(flat_dq[mask], dqflags.pixel['DO_NOT_USE'])
    mask = (wl > image_wl[-1])
    flat_dq[mask] = np.bitwise_or(flat_dq[mask], dqflags.pixel['DO_NOT_USE'])

    # If a pixel is flagged as bad, applying flat_2d should not make any
    # change to the science data.
    flat_bad = np.bitwise_and(flat_dq, dqflags.pixel['DO_NOT_USE'])
    flat_2d[np.where(flat_bad)] = 1.0

    return flat_2d.astype(image_flat.dtype), flat_dq, flat_err


def flat_for_nirspec_ifu(output_model, f_flat_model, s_flat_model, d_flat_model,
                         dispaxis):
    """Create the interpolated flat for NIRSpec IFU

    Parameters
    ----------
    output_model : JWST data model
        Science data model, modified (flat fielded) in-place.

    f_flat_model : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    Returns
    -------
    flat, flat_dq, flat_err, any_updated : numpy.array, numpy.array, numpy.array, bool
        4-tuple of the interpolated flat correction and whether any slice of the IFU
        is actually affected.
    """
    any_updated = False
    exposure_type = output_model.meta.exposure.type
    flat = np.ones_like(output_model.data)
    flat_dq = np.zeros_like(output_model.dq)
    flat_err = np.zeros_like(output_model.data)

    try:
        list_of_wcs = nirspec.nrs_ifu_wcs(output_model)
    except (KeyError, AttributeError):
        if output_model.meta.cal_step.assign_wcs == 'COMPLETE':
            log.error("The input file does not appear to have WCS info.")
            raise RuntimeError("Problem accessing WCS information.")
        else:
            log.error("This mode %s requires WCS information.", exposure_type)
            raise RuntimeError("The assign_wcs step has not been run.")
    for (k, ifu_wcs) in enumerate(list_of_wcs):

        # example:  bounding_box = ((1600.5, 2048.5),   # X
        #                           (1886.5, 1925.5))   # Y
        truncated = False
        try:
            xstart = ifu_wcs.bounding_box[0][0]
            xstop = ifu_wcs.bounding_box[0][1]
            ystart = ifu_wcs.bounding_box[1][0]
            ystop = ifu_wcs.bounding_box[1][1]
            log.debug("Using ifu_wcs.bounding_box.")
        except AttributeError:
            log.info("ifu_wcs.bounding_box not found; using domain instead.")
            xstart = ifu_wcs.domain[0]['lower']
            xstop = ifu_wcs.domain[0]['upper']
            ystart = ifu_wcs.domain[1]['lower']
            ystop = ifu_wcs.domain[1]['upper']

        if xstart < -0.5:
            truncated = True
            log.info("xstart from WCS bounding_box was %g" % xstart)
            xstart = 0.
        if ystart < -0.5:
            truncated = True
            log.info("ystart from WCS bounding_box was %g" % ystart)
            ystart = 0.
        if xstop > 2047.5:
            truncated = True
            log.info("xstop from WCS bounding_box was %g" % xstop)
            xstop = 2047.
        if ystop > 2047.5:
            truncated = True
            log.info("ystop from WCS bounding_box was %g" % ystop)
            ystop = 2047.
        if truncated:
            log.info("WCS bounding_box for stripe %d extended beyond image "
                     "edges, has been truncated to ...", k)
            log.info('  xstart=%g, xstop=%g, ystart=%g, ystop=%g',
                     xstart, xstop, ystart, ystop)

        # Convert these to integers, and add one to the upper limits,
        # because we want to use these as slice limits.
        xstart = int(math.ceil(xstart))
        xstop = int(math.floor(xstop)) + 1
        ystart = int(math.ceil(ystart))
        ystop = int(math.floor(ystop)) + 1

        dx = xstop - xstart
        dy = ystop - ystart
        ind = np.indices((dy, dx))
        x = ind[1] + xstart
        y = ind[0] + ystart
        coords = ifu_wcs(x, y)
        wl = coords[2]
        nan_flag = np.isnan(wl)
        good_flag = np.logical_not(nan_flag)
        if wl[good_flag].max() < MICRONS_100:
            log.warning("Wavelengths in WCS table appear to be in meters")
        # Set NaNs to a relatively harmless value, but don't modify nan_flag.
        wl[nan_flag] = 0.

        flat_2d, flat_dq_2d, flat_err_2d = create_flat_field(
            wl, f_flat_model, s_flat_model, d_flat_model,
            xstart, xstop, ystart, ystop,
            exposure_type, dispaxis, None, None)
        flat_2d[nan_flag] = 1.
        mask = (flat_2d <= 0.)
        nbad = mask.sum(dtype=np.intp)
        if nbad > 0:
            log.debug("%d flat-field values <= 0", nbad)
            flat_2d[mask] = 1.
        del mask

        flat[ystart:ystop, xstart:xstop][good_flag] = flat_2d[good_flag]
        if flat_dq.dtype == flat_dq_2d.dtype:
            flat_dq[ystart:ystop, xstart:xstop] |= flat_dq_2d.copy()
        else:
            log.warning("flat_dq.dtype = {}  flat_dq_2d.dtype = {}"
                        .format(flat_dq.dtype, flat_dq_2d.dtype))
            flat_dq[ystart:ystop, xstart:xstop] |= \
                flat_dq_2d.astype(flat_dq.dtype).copy()
        flat_err[ystart:ystop, xstart:xstop][good_flag] = flat_err_2d[good_flag]
        del nan_flag, good_flag

        any_updated = True

    # That's all folks
    return flat, flat_dq, flat_err, any_updated


def flat_for_nirspec_brightobj(output_model, f_flat_model, s_flat_model, d_flat_model,
                               dispaxis):
    """Create the interpolated flat for NIRSpec IFU

    Parameters
    ----------
    output_model : JWST data model
        Science data model, modified (flat fielded) in-place.

    f_flat_model : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    Returns
    -------
    flat, flat_dq, flat_err, any_updated : numpy.array, numpy.array, numpy.array, bool
        4-tuple of the interpolated flat correction and whether any slice of the IFU
        is actually affected.
    """

    exposure_type = output_model.meta.exposure.type

    got_wcs = (hasattr(output_model.meta, "wcs") and
               output_model.meta.wcs is not None)

    # Create an output model for the interpolated flat fields.
    interpolated_flats = datamodels.ImageModel()
    interpolated_flats.update(output_model, only="PRIMARY")
    if got_wcs:
        interpolated_flats.meta.wcs = output_model.meta.wcs

    slit_name = output_model.name

    # The input may be either 2-D or 3-D; save `shape` for use later.
    shape = output_model.data.shape
    ysize, xsize = shape[-2:]
    # pixels with respect to the original image
    xstart = output_model.meta.subarray.xstart - 1 + output_model.xstart - 1
    ystart = output_model.meta.subarray.ystart - 1 + output_model.ystart - 1
    xstop = xstart + xsize
    ystop = ystart + ysize

    # The wavelength of each pixel in a plane of the data.
    got_wl_attribute = True
    try:
        wl = output_model.wavelength.copy()  # a 2-D array
    except AttributeError:
        got_wl_attribute = False
    if not got_wl_attribute or len(wl) == 0:
        got_wl_attribute = False

    # There must be either a wavelength array or a meta.wcs.
    if not got_wl_attribute or np.nanmin(wl) == 0. and np.nanmax(wl) == 0.:
        log.warning("The wavelength array has not been populated,")
        if got_wcs:
            log.warning("so using wcs instead of the wavelength array.")
            grid = np.indices((ysize, xsize), dtype=np.float64)
            (ra, dec, wl) = output_model.meta.wcs(grid[1], grid[0])
            del ra, dec, grid
        else:
            log.warning("and there is no 'wcs' attribute,")
            if output_model.meta.cal_step.assign_wcs == 'COMPLETE':
                log.warning("assign_wcs has been run, however.")
            else:
                log.warning("likely because assign_wcs has not been run.")
            log.error("Skipping flat_field.")
            output_model.meta.cal_step.flat_field = 'SKIPPED'
            return None
    else:
        log.debug("Wavelengths are from the wavelength array.")

    nan_mask = np.isnan(wl)
    good_mask = np.logical_not(nan_mask)
    sum_nan_mask = nan_mask.sum(dtype=np.intp)
    sum_good_mask = good_mask.sum(dtype=np.intp)
    if sum_nan_mask > 0:
        log.debug("Number of NaNs in wavelength array = %d out of %d",
                  sum_nan_mask, sum_nan_mask + sum_good_mask)
        if sum_good_mask < 1:
            log.warning("(all are NaN)")
        # Replace NaNs with a relatively harmless but out-of-bounds value.
        wl[nan_mask] = 0.

    # Combine the three flat fields.  The same flat will be applied to
    # each plane (integration) in the cube.
    flat_2d, flat_dq_2d, flat_err_2d = create_flat_field(
        wl, f_flat_model, s_flat_model, d_flat_model,
        xstart, xstop, ystart, ystop,
        exposure_type, dispaxis, slit_name, None)
    mask = (flat_2d <= 0.)
    nbad = mask.sum(dtype=np.intp)
    if nbad > 0:
        log.debug("%d flat-field values <= 0", nbad)
        flat_2d[mask] = 1.
    del mask

    flat_dq_2d = flat_dq_2d.astype(output_model.dq.dtype)

    interpolated_flats.data = flat_2d
    interpolated_flats.dq = flat_dq_2d
    interpolated_flats.err = flat_err_2d.astype(output_model.err.dtype)
    interpolated_flats.wavelength = wl

    return interpolated_flats


def flat_for_nirspec_slit(slit, f_flat_model, s_flat_model, d_flat_model,
                          dispaxis, exposure_type, slit_nt, subarray,
                          use_wavecorr):
    """Create the interpolated flat for NIRSpec slit data

    Parameters
    ----------
    slit : SlitModel
        A slit to process

    f_flat_model : ~jwst.datamodels.NirspecFlatModel, ~jwst.datamodels.NirspecQuadFlatModel, or None
        Flat field for the fore optics.

    s_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the spectrograph.

    d_flat_model : ~jwst.datamodels.NirspecFlatModel or None
        Flat field for the detector.

    dispaxis : int
        1 means horizontal dispersion, 2 means vertical dispersion.

    exposure_type : str
        The exposure type

    slit_nt : namedtuple or None
        For MOS data (only), this is used to get the quadrant number and
        the indices of the current shutter in the Y and X directions.

    subarray : DataModel.meta.subarray
        The subarray specification

    use_wavecorr : boolean
        Flag indicating whether or not to use the corrected wavelengths
        provided (upstream) by the wavecorr step.

    Returns
    -------
    flat : SlitModel or None
        The calculated flat. If None, no flat field could be determined
    """
    # Create flat and flat dq arrays with default values
    flat_2d = np.ones_like(slit.data)
    flat_dq_2d = np.zeros_like(slit.dq)
    flat_err_2d = np.zeros_like(slit.err)

    # pixels with respect to the original image
    ysize, xsize = slit.data.shape
    xstart = slit.xstart - 1 + subarray.xstart - 1
    ystart = slit.ystart - 1 + subarray.ystart - 1
    xstop = xstart + xsize
    ystop = ystart + ysize

    got_wcs = hasattr(slit.meta, "wcs") and slit.meta.wcs is not None

    # Get the wavelength at each pixel in the extracted slit data.
    # If the wavelength attribute exists and is populated, use it
    # in preference to the wavelengths returned by the wcs function.
    got_wl_attribute = True
    try:
        wl = slit.wavelength.copy()  # a 2-D array
    except AttributeError:
        got_wl_attribute = False
    if not got_wl_attribute or len(wl) == 0:
        got_wl_attribute = False
    return_dummy = False

    # Has the use_wavecorr param been set?
    if use_wavecorr is not None:
        if use_wavecorr:
            # Need to use the 2D wavelength array, because that's where
            # the corrected wavelengths are stored
            if got_wl_attribute:
                # We've got the "wl" wavelength array we need
                pass
            else:
                # Can't do the computation without the 2D wavelength array
                log.error(f"The wavelength array for slit {slit.name} is not populated")
                log.error("Skipping flat-field correction")
                return_dummy = True
        elif not use_wavecorr:
            # Need to use the WCS object to create an uncorrected 2D wavelength array
            if got_wcs:
                log.info(f"Creating wavelength array from WCS for slit {slit.name}")
                bb = slit.meta.wcs.bounding_box
                grid = grid_from_bounding_box(bb)
                wl = slit.meta.wcs(*grid)[2]
                del grid
            else:
                # Can't create the uncorrected wavelengths without the WCS
                log.error(f"Slit {slit.name} has no WCS object")
                log.error("Skipping flat-field correction")
                return_dummy = True
    else:
        # use_wavecorr was not specified, so use default processing
        if not got_wl_attribute or np.nanmin(wl) == 0. and np.nanmax(wl) == 0.:
            got_wl_attribute = False
            log.warning(f"The wavelength array for slit {slit.name} has not been populated")
            # Try to create it from the WCS
            if got_wcs:
                bb = slit.meta.wcs.bounding_box
                grid = grid_from_bounding_box(bb)
                wl = slit.meta.wcs(*grid)[2]
                del grid
            else:
                log.warning("and this slit does not have a 'wcs' attribute")
                log.warning("likely because assign_wcs has not been run.")
                log.error("skipping ...")
                return_dummy = True
        else:
            log.debug("Wavelengths are from the wavelength array.")

    # Create and return a dummy flat as a placeholder, if necessary
    if return_dummy:
        dummy_flat = datamodels.SlitModel(data=flat_2d, dq=flat_dq_2d, err=flat_err_2d)
        dummy_flat.name = slit.name
        dummy_flat.xstart = slit.xstart
        dummy_flat.xsize = slit.xsize
        dummy_flat.ystart = slit.ystart
        dummy_flat.ysize = slit.ysize
        dummy_flat.wavelength = np.zeros_like(slit.data)

        return dummy_flat

    # We've got everything we need for the rest of processing
    nan_mask = np.isnan(wl)
    good_mask = np.logical_not(nan_mask)
    sum_nan_mask = nan_mask.sum(dtype=np.intp)
    sum_good_mask = good_mask.sum(dtype=np.intp)
    if sum_nan_mask > 0:
        log.debug(f"Number of NaNs in sci wavelength array = {sum_nan_mask} "
                  f"out of {sum_nan_mask + sum_good_mask}")
        if sum_good_mask < 1:
            log.warning("(all are NaN)")
        # Replace NaNs with a relatively harmless but out-of-bounds value.
        wl[nan_mask] = 0.
    max_wavelength = np.nanmax(wl)
    if max_wavelength > 0. and max_wavelength < MICRONS_100:
        log.warning("Wavelengths in science data appear to be in meters.")

    # Combine the three flat fields for the current subarray.
    flat_2d, flat_dq_2d, flat_err_2d = create_flat_field(wl,
                                                         f_flat_model, s_flat_model, d_flat_model,
                                                         xstart, xstop, ystart, ystop,
                                                         exposure_type, dispaxis, slit.name, slit_nt)

    # Mask bad flatfield values
    mask = (flat_2d <= 0.)
    nbad = mask.sum(dtype=np.intp)
    if nbad > 0:
        log.debug("%d flat-field values <= 0", nbad)
        flat_2d[mask] = 1.
    del mask

    # Put the computed flat, flat_dq and flat_err into a datamodel
    new_flat = datamodels.SlitModel(data=flat_2d, dq=flat_dq_2d, err=flat_err_2d)
    new_flat.name = slit.name
    new_flat.xstart = slit.xstart
    new_flat.xsize = slit.xsize
    new_flat.ystart = slit.ystart
    new_flat.ysize = slit.ysize
    new_flat.wavelength = wl.copy()

    # Copy the WCS info from output (same as input).
    if got_wcs:
        new_flat.meta.wcs = slit.meta.wcs

    return new_flat
