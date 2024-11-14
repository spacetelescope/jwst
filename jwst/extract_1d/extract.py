import logging
import json
import numpy as np
from json.decoder import JSONDecodeError

from astropy.modeling import polynomial
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.apcorr import (
    MirLrsApcorrModel, MirMrsApcorrModel, NrcWfssApcorrModel, NrsFsApcorrModel,
    NrsMosApcorrModel, NrsIfuApcorrModel, NisWfssApcorrModel
)

from jwst.datamodels import ModelContainer
from jwst.lib import pipe_utils
from jwst.lib.wcs_utils import get_wavelengths
from jwst.extract_1d import extract1d, spec_wcs, utils
from jwst.extract_1d.apply_apcorr import select_apcorr


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

WFSS_EXPTYPES = ['NIS_WFSS', 'NRC_WFSS', 'NRC_GRISM']
"""Exposure types to be regarded as wide-field slitless spectroscopy."""

ANY = "ANY"
"""Wildcard for slit name.

Extended summary
----------------
For full-frame input data, keyword SLTNAME may not be populated, so the
slit name will be set to this string to indicate that the first slit in
the reference file should be used.
A slit name in the extract1d reference file can also be ANY, in which case the
reference information for that slit will be regarded as matching any slit
name from the input data.
"""

ANY_ORDER = 1000
"""Wildcard for spectral order number in a reference image.

Extended summary
----------------
If the extract1d reference file contains images, keyword SPORDER gives the order
number of the spectrum that would be extracted using a given image in
the reference file.
"""

HORIZONTAL = 1
VERTICAL = 2
"""Dispersion direction, predominantly horizontal or vertical."""

# These values are assigned in get_extract_parameters, using key "match".
# If there was an aperture in the reference file for which the "id" key matched, that's (at least) a partial match.
# If "spectral_order" also matched, that's an exact match.
NO_MATCH = "no match"
PARTIAL = "partial match"
EXACT = "exact match"


class Extract1dError(Exception):
    pass


# Create custom error to pass continue from a function inside of a loop
class ContinueError(Exception):
    pass


def open_extract1d_ref(refname):
    """Open the extract1d reference file.

    Parameters
    ----------
    refname : str
        The name of the extract1d reference file, or 'N/A'.
        If specified, this file is expected to be a JSON file
        giving extraction information.

    Returns
    -------
    ref_dict : dict or None
        If the extract1d reference file is specified, ref_dict will be the
        dictionary returned by json.load().
    """

    refname_type = refname[-4:].lower()
    if refname == "N/A":
        ref_dict = None
    else:
        # If specified, the extract1d reference file can only be in json format
        if refname_type == 'json':
            fd = open(refname)
            try:
                ref_dict = json.load(fd)
                fd.close()
            except (UnicodeDecodeError, JSONDecodeError):
                # Input file does not load correctly as json file.
                # Probably an error in json file
                fd.close()
                log.error("Extract1d json reference file has an error, run a json validator off line and fix the file")
                raise RuntimeError("Invalid json extract 1d reference file, run json validator off line and fix file.")
        else:
            log.error("Invalid Extract 1d reference file, must be json.")
            raise RuntimeError("Invalid Extract 1d reference file, must be json.")

    return ref_dict


def open_apcorr_ref(refname, exptype):
    """Determine the appropriate DataModel class to use when opening the input APCORR reference file.

    Parameters
    ----------
    refname : str
        Path of the APCORR reference file

    exptype : str
        EXPTYPE of the input to the extract_1d step.

    Returns
    -------
    Opened APCORR DataModel.

    Notes
    -----
    This function should be removed after the DATAMODL keyword is required for the APCORR reference file.

    """
    apcorr_model_map = {
        'MIR_LRS-FIXEDSLIT': MirLrsApcorrModel,
        'MIR_LRS-SLITLESS': MirLrsApcorrModel,
        'MIR_MRS': MirMrsApcorrModel,
        'NRC_GRISM': NrcWfssApcorrModel,
        'NRC_WFSS': NrcWfssApcorrModel,
        'NIS_WFSS': NisWfssApcorrModel,
        'NRS_BRIGHTOBJ': NrsFsApcorrModel,
        'NRS_FIXEDSLIT': NrsFsApcorrModel,
        'NRS_IFU': NrsIfuApcorrModel,
        'NRS_MSASPEC': NrsMosApcorrModel
    }

    apcorr_model = apcorr_model_map[exptype]
    return apcorr_model(refname)


def get_extract_parameters(
        ref_dict,
        input_model,
        slitname,
        sp_order,
        meta,
        smoothing_length,
        bkg_fit,
        bkg_order,
        use_source_posn,
        subtract_background
):
    """Get extract1d reference file values.

    Parameters
    ----------
    ref_dict : dict or None
        For an extract1d reference file in JSON format, `ref_dict` will be
        the entire contents of the file.  If there is no extract1d reference
        file, `ref_dict` will be None.

    input_model : data model
        This can be either the input science file or one SlitModel out of
        a list of slits.

    slitname : str
        The name of the slit, or "ANY"

    sp_order : int
        The spectral order number.

    meta : metadata for the actual input model, i.e. not just for the
        current slit.

    smoothing_length : int or None
        Width of a boxcar function for smoothing the background regions.
        If None, the smoothing length will be gotten from `ref_dict`, or
        it will be set to 0 (no background smoothing) if this key is
        not found in `ref_dict`.
        If `smoothing_length` is not None, that means that the user
        explicitly specified the value, so that value will be used.
        This argument is only used if background regions have been
        specified.

    bkg_fit : str
        The type of fit to apply to background values in each
        column (or row, if the dispersion is vertical). The default
        `poly` results in a polynomial fit of order `bkg_order`. Other
        options are `mean` and `median`. If `mean` or `median` is selected,
        `bkg_order` is ignored.

    bkg_order : int or None
        Polynomial order for fitting to each column (or row, if the
        dispersion is vertical) of background.  If None, the polynomial
        order will be gotten from `ref_dict`, or it will be set to 0 if
        not found in `ref_dict`.
        A value of 0 means that a simple average of the background
        regions, column by column (or row by row), will be used.
        If `bkg_order` is not None, that means that the user explicitly
        specified the value, so that value will be used.
        This argument must be positive or zero, and it is only used if
        background regions have been specified.

    use_source_posn : bool or None
        If True, the target and background positions specified in `ref_dict`
        (or a default target position) will be shifted to account for
        the actual source location in the data.
        If None, the value specified in `ref_dict` will be used, or it will
        be set to True if not found in `ref_dict`.

    subtract_background : bool
        If False, all background parameters will be ignored.

    Returns
    -------
    extract_params : dict
        Information copied out of `ref_dict`.  The items will be selected
        based on `slitname` and `sp_order`.  Default values will be
        assigned if `ref_dict` is None.  For a reference image, the key
        'ref_image' gives the (open) image model.
    """

    extract_params = {'match': NO_MATCH}  # initial value
    if ref_dict is None:
        # There is no extract1d reference file; use "reasonable" default values.
        extract_params['match'] = EXACT
        shape = input_model.data.shape
        extract_params['spectral_order'] = sp_order
        extract_params['xstart'] = 0  # first pixel in X
        extract_params['xstop'] = shape[-1] - 1  # last pixel in X
        extract_params['ystart'] = 0  # first pixel in Y
        extract_params['ystop'] = shape[-2] - 1  # last pixel in Y
        extract_params['extract_width'] = None
        extract_params['src_coeff'] = None
        extract_params['bkg_coeff'] = None  # no background sub.
        extract_params['smoothing_length'] = 0  # because no background sub.
        extract_params['bkg_fit'] = None  # because no background sub.
        extract_params['bkg_order'] = 0  # because no background sub.
        extract_params['subtract_background'] = False
        extract_params['extraction_type'] = 'box'

        if use_source_posn is None:
            extract_params['use_source_posn'] = False
        else:
            extract_params['use_source_posn'] = use_source_posn

        extract_params['position_correction'] = 0
        extract_params['independent_var'] = 'pixel'
        # Note that extract_params['dispaxis'] is not assigned.
        # This will be done later, possibly slit by slit.

    else:
        for aper in ref_dict['apertures']:
            if ('id' in aper and aper['id'] != "dummy" and
                    (aper['id'] == slitname or aper['id'] == "ANY" or
                     slitname == "ANY")):
                extract_params['match'] = PARTIAL

                # region_type is retained for backward compatibility; it is
                # not required to be present.
                region_type = aper.get("region_type", "target")
                if region_type != "target":
                    continue

                # spectral_order is a secondary selection criterion.  The
                # default is the expected value, so if the key is not present
                # in the JSON file, the current aperture will be selected.
                # If the current aperture in the JSON file has
                # "spectral_order": "ANY", that aperture will be selected.
                spectral_order = aper.get("spectral_order", sp_order)

                if spectral_order == sp_order or spectral_order == ANY:
                    extract_params['match'] = EXACT
                    extract_params['spectral_order'] = sp_order
                    # Note: extract_params['dispaxis'] is not assigned.
                    # This is done later, possibly slit by slit.

                    # Set default start/stop by shape if not specified
                    shape = input_model.data.shape
                    extract_params['xstart'] = aper.get('xstart', 0)
                    extract_params['xstop'] = aper.get('xstop', shape[-1] - 1)
                    extract_params['ystart'] = aper.get('ystart', 0)
                    extract_params['ystop'] = aper.get('ystop', shape[-2] - 1)

                    extract_params['src_coeff'] = aper.get('src_coeff')
                    extract_params['bkg_coeff'] = aper.get('bkg_coeff')
                    if (extract_params['bkg_coeff'] is not None
                            and subtract_background is not False):
                        extract_params['subtract_background'] = True
                        if bkg_fit is not None:
                            extract_params['bkg_fit'] = bkg_fit
                        else:
                            extract_params['bkg_fit'] = aper.get('bkg_fit', 'poly')
                    else:
                        extract_params['bkg_fit'] = None
                        extract_params['subtract_background'] = False

                    extract_params['independent_var'] = aper.get('independent_var', 'pixel').lower()

                    if bkg_order is None:
                        extract_params['bkg_order'] = aper.get('bkg_order', 0)
                    else:
                        # If the user supplied a value, use that value.
                        extract_params['bkg_order'] = bkg_order

                    # Set use_source_posn based on hierarchy of priorities:
                    # parameter value on the command line is highest precedence,
                    # then parameter value from the extract1d reference file,
                    # and finally a default setting based on exposure type.
                    use_source_posn_aper = aper.get('use_source_posn', None)  # value from the extract1d ref file
                    if use_source_posn is None:  # no value set on command line
                        if use_source_posn_aper is None:  # no value set in ref file
                            # Use a suitable default
                            if meta.exposure.type in ['MIR_LRS-FIXEDSLIT', 'NRS_FIXEDSLIT', 'NRS_MSASPEC']:
                                use_source_posn = True
                                log.info(f"Turning on source position correction for exp_type = {meta.exposure.type}")
                            else:
                                use_source_posn = False
                        else:  # use the value from the ref file
                            use_source_posn = use_source_posn_aper
                    extract_params['use_source_posn'] = use_source_posn

                    extract_params['extract_width'] = aper.get('extract_width')
                    extract_params['position_correction'] = 0  # default value

                    if smoothing_length is None:
                        extract_params['smoothing_length'] = aper.get('smoothing_length', 0)
                    else:
                        # If the user supplied a value, use that value.
                        extract_params['smoothing_length'] = smoothing_length

                    # Default extraction type to box
                    extract_params['extraction_type'] = 'box'

                    break

    return extract_params


def log_initial_parameters(extract_params):
    """Log some of the initial extraction parameters.

    Parameters
    ----------
    extract_params : dict
        Information read from the reference file.
    """
    if "xstart" not in extract_params:
        return

    log.debug("Extraction parameters:")
    log.debug(f"dispaxis = {extract_params['dispaxis']}")
    log.debug(f"spectral order = {extract_params['spectral_order']}")
    log.debug(f"initial xstart = {extract_params['xstart']}")
    log.debug(f"initial xstop = {extract_params['xstop']}")
    log.debug(f"initial ystart = {extract_params['ystart']}")
    log.debug(f"initial ystop = {extract_params['ystop']}")
    log.debug(f"extract_width = {extract_params['extract_width']}")
    log.debug(f"initial src_coeff = {extract_params['src_coeff']}")
    log.debug(f"initial bkg_coeff = {extract_params['bkg_coeff']}")
    log.debug(f"bkg_fit = {extract_params['bkg_fit']}")
    log.debug(f"bkg_order = {extract_params['bkg_order']}")
    log.debug(f"smoothing_length = {extract_params['smoothing_length']}")
    log.debug(f"independent_var = {extract_params['independent_var']}")
    log.debug(f"use_source_posn = {extract_params['use_source_posn']}")
    log.debug(f"extraction_type = {extract_params['extraction_type']}")


def create_poly(coeff):
    """Create a polynomial model from coefficients.

    Parameters
    ----------
    coeff : list of float
        The coefficients of the polynomial, constant term first, highest
        order term last.

    Returns
    -------
    `astropy.modeling.polynomial.Polynomial1D` object, or None if `coeff`
        is empty.
    """
    n = len(coeff)

    if n < 1:
        return None

    coeff_dict = {f'c{i}': coeff[i] for i in range(n)}

    return polynomial.Polynomial1D(degree=n - 1, **coeff_dict)


def run_extract1d(
        input_model,
        extract_ref_name,
        apcorr_ref_name,
        smoothing_length,
        bkg_fit,
        bkg_order,
        log_increment,
        subtract_background,
        use_source_posn
):
    """Extract 1-D spectra.

    Parameters
    ----------
    input_model : data model
        The input science model.

    extract_ref_name : str
        The name of the extract1d reference file, or "N/A".

    apcorr_ref_name : str
        Name of the APCORR reference file. Default is None

    smoothing_length : int or None
        Width of a boxcar function for smoothing the background regions.

    bkg_fit : str
        Type of fitting to apply to background values in each column
        (or row, if the dispersion is vertical).

    bkg_order : int or None
        Polynomial order for fitting to each column (or row, if the
        dispersion is vertical) of background. Only used if `bkg_fit`
        is `poly`.

    log_increment : int
        if `log_increment` is greater than 0 and the input data are
        multi-integration, a message will be written to the log every
        `log_increment` integrations.

    subtract_background : bool or None
        User supplied flag indicating whether the background should be
        subtracted.
        If None, the value in the extract_1d reference file will be used.
        If not None, this parameter overrides the value in the
        extract_1d reference file.

    use_source_posn : bool or None
        If True, the target and background positions specified in the
        reference file (or the default position, if there is no reference
        file) will be shifted to account for source position offset.

    Returns
    -------
    output_model : data model
        A new MultiSpecModel containing the extracted spectra.

    """
    # Set "meta_source" to either the first model in a container,
    # or the individual input model, for convenience
    # of retrieving meta attributes in subsequent statements
    if isinstance(input_model, ModelContainer):
        meta_source = input_model[0]
    else:
        meta_source = input_model

    # Get the exposure type
    exp_type = meta_source.meta.exposure.type

    # Read in the extract1d reference file.
    extract_ref_dict = open_extract1d_ref(extract_ref_name)

    # Read in the aperture correction reference file
    apcorr_ref_model = None
    if apcorr_ref_name is not None and apcorr_ref_name != 'N/A':
        apcorr_ref_model = open_apcorr_ref(apcorr_ref_name, exp_type)

    # Set up the output model
    output_model = datamodels.MultiSpecModel()
    if hasattr(meta_source, "int_times"):
        output_model.int_times = meta_source.int_times.copy()
    output_model.update(meta_source, only='PRIMARY')

    # This will be relevant if we're asked to extract a spectrum
    # and the spectral order is zero.
    # That's only OK if the disperser is a prism.
    prism_mode = is_prism(meta_source)

    # Handle inputs that contain one or more slit models
    if isinstance(input_model, (ModelContainer, datamodels.MultiSlitModel)):

        is_multiple_slits = True
        if isinstance(input_model, ModelContainer):
            slits = input_model
        else:
            slits = input_model.slits

        # Save original use_source_posn value, because it can get
        # toggled within the following loop over slits
        save_use_source_posn = use_source_posn

        for slit in slits:  # Loop over the slits in the input model
            log.info(f'Working on slit {slit.name}')
            log.debug(f'Slit is of type {type(slit)}')

            slitname = slit.name
            use_source_posn = save_use_source_posn  # restore original value

            if np.size(slit.data) <= 0:
                log.info(f'No data for slit {slit.name}, skipping ...')
                continue

            sp_order = get_spectral_order(slit)
            if sp_order == 0 and not prism_mode:
                log.info("Spectral order 0 is a direct image, skipping ...")
                continue

            try:
                output_model = create_extraction(
                    extract_ref_dict, slit, slitname, sp_order,
                    smoothing_length, bkg_fit, bkg_order, use_source_posn,
                    exp_type, subtract_background, meta_source,
                    output_model, apcorr_ref_model, log_increment,
                    is_multiple_slits
                )
            except ContinueError:
                continue

    else:
        # Define source of metadata
        slit = None
        is_multiple_slits = False

        # These default values for slitname are not really slit names,
        # and slitname may be assigned a better value below, in the
        # sections for input_model being an ImageModel or a SlitModel.
        slitname = exp_type
        if slitname is None:
            slitname = ANY

        if isinstance(input_model, datamodels.ImageModel):
            if hasattr(input_model, "name"):
                slitname = input_model.name

            sp_order = get_spectral_order(input_model)
            if sp_order == 0 and not prism_mode:
                log.info("Spectral order 0 is a direct image, skipping ...")
            else:
                log.info(f'Processing spectral order {sp_order}')
                try:
                    output_model = create_extraction(
                        extract_ref_dict, slit, slitname, sp_order,
                        smoothing_length, bkg_fit, bkg_order, use_source_posn,
                        exp_type, subtract_background, input_model,
                        output_model, apcorr_ref_model, log_increment,
                        is_multiple_slits
                    )
                except ContinueError:
                    pass

        elif isinstance(input_model, (datamodels.CubeModel, datamodels.SlitModel)):
            # This branch will be invoked for inputs that are a CubeModel, which typically includes
            # NIRSpec BrightObj (fixed slit) mode, as well as inputs that are a
            # single SlitModel, which typically includes data from a single resampled/combined slit
            # instance from level-3 processing of NIRSpec fixed slits and MOS modes.

            # Replace the default value for slitname with a more accurate value, if possible.
            # For NRS_BRIGHTOBJ, the slit name comes from the slit model info
            if exp_type == 'NRS_BRIGHTOBJ' and hasattr(input_model, "name"):
                slitname = input_model.name

            # For NRS_FIXEDSLIT, the slit name comes from the FXD_SLIT keyword
            # in the model meta if not present in the input model
            if exp_type == 'NRS_FIXEDSLIT':
                if hasattr(input_model, "name") and input_model.name is not None:
                    slitname = input_model.name
                else:
                    slitname = input_model.meta.instrument.fixed_slit

            sp_order = get_spectral_order(input_model)
            if sp_order == 0 and not prism_mode:
                log.info("Spectral order 0 is a direct image, skipping ...")
            else:
                log.info(f'Processing spectral order {sp_order}')

                try:
                    output_model = create_extraction(
                        extract_ref_dict, slit, slitname, sp_order,
                        smoothing_length, bkg_fit, bkg_order, use_source_posn,
                        exp_type, subtract_background, input_model,
                        output_model, apcorr_ref_model, log_increment,
                        is_multiple_slits
                    )
                except ContinueError:
                    pass

        else:
            log.error("The input file is not supported for this step.")
            raise RuntimeError("Can't extract a spectrum from this file.")

    # Copy the integration time information from the INT_TIMES table to keywords in the output file.
    if pipe_utils.is_tso(input_model):
        populate_time_keywords(input_model, output_model)
    else:
        log.debug("Not copying from the INT_TIMES table because this is not a TSO exposure.")
        if hasattr(output_model, "int_times"):
            del output_model.int_times

    output_model.meta.wcs = None  # See output_model.spec[i].meta.wcs instead.

    if apcorr_ref_model is not None:
        apcorr_ref_model.close()

    # Remove target.source_type from the output model, so that it
    # doesn't force creation of an empty SCI extension in the output
    # x1d product just to hold this keyword.
    output_model.meta.target.source_type = None

    return output_model


def populate_time_keywords(input_model, output_model):
    """Copy the integration times keywords to header keywords.

    Parameters
    ----------
    input_model : data model
        The input science model.

    output_model : data model
        The output science model.  This may be modified in-place.

    """
    nints = input_model.meta.exposure.nints
    int_start = input_model.meta.exposure.integration_start

    if hasattr(input_model, 'data'):
        shape = input_model.data.shape

        if len(shape) == 2:
            num_integ = 1
        else:  # len(shape) == 3
            num_integ = shape[0]
    else:  # e.g. MultiSlit data
        num_integ = 1

    # This assumes that the spec attribute of output_model has already been created,
    # and spectra have been appended.
    n_output_spec = len(output_model.spec)

    # num_j is the number of spectra per integration, e.g. the number
    # of fixed-slit spectra, MSA spectra, or different
    # spectral orders; num_integ is the number of integrations.
    # The total number of output spectra is n_output_spec = num_integ * num_j
    num_j = n_output_spec // num_integ

    if n_output_spec != num_j * num_integ:  # sanity check
        log.warning(
            f"populate_time_keywords:  Don't understand n_output_spec = {n_output_spec}, num_j = {num_j}, num_integ = "
            f"{num_integ}"
        )
    else:
        log.debug(
            f"Number of output spectra = {n_output_spec}; number of spectra for each integration = {num_j}; "
            f"number of integrations = {num_integ}"
        )

    if int_start is None:
        log.warning("INTSTART not found; assuming a value of 1.")
        int_start = 1

    int_start -= 1  # zero indexed
    int_end = input_model.meta.exposure.integration_end

    if int_end is None:
        log.warning(f"INTEND not found; assuming a value of {nints}.")
        int_end = nints

    int_end -= 1  # zero indexed

    if nints > 1:
        num_integrations = int_end - int_start + 1
    else:
        num_integrations = 1

    if hasattr(input_model, 'int_times') and input_model.int_times is not None:
        nrows = len(input_model.int_times)
    else:
        nrows = 0

    if nrows < 1:
        log.warning("There is no INT_TIMES table in the input file - "
                    "Making best guess on integration numbers.")
        for j in range(num_j):  # for each spectrum or order
            for k in range(num_integ):  # for each integration
                output_model.spec[(j * num_integ) + k].int_num = k + 1  # set int_num to (k+1) - 1-indexed integration
        return

    # If we have a single plane (e.g. ImageModel or MultiSlitModel),
    # we will only populate the keywords if the corresponding uncal file
    # had one integration.
    # If the data were or might have been segmented, we use the first and
    # last integration numbers to determine whether the data were in fact
    # averaged over integrations, and if so, we should not populate the
    # int_times-related header keywords.
    skip = False  # initial value

    if isinstance(input_model, (datamodels.MultiSlitModel, datamodels.ImageModel)):
        if num_integrations > 1:
            log.warning("Not using INT_TIMES table because the data have been averaged over integrations.")
            skip = True
    elif isinstance(input_model, (datamodels.CubeModel, datamodels.SlitModel)):
        shape = input_model.data.shape

        if len(shape) == 2 and num_integrations > 1:
            log.warning("Not using INT_TIMES table because the data have been averaged over integrations.")
            skip = True
        elif len(shape) != 3 or shape[0] > nrows:
            # Later, we'll check that the integration_number column actually
            # has a row corresponding to every integration in the input.
            log.warning(
                "Not using INT_TIMES table because the data shape is not consistent with the number of table rows."
            )
            skip = True
    elif isinstance(input_model, datamodels.IFUCubeModel):
        log.warning("The INT_TIMES table will be ignored for IFU data.")
        skip = True

    if skip:
        return

    int_num = input_model.int_times['integration_number']
    start_time_mjd = input_model.int_times['int_start_MJD_UTC']
    mid_time_mjd = input_model.int_times['int_mid_MJD_UTC']
    end_time_mjd = input_model.int_times['int_end_MJD_UTC']
    start_tdb = input_model.int_times['int_start_BJD_TDB']
    mid_tdb = input_model.int_times['int_mid_BJD_TDB']
    end_tdb = input_model.int_times['int_end_BJD_TDB']

    data_range = (int_start, int_end)  # Inclusive range of integration numbers in the input data, zero indexed.

    # Inclusive range of integration numbers in the INT_TIMES table, zero indexed.
    table_range = (int_num[0] - 1, int_num[-1] - 1)
    offset = data_range[0] - table_range[0]

    if data_range[0] < table_range[0] or data_range[1] > table_range[1]:
        log.warning("Not using the INT_TIMES table because it does not include rows for all integrations in the data.")
        return

    log.debug("TSO data, so copying times from the INT_TIMES table.")

    n = 0  # Counter for spectra in output_model.

    for k in range(num_integ):  # for each spectrum or order
        for j in range(num_j):  # for each integration
            row = k + offset
            spec = output_model.spec[n]  # n is incremented below
            spec.int_num = int_num[row]
            spec.time_scale = "UTC"
            spec.start_time_mjd = start_time_mjd[row]
            spec.mid_time_mjd = mid_time_mjd[row]
            spec.end_time_mjd = end_time_mjd[row]
            spec.start_tdb = start_tdb[row]
            spec.mid_tdb = mid_tdb[row]
            spec.end_tdb = end_tdb[row]
            n += 1


def get_spectral_order(slit):
    """Get the spectral order number.

    Parameters
    ----------
    slit : SlitModel object
        One slit from an input MultiSlitModel or similar.

    Returns
    -------
    int
        Spectral order number for `slit`.  If no information about spectral
        order is available in `wcsinfo`, a default value of 1 will be
        returned.

    """
    if hasattr(slit.meta, 'wcsinfo'):
        sp_order = slit.meta.wcsinfo.spectral_order

        if sp_order is None:
            log.warning("spectral_order is None; using 1")
            sp_order = 1
    else:
        log.warning("slit.meta doesn't have attribute wcsinfo; setting spectral order to 1")
        sp_order = 1

    return sp_order


def is_prism(input_model):
    """Determine whether the current observing mode used a prism.

    Extended summary
    ----------------
    The reason for this test is so we can skip spectral extraction if the
    spectral order is zero and the exposure was not made using a prism.
    In this context, therefore, a grism is not considered to be a prism.

    Parameters
    ----------
    input_model : data model
        The input science model.

    Returns
    -------
    bool
        True if the exposure used a prism; False otherwise.

    """
    instrument = input_model.meta.instrument.name

    if instrument is None:
        return False

    instrument_filter = input_model.meta.instrument.filter

    if instrument_filter is None:
        instrument_filter = "NONE"
    else:
        instrument_filter = instrument_filter.upper()

    grating = input_model.meta.instrument.grating

    if grating is None:
        grating = "NONE"
    else:
        grating = grating.upper()

    prism_mode = False

    if ((instrument == "MIRI" and instrument_filter.find("P750L") >= 0) or
            (instrument == "NIRSPEC" and grating.find("PRISM") >= 0)):
        prism_mode = True

    return prism_mode


def copy_keyword_info(slit, slitname, spec):
    """Copy metadata from the input to the output spectrum.

    Parameters
    ----------
    slit : A SlitModel object
        Metadata will be copied from the input `slit` to output `spec`.

    slitname : str or None
        The name of the slit.

    spec : One element of MultiSpecModel.spec
        Metadata attributes will be updated in-place.

    """
    if slitname is not None and slitname != "ANY":
        spec.name = slitname

    if hasattr(slit, "slitlet_id"):
        spec.slitlet_id = slit.slitlet_id

    if hasattr(slit, "source_id"):
        spec.source_id = slit.source_id

    if hasattr(slit, "source_name") and slit.source_name is not None:
        spec.source_name = slit.source_name

    if hasattr(slit, "source_alias") and slit.source_alias is not None:
        spec.source_alias = slit.source_alias

    if hasattr(slit, "source_type") and slit.source_type is not None:
        spec.source_type = slit.source_type

    if hasattr(slit, "stellarity") and slit.stellarity is not None:
        spec.stellarity = slit.stellarity

    if hasattr(slit, "source_xpos"):
        spec.source_xpos = slit.source_xpos

    if hasattr(slit, "source_ypos"):
        spec.source_ypos = slit.source_ypos

    if hasattr(slit, "source_ra"):
        spec.source_ra = slit.source_ra

    if hasattr(slit, "source_dec"):
        spec.source_dec = slit.source_dec

    if hasattr(slit, "shutter_state"):
        spec.shutter_state = slit.shutter_state


def _set_weight_from_limits(profile, idx, lower_limit, upper_limit, allow_partial=True):
    # Both limits are inclusive
    profile[(idx >= lower_limit) & (idx <= upper_limit)] = 1.0

    if allow_partial:
        for partial_pixel_weight in [lower_limit - idx, idx - upper_limit]:
            test = (partial_pixel_weight > 0) & (partial_pixel_weight < 1)
            profile[test] = 1 - partial_pixel_weight[test]


def box_profile(shape, extract_params, wl_array, coefficients='src_coeff',
                label='aperture', return_limits=False):
    # Get pixel index values for the array
    yidx, xidx = np.mgrid[:shape[0], :shape[1]]
    if extract_params['dispaxis'] == HORIZONTAL:
        dval = yidx
    else:
        dval = xidx

    # Get start/stop values from parameters if present,
    # or default to data shape
    xstart = extract_params.get('xstart', 0)
    xstop = extract_params.get('xstop', shape[1] - 1)
    ystart = extract_params.get('ystart', 0)
    ystop = extract_params.get('ystop', shape[0] - 1)

    # Check if the profile should contain partial pixel weights
    if coefficients == 'bkg_coeff' and extract_params['bkg_fit'] == 'median':
        allow_partial = False
    else:
        allow_partial = True

    # Set aperture region, in this priority order:
    # 1. src_coeff upper and lower limits (or bkg_coeff, for background profile)
    # 2. center of start/stop values +/- extraction width
    # 3. start/stop values
    profile = np.full(shape, 0.0)
    if extract_params[coefficients] is not None:
        # Limits from source coefficients: ignore ystart/stop/width
        if extract_params['independent_var'].startswith("wavelength"):
            ival = wl_array
        elif extract_params['dispaxis'] == HORIZONTAL:
            ival = xidx.astype(np.float32)
        else:
            ival = yidx.astype(np.float32)

        # The source extraction can include more than one region,
        # but must contain pairs of lower and upper limits.
        n_src_coeff = len(extract_params[coefficients])
        if n_src_coeff % 2 != 0:
            raise RuntimeError(f"{coefficients} must contain alternating "
                               f"lists of lower and upper limits.")

        lower = None
        lower_limit = None
        upper_limit = None
        for i, coeff_list in enumerate(extract_params[coefficients]):
            if i % 2 == 0:
                lower = create_poly(coeff_list)
            else:
                upper = create_poly(coeff_list)

                # NOTE: source coefficients currently have a different
                # definition for pixel inclusion than the start/stop input.
                # Source coefficients define 0 at the center of the pixel,
                # start/stop defines 0 at the start of the lower pixel and
                # includes the upper pixel. Here, we are setting limits
                # in the spatial profile according to the more commonly used
                # start/stop definition, so we need to modify the lower and
                # upper limits from polynomial coefficients to match.
                lower_limit_region = lower(ival) + 0.5
                upper_limit_region = upper(ival) - 0.5

                _set_weight_from_limits(profile, dval, lower_limit_region,
                                        upper_limit_region,
                                        allow_partial=allow_partial)
                mean_lower = np.mean(lower_limit_region)
                mean_upper = np.mean(upper_limit_region)
                log.info(f'Mean {label} start/stop from {coefficients}: '
                         f'{mean_lower:.2f} -> {mean_upper:.2f} (inclusive)')

                if lower_limit is None:
                    lower_limit = mean_lower
                    upper_limit = mean_upper
                else:
                    if mean_lower < lower_limit:
                        lower_limit = mean_lower
                    if mean_upper > upper_limit:
                        upper_limit = mean_upper

    elif extract_params['extract_width'] is not None:
        # Limits from extraction width at center of ystart/stop if present,
        # center of array if not
        if extract_params['dispaxis'] == HORIZONTAL:
            nominal_middle = (ystart + ystop) / 2.0
        else:
            nominal_middle = (xstart + xstop) / 2.0

        width = extract_params['extract_width']
        lower_limit = nominal_middle - (width - 1.0) / 2.0
        upper_limit = lower_limit + width - 1

        _set_weight_from_limits(profile, dval, lower_limit, upper_limit)
        log.info(f'{label.capitalize()} start/stop: '
                 f'{lower_limit:.2f} -> {upper_limit:.2f} (inclusive)')

    else:
        # Limits from start/stop only
        if extract_params['dispaxis'] == HORIZONTAL:
            lower_limit = ystart
            upper_limit = ystop
        else:
            lower_limit = xstart
            upper_limit = xstop

        _set_weight_from_limits(profile, dval, lower_limit, upper_limit)
        log.info(f'{label.capitalize()} start/stop: '
                 f'{lower_limit:.2f} -> {upper_limit:.2f} (inclusive)')

    # Set weights to zero outside left and right limits
    if extract_params['dispaxis'] == HORIZONTAL:
        profile[:, :int(round(xstart))] = 0
        profile[:, int(round(xstop)) + 1:] = 0
    else:
        profile[:int(round(ystart)), :] = 0
        profile[int(round(ystop)) + 1:, :] = 0

    if return_limits:
        return profile, lower_limit, upper_limit
    else:
        return profile


def aperture_center(profile, dispaxis=1, middle_pix=None):
    if middle_pix is not None and np.sum(profile) > 0:
        spec_center = middle_pix
        if dispaxis == HORIZONTAL:
            slit_center = np.average(np.arange(profile.shape[0]),
                                     weights=profile[:, middle_pix])
        else:
            slit_center = np.average(np.arange(profile.shape[1]),
                                     weights=profile[middle_pix, :])
    else:
        yidx, xidx = np.mgrid[:profile.shape[0], :profile.shape[1]]
        if np.sum(profile) > 0:
            center_y = np.average(yidx, weights=profile)
            center_x = np.average(xidx, weights=profile)
        else:
            center_y = profile.shape[0] // 2
            center_x = profile.shape[1] // 2
        if dispaxis == HORIZONTAL:
            spec_center = center_y
            slit_center = center_x
        else:
            spec_center = center_x
            slit_center = center_y

    # if dispaxis == 1 (default), this returns center_x, center_y
    return slit_center, spec_center


def shift_by_source_location(location, nominal_location, extract_params):

    # Get the center of the nominal aperture
    offset = location - nominal_location
    log.info(f"Nominal location is {nominal_location:.2f}, "
             f"so offset is {offset:.2f} pixels")

    # Shift aperture limits by the difference between the
    # source location and the nominal center
    coeff_params = ['src_coeff', 'bkg_coeff']
    for params in coeff_params:
        if extract_params[params] is not None:
            for coeff_list in extract_params[params]:
                coeff_list[0] += offset
    if extract_params['dispaxis'] == HORIZONTAL:
        start_stop_params = ['ystart', 'ystop']
    else:
        start_stop_params = ['xstart', 'xstop']
    for params in start_stop_params:
        if extract_params[params] is not None:
            extract_params[params] += offset


def define_aperture(input_model, slit, extract_params, exp_type):
    if slit is None:
        data_model = input_model
    else:
        data_model = slit
    data_shape = data_model.data.shape[-2:]

    # Get a wavelength array for the data
    wl_array = get_wavelengths(data_model, exp_type, extract_params['spectral_order'])

    # Shift aperture definitions by source position if needed
    # Extract parameters are updated in place
    if extract_params['use_source_posn']:
        # Source location from WCS
        targ_ra, targ_dec = utils.get_target_coordinates(input_model, slit)
        middle_pix, middle_wl, location = utils.locn_from_wcs(input_model, slit, targ_ra, targ_dec)

        if location is not None:
            log.info(f"Computed source location is {location:.2f}, "
                     f"at pixel {middle_pix}, wavelength {middle_wl:.2f}")

            # Nominal location from extract params + located middle
            nominal_profile = box_profile(data_shape, extract_params, wl_array,
                                          label='nominal aperture')
            nominal_location, _ = aperture_center(
                nominal_profile, extract_params['dispaxis'], middle_pix=middle_pix)

            # Offet extract parameters by location - nominal
            shift_by_source_location(location, nominal_location, extract_params)
    else:
        middle_pix, middle_wl, location = None, None, None

    # Make a spatial profile, including source shifts if necessary
    profile, lower_limit, upper_limit = box_profile(data_shape, extract_params, wl_array,
                                                    return_limits=True)

    # Make sure profile weights are zero where wavelengths are invalid
    profile[~np.isfinite(wl_array)] = 0.0

    # Get the effective left and right limits from the profile weights
    nonzero_weight = np.where(np.sum(profile, axis=extract_params['dispaxis'] - 1) > 0)[0]
    if len(nonzero_weight) > 0:
        left_limit = nonzero_weight[0]
        right_limit = nonzero_weight[-1]
    else:
        left_limit = None
        right_limit = None

    # Make a background profile if necessary
    # (will also include source shifts)
    if (extract_params['subtract_background']
            and extract_params['bkg_coeff'] is not None):
        bg_profile = box_profile(data_shape, extract_params, wl_array,
                                 coefficients='bkg_coeff')
    else:
        bg_profile = None

    # Get 1D wavelength corresponding to the spatial profile
    mask = np.isnan(wl_array) | (profile == 0)
    masked_wl = np.ma.masked_array(wl_array, mask=mask)
    masked_weights = np.ma.masked_array(profile, mask=mask)
    if extract_params['dispaxis'] == HORIZONTAL:
        wavelength = np.average(masked_wl, weights=masked_weights, axis=0).filled(np.nan)
    else:
        wavelength = np.average(masked_wl, weights=masked_weights, axis=1).filled(np.nan)

    # Get RA and Dec corresponding to the center of the array,
    # weighted by the spatial profile
    center_x, center_y = aperture_center(profile, 1)
    coords = data_model.meta.wcs(center_x, center_y)
    ra = float(coords[0])
    dec = float(coords[1])

    # Return limits as a tuple with 4 elements: lower, upper, left, right
    limits = (lower_limit, upper_limit, left_limit, right_limit)

    return ra, dec, wavelength, profile, bg_profile, limits


def extract_one_slit(data_model, integ, profile, bg_profile, extract_params):
    """Extract data for one slit, or spectral order, or integration.

    Parameters
    ----------
    data_model : data model
        The input science model. May be a single slit from a MultiSlitModel
        (or similar), or a single data type, like an ImageModel, SlitModel,
        or CubeModel.

    integ : int
        For the case that data_model is a SlitModel or a CubeModel,
        `integ` is the integration number.  If the integration number is
        not relevant (i.e. the data array is 2-D), `integ` should be -1.

    profile : ndarray of float
        Spatial profile indicating the aperture location. Data is a
        2D image matching the input, with floating point values between 0
        and 1 assigning a weight to each pixel.  0 means the pixel is not used,
        1 means the pixel is fully included in the aperture.

    bg_profile : ndarray of float or None
        Background profile indicating any background regions to use, following
        the same format as the spatial profile. Ignored if
        extract_params['subtract_background'] is False.

    extract_params : dict
        Parameters read from the extract1d reference file.

    Returns
    -------
    sum_flux : ndarray, 1-D, float64
        The sum of the data values in the extraction region minus the sum
        of the data values in the background regions (scaled by the ratio
        of the numbers of pixels), for each pixel.
        The data values are usually in units of surface brightness,
        so this value isn't the flux, it's an intermediate value.
        Multiply `sum_flux` by the solid angle of a pixel to get the flux for a
        point source (column "flux").  Divide `sum_flux` by `npixels` (to
        compute the average) to get the array for the "surf_bright"
        (surface brightness) output column.

    f_var_poisson : ndarray, 1-D
        The extracted poisson variance values to go along with the
        sum_flux array.

    f_var_rnoise : ndarray, 1-D
        The extracted read noise variance values to go along with the
        sum_flux array.

    f_var_flat : ndarray, 1-D
        The extracted flat field variance values to go along with the
        sum_flux array.

    background : ndarray, 1-D
        The background count rate that was subtracted from the sum of
        the source data values to get `sum_flux`.

    b_var_poisson : ndarray, 1-D
        The extracted poisson variance values to go along with the
        background array.

    b_var_rnoise : ndarray, 1-D
        The extracted read noise variance values to go along with the
        background array.

    b_var_flat : ndarray, 1-D
        The extracted flat field variance values to go along with the
        background array.

    npixels : ndarray, 1-D, float64
        The number of pixels that were added together to get `sum_flux`,
        including any fractional pixels included via non-integer weights
        in the input profile.

    """
    # Get the data and variance arrays
    if integ > -1:
        log.info(f"Extracting integration {integ + 1}")
        data = data_model.data[integ]
        var_poisson = data_model.var_poisson[integ]
        var_rnoise = data_model.var_rnoise[integ]
        var_flat = data_model.var_flat[integ]
    else:
        data = data_model.data
        var_poisson = data_model.var_poisson
        var_rnoise = data_model.var_rnoise
        var_flat = data_model.var_flat

    # Transpose data for extraction
    if extract_params['dispaxis'] == HORIZONTAL:
        profile_view = profile
        bg_profile_view = bg_profile
    else:
        data = data.T
        profile_view = profile.T
        var_rnoise = var_rnoise.T
        var_poisson = var_poisson.T
        var_flat = var_flat.T
        if bg_profile is not None:
            bg_profile_view = bg_profile.T
        else:
            bg_profile_view = None

    # Extract spectra from the data
    result = extract1d.extract1d(data, [profile_view], var_rnoise, var_poisson, var_flat,
                                 profile_bg=bg_profile_view,
                                 bg_smooth_length=extract_params['smoothing_length'],
                                 fit_bkg=extract_params['subtract_background'],
                                 bkg_fit_type=extract_params['bkg_fit'],
                                 bkg_order=extract_params['bkg_order'],
                                 extraction_type=extract_params['extraction_type'])

    # Extraction routine can return multiple spectra;
    # here, we just want the first result
    first_result = []
    for r in result:
        first_result.append(r[0])
    return first_result


def create_extraction(
        extract_ref_dict,
        slit,
        slitname,
        sp_order,
        smoothing_length,
        bkg_fit,
        bkg_order,
        use_source_posn,
        exp_type,
        subtract_background,
        input_model,
        output_model,
        apcorr_ref_model,
        log_increment,
        is_multiple_slits
):
    if slit is None:
        data_model = input_model
    else:
        data_model = slit

    # Make sure NaNs and DQ flags match up in input
    pipe_utils.match_nans_and_flags(data_model)

    if exp_type in WFSS_EXPTYPES:
        instrument = input_model.meta.instrument.name
    else:
        instrument = data_model.meta.instrument.name
    if instrument is not None:
        instrument = instrument.upper()

    # We need a flag to indicate whether the photom step has been run.
    # If it hasn't, we'll copy the count rate to the flux column.
    try:
        s_photom = input_model.meta.cal_step.photom
    except AttributeError:
        s_photom = None

    if s_photom is not None and s_photom.upper() == 'COMPLETE':
        photom_has_been_run = True
        flux_units = 'Jy'
        f_var_units = 'Jy^2'
        sb_units = 'MJy/sr'
        sb_var_units = 'MJy^2 / sr^2'
    else:
        photom_has_been_run = False
        flux_units = 'DN/s'
        f_var_units = 'DN^2 / s^2'
        sb_units = 'DN/s'
        sb_var_units = 'DN^2 / s^2'
        log.warning("The photom step has not been run.")

    # Get the source type for the data
    if is_multiple_slits:
        source_type = data_model.source_type
    else:
        if isinstance(input_model, datamodels.SlitModel):
            source_type = input_model.source_type
            if source_type is None:
                source_type = input_model.meta.target.source_type
                input_model.source_type = source_type
        else:
            source_type = input_model.meta.target.source_type

    # Turn off use_source_posn if the source is not POINT
    if source_type != 'POINT' or exp_type in WFSS_EXPTYPES:
        use_source_posn = False
        log.info(f"Setting use_source_posn to False for exposure type {exp_type}, "
                 f"source type {source_type}")

    if photom_has_been_run:
        pixel_solid_angle = data_model.meta.photometry.pixelarea_steradians
        if pixel_solid_angle is None:
            pixel_solid_angle = 1.
            log.warning("Pixel area (solid angle) is not populated; the flux will not be correct.")
    else:
        pixel_solid_angle = 1.  # not needed

    extract_params = get_extract_parameters(
        extract_ref_dict, data_model, slitname, sp_order, input_model.meta,
        smoothing_length, bkg_fit,bkg_order, use_source_posn,
        subtract_background
    )

    if extract_params['match'] == NO_MATCH:
        log.critical('Missing extraction parameters.')
        raise ValueError('Missing extraction parameters.')
    elif extract_params['match'] == PARTIAL:
        log.info(f'Spectral order {sp_order} not found, skipping ...')
        raise ContinueError()

    extract_params['dispaxis'] = data_model.meta.wcsinfo.dispersion_direction
    if extract_params['dispaxis'] is None:
        log.warning("The dispersion direction information is missing, so skipping ...")
        raise ContinueError()

    # Set up spatial profiles and wavelength array,
    # to be used for every integration
    (ra, dec, wavelength, profile, bg_profile, limits) = define_aperture(
        input_model, slit, extract_params, exp_type)

    valid = ~np.isnan(wavelength)
    wavelength = wavelength[valid]
    if np.sum(valid) == 0:
        log.error("Spectrum is empty; no valid data.")
        raise ContinueError()

    # Set up aperture correction, to be used for every integration
    apcorr_available = False
    if source_type is not None and source_type.upper() == 'POINT' and apcorr_ref_model is not None:
        # NIRSpec needs to use a wavelength in the middle of the
        # range rather than the beginning of the range
        # for calculating the pixel scale since some wavelengths at the
        # edges of the range won't map to the sky
        if instrument == 'NIRSPEC':
            wl = np.median(wavelength)
        else:
            wl = wavelength.min()

        kwargs = {'location': (ra, dec, wl)}
        if not isinstance(input_model, datamodels.ImageModel):
            kwargs['slit_name'] = slitname
            if exp_type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
                kwargs['slit'] = slitname

        apcorr_model = select_apcorr(input_model)
        apcorr = apcorr_model(input_model, apcorr_ref_model.apcorr_table,
                              apcorr_ref_model.sizeunit, **kwargs)
    else:
        apcorr = None

    # Log the parameters before extracting
    log_initial_parameters(extract_params)

    # Set up integration iterations and progress messages
    progress_msg_printed = True
    shape = data_model.data.shape
    if len(shape) == 3 and shape[0] == 1:
        integrations = [0]
    elif len(shape) == 2:
        integrations = [-1]
    else:
        log.info(f"Beginning loop over {shape[0]} integrations ...")
        integrations = range(shape[0])
        progress_msg_printed = False

    # Extract each integration
    for integ in integrations:
        (sum_flux, f_var_rnoise, f_var_poisson,
            f_var_flat, background, b_var_rnoise, b_var_poisson,
            b_var_flat, npixels, flux_model) = extract_one_slit(
                data_model,
                integ,
                profile,
                bg_profile,
                extract_params
        )

        # Convert the sum to an average, for surface brightness.
        npixels_temp = np.where(npixels > 0., npixels, 1.)
        surf_bright = sum_flux / npixels_temp  # may be reset below
        sb_var_poisson = f_var_poisson / npixels_temp / npixels_temp
        sb_var_rnoise = f_var_rnoise / npixels_temp / npixels_temp
        sb_var_flat = f_var_flat / npixels_temp / npixels_temp
        background /= npixels_temp
        b_var_poisson = b_var_poisson / npixels_temp / npixels_temp
        b_var_rnoise = b_var_rnoise / npixels_temp / npixels_temp
        b_var_flat = b_var_flat / npixels_temp / npixels_temp

        del npixels_temp

        # Convert to flux density.
        # The input units will normally be MJy / sr, but for NIRSpec
        # point-source spectra the units will be MJy.
        input_units_are_megajanskys = (
            photom_has_been_run
            and source_type == 'POINT'
            and instrument == 'NIRSPEC'
        )

        if photom_has_been_run:
            # for NIRSpec point sources
            if input_units_are_megajanskys:
                flux = sum_flux * 1.e6  # MJy --> Jy
                f_var_poisson *= 1.e12  # MJy**2 --> Jy**2
                f_var_rnoise *= 1.e12  # MJy**2 --> Jy**2
                f_var_flat *= 1.e12  # MJy**2 --> Jy**2
                surf_bright[:] = 0.
                sb_var_poisson[:] = 0.
                sb_var_rnoise[:] = 0.
                sb_var_flat[:] = 0.
                background[:] /= pixel_solid_angle  # MJy / sr
                b_var_poisson = b_var_poisson / pixel_solid_angle / pixel_solid_angle
                b_var_rnoise = b_var_rnoise / pixel_solid_angle / pixel_solid_angle
                b_var_flat = b_var_flat / pixel_solid_angle / pixel_solid_angle
            else:
                flux = sum_flux * pixel_solid_angle * 1.e6  # MJy / steradian --> Jy
                f_var_poisson *= (pixel_solid_angle ** 2 * 1.e12)  # (MJy / sr)**2 --> Jy**2
                f_var_rnoise *= (pixel_solid_angle ** 2 * 1.e12)  # (MJy / sr)**2 --> Jy**2
                f_var_flat *= (pixel_solid_angle ** 2 * 1.e12)  # (MJy / sr)**2 --> Jy**2
        else:
            flux = sum_flux  # count rate

        del sum_flux

        error = np.sqrt(f_var_poisson + f_var_rnoise + f_var_flat)
        sb_error = np.sqrt(sb_var_poisson + sb_var_rnoise + sb_var_flat)
        berror = np.sqrt(b_var_poisson + b_var_rnoise + b_var_flat)

        # Set DQ from the flux value
        dq = np.zeros(flux.shape, dtype=np.uint32)
        dq[np.isnan(flux)] = datamodels.dqflags.pixel['DO_NOT_USE']

        # Make a table of the values, trimming to points with valid wavelengths only
        otab = np.array(
            list(
                zip(wavelength, flux[valid], error[valid],
                    f_var_poisson[valid], f_var_rnoise[valid], f_var_flat[valid],
                    surf_bright[valid], sb_error[valid], sb_var_poisson[valid],
                    sb_var_rnoise[valid], sb_var_flat[valid],
                    dq[valid], background[valid], berror[valid],
                    b_var_poisson[valid], b_var_rnoise[valid], b_var_flat[valid],
                    npixels[valid])
            ),
            dtype=datamodels.SpecModel().spec_table.dtype
        )

        spec = datamodels.SpecModel(spec_table=otab)
        spec.meta.wcs = spec_wcs.create_spectral_wcs(ra, dec, wavelength)
        spec.spec_table.columns['wavelength'].unit = 'um'
        spec.spec_table.columns['flux'].unit = flux_units
        spec.spec_table.columns['flux_error'].unit = flux_units
        spec.spec_table.columns['flux_var_poisson'].unit = f_var_units
        spec.spec_table.columns['flux_var_rnoise'].unit = f_var_units
        spec.spec_table.columns['flux_var_flat'].unit = f_var_units
        spec.spec_table.columns['surf_bright'].unit = sb_units
        spec.spec_table.columns['sb_error'].unit = sb_units
        spec.spec_table.columns['sb_var_poisson'].unit = sb_var_units
        spec.spec_table.columns['sb_var_rnoise'].unit = sb_var_units
        spec.spec_table.columns['sb_var_flat'].unit = sb_var_units
        spec.spec_table.columns['background'].unit = sb_units
        spec.spec_table.columns['bkgd_error'].unit = sb_units
        spec.spec_table.columns['bkgd_var_poisson'].unit = sb_var_units
        spec.spec_table.columns['bkgd_var_rnoise'].unit = sb_var_units
        spec.spec_table.columns['bkgd_var_flat'].unit = sb_var_units
        spec.slit_ra = ra
        spec.slit_dec = dec
        spec.spectral_order = sp_order
        spec.dispersion_direction = extract_params['dispaxis']

        # Record aperture limits as x/y start/stop values
        lower_limit, upper_limit, left_limit, right_limit = limits
        if spec.dispersion_direction == HORIZONTAL:
            spec.extraction_xstart = left_limit + 1
            spec.extraction_xstop = right_limit + 1
            spec.extraction_ystart = lower_limit + 1
            spec.extraction_ystop = upper_limit + 1
        else:
            spec.extraction_xstart = lower_limit + 1
            spec.extraction_xstop = upper_limit + 1
            spec.extraction_ystart = left_limit + 1
            spec.extraction_ystop = right_limit + 1

        copy_keyword_info(data_model, slitname, spec)

        if apcorr is not None:
            log.info('Applying Aperture correction.')
            if hasattr(apcorr, 'tabulated_correction'):
                if apcorr.tabulated_correction is not None:
                    apcorr_available = True

            # See whether we can reuse the previous aperture correction
            # object.  If so, just apply the pre-computed correction to
            # save a ton of time.
            if apcorr_available:
                # re-use the last aperture correction
                apcorr.apply(spec.spec_table, use_tabulated=True)
            else:
                # Attempt to tabulate the aperture correction for later use.
                # If this fails, fall back on the old method.
                try:
                    apcorr.tabulate_correction(spec.spec_table)
                    apcorr.apply(spec.spec_table, use_tabulated=True)
                    log.info("Tabulating aperture correction for use in multiple integrations.")
                except AttributeError:
                    log.info("Computing aperture correction.")
                    apcorr.apply(spec.spec_table)

        output_model.spec.append(spec)

        if log_increment > 0 and (integ + 1) % log_increment == 0:
            if integ == -1:
                pass
            elif integ == 0:
                if input_model.data.shape[0] == 1:
                    log.info("1 integration done")
                    progress_msg_printed = True
                else:
                    log.info("... 1 integration done")
            elif integ == input_model.data.shape[0] - 1:
                log.info(f"All {input_model.data.shape[0]} integrations done")
                progress_msg_printed = True
            else:
                log.info(f"... {integ + 1} integrations done")

    if not progress_msg_printed:
        if input_model.data.shape[0] == 1:
            log.info("1 integration done")
        else:
            log.info(f"All {input_model.data.shape[0]} integrations done")

    return output_model
