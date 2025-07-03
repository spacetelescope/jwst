import logging
import json
from json.decoder import JSONDecodeError
from pathlib import Path

from astropy.modeling import polynomial
import numpy as np
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.apcorr import (
    MirLrsApcorrModel,
    MirMrsApcorrModel,
    NrcWfssApcorrModel,
    NrsFsApcorrModel,
    NrsMosApcorrModel,
    NrsIfuApcorrModel,
    NisWfssApcorrModel,
)

from jwst.datamodels import ModelContainer
from jwst.datamodels.utils import attrs_to_group_id
from jwst.datamodels.utils.tso_multispec import make_tso_specmodel
from jwst.lib import pipe_utils
from jwst.lib.wcs_utils import get_wavelengths
from jwst.extract_1d import extract1d, spec_wcs
from jwst.extract_1d.apply_apcorr import select_apcorr
from jwst.extract_1d.psf_profile import psf_profile
from jwst.extract_1d.source_location import location_from_wcs

__all__ = [
    "run_extract1d",
    "read_extract1d_ref",
    "read_apcorr_ref",
    "get_extract_parameters",
    "box_profile",
    "aperture_center",
    "shift_by_offset",
    "define_aperture",
    "extract_one_slit",
    "create_extraction",
]


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

WFSS_EXPTYPES = ["NIS_WFSS", "NRC_WFSS", "NRC_GRISM"]
"""Exposure types to be regarded as wide-field slitless spectroscopy."""

SRCPOS_EXPTYPES = ["MIR_LRS-FIXEDSLIT", "NRS_FIXEDSLIT", "NRS_MSASPEC", "NRS_BRIGHTOBJ"]
"""Exposure types for which source position can be estimated."""

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

HORIZONTAL = 1
"""Horizontal dispersion axis."""
VERTICAL = 2
"""Vertical dispersion axis."""

# These values are assigned in get_extract_parameters, using key "match".
# If there was an aperture in the reference file for which the "id" key matched,
# that's (at least) a partial match.
# If "spectral_order" also matched, that's an exact match.
NO_MATCH = "no match"
PARTIAL = "partial match"
EXACT = "exact match"


class ContinueError(Exception):
    """Custom error to pass continue from a function inside a loop."""

    pass


def read_extract1d_ref(refname):
    """
    Read the extract1d reference file.

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
        if refname_type == "json":
            fd = Path(refname).open()
            try:
                ref_dict = json.load(fd)
                fd.close()
            except (UnicodeDecodeError, JSONDecodeError):
                # Input file does not load correctly as json file.
                # Probably an error in json file
                fd.close()
                log.error(
                    "Extract1d JSON reference file has an error. "
                    "Run a json validator off line and fix the file."
                )
                raise RuntimeError("Invalid JSON extract1d reference file.") from None
        else:
            log.error("Invalid Extract1d reference file: must be JSON.")
            raise RuntimeError("Invalid extract1d reference file: must be JSON.")

    return ref_dict


def read_apcorr_ref(refname, exptype):
    """
    Read the apcorr reference file.

    Determine the appropriate DataModel class to use for an APCORR reference file
    and read the file into it.

    Parameters
    ----------
    refname : str
        Path to the APCORR reference file.
    exptype : str
        EXP_TYPE of the input to the extract_1d step.

    Returns
    -------
    DataModel
        A datamodel containing the reference file input.
    """
    apcorr_model_map = {
        "MIR_LRS-FIXEDSLIT": MirLrsApcorrModel,
        "MIR_LRS-SLITLESS": MirLrsApcorrModel,
        "MIR_MRS": MirMrsApcorrModel,
        "NRC_GRISM": NrcWfssApcorrModel,
        "NRC_WFSS": NrcWfssApcorrModel,
        "NIS_WFSS": NisWfssApcorrModel,
        "NRS_BRIGHTOBJ": NrsFsApcorrModel,
        "NRS_FIXEDSLIT": NrsFsApcorrModel,
        "NRS_IFU": NrsIfuApcorrModel,
        "NRS_MSASPEC": NrsMosApcorrModel,
    }

    apcorr_model = apcorr_model_map[exptype]
    return apcorr_model(refname)


def get_extract_parameters(
    ref_dict,
    input_model,
    slitname,
    sp_order,
    meta,
    smoothing_length=None,
    bkg_fit=None,
    bkg_order=None,
    subtract_background=None,
    use_source_posn=None,
    position_offset=0.0,
    model_nod_pair=False,
    optimize_psf_location=False,
    extraction_type="box",
    psf_ref_name="N/A",
):
    """
    Get extraction parameter values.

    Parameters
    ----------
    ref_dict : dict or None
        For an extract1d reference file in JSON format, ``ref_dict`` will be
        the entire contents of the file.  If there is no extract1d reference
        file, ``ref_dict`` will be None.
    input_model : JWSTDataModel
        This can be either the input science file or one SlitModel out of
        a list of slits.
    slitname : str
        The name of the slit, or "ANY".
    sp_order : int
        The spectral order number.
    meta : ObjectNode
        The metadata for the actual input model, i.e. not just for the
        current slit, from input_model.meta.
    smoothing_length : int or None, optional
        Width of a boxcar function for smoothing the background regions.
        If None, the smoothing length will be retrieved from ``ref_dict``, or
        it will be set to 0 (no background smoothing) if this key is
        not found in ``ref_dict``.
        If ``smoothing_length`` is not None, that means that the user
        explicitly specified the value, so that value will be used.
        This argument is only used if background regions have been
        specified.
    bkg_fit : str or None, optional
        The type of fit to apply to background values in each
        column (or row, if the dispersion is vertical). The default
        'poly' results in a polynomial fit of order ``bkg_order``. Other
        options are 'mean' and 'median'. If 'mean' or 'median' is selected,
        ``bkg_order`` is ignored.
    bkg_order : int or None, optional
        Polynomial order for fitting to each column (or row, if the
        dispersion is vertical) of background.  If None, the polynomial
        order will be retrieved from ``ref_dict``, or it will be set to 0 if
        not found in ``ref_dict``.
        A value of 0 means that a simple average of the background
        regions, column by column (or row by row), will be used.
        If ``bkg_order`` is not None, that means that the user explicitly
        specified the value, so that value will be used.
        This argument must be positive or zero, and it is only used if
        background regions have been specified.
    subtract_background : bool or None, optional
        If False, all background parameters will be ignored.
    use_source_posn : bool or None, optional
        If True, the target and background positions specified in ``ref_dict``
        (or a default target position) will be shifted to account for
        the actual source location in the data.
        If None, a default value will be set, based on the exposure type.
    position_offset : float or None, optional
        Pixel offset to apply to the nominal source location.
        If None, the value specified in ``ref_dict`` will be used or it
        will default to 0.
    model_nod_pair : bool, optional
        If True, and if ``extraction_type`` is 'optimal', then a negative
        trace from nod subtraction will be modeled alongside the positive
        source, if possible.
    optimize_psf_location : bool
        If True, and if ``extraction_type`` is 'optimal', then the source
        location will be optimized, via iterative comparisons of the scene
        model with the input data.
    extraction_type : str, optional
        Extraction type ('box' or 'optimal').  Optimal extraction is
        only available if ``psf_ref_name`` is not 'N/A'.  If set to None,
        optimal extraction will be used if ``use_source_posn`` is True.
    psf_ref_name : str, optional
        The name of the PSF reference file, or "N/A".

    Returns
    -------
    extract_params : dict
        Information copied from ``ref_dict``.  The items will be selected
        based on ``slitname`` and ``sp_order``.  Default values will be
        assigned if ``ref_dict`` is None.
    """
    extract_params = {"match": NO_MATCH}  # initial value
    if ref_dict is None:
        # There is no extract1d reference file; use "reasonable" default values.
        extract_params["match"] = EXACT
        shape = input_model.data.shape
        extract_params["spectral_order"] = sp_order
        extract_params["xstart"] = 0  # first pixel in X
        extract_params["xstop"] = shape[-1] - 1  # last pixel in X
        extract_params["ystart"] = 0  # first pixel in Y
        extract_params["ystop"] = shape[-2] - 1  # last pixel in Y
        extract_params["extract_width"] = None
        extract_params["src_coeff"] = None
        extract_params["bkg_coeff"] = None  # no background sub.
        extract_params["smoothing_length"] = 0  # because no background sub.
        extract_params["bkg_fit"] = None  # because no background sub.
        extract_params["bkg_order"] = 0  # because no background sub.
        extract_params["subtract_background"] = False
        extract_params["extraction_type"] = "box"
        extract_params["use_source_posn"] = False  # no source position correction
        extract_params["position_offset"] = 0.0
        extract_params["model_nod_pair"] = False
        extract_params["optimize_psf_location"] = False
        extract_params["psf"] = "N/A"
        extract_params["independent_var"] = "pixel"
        extract_params["trace"] = None
        # Note that extract_params['dispaxis'] is not assigned.
        # This will be done later, possibly slit by slit.

    else:
        for aper in ref_dict["apertures"]:
            if (
                "id" in aper
                and aper["id"] != "dummy"
                and (aper["id"] == slitname or aper["id"] == ANY or slitname == ANY)
            ):
                # region_type is retained for backward compatibility; it is
                # not required to be present, but if it is and is not set
                # to target, then the aperture is not matched.
                region_type = aper.get("region_type", "target")
                if region_type != "target":
                    continue

                extract_params["match"] = PARTIAL

                # spectral_order is a secondary selection criterion.  The
                # default is the expected value, so if the key is not present
                # in the JSON file, the current aperture will be selected.
                # If the current aperture in the JSON file has
                # "spectral_order": "ANY", that aperture will be selected.
                spectral_order = aper.get("spectral_order", sp_order)

                if spectral_order == sp_order or spectral_order == ANY:
                    extract_params["match"] = EXACT
                    extract_params["spectral_order"] = sp_order
                    # Note: extract_params['dispaxis'] is not assigned.
                    # This is done later, possibly slit by slit.

                    # Set default start/stop by shape if not specified
                    shape = input_model.data.shape
                    extract_params["xstart"] = aper.get("xstart", 0)
                    extract_params["xstop"] = aper.get("xstop", shape[-1] - 1)
                    extract_params["ystart"] = aper.get("ystart", 0)
                    extract_params["ystop"] = aper.get("ystop", shape[-2] - 1)

                    extract_params["src_coeff"] = aper.get("src_coeff")
                    extract_params["bkg_coeff"] = aper.get("bkg_coeff")

                    if extract_params["bkg_coeff"] is not None and subtract_background is None:
                        subtract_background = True

                    if subtract_background:
                        extract_params["subtract_background"] = True
                        if bkg_fit is not None:
                            # Mean value for background fitting is equivalent
                            # to a polynomial fit with order 0.
                            if bkg_fit == "mean":
                                bkg_fit = "poly"
                                bkg_order = 0
                            extract_params["bkg_fit"] = bkg_fit
                        else:
                            extract_params["bkg_fit"] = aper.get("bkg_fit", "poly")
                    else:
                        extract_params["subtract_background"] = False
                        extract_params["bkg_fit"] = None

                    extract_params["independent_var"] = aper.get("independent_var", "pixel").lower()

                    if bkg_order is None:
                        extract_params["bkg_order"] = aper.get("bkg_order", 0)
                    else:
                        # If the user supplied a value, use that value.
                        extract_params["bkg_order"] = bkg_order

                    # Set use_source_posn based on hierarchy of priorities:
                    # parameter value on the command line is highest precedence,
                    # then parameter value from the extract1d reference file,
                    # and finally a default setting based on exposure type.

                    # value from the ref file
                    use_source_posn_aper = aper.get("use_source_posn", None)
                    if use_source_posn is None:
                        # no value set on command line
                        if use_source_posn_aper is None:
                            # no value set in ref file - use a suitable default
                            if meta.exposure.type in SRCPOS_EXPTYPES:
                                use_source_posn = True
                                log.info(
                                    f"Turning on source position correction "
                                    f"for exp_type = {meta.exposure.type}"
                                )
                            else:
                                use_source_posn = False
                        else:
                            # use the value from the ref file
                            use_source_posn = use_source_posn_aper
                    extract_params["use_source_posn"] = use_source_posn
                    extract_params["position_offset"] = position_offset
                    extract_params["trace"] = None

                    extract_params["extract_width"] = aper.get("extract_width")

                    if smoothing_length is None:
                        extract_params["smoothing_length"] = aper.get("smoothing_length", 0)
                    else:
                        # If the user supplied a value, use that value.
                        extract_params["smoothing_length"] = smoothing_length

                    # Check that the smoothing length has a reasonable value
                    sm_length = extract_params["smoothing_length"]
                    if sm_length != int(sm_length):
                        sm_length = int(np.round(sm_length))
                        log.warning(f"Smoothing length must be an integer. Rounding to {sm_length}")
                    if sm_length != 0:
                        if sm_length < 3:
                            log.warning(
                                f"Smoothing length {sm_length} is not allowed. Setting it to 0."
                            )
                            sm_length = 0
                        elif sm_length % 2 == 0:
                            sm_length -= 1
                            log.warning(
                                "Even smoothing lengths are not supported. "
                                "Rounding the smoothing length down to "
                                f"{sm_length}."
                            )
                    extract_params["smoothing_length"] = sm_length

                    extract_params["psf"] = psf_ref_name
                    extract_params["optimize_psf_location"] = optimize_psf_location
                    extract_params["model_nod_pair"] = model_nod_pair

                    # Check for a valid PSF file
                    extraction_type = str(extraction_type).lower()
                    if extract_params["psf"] == "N/A" and extraction_type != "box":
                        log.warning("No PSF file available. Setting extraction type to 'box'.")
                        extraction_type = "box"

                    # Set the extraction type to 'box' or 'optimal'
                    if extraction_type == "none":
                        if extract_params["use_source_posn"]:
                            extract_params["extraction_type"] = "optimal"
                        else:
                            extract_params["extraction_type"] = "box"
                        log.info(
                            f"Using extraction type '{extract_params['extraction_type']}' "
                            f"for use_source_posn = {extract_params['use_source_posn']}"
                        )
                    else:
                        extract_params["extraction_type"] = extraction_type

                    break

    return extract_params


def log_initial_parameters(extract_params):
    """
    Log the initial extraction parameters.

    Parameters
    ----------
    extract_params : dict
        Information read from the reference file.
    """
    if "dispaxis" not in extract_params:
        return

    log.debug("Extraction parameters:")
    skip_keys = {"match", "trace"}
    for key, value in extract_params.items():
        if key not in skip_keys:
            log.debug(f"  {key} = {value}")


def create_poly(coeff):
    """
    Create a polynomial model from coefficients.

    Parameters
    ----------
    coeff : list of float
        The coefficients of the polynomial, constant term first, highest
        order term last.

    Returns
    -------
    `astropy.modeling.polynomial.Polynomial1D` or None
        None is returned if ``coeff`` is empty.
    """
    n = len(coeff)
    if n < 1:
        return None

    coeff_dict = {f"c{i}": coeff[i] for i in range(n)}

    return polynomial.Polynomial1D(degree=n - 1, **coeff_dict)


def populate_time_keywords(input_model, output_model):
    """
    Copy the integration time keywords to header keywords.

    Parameters
    ----------
    input_model : TSOMultiSpecModel
        The input science model.
    output_model : TSOMultiSpecModel
        The output science model.  This may be modified in-place.
    """
    nints = input_model.meta.exposure.nints
    int_start = input_model.meta.exposure.integration_start

    if hasattr(input_model, "data"):
        shape = input_model.data.shape
        if len(shape) == 2:
            num_integ = 1
        else:
            # len(shape) == 3
            num_integ = shape[0]
    else:
        # e.g. MultiSlit data
        num_integ = 1

    # This assumes that the spec attribute of output_model has already been created,
    # and spectra have been appended to the spec_table
    n_output_spec = 0
    for spec in output_model.spec:
        n_output_spec += len(spec.spec_table)

    # num_j is the number of spectral tables, e.g. the number
    # of fixed-slit spectra, MSA spectra, or different
    # spectral orders; num_integ is the number of integrations.
    # The total number of output spectra is n_output_spec = num_integ * num_j
    num_j = len(output_model.spec)

    # check that the number of spectra present matches expectations
    if n_output_spec != num_j * num_integ:
        log.warning(
            f"populate_time_keywords:  Don't understand n_output_spec = {n_output_spec}, "
            f"num_j = {num_j}, num_integ = {num_integ}"
        )
    else:
        log.debug(
            f"Number of output spectra = {n_output_spec}; "
            f"number of spectra for each integration = {num_j}; "
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

    if hasattr(input_model, "int_times") and input_model.int_times is not None:
        nrows = len(input_model.int_times)
    else:
        nrows = 0

    if nrows < 1:
        log.warning(
            "There is no INT_TIMES table in the input file - "
            "Making best guess on integration numbers."
        )

        # Set a simple integration index
        for spec in output_model.spec:
            spec.spec_table["INT_NUM"] = np.arange(1, len(spec.spec_table) + 1)
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
            log.warning(
                "Not using INT_TIMES table because the data have been averaged over integrations."
            )
            skip = True
    elif isinstance(input_model, (datamodels.CubeModel, datamodels.SlitModel)):
        shape = input_model.data.shape

        if len(shape) == 2 and num_integrations > 1:
            log.warning(
                "Not using INT_TIMES table because the data have been averaged over integrations."
            )
            skip = True
        elif len(shape) != 3 or shape[0] > nrows:
            # Later, we'll check that the integration_number column actually
            # has a row corresponding to every integration in the input.
            log.warning(
                "Not using INT_TIMES table because the data shape is "
                "not consistent with the number of table rows."
            )
            skip = True
    elif isinstance(input_model, datamodels.IFUCubeModel):
        log.warning("The INT_TIMES table will be ignored for IFU data.")
        skip = True

    if skip:
        return

    int_num = input_model.int_times["integration_number"]
    start_time_mjd = input_model.int_times["int_start_MJD_UTC"]
    mid_time_mjd = input_model.int_times["int_mid_MJD_UTC"]
    end_time_mjd = input_model.int_times["int_end_MJD_UTC"]
    start_tdb = input_model.int_times["int_start_BJD_TDB"]
    mid_tdb = input_model.int_times["int_mid_BJD_TDB"]
    end_tdb = input_model.int_times["int_end_BJD_TDB"]

    # Inclusive range of integration numbers in the input data, zero indexed.
    data_range = (int_start, int_end)

    # Inclusive range of integration numbers in the INT_TIMES table, zero indexed.
    table_range = (int_num[0] - 1, int_num[-1] - 1)
    offset = data_range[0] - table_range[0]
    data_range_test = data_range[0] < table_range[0] or data_range[1] > table_range[1]
    if data_range_test or offset > table_range[1]:
        log.warning(
            "Not using the INT_TIMES table because it does not include "
            "rows for all integrations in the data."
        )
        return

    log.debug("TSO data, so copying times from the INT_TIMES table.")

    for j in range(num_j):  # for each slit or order (EXTRACT1D extension)
        spec = output_model.spec[j]
        spec.time_scale = "UTC"
        for k in range(len(spec.spec_table)):  # for each integration (table row)
            row = k + offset
            spec.spec_table["INT_NUM"][k] = int_num[row]
            spec.spec_table["MJD-BEG"][k] = start_time_mjd[row]
            spec.spec_table["MJD-AVG"][k] = mid_time_mjd[row]
            spec.spec_table["MJD-END"][k] = end_time_mjd[row]
            spec.spec_table["TDB-BEG"][k] = start_tdb[row]
            spec.spec_table["TDB-MID"][k] = mid_tdb[row]
            spec.spec_table["TDB-END"][k] = end_tdb[row]


def get_spectral_order(slit):
    """
    Get the spectral order number.

    Parameters
    ----------
    slit : SlitModel
        One slit from an input MultiSlitModel or similar.

    Returns
    -------
    int
        Spectral order number for ``slit``.  If no information about spectral
        order is available in ``wcsinfo``, a default value of 1 will be
        returned.
    """
    sp_order = slit.meta.wcsinfo.spectral_order

    if sp_order is None:
        log.warning("spectral_order is None; using 1")
        sp_order = 1

    return sp_order


def is_prism(input_model):
    """
    Determine whether the current observing mode used a prism.

    The reason for this test is so we can skip spectral extraction if the
    spectral order is zero and the exposure was not made using a prism.
    In this context, therefore, a grism is not considered to be a prism.

    Parameters
    ----------
    input_model : JWSTDataModel
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

    if (instrument == "MIRI" and instrument_filter.find("P750L") >= 0) or (
        instrument == "NIRSPEC" and grating.find("PRISM") >= 0
    ):
        prism_mode = True

    return prism_mode


def copy_keyword_info(slit, slitname, spec):
    """
    Copy metadata from the input to the output spectrum.

    Parameters
    ----------
    slit : A SlitModel object
        Metadata will be copied from the input ``slit`` to output ``spec``.
    slitname : str or None
        The name of the slit.
    spec : One element of MultiSpecModel.spec
        Metadata attributes will be updated in-place.
    """
    if slitname is not None and slitname != "ANY":
        spec.name = slitname

    # Copy over some attributes if present, even if they are None
    copy_attributes = [
        "slitlet_id",
        "source_id",
        "source_xpos",
        "source_ypos",
        "source_ra",
        "source_dec",
        "shutter_state",
    ]
    for key in copy_attributes:
        if hasattr(slit, key):
            setattr(spec, key, getattr(slit, key))

    # Copy over some attributes only if they are present and not None
    copy_populated_attributes = [
        "source_name",
        "source_alias",
        "source_type",
        "stellarity",
        "quadrant",
        "slit_xscale",
        "slit_yscale",
        "wavelength_corrected",
        "pathloss_correction_type",
        "barshadow_corrected",
    ]
    for key in copy_populated_attributes:
        if getattr(slit, key, None) is not None:
            setattr(spec, key, getattr(slit, key))


def _set_weight_from_limits(profile, idx, lower_limit, upper_limit, allow_partial=True):
    """
    Set profile pixel weighting from a lower and upper limit.

    Pixels to be fully included in the aperture are set to 1.0. Pixels partially
    included are set to fractional values.

    The profile is updated in place, so aperture setting is cumulative.  If
    there are overlapping apertures specified, later ones will overwrite
    earlier values.

    Parameters
    ----------
    profile : ndarray of float
        The spatial profile to update.
    idx : ndarray of int
        Index values for the profile array, corresponding to the cross-dispersion
        axis. Dimensions must match ``profile`` shape.
    lower_limit : float or ndarray of float
        Lower limit for the aperture. If not a single value, dimensions must
        match ``profile`` shape.
    upper_limit : float or ndarray of float
        Upper limit for the aperture. If not a single value, dimensions must
        match ``profile`` shape.
    allow_partial : bool, optional
        If True, partial pixel weights are set where the pixel index intersects
        the limit values.  If False, only whole integer weights are set.
    """
    # Both limits are inclusive
    profile[(idx >= lower_limit) & (idx <= upper_limit)] = 1.0

    if allow_partial:
        # For each pixel, get the distance from the lower and upper limits,
        # to set partial pixel weights to the fraction of the pixel included.
        # For example, if the lower limit is 1.7, then the weight of pixel 1
        # should be set to 0.3. If the upper limit is 10.7, then the weight of
        # pixel 11 should be set to 0.7.
        lower_weight = idx - lower_limit + 1
        upper_weight = upper_limit - idx + 1

        for partial_pixel_weight in [lower_weight, upper_weight]:
            # Check for overlap values that are between 0 and 1, for which
            # the profile does not already contain a higher weight
            test = (
                (partial_pixel_weight > 0)
                & (partial_pixel_weight < 1)
                & (profile < partial_pixel_weight)
            )

            # Set these values to the partial pixel weight
            profile[test] = partial_pixel_weight[test]


def box_profile(
    shape, extract_params, wl_array, coefficients="src_coeff", label="aperture", return_limits=False
):
    """
    Create a spatial profile for box extraction.

    The output profile is an image matching the input data shape,
    containing weights for each pixel in the image.  Pixels to be fully
    included in an extraction aperture are set to 1.0.  Pixels partially
    included are set to values between 0 and 1, depending on the overlap
    with the aperture limits.  Pixels not included in the aperture are
    set to 0.0.

    Upper and lower limits for the aperture are determined from the
    ``extract_params``, in this priority order:

       1. src_coeff upper and lower limits (or bkg_coeff, for a background profile)
       2. trace +/- extraction width / 2
       3. center of start/stop values +/- extraction width / 2
       4. cross-dispersion start/stop values
       5. array limits.

    Left and right limits are set from start/stop values only.

    Only a single aperture is supported at this time.  The 'src_coeff'
    and 'bkg_coeff' specifications allow multiple regions to be specified,
    but all regions are considered part of the same aperture and are
    added to the same profile, to be extracted at once.

    Parameters
    ----------
    shape : tuple of int
        Data shape for the output profile, to match the spectral image.
    extract_params : dict
        Extraction parameters, as returned from `get_extract_parameters`.
    wl_array : ndarray
        Array of wavelength values, matching ``shape``, for each pixel in
        the array.
    coefficients : {'src_coeff', 'bkg_coeff'}, optional
        The polynomial coefficients to look for in the ``extract_params``
        dictionary. If 'bkg_coeff', the output aperture contains background
        regions; otherwise, it contains target source regions.
    label : str, optional
        A label to use for the aperture, while logging limit values.
    return_limits : bool, optional
        If True, an upper and lower limit value for the aperture are
        returned along with the spatial profile.  These are used for
        recording the aperture extent in output metadata.  For
        apertures set from polynomial coefficients, the returned values
        are averages of the true upper and lower limits, which may vary
        by dispersion element.

    Returns
    -------
    profile : ndarray of float
        Aperture weights to use in box extraction from a spectral image.
    lower_limit : float, optional
        Average lower limit for the aperture. Returned only if ``return_limits``
        is set.
    upper_limit : float, optional
        Average upper limit for the aperture. Returned only if ``return_limits``
        is set.
    """
    # Get pixel index values for the array
    yidx, xidx = np.mgrid[: shape[0], : shape[1]]
    if extract_params["dispaxis"] == HORIZONTAL:
        dval = yidx
    else:
        dval = xidx

    # Get start/stop values from parameters if present,
    # or default to data shape
    xstart = extract_params.get("xstart", 0)
    xstop = extract_params.get("xstop", shape[1] - 1)
    ystart = extract_params.get("ystart", 0)
    ystop = extract_params.get("ystop", shape[0] - 1)

    # Check if the profile should contain partial pixel weights
    if coefficients == "bkg_coeff" and extract_params["bkg_fit"] == "median":
        allow_partial = False
    else:
        allow_partial = True

    # Set aperture region, in this priority order:
    # 1. src_coeff upper and lower limits (or bkg_coeff, for background profile)
    # 2. trace +/- extraction width / 2
    # 3. center of start/stop values +/- extraction width / 2
    # 4. start/stop values
    profile = np.full(shape, 0.0)
    if extract_params[coefficients] is not None:
        # Limits from source coefficients: ignore ystart/stop/width
        if extract_params["independent_var"].startswith("wavelength"):
            ival = wl_array
        elif extract_params["dispaxis"] == HORIZONTAL:
            ival = xidx.astype(np.float32)
        else:
            ival = yidx.astype(np.float32)

        # The source extraction can include more than one region,
        # but must contain pairs of lower and upper limits.
        n_src_coeff = len(extract_params[coefficients])
        if n_src_coeff % 2 != 0:
            raise RuntimeError(
                f"{coefficients} must contain alternating lists of lower and upper limits."
            )

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

                _set_weight_from_limits(
                    profile,
                    dval,
                    lower_limit_region,
                    upper_limit_region,
                    allow_partial=allow_partial,
                )
                mean_lower = np.mean(lower_limit_region)
                mean_upper = np.mean(upper_limit_region)
                log.info(
                    f"Mean {label} start/stop from {coefficients}: "
                    f"{mean_lower:.2f} -> {mean_upper:.2f} (inclusive)"
                )

                if lower_limit is None:
                    lower_limit = mean_lower
                    upper_limit = mean_upper
                else:
                    if mean_lower < lower_limit:
                        lower_limit = mean_lower
                    if mean_upper > upper_limit:
                        upper_limit = mean_upper

    elif extract_params["extract_width"] is not None and extract_params["trace"] is not None:
        width = extract_params["extract_width"]
        trace = extract_params["trace"]

        if extract_params["dispaxis"] != HORIZONTAL:
            trace = np.tile(trace, (shape[1], 1)).T

        lower_limit_region = trace - (width - 1.0) / 2.0
        upper_limit_region = lower_limit_region + width - 1

        _set_weight_from_limits(profile, dval, lower_limit_region, upper_limit_region)

        lower_limit = np.nanmean(lower_limit_region)
        upper_limit = np.nanmean(upper_limit_region)
        log.info(
            f"Mean {label} start/stop from trace: "
            f"{lower_limit:.2f} -> {upper_limit:.2f} (inclusive)"
        )

    elif extract_params["extract_width"] is not None:
        # Limits from extraction width at center of ystart/stop if present,
        # center of array if not
        if extract_params["dispaxis"] == HORIZONTAL:
            nominal_middle = (ystart + ystop) / 2.0
        else:
            nominal_middle = (xstart + xstop) / 2.0

        width = extract_params["extract_width"]
        lower_limit = nominal_middle - (width - 1.0) / 2.0
        upper_limit = lower_limit + width - 1

        _set_weight_from_limits(profile, dval, lower_limit, upper_limit)
        log.info(
            f"{label.capitalize()} start/stop: {lower_limit:.2f} -> {upper_limit:.2f} (inclusive)"
        )

    else:
        # Limits from start/stop only, defaulting to the full array
        # if not specified
        if extract_params["dispaxis"] == HORIZONTAL:
            lower_limit = ystart
            upper_limit = ystop
        else:
            lower_limit = xstart
            upper_limit = xstop

        _set_weight_from_limits(profile, dval, lower_limit, upper_limit)
        log.info(
            f"{label.capitalize()} start/stop: {lower_limit:.2f} -> {upper_limit:.2f} (inclusive)"
        )

    # Set weights to zero outside left and right limits
    if extract_params["dispaxis"] == HORIZONTAL:
        profile[:, : int(np.ceil(xstart))] = 0
        profile[:, int(np.floor(xstop)) + 1 :] = 0
    else:
        profile[: int(np.ceil(ystart)), :] = 0
        profile[int(np.floor(ystop)) + 1 :, :] = 0

    if return_limits:
        return profile, lower_limit, upper_limit
    return profile


def aperture_center(profile, dispaxis=1, middle_pix=None):
    """
    Determine the nominal center of an aperture.

    The center is determined from a weighted average of the pixel
    coordinates, where the weights are set by the profile image.

    If ``middle_pix`` is specified, it is expected to be the
    dispersion element at which the cross-dispersion center
    should be determined.  In this case, the profile must contain
    some non-zero elements at that dispersion coordinate.

    If ``dispaxis`` is 1 (the default), the return values are in (y,x)
    order. Otherwise, the return values are in (x,y) order.

    Parameters
    ----------
    profile : ndarray
        Pixel weights for the aperture.
    dispaxis : int, optional
        If 1, the dispersion axis is horizontal.  Otherwise,
        the dispersion axis is vertical.
    middle_pix : int or None
        Index value for the center of the slit along the
        dispersion axis.  If specified, it is returned as
        ``spec_center``.

    Returns
    -------
    slit_center : float
        Index value for the center of the aperture along the
        cross-dispersion axis.
    spec_center : float
        Index value for the center of the aperture along the
        dispersion axis.
    """
    weights = profile.copy()
    weights[weights <= 0] = 0.0

    yidx, xidx = np.mgrid[: profile.shape[0], : profile.shape[1]]
    if dispaxis == HORIZONTAL:
        spec_idx = xidx
        if middle_pix is not None:
            slit_idx = yidx[:, middle_pix]
            weights = weights[:, middle_pix]
        else:
            slit_idx = yidx
    else:
        spec_idx = yidx
        if middle_pix is not None:
            slit_idx = xidx[middle_pix, :]
            weights = weights[middle_pix, :]
        else:
            slit_idx = xidx

    if np.sum(weights) == 0:
        weights[:] = 1

    slit_center = np.average(slit_idx, weights=weights)
    if middle_pix is not None:
        spec_center = middle_pix
    else:
        spec_center = np.average(spec_idx, weights=weights)

    # if dispaxis == 1 (default), this returns center_x, center_y
    return slit_center, spec_center


def shift_by_offset(offset, extract_params, update_trace=True):
    """
    Shift the nominal extraction parameters by a pixel offset.

    Start, stop, and polynomial coefficient values for source and
    background are updated in place in the ``extract_params`` dictionary.
    The source trace value, if present, is also updated if desired.

    Parameters
    ----------
    offset : float
        Cross-dispersion offset to apply, in pixels.
    extract_params : dict
        Extraction parameters to update, as created by
        `get_extraction_parameters`.
    update_trace : bool
        If True, the trace in ``extract_params['trace']`` is also updated
        if present.
    """
    # Shift polynomial coefficients
    coeff_params = ["src_coeff", "bkg_coeff"]
    for params in coeff_params:
        if extract_params[params] is not None:
            for coeff_list in extract_params[params]:
                coeff_list[0] += offset

    # Shift start/stop values
    if extract_params["dispaxis"] == HORIZONTAL:
        start_stop_params = ["ystart", "ystop"]
    else:
        start_stop_params = ["xstart", "xstop"]
    for params in start_stop_params:
        if extract_params[params] is not None:
            extract_params[params] += offset

    # Shift the full trace
    if update_trace and extract_params["trace"] is not None:
        extract_params["trace"] += offset


def define_aperture(input_model, slit, extract_params, exp_type):
    """
    Define an extraction aperture from input parameters.

    Parameters
    ----------
    input_model : DataModel
        The input science model containing metadata information.
    slit : DataModel or None
        One slit from a MultiSlitModel (or similar), or None.
        The spectral image and WCS information will be retrieved from ``slit``
        unless ``slit`` is None. In that case, they will be retrieved
        from ``input_model``.
    extract_params : dict
        Extraction parameters, as created by `get_extraction_parameters`.
    exp_type : str
        Exposure type for the input data.

    Returns
    -------
    ra : float
        A representative RA value for the source centered in the
        aperture.
    dec : float
        A representative Dec value for the source centered in the
        aperture.
    wavelength : ndarray of float
        The 1D wavelength array, matching the dispersion dimension of
        the input spectral image, corresponding to the aperture. May
        contain NaNs for invalid dispersion elements.
    profile : ndarray of float
        A 2D image containing pixel weights for the extraction aperture,
        matching the dimensions of the input spectral image.  Values
        are between 0.0 (pixel not included in the extraction aperture)
        and 1.0 (pixel fully included in the aperture).
    bg_profile : ndarray of float or None
        If background regions are specified in ``extract_params['bkg_coeff']``,
        and ``extract_params['subtract_background']`` is True, then
        ``bg_profile`` is a 2D image containing pixel weights for background
        regions, to be fit during extraction.  Otherwise, ``bg_profile`` is
        None.
    nod_profile : ndarray of float or None
        For optimal extraction, if nod subtraction was performed, a
        second spatial profile is generated, modeling the negative source
        in the slit. This second spatial profile is returned in ``nod_profile``
        if generated.  Otherwise, ``nod_profile`` is None.
    limits : tuple of float
        Index limit values for the aperture, returned as (lower_limit, upper_limit,
        left_limit, right_limit).  Upper/lower limits are along the
        cross-dispersion axis.  Left/right limits are along the dispersion axis.
        All limits are inclusive and start at zero index value.
    """
    if slit is None:
        data_model = input_model
    else:
        data_model = slit
    data_shape = data_model.data.shape[-2:]

    # Get a wavelength array for the data
    wl_array = get_wavelengths(data_model, exp_type, extract_params["spectral_order"])

    # Shift aperture definitions by source position if needed
    # Extract parameters are updated in place
    if extract_params["use_source_posn"]:
        # Source location from WCS
        middle_pix, middle_wl, location, trace = location_from_wcs(input_model, slit)
        if location is not None:
            log.info(
                f"Computed source location is {location:.2f}, "
                f"at pixel {middle_pix}, wavelength {middle_wl:.2f}"
            )

            # Nominal location from extract params + located middle
            nominal_profile = box_profile(
                data_shape, extract_params, wl_array, label="nominal aperture"
            )
            nominal_location, _ = aperture_center(
                nominal_profile, extract_params["dispaxis"], middle_pix=middle_pix
            )

            # Offset extract parameters by location - nominal
            offset = location - nominal_location
            log.info(
                f"Nominal location is {nominal_location:.2f}, so offset is {offset:.2f} pixels"
            )
            shift_by_offset(offset, extract_params, update_trace=False)
    else:
        middle_pix, middle_wl, location, trace = None, None, None, None

    # Store the trace, if computed
    extract_params["trace"] = trace

    # Add an extra position offset if desired, from extract_params['position_offset']
    offset = extract_params.get("position_offset", 0.0)
    if offset != 0.0:
        log.info(f"Applying additional cross-dispersion offset {offset:.2f} pixels")
        shift_by_offset(offset, extract_params, update_trace=True)

    # Make a spatial profile, including source shifts if necessary
    nod_profile = None
    if extract_params["extraction_type"] == "optimal":
        profiles, lower_limit, upper_limit = psf_profile(
            data_model,
            extract_params["trace"],
            wl_array,
            extract_params["psf"],
            model_nod_pair=extract_params["model_nod_pair"],
            optimize_shifts=extract_params["optimize_psf_location"],
        )
        if len(profiles) > 1:
            profile, nod_profile = profiles
        else:
            profile = profiles[0]
    else:
        profile, lower_limit, upper_limit = box_profile(
            data_shape, extract_params, wl_array, return_limits=True
        )

    # Make sure profile weights are zero where wavelengths are invalid
    profile[~np.isfinite(wl_array)] = 0.0
    if nod_profile is not None:
        nod_profile[~np.isfinite(wl_array)] = 0.0

    # Get the effective left and right limits from the profile weights
    nonzero_weight = np.where(np.sum(profile, axis=extract_params["dispaxis"] - 1) > 0)[0]
    if len(nonzero_weight) > 0:
        left_limit = nonzero_weight[0]
        right_limit = nonzero_weight[-1]
    else:
        left_limit = None
        right_limit = None

    # Make a background profile if necessary
    # (will also include source shifts)
    if extract_params["subtract_background"] and extract_params["bkg_coeff"] is not None:
        bg_profile = box_profile(data_shape, extract_params, wl_array, coefficients="bkg_coeff")
    else:
        bg_profile = None

    # Get 1D wavelength corresponding to the spatial profile
    mask = np.isnan(wl_array) | (profile <= 0)
    masked_wl = np.ma.masked_array(wl_array, mask=mask)
    masked_weights = np.ma.masked_array(profile, mask=mask)
    if extract_params["dispaxis"] == HORIZONTAL:
        wavelength = np.average(masked_wl, weights=masked_weights, axis=0).filled(np.nan)
    else:
        wavelength = np.average(masked_wl, weights=masked_weights, axis=1).filled(np.nan)

    # Get RA and Dec corresponding to the center of the array,
    # weighted by the spatial profile
    center_y, center_x = aperture_center(profile, 1)
    coords = data_model.meta.wcs(center_x, center_y)
    if np.any(np.isnan(coords)):
        ra = None
        dec = None
    else:
        ra = float(coords[0])
        dec = float(coords[1])

    # Return limits as a tuple with 4 elements: lower, upper, left, right
    limits = (lower_limit, upper_limit, left_limit, right_limit)

    return ra, dec, wavelength, profile, bg_profile, nod_profile, limits


def extract_one_slit(data_model, integration, profile, bg_profile, nod_profile, extract_params):
    """
    Extract data for one slit, or spectral order, or integration.

    Parameters
    ----------
    data_model : JWSTDataModel
        The input science model. May be a single slit from a MultiSlitModel
        (or similar), or a single data type, like an ImageModel, SlitModel,
        or CubeModel.
    integration : int
        For the case that data_model is a SlitModel or a CubeModel,
        ``integration`` is the integration number.  If the integration number is
        not relevant (i.e. the data array is 2-D), ``integration`` should be -1.
    profile : ndarray of float
        Spatial profile indicating the aperture location. Must be a
        2D image matching the input, with floating point values between 0
        and 1 assigning a weight to each pixel.  0 means the pixel is not used,
        1 means the pixel is fully included in the aperture.
    bg_profile : ndarray of float or None
        Background profile indicating any background regions to use, following
        the same format as the spatial profile. Ignored if
        extract_params['subtract_background'] is False.
    nod_profile : ndarray of float or None
        For optimal extraction, if nod subtraction was performed, a
        second spatial profile is generated, modeling the negative source
        in the slit. This second spatial profile may be passed in ``nod_profile``
        for simultaneous fitting with the primary source in ``profile``.
        Otherwise, ``nod_profile`` should be None.
    extract_params : dict
        Parameters read from the extract1d reference file, as returned by
        `get_extract_parameters`.

    Returns
    -------
    sum_flux : ndarray, 1-D, float64
        The sum of the data values in the extraction region minus the sum
        of the data values in the background regions (scaled by the ratio
        of the numbers of pixels), for each pixel.
        The data values are usually in units of surface brightness,
        so this value isn't the flux, it's an intermediate value.
        Multiply ``sum_flux`` by the solid angle of a pixel to get the flux for a
        point source (column "flux").  Divide ``sum_flux`` by ``npixels`` (to
        compute the average) to get the array for the "surf_bright"
        (surface brightness) output column.
    f_var_rnoise : ndarray, 1-D
        The extracted read noise variance values to go along with the
        sum_flux array.
    f_var_poisson : ndarray, 1-D
        The extracted poisson variance values to go along with the
        sum_flux array.
    f_var_flat : ndarray, 1-D
        The extracted flat field variance values to go along with the
        sum_flux array.
    background : ndarray, 1-D
        The background count rate that was subtracted from the sum of
        the source data values to get ``sum_flux``.
    b_var_rnoise : ndarray, 1-D
        The extracted read noise variance values to go along with the
        background array.
    b_var_poisson : ndarray, 1-D
        The extracted poisson variance values to go along with the
        background array.
    b_var_flat : ndarray, 1-D
        The extracted flat field variance values to go along with the
        background array.
    npixels : ndarray, 1-D, float64
        The number of pixels that were added together to get ``sum_flux``,
        including any fractional pixels included via non-integer weights
        in the input profile.
    scene_model : ndarray, 2-D, float64
        A 2D model of the flux in the spectral image, corresponding to
        the extracted aperture.
    residual : ndarray, 2-D, float64
        Residual image from the input minus the scene model.
    """
    # Get the data and variance arrays
    if integration > -1:
        log.debug(f"Extracting integration {integration + 1}")
        data = data_model.data[integration]
        var_rnoise = data_model.var_rnoise[integration]
        var_poisson = data_model.var_poisson[integration]
        var_flat = data_model.var_flat[integration]
    else:
        data = data_model.data
        var_rnoise = data_model.var_rnoise
        var_poisson = data_model.var_poisson
        var_flat = data_model.var_flat

    # Make sure variances match data
    if var_rnoise.shape != data.shape:
        var_rnoise = np.zeros_like(data)
    if var_poisson.shape != data.shape:
        var_poisson = np.zeros_like(data)
    if var_flat.shape != data.shape:
        var_flat = np.zeros_like(data)

    # Transpose data for extraction
    if extract_params["dispaxis"] == HORIZONTAL:
        bg_profile_view = bg_profile
        if nod_profile is not None:
            profiles = [profile, nod_profile]
        else:
            profiles = [profile]
    else:
        data = data.T
        var_rnoise = var_rnoise.T
        var_poisson = var_poisson.T
        var_flat = var_flat.T
        if bg_profile is not None:
            bg_profile_view = bg_profile.T
        else:
            bg_profile_view = None
        if nod_profile is not None:
            profiles = [profile.T, nod_profile.T]
        else:
            profiles = [profile.T]

    # Extract spectra from the data
    result = extract1d.extract1d(
        data,
        profiles,
        var_rnoise,
        var_poisson,
        var_flat,
        profile_bg=bg_profile_view,
        bg_smooth_length=extract_params["smoothing_length"],
        fit_bkg=extract_params["subtract_background"],
        bkg_fit_type=extract_params["bkg_fit"],
        bkg_order=extract_params["bkg_order"],
        extraction_type=extract_params["extraction_type"],
    )

    # Extraction routine can return multiple spectra;
    # here, we just want the first result
    first_result = []
    for r in result[:-1]:
        first_result.append(r[0])

    # The last return value is the 2D model - there is only one, regardless
    # of the number of input profiles. It may need to be transposed to match
    # the input data.
    scene_model = result[-1]
    residual = data - scene_model
    if extract_params["dispaxis"] == HORIZONTAL:
        first_result.append(scene_model)
        first_result.append(residual)
    else:
        first_result.append(scene_model.T)
        first_result.append(residual.T)
    return first_result


def create_extraction(
    input_model,
    slit,
    output_model,
    extract_ref_dict,
    slitname,
    sp_order,
    exp_type,
    apcorr_ref_model=None,
    log_increment=50,
    save_profile=False,
    save_scene_model=False,
    save_residual_image=False,
    **kwargs,
):
    """
    Extract spectra from an input model and append to an output model.

    Input data, specified in the ``slit`` or ``input_model``, should contain data
    with consistent wavelengths, target positions, and spectral image cutouts,
    suitable for extraction with a shared aperture.  Extraction parameters
    are determined collectively, then multiple integrations, if present,
    are each extracted separately.

    The output model must be a `MultiSpecModel` or `TSOMultiSpecModel`,
    created before calling this function, and passed as ``output_model``.
    It is updated in place, with new spectral tables appended as they are created.

    The process is:

       1. Retrieve extraction parameters from ``extract_ref_dict`` and
          set defaults for any missing values as needed.
       2. Define an extraction aperture from the input parameters, as
          well as the nominal RA, Dec, and 1D wavelength arrays for
          the output spectrum.
       3. Set up an aperture correction to apply to each spectrum,
          if ``apcorr_ref_model`` is not None.
       4. Loop over integrations to extract all spectra.

    For each integration, the extraction process is:

       1. Extract summed flux and variance values from the aperture.
       2. Compute an average value (from flux / npixels), to be stored
          as the surface brightness.
       3. Convert the summed flux to flux density (Jy).
       4. Compute an error spectrum from the square root of the sum
          of the variance components.
       5. Set a DQ array, with DO_NOT_USE flags set where the
          flux is NaN.
       6. Create a spectral table to contain all extracted values
          and store it in a `SpecModel`.
       7. Apply the aperture correction to the spectral table, if
          available.
       8. Append the new SpecModel to the output_model.

    Parameters
    ----------
    input_model : DataModel
        Top-level datamodel containing metadata for the input data.
        If slit is not specified, ``input_model`` must also contain the
        spectral image(s) to extract, in the data attribute.
    slit : SlitModel or None
        One slit from an input MultiSlitModel or similar. If not None,
        ``slit`` must contain the spectral image(s) to extract, in the data
        attribute, along with appropriate WCS metadata.
    output_model : MultiSpecModel, TSOMultiSpecModel
        The output model to append spectral tables to.
    extract_ref_dict : dict or None
        Extraction parameters read in from the extract_1d reference file,
        or None, if there was no reference file.
    slitname : str
        Slit name for the input data.
    sp_order : int
        Spectral order for the input data.
    exp_type : str
        Exposure type for the input data.
    apcorr_ref_model : DataModel or None, optional
        The aperture correction reference datamodel, containing the
        APCORR reference file data.
    log_increment : int, optional
        If greater than 0 and the input data are multi-integration, a message
        will be written to the log every ``log_increment`` integrations.
    save_profile : bool, optional
        If True, the spatial profile created for the aperture will be returned
        as an ImageModel.  If False, the return value is None.
    save_scene_model : bool, optional
        If True, the flux model created during extraction will be returned
        as an ImageModel or CubeModel.  If False, the return value is None.
    save_residual_image : bool, optional
        If True, the residual image (from input minus scene model) will be returned
        as an ImageModel or CubeModel.  If False, the return value is None.
    **kwargs
        Additional options to pass to `get_extract_parameters`.

    Returns
    -------
    profile_model : ImageModel or None
        If ``save_profile`` is True, the return value is an ImageModel containing
        the spatial profile with aperture weights, used in extracting all
        integrations.
    scene_model : ImageModel, CubeModel, or None
        If ``save_scene_model`` is True, the return value is an ImageModel or CubeModel
        matching the input data, containing the flux model generated during
        extraction.
    residual : ImageModel, CubeModel, or None
        If ``save_residual_image`` is True, the return value is an ImageModel or CubeModel
        matching the input data, containing the residual image.
    """
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
    s_photom = input_model.meta.cal_step.photom
    if s_photom is not None and s_photom.upper() == "COMPLETE":
        photom_has_been_run = True
        flux_units = "Jy"
        f_var_units = "Jy^2"
        sb_units = "MJy/sr"
        sb_var_units = "MJy^2 / sr^2"
    else:
        photom_has_been_run = False
        flux_units = "DN/s"
        f_var_units = "DN^2 / s^2"
        sb_units = "DN/s"
        sb_var_units = "DN^2 / s^2"
        log.warning("The photom step has not been run.")

    # Get the source type for the data
    if slit is not None:
        source_type = slit.source_type
    else:
        if isinstance(input_model, datamodels.SlitModel):
            source_type = input_model.source_type
            if source_type is None:
                source_type = input_model.meta.target.source_type
                input_model.source_type = source_type
        else:
            source_type = input_model.meta.target.source_type

    # Turn off use_source_posn if the source is not POINT
    if source_type != "POINT" or exp_type in WFSS_EXPTYPES:
        if kwargs.get("use_source_posn") is None:
            kwargs["use_source_posn"] = False
            log.info(
                f"Setting use_source_posn to False for exposure type {exp_type}, "
                f"source type {source_type}"
            )

    if photom_has_been_run:
        pixel_solid_angle = data_model.meta.photometry.pixelarea_steradians
        if pixel_solid_angle is None:
            pixel_solid_angle = 1.0
            log.warning("Pixel area (solid angle) is not populated; the flux will not be correct.")
    else:
        pixel_solid_angle = 1.0  # not needed

    extract_params = get_extract_parameters(
        extract_ref_dict, data_model, slitname, sp_order, input_model.meta, **kwargs
    )

    if extract_params["match"] == NO_MATCH:
        log.critical("Missing extraction parameters.")
        raise ValueError("Missing extraction parameters.")
    elif extract_params["match"] == PARTIAL:
        log.info(f"Spectral order {sp_order} not found, skipping ...")
        raise ContinueError()

    extract_params["dispaxis"] = data_model.meta.wcsinfo.dispersion_direction
    if extract_params["dispaxis"] is None:
        log.warning("The dispersion direction information is missing, so skipping ...")
        raise ContinueError()

    # Set up spatial profiles and wavelength array,
    # to be used for every integration
    (ra, dec, wavelength, profile, bg_profile, nod_profile, limits) = define_aperture(
        input_model, slit, extract_params, exp_type
    )

    valid = np.isfinite(wavelength)
    wavelength = wavelength[valid]
    if np.sum(valid) == 0:
        log.error("Spectrum is empty; no valid data.")
        raise ContinueError()

    # Save the profile if desired
    if save_profile:
        if nod_profile is not None:
            profile_model = datamodels.ImageModel(profile + nod_profile)
        else:
            profile_model = datamodels.ImageModel(profile)
        profile_model.update(input_model, only="PRIMARY")
        profile_model.name = slitname
    else:
        profile_model = None

    # Set up aperture correction, to be used for every integration
    apcorr_available = False
    if source_type is not None and source_type.upper() == "POINT" and apcorr_ref_model is not None:
        log.info("Creating aperture correction.")
        # NIRSpec needs to use a wavelength in the middle of the
        # range rather than the beginning of the range
        # for calculating the pixel scale since some wavelengths at the
        # edges of the range won't map to the sky
        if instrument == "NIRSPEC":
            wl = np.median(wavelength)
        else:
            wl = wavelength.min()

        kwargs = {"location": (ra, dec, wl)}
        if not isinstance(input_model, datamodels.ImageModel):
            kwargs["slit_name"] = slitname
            if exp_type in ["NRS_FIXEDSLIT", "NRS_BRIGHTOBJ"]:
                kwargs["slit"] = slitname

        apcorr_model = select_apcorr(input_model)
        apcorr = apcorr_model(
            input_model, apcorr_ref_model.apcorr_table, apcorr_ref_model.sizeunit, **kwargs
        )
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

    # Set up a scene model and residual image to update if desired
    if save_scene_model:
        if len(integrations) > 1:
            scene_model = datamodels.CubeModel(shape)
        else:
            scene_model = datamodels.ImageModel()
        scene_model.update(input_model, only="PRIMARY")
        scene_model.name = slitname
    else:
        scene_model = None
    if save_residual_image:
        if len(integrations) > 1:
            residual = datamodels.CubeModel(shape)
        else:
            residual = datamodels.ImageModel()
        residual.update(input_model, only="PRIMARY")
        residual.name = slitname
    else:
        residual = None

    # Extract each integration
    spec_list = []
    for integ in integrations:
        (
            sum_flux,
            f_var_rnoise,
            f_var_poisson,
            f_var_flat,
            background,
            b_var_rnoise,
            b_var_poisson,
            b_var_flat,
            npixels,
            scene_model_2d,
            residual_2d,
        ) = extract_one_slit(data_model, integ, profile, bg_profile, nod_profile, extract_params)

        # Save the scene model and residual
        if save_scene_model:
            if isinstance(scene_model, datamodels.CubeModel):
                scene_model.data[integ] = scene_model_2d
            else:
                scene_model.data = scene_model_2d
        if save_residual_image:
            if isinstance(residual, datamodels.CubeModel):
                residual.data[integ] = residual_2d
            else:
                residual.data = residual_2d

        # Convert the sum to an average, for surface brightness.
        npixels_temp = np.where(npixels > 0.0, npixels, 1.0)
        npixels_squared = npixels_temp**2
        if extract_params["extraction_type"] == "optimal":
            # surface brightness makes no sense for an optimal extraction
            surf_bright = np.zeros_like(sum_flux)
            sb_var_poisson = np.zeros_like(sum_flux)
            sb_var_rnoise = np.zeros_like(sum_flux)
            sb_var_flat = np.zeros_like(sum_flux)
        else:
            surf_bright = sum_flux / npixels_temp  # may be reset below
            sb_var_poisson = f_var_poisson / npixels_squared
            sb_var_rnoise = f_var_rnoise / npixels_squared
            sb_var_flat = f_var_flat / npixels_squared
        background /= npixels_temp
        b_var_poisson = b_var_poisson / npixels_squared
        b_var_rnoise = b_var_rnoise / npixels_squared
        b_var_flat = b_var_flat / npixels_squared

        del npixels_temp

        # Convert to flux density.
        # The input units will normally be MJy / sr, but for NIRSpec
        # point-source spectra the units will be MJy.
        input_units_are_megajanskys = (
            photom_has_been_run and source_type == "POINT" and instrument == "NIRSPEC"
        )

        if photom_has_been_run:
            # for NIRSpec point sources
            if input_units_are_megajanskys:
                flux = sum_flux * 1.0e6  # MJy --> Jy
                f_var_poisson *= 1.0e12  # MJy**2 --> Jy**2
                f_var_rnoise *= 1.0e12  # MJy**2 --> Jy**2
                f_var_flat *= 1.0e12  # MJy**2 --> Jy**2
                surf_bright[:] = 0.0
                sb_var_poisson[:] = 0.0
                sb_var_rnoise[:] = 0.0
                sb_var_flat[:] = 0.0
                background[:] /= pixel_solid_angle  # MJy / sr
                b_var_poisson = b_var_poisson / pixel_solid_angle**2
                b_var_rnoise = b_var_rnoise / pixel_solid_angle**2
                b_var_flat = b_var_flat / pixel_solid_angle**2
            else:
                flux = sum_flux * pixel_solid_angle * 1.0e6  # MJy / steradian --> Jy
                f_var_poisson *= pixel_solid_angle**2 * 1.0e12  # (MJy / sr)**2 --> Jy**2
                f_var_rnoise *= pixel_solid_angle**2 * 1.0e12  # (MJy / sr)**2 --> Jy**2
                f_var_flat *= pixel_solid_angle**2 * 1.0e12  # (MJy / sr)**2 --> Jy**2
        else:
            flux = sum_flux  # count rate

        del sum_flux

        error = np.sqrt(f_var_poisson + f_var_rnoise + f_var_flat)
        sb_error = np.sqrt(sb_var_poisson + sb_var_rnoise + sb_var_flat)
        berror = np.sqrt(b_var_poisson + b_var_rnoise + b_var_flat)

        # Set DQ from the flux value
        dq = np.zeros(flux.shape, dtype=np.uint32)
        dq[np.isnan(flux)] = datamodels.dqflags.pixel["DO_NOT_USE"]

        # Make a table of the values, trimming to points with valid wavelengths only
        otab = np.array(
            list(
                zip(
                    wavelength,
                    flux[valid],
                    error[valid],
                    f_var_poisson[valid],
                    f_var_rnoise[valid],
                    f_var_flat[valid],
                    surf_bright[valid],
                    sb_error[valid],
                    sb_var_poisson[valid],
                    sb_var_rnoise[valid],
                    sb_var_flat[valid],
                    dq[valid],
                    background[valid],
                    berror[valid],
                    b_var_poisson[valid],
                    b_var_rnoise[valid],
                    b_var_flat[valid],
                    npixels[valid],
                    strict=False,
                )
            ),
            dtype=datamodels.SpecModel().spec_table.dtype,
        )

        spec = datamodels.SpecModel(spec_table=otab)
        spec.meta.wcs = spec_wcs.create_spectral_wcs(ra, dec, wavelength)
        spec.spec_table.columns["wavelength"].unit = "um"
        spec.spec_table.columns["flux"].unit = flux_units
        spec.spec_table.columns["flux_error"].unit = flux_units
        spec.spec_table.columns["flux_var_poisson"].unit = f_var_units
        spec.spec_table.columns["flux_var_rnoise"].unit = f_var_units
        spec.spec_table.columns["flux_var_flat"].unit = f_var_units
        spec.spec_table.columns["surf_bright"].unit = sb_units
        spec.spec_table.columns["sb_error"].unit = sb_units
        spec.spec_table.columns["sb_var_poisson"].unit = sb_var_units
        spec.spec_table.columns["sb_var_rnoise"].unit = sb_var_units
        spec.spec_table.columns["sb_var_flat"].unit = sb_var_units
        spec.spec_table.columns["background"].unit = sb_units
        spec.spec_table.columns["bkgd_error"].unit = sb_units
        spec.spec_table.columns["bkgd_var_poisson"].unit = sb_var_units
        spec.spec_table.columns["bkgd_var_rnoise"].unit = sb_var_units
        spec.spec_table.columns["bkgd_var_flat"].unit = sb_var_units
        spec.slit_ra = ra
        spec.slit_dec = dec
        spec.spectral_order = sp_order
        spec.dispersion_direction = extract_params["dispaxis"]
        spec.detector = input_model.meta.instrument.detector

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

        if exp_type in WFSS_EXPTYPES:
            spectral_order = data_model.meta.wcsinfo.spectral_order
            if hasattr(data_model.meta, "filename"):
                # calwebb_spec3 case: no separate slit input to function
                spec.meta.filename = data_model.meta.filename
                spec.meta.group_id = _make_group_id(data_model, spectral_order)
            else:
                # calwebb_spec2 case: data_model is a slit so need to get this meta from input_model
                # In this case, group_id is expected to be the same for all slits
                # because spec list corresponds to different sources in the same exposure
                # The exception is spectral_order so we get that from data_model and it
                # is handled separately
                spec.meta.filename = getattr(input_model.meta, "filename", None)
                spec.meta.group_id = _make_group_id(input_model, spectral_order)
            spec.extract2d_xstart = data_model.xstart
            spec.extract2d_ystart = data_model.ystart
            if data_model.xstart is None or data_model.ystart is None:
                spec.extract2d_xstop = None
                spec.extract2d_ystop = None
            else:
                spec.extract2d_xstop = data_model.xstart + data_model.xsize
                spec.extract2d_ystop = data_model.ystart + data_model.ysize

        if apcorr is not None:
            if hasattr(apcorr, "tabulated_correction"):
                if apcorr.tabulated_correction is not None:
                    apcorr_available = True

            # See whether we can reuse the previous aperture correction
            # object.  If so, just apply the pre-computed correction to
            # save a ton of time.
            if apcorr_available:
                # reuse the last aperture correction
                apcorr.apply(spec.spec_table, use_tabulated=True)
            else:
                # Attempt to tabulate the aperture correction for later use.
                # If this fails, fall back on the old method.
                try:
                    apcorr.tabulate_correction(spec.spec_table)
                    apcorr.apply(spec.spec_table, use_tabulated=True)
                    log.debug("Tabulating aperture correction for use in multiple integrations.")
                except AttributeError:
                    log.debug("Computing aperture correction.")
                    apcorr.apply(spec.spec_table)

        spec_list.append(spec)

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
        log.info(f"All {input_model.data.shape[0]} integrations done")

    if len(spec_list) > 1:
        # For multi-int data, assemble a single TSOSpecModel from the list of spectra
        tso_spec = make_tso_specmodel(spec_list, segment=input_model.meta.exposure.segment_number)

        # Add to the output model
        output_model.spec.append(tso_spec)
    else:
        # Nothing further needed, just append the only spectrum
        # created to the output model
        output_model.spec.append(spec_list[0])

    return profile_model, scene_model, residual


def _make_group_id(model, spectral_order):
    """
    Generate a unique ID for an exposure group, including the spectral order.

    Parameters
    ----------
    model : DataModel
        The input data model.
    spectral_order : int
        The spectral order for the exposure.

    Returns
    -------
    str
        The group ID.
    """
    group_id = attrs_to_group_id(model.meta.observation)
    return group_id + f"_{spectral_order}"


def _make_output_model(data_model, meta_source):
    """
    Set up an output model matching the input.

    Parameters
    ----------
    data_model : DataModel
        A data model containing a "data" attribute.
    meta_source : DataModel
        A data model containing top-level metadata to copy to the output.

    Returns
    -------
    MultiSpecModel or TSOMultiSpecModel
        If the input data is multi-integration, a TSOMultiSpecModel is
        returned.  Otherwise, a MultiSpecModel is returned.
    """
    multi_int = (data_model.data.ndim == 3) and (data_model.data.shape[0] > 1)
    if multi_int:
        output_model = datamodels.TSOMultiSpecModel()
    else:
        output_model = datamodels.MultiSpecModel()
    if hasattr(meta_source, "int_times"):
        output_model.int_times = meta_source.int_times.copy()
    output_model.update(meta_source, only="PRIMARY")
    return output_model


def run_extract1d(
    input_model,
    extract_ref_name="N/A",
    apcorr_ref_name=None,
    psf_ref_name="N/A",
    extraction_type="box",
    smoothing_length=None,
    bkg_fit=None,
    bkg_order=None,
    log_increment=50,
    subtract_background=None,
    use_source_posn=None,
    position_offset=0.0,
    model_nod_pair=False,
    optimize_psf_location=True,
    save_profile=False,
    save_scene_model=False,
    save_residual_image=False,
):
    """
    Extract all 1-D spectra from an input model.

    Parameters
    ----------
    input_model : JWSTDataModel
        The input science model.
    extract_ref_name : str
        The name of the extract1d reference file, or "N/A".
    apcorr_ref_name : str or None
        Name of the APCORR reference file. Default is None
    psf_ref_name : str
        The name of the PSF reference file, or "N/A".
    extraction_type : str
        Extraction type ('box' or 'optimal').  Optimal extraction is
        only available if ``psf_ref_name`` is not "N/A".
    smoothing_length : int or None
        Width of a boxcar function for smoothing the background regions.
    bkg_fit : str or None
        Type of fitting to apply to background values in each column
        (or row, if the dispersion is vertical).  Allowed values are
        'mean', 'median', 'poly', or None.
    bkg_order : int or None
        Polynomial order for fitting to each column (or row, if the
        dispersion is vertical) of background. Only used if ``bkg_fit``
        is 'poly'.  Allowed values are >= 0.
    log_increment : int
        If ``log_increment`` is greater than 0 and the input data are
        multi-integration, a message will be written to the log every
        ``log_increment`` integrations.
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
    position_offset : float
        Number of pixels to shift the nominal source position in the
        cross-dispersion direction.
    model_nod_pair : bool
        If True, and if ``extraction_type`` is 'optimal', then a negative trace
        from nod subtraction is modeled alongside the positive source during
        extraction.  Even if set to True, this will be attempted only if the
        input data has been background subtracted and the dither pattern
        indicates that only 2 nods were used.
    optimize_psf_location : bool
        If True, and if ``extraction_type`` is 'optimal', then the source
        location will be optimized, via iterative comparisons of the scene
        model with the input data.
    save_profile : bool
        If True, the spatial profiles created for the input model will be returned
        as ImageModels. If False, the return value is None.
    save_scene_model : bool
        If True, a model of the 2D flux as defined by the extraction aperture
        is returned as an ImageModel or CubeModel.  If False, the return value
        is None.
    save_residual_image : bool
        If True, the residual image (from the input minus the scene model)
        is returned as an ImageModel or CubeModel.  If False, the return value
        is None.

    Returns
    -------
    output_model : MultiSpecModel or TSOMultiSpecModel
        A new data model containing the extracted spectra.
    profile_model : ModelContainer, ImageModel, or None
        If ``save_profile`` is True, the return value is an ImageModel containing
        the spatial profile with aperture weights, used in extracting a single
        slit, or else a container of ImageModels, one for each slit extracted.
        Otherwise, the return value is None.
    scene_model : ModelContainer, ImageModel, CubeModel, or None
        If ``save_scene_model`` is True, the return value is an ImageModel or CubeModel
        matching the input data, containing a model of the flux as defined by the
        aperture, created during extraction. Otherwise, the return value is None.
    residual : ModelContainer, ImageModel, CubeModel, or None
        If ``save_residual_image`` is True, the return value is an ImageModel or CubeModel
        matching the input data, containing the residual image (from the input minus
        the scene model). Otherwise, the return value is None.
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
    extract_ref_dict = read_extract1d_ref(extract_ref_name)

    # Check for non-null PSF reference file
    if psf_ref_name == "N/A":
        if extraction_type != "box":
            log.warning(f"Optimal extraction is not available for EXP_TYPE {exp_type}")
            log.warning("Defaulting to box extraction.")
            extraction_type = "box"

    # Read in the aperture correction reference file
    apcorr_ref_model = None
    if apcorr_ref_name is not None and apcorr_ref_name != "N/A":
        if extraction_type == "optimal":
            log.warning("Turning off aperture correction for optimal extraction")
        else:
            apcorr_ref_model = read_apcorr_ref(apcorr_ref_name, exp_type)

    # This will be relevant if we're asked to extract a spectrum
    # and the spectral order is zero.
    # That's only OK if the disperser is a prism.
    prism_mode = is_prism(meta_source)

    # Handle inputs that contain one or more slit models
    profile_model = None
    scene_model = None
    residual = None
    if isinstance(input_model, (ModelContainer, datamodels.MultiSlitModel)):
        if isinstance(input_model, ModelContainer):
            slits = input_model
        else:
            slits = input_model.slits

        # Save original use_source_posn value, because it can get
        # toggled within the following loop over slits
        save_use_source_posn = use_source_posn

        # Make a container for the profile models, if needed
        if save_profile:
            profile_model = ModelContainer()
        if save_scene_model:
            scene_model = ModelContainer()
        if save_residual_image:
            residual = ModelContainer()

        # Set up the output model
        output_model = _make_output_model(slits[0], meta_source)

        for slit in slits:  # Loop over the slits in the input model
            log.info(f"Working on slit {slit.name}")
            log.debug(f"Slit is of type {type(slit)}")

            slitname = slit.name
            use_source_posn = save_use_source_posn  # restore original value

            if np.size(slit.data) <= 0:
                log.info(f"No data for slit {slit.name}, skipping ...")
                continue

            sp_order = get_spectral_order(slit)
            if sp_order == 0 and not prism_mode:
                log.info("Spectral order 0 is a direct image, skipping ...")
                continue

            try:
                profile, slit_scene_model, slit_residual = create_extraction(
                    meta_source,
                    slit,
                    output_model,
                    extract_ref_dict,
                    slitname,
                    sp_order,
                    exp_type,
                    apcorr_ref_model=apcorr_ref_model,
                    log_increment=log_increment,
                    save_profile=save_profile,
                    save_scene_model=save_scene_model,
                    save_residual_image=save_residual_image,
                    psf_ref_name=psf_ref_name,
                    extraction_type=extraction_type,
                    smoothing_length=smoothing_length,
                    bkg_fit=bkg_fit,
                    bkg_order=bkg_order,
                    subtract_background=subtract_background,
                    use_source_posn=use_source_posn,
                    position_offset=position_offset,
                    model_nod_pair=model_nod_pair,
                    optimize_psf_location=optimize_psf_location,
                )
            except ContinueError:
                continue

            if save_profile:
                profile_model.append(profile)
            if save_scene_model:
                scene_model.append(slit_scene_model)
            if save_residual_image:
                residual.append(slit_residual)

    else:
        # Define source of metadata
        slit = None

        # Default the slitname to the exp_type, in case there's no
        # better value available
        slitname = exp_type
        if isinstance(input_model, datamodels.ImageModel):
            # Get the slit name from the input model
            if hasattr(input_model, "name") and input_model.name is not None:
                slitname = input_model.name

        elif isinstance(input_model, (datamodels.CubeModel, datamodels.SlitModel)):
            # For NRS_BRIGHTOBJ, the slit name comes from the slit model info
            if exp_type == "NRS_BRIGHTOBJ" and hasattr(input_model, "name"):
                slitname = input_model.name

            # For NRS_FIXEDSLIT, the slit name comes from the FXD_SLIT keyword
            # in the model meta if not present in the input model
            if exp_type == "NRS_FIXEDSLIT":
                if hasattr(input_model, "name") and input_model.name is not None:
                    slitname = input_model.name
                else:
                    slitname = input_model.meta.instrument.fixed_slit

        else:
            log.error("The input file is not supported for this step.")
            raise TypeError("Can't extract a spectrum from this file.")

        # Set up the output model
        output_model = _make_output_model(input_model, meta_source)

        sp_order = get_spectral_order(input_model)
        if sp_order == 0 and not prism_mode:
            log.info("Spectral order 0 is a direct image, skipping ...")
        else:
            log.info(f"Processing spectral order {sp_order}")

            try:
                profile_model, scene_model, residual = create_extraction(
                    input_model,
                    slit,
                    output_model,
                    extract_ref_dict,
                    slitname,
                    sp_order,
                    exp_type,
                    apcorr_ref_model=apcorr_ref_model,
                    log_increment=log_increment,
                    save_profile=save_profile,
                    save_scene_model=save_scene_model,
                    save_residual_image=save_residual_image,
                    psf_ref_name=psf_ref_name,
                    extraction_type=extraction_type,
                    smoothing_length=smoothing_length,
                    bkg_fit=bkg_fit,
                    bkg_order=bkg_order,
                    subtract_background=subtract_background,
                    use_source_posn=use_source_posn,
                    position_offset=position_offset,
                    model_nod_pair=model_nod_pair,
                    optimize_psf_location=optimize_psf_location,
                )
            except ContinueError:
                pass

    # Copy the integration time information from the INT_TIMES table to keywords in the output file.
    if pipe_utils.is_tso(input_model):
        # int_times extension can stay for all TSO observations, but time
        # keywords are only populated for multi-int spectra
        if isinstance(output_model, datamodels.TSOMultiSpecModel):
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

    return output_model, profile_model, scene_model, residual
