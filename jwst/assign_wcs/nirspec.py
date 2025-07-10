"""
Tools to create the WCS pipeline NIRSPEC modes.

Calls create_pipeline() which redirects based on EXP_TYPE.
"""

import copy
import logging
import warnings

import gwcs
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits
from astropy.modeling import bounding_box as mbbox
from astropy.modeling import fix_inputs, models, bind_bounding_box
from astropy.modeling.models import Mapping, Identity, Const1D, Scale, Tabular1D
from gwcs import coordinate_frames as cf
from gwcs import selector
from gwcs.wcstools import grid_from_bounding_box

from stdatamodels.jwst.datamodels import (
    CollimatorModel,
    CameraModel,
    DisperserModel,
    FOREModel,
    IFUFOREModel,
    MSAModel,
    OTEModel,
    IFUPostModel,
    IFUSlicerModel,
    WavelengthrangeModel,
    FPAModel,
)
from stdatamodels.jwst.transforms.models import (
    Rotation3DToGWA,
    DirCos2Unitless,
    Slit2Msa,
    Slit2MsaLegacy,
    AngleFromGratingEquation,
    WavelengthFromGratingEquation,
    Gwa2Slit,
    Unitless2DirCos,
    Logical,
    Slit,
    Snell,
    RefractionIndexFromPrism,
)

from jwst.lib.exposure_types import is_nrs_ifu_lamp
from .util import MSAFileError, NoDataOnDetectorError, not_implemented_mode, velocity_correction
from . import pointing

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

FIXED_SLIT_NUMS = {"NONE": 0, "S200A1": 1, "S200A2": 2, "S400A1": 3, "S1600A1": 4, "S200B1": 5}

# Approximate fallback values for MSA slit scaling
MSA_SLIT_SCALES = (1.35, 1.15)

__all__ = [
    "create_pipeline",
    "imaging",
    "ifu",
    "slits_wcs",
    "get_open_slits",
    "generate_compound_bbox",
    "nrs_fs_slit_id",
    "nrs_fs_slit_name",
    "nrs_wcs_set_input",
    "nrs_ifu_wcs",
    "get_spectral_order_wrange",
]


def create_pipeline(input_model, reference_files, slit_y_range):
    """
    Create a pipeline list based on EXP_TYPE.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    exp_type = input_model.meta.exposure.type.lower()
    if input_model.meta.instrument.grating.lower() == "mirror":
        pipeline = imaging(input_model, reference_files)
    else:
        pipeline = exp_type2transform[exp_type](
            input_model, reference_files, slit_y_range=slit_y_range
        )
    if pipeline:
        log.info(f"Created a NIRSPEC {exp_type} pipeline with references {reference_files}")
    return pipeline


def imaging(input_model, reference_files):
    """
    Create the WCS pipeline for NIRSpec imaging data.

    It has the following coordinate frames:
    "detector" : the science frame
    "sca" : frame associated with the SCA
    "gwa" " just before the GWA going from detector to sky
    "msa_frame" : at the MSA
    "oteip" : after the FWA
    "v2v3" and "world"

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'disperser', 'collimator', 'wavelengthrange', 'fpa', 'camera',
        'ote', and 'fore' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    # Get the corrected disperser model
    disperser = get_disperser(input_model, reference_files["disperser"])

    # DMS to SCA transform
    dms2detector = dms_to_sca(input_model)[:-1]
    dms2detector.name = "dms2sca"
    dms2detector.inputs = ("x", "y")
    dms2detector.outputs = ("x", "y")

    # DETECTOR to GWA transform
    det2gwa = detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)[:-2]
    det2gwa.name = "det2gwa"
    det2gwa.inputs = ("x", "y")
    det2gwa.outputs = ("angle1", "angle2", "angle3")

    gwa_through = Const1D(-1) * Identity(1) & Const1D(-1) * Identity(1) & Identity(1)

    angles = [disperser["theta_x"], disperser["theta_y"], disperser["theta_z"], disperser["tilt_y"]]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name="rotation").inverse
    dircos2unitless = DirCos2Unitless(name="directional_cosines2unitless")

    col_model = CollimatorModel(reference_files["collimator"])
    col = col_model.model
    col_model.close()

    # Get the default spectral order and wavelength range and record them in the model.
    sporder, wrange = get_spectral_order_wrange(input_model, reference_files["wavelengthrange"])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = sporder

    lam = wrange[0] + (wrange[1] - wrange[0]) * 0.5

    # Scale wavelengths to microns if msa coordinates are terminal
    if input_model.meta.instrument.filter == "OPAQUE":
        lam *= 1e6

    lam_model = Mapping((0, 1, 1)) | Identity(2) & Const1D(lam)

    gwa2msa = gwa_through | rotation | dircos2unitless | col | lam_model
    gwa2msa.inverse = col.inverse | dircos2unitless.inverse | rotation.inverse | gwa_through
    gwa2msa.name = "gwa2msa"
    gwa2msa.inputs = ("angle1", "angle2", "angle3")
    gwa2msa.outputs = ("x_msa", "y_msa", "lam")

    # Create coordinate frames in the NIRSPEC WCS pipeline
    # "detector", "gwa", "msa", "oteip", "v2v3", "v2v3vacorr", "world"
    det, sca, gwa, msa_frame, oteip, v2v3, v2v3vacorr, world = create_imaging_frames()
    if input_model.meta.instrument.filter == "OPAQUE":
        # Pipeline ends with MSA coordinates
        imaging_pipeline = [(det, dms2detector), (sca, det2gwa), (gwa, gwa2msa), (msa_frame, None)]
        return imaging_pipeline

    # MSA to OTEIP transform
    msa2ote = msa_to_oteip(reference_files)
    msa2oteip = msa2ote | Mapping((0, 1), n_inputs=3)
    map1 = Mapping((0, 1, 0, 1))
    minv = msa2ote.inverse
    del minv.inverse
    msa2oteip.inverse = map1 | minv | Mapping((0, 1), n_inputs=3)
    msa2oteip.name = "msa2oteip"
    msa2oteip.inputs = ("x_msa", "y_msa", "lam")
    msa2oteip.outputs = ("x_ote", "y_ote")

    # OTEIP to V2,V3 transform
    with OTEModel(reference_files["ote"]) as f:
        oteip2v23 = f.model
    oteip2v23.name = "oteip2v23"
    oteip2v23.inputs = ("x_ote", "y_ote")
    oteip2v23.outputs = ("v2", "v3")

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    )
    va_corr.name = "va_corr"
    va_corr.inputs = ("v2", "v3")
    va_corr.outputs = ("v2", "v3")

    # V2, V3 to world (RA, DEC) transform
    tel2sky = pointing.v23tosky(input_model)
    tel2sky.name = "v2v3_to_sky"
    tel2sky.inputs = ("v2", "v3")
    tel2sky.outputs = ("ra", "dec")

    imaging_pipeline = [
        (det, dms2detector),
        (sca, det2gwa),
        (gwa, gwa2msa),
        (msa_frame, msa2oteip),
        (oteip, oteip2v23),
        (v2v3, va_corr),
        (v2v3vacorr, tel2sky),
        (world, None),
    ]

    return imaging_pipeline


def ifu(input_model, reference_files, slit_y_range=(-0.55, 0.55)):
    """
    Create the WCS pipeline for NIRSpec IFU data.

    The coordinate frames are:
    "detector" : the science frame
    "sca" : frame associated with the SCA
    "gwa" " just before the GWA going from detector to sky
    "slit_frame" : frame associated with the virtual slit
    "slicer' : frame associated with the slicer
    "msa_frame" : at the MSA
    "oteip" : after the FWA
    "v2v3" and "world"

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'ifufore', 'ifuslicer', 'ifupost', 'disperser', 'wavelengthrange',
        'fpa', 'camera', 'collimator', 'fore', and 'ote' reference files.
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    detector = input_model.meta.instrument.detector
    grating = input_model.meta.instrument.grating
    filt = input_model.meta.instrument.filter
    # Check for ifu reference files
    if (
        reference_files["ifufore"] is None
        and reference_files["ifuslicer"] is None
        and reference_files["ifupost"] is None
    ):
        # No ifu reference files, won't be able to create pipeline
        log_message = "No ifufore, ifuslicer or ifupost reference files"
        log.critical(log_message)
        raise RuntimeError(log_message)
    # Check for data actually being present on NRS2
    log_message = f"No IFU slices fall on detector {detector}"
    if detector == "NRS2" and grating.endswith("M"):
        # Mid-resolution gratings do not project on NRS2.
        log.critical(log_message)
        raise NoDataOnDetectorError(log_message)
    if detector == "NRS2" and grating == "G140H" and filt == "F070LP":
        # This combination of grating and filter does not project on NRS2.
        log.critical(log_message)
        raise NoDataOnDetectorError(log_message)

    slits = np.arange(30)
    # Get the corrected disperser model
    disperser = get_disperser(input_model, reference_files["disperser"])

    # Get the default spectral order and wavelength range and record them in the model.
    sporder, wrange = get_spectral_order_wrange(input_model, reference_files["wavelengthrange"])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = sporder

    # DMS to SCA transform
    dms2detector = dms_to_sca(input_model)
    dms2detector.name = "dms2sca"

    # DETECTOR to GWA transform
    det2gwa = detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)
    det2gwa.name = "det2gwa"

    # GWA to SLIT
    gwa2slit = gwa_to_ifuslit(slits, input_model, disperser, reference_files, slit_y_range)
    gwa2slit.name = "gwa2slit"

    # SLIT to MSA transform
    slit2slicer = ifuslit_to_slicer(slits, reference_files)
    slit2slicer.name = "slit2slicer"

    # SLICER to MSA Entrance
    slicer2msa = slicer_to_msa(reference_files)
    slicer2msa.name = "slicer2msa"

    det, sca, gwa, slit_frame, slicer_frame, msa_frame, oteip, v2v3, v2v3vacorr, world = (
        create_frames()
    )

    exp_type = input_model.meta.exposure.type.upper()

    is_lamp_exposure = exp_type in ["NRS_LAMP", "NRS_AUTOWAVE", "NRS_AUTOFLAT"]

    if input_model.meta.instrument.filter == "OPAQUE" or is_lamp_exposure:
        # If filter is "OPAQUE" or if internal lamp exposure,
        # the NIRSPEC WCS pipeline stops at the MSA.
        pipeline = [
            (det, dms2detector),
            (sca, det2gwa),
            (gwa, gwa2slit),
            (slit_frame, slit2slicer),
            (slicer_frame, slicer2msa),
            (msa_frame, None),
        ]
        return pipeline

    # MSA to OTEIP transform
    msa2oteip = ifu_msa_to_oteip(reference_files) & Identity(1)  # for slit
    msa2oteip.name = "msa2oteip"
    msa2oteip.inputs = ("x_msa", "y_msa", "lam", "name")
    msa2oteip.outputs = ("x_ote", "y_ote", "lam", "name")

    # OTEIP to V2,V3 transform
    # This includes a wavelength unit conversion from meters to microns.
    oteip2v23 = oteip_to_v23(reference_files) & Identity(1)
    oteip2v23.name = "oteip2v23"
    oteip2v23.inputs = ("x_ote", "y_ote", "lam", "name")
    oteip2v23.outputs = ("v2", "v3", "lam", "name")

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    ) & Identity(2)
    va_corr.name = "va_corr"
    va_corr.inputs = ("v2", "v3", "lam", "name")
    va_corr.outputs = ("v2", "v3", "lam", "name")

    # V2, V3 to sky
    tel2sky = pointing.v23tosky(input_model) & Identity(2)
    tel2sky.name = "v2v3_to_sky"
    tel2sky.inputs = ("v2", "v3", "lam", "name")
    tel2sky.outputs = ("ra", "dec", "lam", "name")

    # Create coordinate frames in the NIRSPEC WCS pipeline"
    #
    # The oteip2v2v3 transform converts the wavelength from meters (which is assumed
    # in the whole pipeline) to microns (which is the expected output)
    #
    # "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world"

    pipeline = [
        (det, dms2detector),
        (sca, det2gwa),
        (gwa, gwa2slit),
        (slit_frame, slit2slicer),
        (slicer_frame, slicer2msa),
        (msa_frame, msa2oteip),
        (oteip, oteip2v23),
        (v2v3, va_corr),
        (v2v3vacorr, tel2sky),
        (world, None),
    ]

    return pipeline


def slits_wcs(input_model, reference_files, slit_y_range):
    """
    Create the WCS pipeline for MOS and fixed slits.

    The coordinate frames are:
    "detector" : the science frame
    "sca" : frame associated with the SCA
    "gwa" " just before the GWA going from detector to sky
    "slit_frame" : frame associated with the virtual slit
    "msa_frame" : at the MSA
    "oteip" : after the FWA
    "v2v3" : at V2V3
    "world" : sky and spectral

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'msa', 'msametafile', 'disperser', 'wavelengthrange',
        'fpa', 'camera', 'collimator', 'fore', and 'ote' reference files.
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    open_slits_id = get_open_slits(input_model, reference_files, slit_y_range)
    if not open_slits_id:
        return None
    n_slits = len(open_slits_id)
    log.info(f"Computing WCS for {n_slits} open slitlets")

    msa_pipeline = slitlets_wcs(input_model, reference_files, open_slits_id)

    return msa_pipeline


def slitlets_wcs(input_model, reference_files, open_slits_id):
    """
    Create WCS pipeline for MOS and Fixed slits for the specific open shutters/slits.

    ``slit_y_range`` is taken from ``slit.ymin`` and ``slit.ymax``.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'msa', 'disperser', 'wavelengthrange',
        'fpa', 'camera', 'collimator', 'fore', and 'ote' reference files.
    open_slits_id : list
        A list of slit IDs that are open.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.

    Notes
    -----
    This function is also used by the ``msaflagopen`` step.
    """
    # Get the corrected disperser model
    disperser = get_disperser(input_model, reference_files["disperser"])

    # Get the default spectral order and wavelength range and record them in the model.
    sporder, wrange = get_spectral_order_wrange(input_model, reference_files["wavelengthrange"])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    log.info(f"SPORDER= {sporder}, wrange={wrange}")
    input_model.meta.wcsinfo.spectral_order = sporder

    # DMS to SCA transform
    dms2detector = dms_to_sca(input_model)
    dms2detector.name = "dms2sca"

    # DETECTOR to GWA transform
    det2gwa = detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)
    det2gwa.name = "det2gwa"

    # GWA to SLIT
    gwa2slit = gwa_to_slit(open_slits_id, input_model, disperser, reference_files)
    gwa2slit.name = "gwa2slit"

    # SLIT to MSA transform
    slit2msa = slit_to_msa(open_slits_id, reference_files["msa"])
    slit2msa.name = "slit2msa"

    # Create coordinate frames in the NIRSPEC WCS pipeline"
    # "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "v2v3vacorr", "world"
    # _ would be the slicer_frame that is not used
    det, sca, gwa, slit_frame, _, msa_frame, oteip, v2v3, v2v3vacorr, world = create_frames()

    exp_type = input_model.meta.exposure.type.upper()

    is_lamp_exposure = exp_type in ["NRS_LAMP", "NRS_AUTOWAVE", "NRS_AUTOFLAT"]

    if input_model.meta.instrument.filter == "OPAQUE" or is_lamp_exposure:
        # convert to microns if the pipeline ends earlier
        msa_pipeline = [
            (det, dms2detector),
            (sca, det2gwa),
            (gwa, gwa2slit),
            (slit_frame, slit2msa),
            (msa_frame, None),
        ]
        return msa_pipeline

    # MSA to OTEIP transform
    msa2oteip = msa_to_oteip(reference_files) & Identity(1)
    msa2oteip.name = "msa2oteip"
    msa2oteip.inputs = ("x_msa", "y_msa", "lam", "name")
    msa2oteip.outputs = ("x_ote", "y_ote", "lam", "name")

    # OTEIP to V2,V3 transform
    # This includes a wavelength unit conversion from meters to microns.
    oteip2v23 = oteip_to_v23(reference_files) & Identity(1)
    oteip2v23.name = "oteip2v23"
    oteip2v23.inputs = ("x_ote", "y_ote", "lam", "name")
    oteip2v23.outputs = ("v2", "v3", "lam", "name")

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    ) & Identity(2)
    va_corr.name = "va_corr"
    va_corr.inputs = ("v2", "v3", "lam", "name")
    va_corr.outputs = ("v2", "v3", "lam", "name")

    # V2, V3 to sky
    tel2sky = pointing.v23tosky(input_model) & Identity(2)
    tel2sky.name = "v2v3_to_sky"
    tel2sky.inputs = ("v2", "v3", "lam", "name")
    tel2sky.outputs = ("ra", "dec", "lam", "name")

    msa_pipeline = [
        (det, dms2detector),
        (sca, det2gwa),
        (gwa, gwa2slit),
        (slit_frame, slit2msa),
        (msa_frame, msa2oteip),
        (oteip, oteip2v23),
        (v2v3, va_corr),
        (v2v3vacorr, tel2sky),
        (world, None),
    ]

    return msa_pipeline


def get_open_slits(input_model, reference_files=None, slit_y_range=(-0.55, 0.55)):
    """
    Return the open slits/shutters in a MOS or Fixed Slits exposure.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'msa', 'msametafile', 'disperser', 'wavelengthrange',
        'fpa', 'camera', and 'collimator' reference files.
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.

    Returns
    -------
    slits : list[Slit]
        A list of `~stdatamodels.jwst.transforms.models.Slit` objects.
    """
    exp_type = input_model.meta.exposure.type.lower()
    lamp_mode = input_model.meta.instrument.lamp_mode
    if isinstance(lamp_mode, str):
        lamp_mode = lamp_mode.lower()
    else:
        lamp_mode = "none"

    # MOS/MSA exposure requiring MSA metadata file
    if exp_type in ["nrs_msaspec", "nrs_autoflat"] or (
        (exp_type in ["nrs_lamp", "nrs_autowave"]) and (lamp_mode == "msaspec")
    ):
        prog_id = input_model.meta.observation.program_number.lstrip("0")
        msa_metadata_file, msa_metadata_id, dither_point = get_msa_metadata(
            input_model, reference_files
        )
        if reference_files is not None and "msa" in reference_files:
            slit_scales = get_msa_slit_scales(reference_files["msa"])
        else:
            slit_scales = None
        slits = get_open_msa_slits(
            prog_id, msa_metadata_file, msa_metadata_id, dither_point, slit_y_range, slit_scales
        )

    # Fixed slits exposure (non-TSO)
    elif exp_type == "nrs_fixedslit":
        slits = get_open_fixed_slits(input_model, slit_y_range)

    # Bright object (TSO) exposure in S1600A1 fixed slit
    elif exp_type == "nrs_brightobj":
        slits = [
            Slit(
                "S1600A1",
                3,
                0,
                0,
                0,
                slit_y_range[0],
                slit_y_range[1],
                5,
                1,
                slit_id=nrs_fs_slit_id("S1600A1"),
            )
        ]

    # Lamp exposure using fixed slits
    elif exp_type in ["nrs_lamp", "nrs_autowave"]:
        if lamp_mode in ["fixedslit", "brightobj"]:
            slits = get_open_fixed_slits(input_model, slit_y_range)
    else:
        raise ValueError(f"EXP_TYPE {exp_type.upper()} is not supported")

    if reference_files is not None and slits:
        slits = validate_open_slits(input_model, slits, reference_files)
        log.info(
            f"Slits projected on detector {input_model.meta.instrument.detector}: "
            f"{[str(sl.name) for sl in slits]}"
        )
    if not slits:
        log_message = f"No open slits fall on detector {input_model.meta.instrument.detector}."
        log.critical(log_message)
        raise NoDataOnDetectorError(log_message)
    return slits


def get_open_fixed_slits(input_model, slit_y_range=(-0.55, 0.55)):
    """
    Return the open fixed slits.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.

    Returns
    -------
    slits : list[Slit]
        A list of `~stdatamodels.jwst.transforms.models.Slit` objects.
    """
    if input_model.meta.subarray.name is None:
        raise ValueError("Input file is missing SUBARRAY value/keyword.")
    if input_model.meta.instrument.fixed_slit is None:
        input_model.meta.instrument.fixed_slit = "NONE"

    primary_slit = input_model.meta.instrument.fixed_slit
    ylow, yhigh = slit_y_range

    # Slits are defined with hardwired source ID's, based on the assignments
    # in the "FIXED_SLIT_NUMS" dictionary. Exact assignments depend on whether the
    # slit is the "primary" and hence contains the target of interest. The
    # source_id for the primary slit is always 1, while source_ids for secondary
    # slits is a two-digit value, where the first (tens) digit corresponds to
    # primary slit in use and the second (ones) digit corresponds to the
    # secondary slit number. All fixed slits are assigned to the virtual MSA
    # quadrant 5.
    #
    # Slit(Name, ShutterID, DitherPos, Xcen, Ycen, Ymin, Ymax, Quad, SourceID)
    s2a1 = Slit(
        "S200A1",
        0,
        0,
        0,
        0,
        ylow,
        yhigh,
        5,
        1 if primary_slit == "S200A1" else 10 * FIXED_SLIT_NUMS[primary_slit] + 1,
        slit_id=nrs_fs_slit_id("S200A1"),
    )
    s2a2 = Slit(
        "S200A2",
        1,
        0,
        0,
        0,
        ylow,
        yhigh,
        5,
        1 if primary_slit == "S200A2" else 10 * FIXED_SLIT_NUMS[primary_slit] + 2,
        slit_id=nrs_fs_slit_id("S200A2"),
    )
    s4a1 = Slit(
        "S400A1",
        2,
        0,
        0,
        0,
        ylow,
        yhigh,
        5,
        1 if primary_slit == "S400A1" else 10 * FIXED_SLIT_NUMS[primary_slit] + 3,
        slit_id=nrs_fs_slit_id("S400A1"),
    )
    s16a1 = Slit(
        "S1600A1",
        3,
        0,
        0,
        0,
        ylow,
        yhigh,
        5,
        1 if primary_slit == "S1600A1" else 10 * FIXED_SLIT_NUMS[primary_slit] + 4,
        slit_id=nrs_fs_slit_id("S1600A1"),
    )
    s2b1 = Slit(
        "S200B1",
        4,
        0,
        0,
        0,
        ylow,
        yhigh,
        5,
        1 if primary_slit == "S200B1" else 10 * FIXED_SLIT_NUMS[primary_slit] + 5,
        slit_id=nrs_fs_slit_id("S200B1"),
    )

    # Decide which slits need to be added to this exposure
    subarray = input_model.meta.subarray.name.upper()
    slits = []
    if subarray == "SUBS200A1":
        slits.append(s2a1)
    elif subarray == "SUBS200A2":
        slits.append(s2a2)
    elif subarray == "SUBS400A1":
        slits.append(s4a1)
    elif subarray in ("SUB2048", "SUB512", "SUB512S", "SUB1024A", "SUB1024B"):
        slits.append(s16a1)
    elif subarray == "SUBS200B1":
        slits.append(s2b1)
    else:
        slits.extend([s2a1, s2a2, s4a1, s16a1, s2b1])

    return slits


def get_msa_metadata(input_model, reference_files):
    """
    Get the MSA metadata file (MSAMTFL) and the msa metadata ID (MSAMETID).

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'msametafile' reference file.

    Returns
    -------
    slits : list[Slit]
        A list of `~stdatamodels.jwst.transforms.models.Slit` objects.
    """
    try:
        msa_config = reference_files["msametafile"]
    except (KeyError, TypeError):
        log.info("MSA metadata file not in reference files dict")
        log.info("Getting MSA metadata file from MSAMETFL keyword")
        msa_config = input_model.meta.instrument.msa_metadata_file
        if msa_config is None:
            message = "msa_metadata_file is None."
            log.critical(message)
            raise MSAFileError(message) from None
    msa_metadata_id = input_model.meta.instrument.msa_metadata_id
    if msa_metadata_id is None:
        message = "Missing msa_metadata_id (keyword MSAMETID)."
        log.critical(message)
        raise MSAFileError(message)
    dither_position = input_model.meta.dither.position_number
    if dither_position is None:
        message = "Missing dither pattern number (keyword PATT_NUM)."
        log.critical(message)
        raise MSAFileError(message)
    return msa_config, msa_metadata_id, dither_position


def get_msa_slit_scales(msa_ref_file):
    """
    Get slit area scaling factors for MSA shutters.

    Parameters
    ----------
    msa_ref_file : str
        The name of the MSA reference file (not metadata file).

    Returns
    -------
    scales : dict
        Keys are integer quadrant values (one-indexed).  Values
        are 2-tuples of float values (scale_x, scale_y).
    """
    msa = MSAModel(msa_ref_file)
    scales = {}
    for quadrant in range(1, 5):
        msa_quadrant = getattr(msa, f"Q{quadrant}")

        msa_data = msa_quadrant.data
        scale_x = (msa_data["XC"][1] - msa_data["XC"][0]) / msa_data["SIZEX"][0]
        scale_y = msa_data["YC"][365] / msa_data["SIZEY"][0]
        scales[quadrant] = (scale_x, scale_y)
    return scales


def get_open_msa_slits(
    prog_id,
    msa_file,
    msa_metadata_id,
    dither_position,
    slit_y_range=(-0.55, 0.55),
    slit_scales=None,
):
    """
    Return the open MOS slitlets.

    Computes (ymin, ymax) for each open slitlet.

    The msa_file is expected to contain data (tuples) with the following fields:

        ('slitlet_id', '>i2'),
        ('msa_metadata_id', '>i2'),
        ('shutter_quadrant', '>i2'),
        ('shutter_row', '>i2'),
        ('shutter_column', '>i2'),
        ('source_id', '>i2'),
        ('background', 'S1'),
        ('shutter_state', 'S6'),
        ('estimated_source_in_shutter_x', '>f4'),
        ('estimated_source_in_shutter_y', '>f4'),
        ('dither_point_index', '>i2'),
        ('primary_source', 'S1')

    For example, something like:
        (12, 2, 4, 251, 22, 1, 'Y', 'OPEN', nan, nan, 1, 'N'),

    Parameters
    ----------
    prog_id : str
        The program number
    msa_file : str
        MSA meta data file name, FITS keyword ``MSAMETFL``.
    msa_metadata_id : int
        The MSA meta id for the science file, FITS keyword ``MSAMETID``.
    dither_position : int
        The index in the dither pattern, FITS keyword ``PATT_NUM``.
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.
    slit_scales : dict or None
        A dictionary of scaling factors for MSA shutters.  Keys are integer quadrant values
        (one-indexed).  Values are 2-tuples of float values (scale_x, scale_y).
        If not provided, default values from `MSA_SLIT_SCALES` are used.

    Returns
    -------
    slitlets : list
        A list of `~stdatamodels.jwst.transforms.models.Slit` objects. Each slitlet is a tuple with
        ("name", "shutter_id", "xcen", "ycen", "ymin", "ymax",
        "quadrant", "source_id", "shutter_state", "source_name", "source_alias", "stellarity",
        "source_xpos", "source_ypos", "source_ra", "source_dec")
    """
    slitlets = []
    ylow, yhigh = slit_y_range
    # If they passed in a string then we shall assume it is the filename
    # of the configuration file.
    try:
        msa_file = fits.open(msa_file, memmap=False)
    except FileNotFoundError:
        message = f"Missing MSA meta (MSAMETFL) file {msa_file}"
        log.error(message)
        raise MSAFileError(message) from None
    except OSError:
        message = f"Unable to read MSA FITS file (MSAMETFL) {msa_file}"
        log.error(message)
        raise MSAFileError(message) from None
    except Exception:
        message = f"Problem reading MSA metafile (MSAMETFL) {msa_file}"
        log.error(message)
        raise MSAFileError(message) from None

    # Set an empty dictionary for slit_scales if not provided
    if slit_scales is None:
        slit_scales = {}

    # Get the shutter and source info tables from the _msa.fits file.
    msa_conf = msa_file[("SHUTTER_INFO", 1)]  # EXTNAME = 'SHUTTER_INFO'
    msa_source = msa_file[("SOURCE_INFO", 1)].data  # EXTNAME = 'SOURCE_INFO'

    # First we are going to filter the msa_file data on the msa_metadata_id
    # and dither_point_index.
    msa_data = [
        x
        for x in msa_conf.data
        if x["msa_metadata_id"] == msa_metadata_id and x["dither_point_index"] == dither_position
    ]
    log.info(
        f"Retrieving open MSA slitlets for msa_metadata_id = {msa_metadata_id} "
        f"and dither_index = {dither_position}"
    )

    # Sort the MSA rows by slitlet_id
    slitlet_sets = {}
    for row in msa_data:
        # Check for fixed slit: if set, then slitlet_id is null
        is_fs = False
        try:
            fixed_slit = row["fixed_slit"]
            if fixed_slit in FIXED_SLIT_NUMS.keys() and fixed_slit != "NONE":
                is_fs = True
        except (IndexError, ValueError, KeyError):
            # May be old-style MSA file without a fixed_slit column
            fixed_slit = None

        if is_fs:
            # Fixed slit - use the slit name as the ID
            slitlet_id = fixed_slit
        else:
            # MSA - use the slitlet ID
            slitlet_id = row["slitlet_id"]

        # Append the row for the slitlet
        if slitlet_id in slitlet_sets:
            slitlet_sets[slitlet_id].append(row)
        else:
            slitlet_sets[slitlet_id] = [row]

    # Add a margin to the slit y limits
    margin = 0.5

    # Now let's look at each unique slitlet id
    for slitlet_id, slitlet_rows in slitlet_sets.items():
        # Get the open shutter information from the slitlet rows
        open_shutters = [x["shutter_column"] for x in slitlet_rows]

        # How many shutters in the slitlet are labeled as "main" or "primary"?
        n_main_shutter = len([s for s in slitlet_rows if s["primary_source"] == "Y"])

        # Check for fixed slit sources defined in the MSA file
        is_fs = [False] * len(slitlet_rows)
        for i, slitlet in enumerate(slitlet_rows):
            try:
                if (
                    slitlet["fixed_slit"] in FIXED_SLIT_NUMS.keys()
                    and slitlet["fixed_slit"] != "NONE"
                ):
                    is_fs[i] = True
            except (IndexError, ValueError, KeyError):
                # May be old-style MSA file without a fixed_slit column
                pass

        # In the next part we need to calculate, find, or determine 5 things for each slit:
        #    quadrant, xcen, ycen, ymin, ymax

        # First, check for a fixed slit
        slit_id_number = None
        if all(is_fs) and len(slitlet_rows) == 1:
            # One fixed slit open for the source
            slitlet = slitlet_rows[0]

            # Use a standard number for fixed slit shutter id
            shutter_id = FIXED_SLIT_NUMS[slitlet_id] - 1
            xcen = ycen = 0
            quadrant = 5

            # Get a slit ID number for WCS propagation purposes
            slit_id_number = nrs_fs_slit_id(slitlet_id)

            # No additional margin for fixed slit bounding boxes
            ymin = ylow
            ymax = yhigh

            # Source position and id
            if n_main_shutter == 1:
                # Source is marked primary
                source_id = slitlet["source_id"]
                source_xpos = np.nan_to_num(slitlet["estimated_source_in_shutter_x"], nan=0.5)
                source_ypos = np.nan_to_num(slitlet["estimated_source_in_shutter_y"], nan=0.5)

                log.info(f"Found fixed slit {slitlet_id} with source_id = {source_id}.")

                # Get source info for this slitlet:
                # note that slits with a real source assigned have source_id > 0,
                # while slits with source_id < 0 contain "virtual" sources
                try:
                    source_name, source_alias, stellarity, source_ra, source_dec = [
                        (s["source_name"], s["alias"], s["stellarity"], s["ra"], s["dec"])
                        for s in msa_source
                        if s["source_id"] == source_id
                    ][0]
                except IndexError:
                    # Missing source information: assign a virtual source name
                    log.warning("Could not retrieve source info from MSA file")
                    source_name = f"{prog_id}_VRT{slitlet_id}"
                    source_alias = f"VRT{slitlet_id}"
                    stellarity = 0.0
                    source_ra = 0.0
                    source_dec = 0.0

            else:
                log.warning(f"Fixed slit {slitlet_id} is not a primary source; skipping it.")
                continue

        elif any(is_fs):
            # Unsupported fixed slit configuration
            message = (
                f"For slitlet_id = {slitlet_id}, metadata_id = {msa_metadata_id}, "
                f"dither_index = {dither_position}"
            )
            log.warning(message)
            message = "MSA configuration file has an unsupported fixed slit configuration."
            log.warning(message)
            msa_file.close()
            raise MSAFileError(message)

        # Now check for regular MSA slitlets
        elif n_main_shutter == 0:
            # There are no main shutters: all are background
            if len(open_shutters) == 1:
                jmin = jmax = j = open_shutters[0]
            else:
                jmin = min([s["shutter_column"] for s in slitlet_rows])
                jmax = max([s["shutter_column"] for s in slitlet_rows])
                j = jmin + (jmax - jmin) // 2
            ymax = yhigh + margin + (jmax - j) * 1.15
            ymin = -(-ylow + margin) + (jmin - j) * 1.15
            quadrant = slitlet_rows[0]["shutter_quadrant"]
            ycen = j
            xcen = slitlet_rows[0]["shutter_row"]  # grab the first as they are all the same
            # shutter numbers in MSA file are 1-indexed
            shutter_id = np.int64(xcen) + (np.int64(ycen) - 1) * 365
            # Background slits all have source_id=0 in the msa_file,
            # so assign a unique id based on the slitlet_id
            source_id = slitlet_id

            # Hardwire the source info for background slits, because there's
            # no source info for them in the msa_file
            source_xpos = 0.5
            source_ypos = 0.5
            source_name = f"{prog_id}_BKG{slitlet_id}"
            source_alias = f"BKG{slitlet_id}"
            stellarity = 0.0
            source_ra = 0.0
            source_dec = 0.0
            log.info(f"Slitlet {slitlet_id} is background only; assigned source_id={source_id}")

        # There is 1 main shutter: this is a slit containing either a real or virtual source
        elif n_main_shutter == 1:
            xcen, ycen, quadrant, source_xpos, source_ypos = [
                (
                    s["shutter_row"],
                    s["shutter_column"],
                    s["shutter_quadrant"],
                    np.nan_to_num(s["estimated_source_in_shutter_x"], nan=0.5),
                    np.nan_to_num(s["estimated_source_in_shutter_y"], nan=0.5),
                )
                for s in slitlet_rows
                if s["background"] == "N"
            ][0]
            # shutter numbers in MSA file are 1-indexed
            shutter_id = np.int64(xcen) + (np.int64(ycen) - 1) * 365

            # y-size
            jmin = min([s["shutter_column"] for s in slitlet_rows])
            jmax = max([s["shutter_column"] for s in slitlet_rows])
            j = ycen
            ymax = yhigh + margin + (jmax - j) * 1.15
            ymin = -(-ylow + margin) + (jmin - j) * 1.15

            # Get the source_id from the primary shutter entry
            source_id = None
            for i in range(len(slitlet_rows)):
                if slitlet_rows[i]["primary_source"] == "Y":
                    source_id = slitlet_rows[i]["source_id"]

            # Get source info for this slitlet;
            # note that slits with a real source assigned have source_id > 0,
            # while slits with source_id < 0 contain "virtual" sources
            try:
                source_name, source_alias, stellarity, source_ra, source_dec = [
                    (s["source_name"], s["alias"], s["stellarity"], s["ra"], s["dec"])
                    for s in msa_source
                    if s["source_id"] == source_id
                ][0]
            except IndexError:
                source_name = f"{prog_id}_VRT{slitlet_id}"
                source_alias = f"VRT{slitlet_id}"
                stellarity = 0.0
                source_ra = 0.0
                source_dec = 0.0
                log.warning(
                    f"Could not retrieve source info from MSA file; "
                    f"assigning virtual source_name={source_name}"
                )

            if source_id < 0:
                log.info(
                    f"Slitlet {slitlet_id} contains virtual source, with source_id={source_id}"
                )

        # More than 1 main shutter: Not allowed!
        else:
            message = (
                f"For slitlet_id = {slitlet_id}, metadata_id = {msa_metadata_id}, "
                f"and dither_index = {dither_position}"
            )
            log.warning(message)
            message = "MSA configuration file has more than 1 shutter with primary source"
            log.warning(message)
            msa_file.close()
            raise MSAFileError(message)

        # Convert source positions from PPS to Model coordinate frame.
        # The source x,y position in the shutter is given in the msa configuration file,
        # columns "estimated_source_in_shutter_x" and "estimated_source_in_shutter_y".
        # The source position is in a coordinate system associated with each shutter whose
        # origin is the lower left corner of the shutter, positive x is to the right
        # and positive y is upwards.
        source_xpos -= 0.5
        source_ypos -= 0.5

        # Create the shutter_state string
        all_shutters = _shutter_id_to_str(open_shutters, ycen)

        # Get the slit scale by quadrant, or return approximate
        # fallback values if not available
        scale_x, scale_y = slit_scales.get(quadrant, MSA_SLIT_SCALES)

        # Set the slit ID number for WCS propagation purposes to the slit name
        # if not already set
        if slit_id_number is None:
            slit_id_number = slitlet_id

        # Create the output list of tuples that contain the required
        # data for further computations
        slit_parameters = (
            slitlet_id,
            shutter_id,
            dither_position,
            xcen,
            ycen,
            ymin,
            ymax,
            quadrant,
            source_id,
            all_shutters,
            source_name,
            source_alias,
            stellarity,
            source_xpos,
            source_ypos,
            source_ra,
            source_dec,
            scale_x,
            scale_y,
            slit_id_number,
        )
        log.debug(f"Appending slit: {[str(s) for s in slit_parameters]}")
        slitlets.append(Slit(*slit_parameters))

    msa_file.close()
    return slitlets


def _shutter_id_to_str(open_shutters, ycen):
    """
    Return a string representing the open and closed shutters in a slitlet.

    Parameters
    ----------
    open_shutters : list
        List of IDs (shutter_id) of open shutters.
    ycen : int
        Y coordinate of main shutter.

    Returns
    -------
    all_shutters : str
        String representing the state of the shutters.
        "1" indicates an open shutter, "0" - a closed one, and
        "x" - the main shutter.
    """
    all_shutters = np.array(range(min(open_shutters), max(open_shutters) + 1))
    cen_ind = (all_shutters == ycen).nonzero()[0].item()
    for i in open_shutters:
        all_shutters[all_shutters == i] = 1
    all_shutters[all_shutters != 1] = 0
    all_shutters = all_shutters.astype(str)
    all_shutters[cen_ind] = "x"
    return "".join(all_shutters)


def get_spectral_order_wrange(input_model, wavelengthrange_file):
    """
    Read the spectral order and wavelength range from the reference file.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    wavelengthrange_file : str
        Reference file of type "wavelengthrange".

    Returns
    -------
    order : int
        The spectral order.
    wrange : list
        The wavelength range.
    """
    # Nirspec full spectral range
    full_range = [0.6e-6, 5.3e-6]

    filt = input_model.meta.instrument.filter
    lamp = input_model.meta.instrument.lamp_state
    grating = input_model.meta.instrument.grating
    exp_type = input_model.meta.exposure.type

    is_lamp_exposure = exp_type in ["NRS_LAMP", "NRS_AUTOWAVE", "NRS_AUTOFLAT"]

    wave_range_model = WavelengthrangeModel(wavelengthrange_file)
    wrange_selector = wave_range_model.waverange_selector
    if filt == "OPAQUE" or is_lamp_exposure:
        keyword = lamp + "_" + grating
    else:
        keyword = filt + "_" + grating
    try:
        index = wrange_selector.index(keyword)
    except (KeyError, ValueError):
        # Combination of filter_grating is not in wavelengthrange file.
        gratings = [s.split("_")[1] for s in wrange_selector]
        try:
            index = gratings.index(grating)
        except ValueError:  # grating not in list
            order = -1
            wrange = full_range
        else:
            order = wave_range_model.order[index]
            wrange = wave_range_model.wavelengthrange[index]
        log.info(
            f"Combination {keyword} missing in wavelengthrange file, setting "
            f"order to {order} and range to {wrange}."
        )
    else:
        # Combination of filter_grating is found in wavelengthrange file.
        order = wave_range_model.order[index]
        wrange = wave_range_model.wavelengthrange[index]

    wave_range_model.close()
    return order, wrange


def ifuslit_to_slicer(slits, reference_files):
    """
    Create the transform from ``slit_frame`` to ``slicer`` frame.

    Parameters
    ----------
    slits : list
        A list of slit IDs for all slices.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'ifuslicer' reference file.

    Returns
    -------
    model : `~astropy.modeling.Model`.
        Transform from ``slit_frame`` to ``slicer`` frame.
    """
    ifuslicer = IFUSlicerModel(reference_files["ifuslicer"])
    models = []
    ifuslicer_model = ifuslicer.model
    for slit in slits:
        slitdata = ifuslicer.data[slit]
        slitdata_model = (get_slit_location_model(slitdata)).rename("slitdata_model")
        slicer_model = slitdata_model | ifuslicer_model

        msa_transform = slicer_model
        models.append(msa_transform)
    ifuslicer.close()
    s2m = Slit2Msa(slits, models)

    # Identity is for passing the computed wavelength
    mapping = Mapping((0, 1, 3, 2))
    mapping.inverse = Mapping((0, 1, 3, 2))
    slit2slicer = s2m & Identity(1) | mapping
    slit2slicer.inputs = ("name", "x_slit", "y_slit", "lam")
    slit2slicer.outputs = ("x_slice", "y_slice", "lam", "name")
    return slit2slicer


def slicer_to_msa(reference_files):
    """
    Transform from slicer coordinates to MSA entrance (the IFUFORE transform).

    Parameters
    ----------
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'ifufore' reference file.

    Returns
    -------
    model : `~astropy.modeling.Model`
        Transform from ``slicer`` frame to ``msa_frame``.
    """
    with IFUFOREModel(reference_files["ifufore"]) as f:
        ifufore = f.model
    slicer2fore_mapping = Mapping((0, 1, 2, 2, 3))
    slicer2fore_mapping.inverse = Identity(4)
    ifu_fore_transform = slicer2fore_mapping | ifufore & Identity(2)
    ifu_fore_transform.inputs = ("x_slice", "y_slice", "lam", "name")
    ifu_fore_transform.outputs = ("x_msa", "y_msa", "lam", "name")
    return ifu_fore_transform


def slit_to_msa(open_slits, msafile):
    """
    Create the transform from ``slit_frame`` to ``msa_frame``.

    Parameters
    ----------
    open_slits : list
        A list of slit IDs for all open shutters/slitlets.
    msafile : str
        The name of the msa reference file.

    Returns
    -------
    model : `~stdatamodels.jwst.transforms.Slit2Msa` model.
        Transform from ``slit_frame`` to ``msa_frame``.
    """
    msa = MSAModel(msafile)
    models = []
    slits = []
    for quadrant in range(1, 6):
        slits_in_quadrant = [s for s in open_slits if s.quadrant == quadrant]
        msa_quadrant = getattr(msa, f"Q{quadrant}")

        if any(slits_in_quadrant):
            msa_data = msa_quadrant.data
            msa_model = msa_quadrant.model

            for slit in slits_in_quadrant:
                slit_id = slit.shutter_id
                # Shutters are numbered starting from 1.
                # Fixed slits (Quadrant 5) are mapped starting from 0.
                if quadrant != 5:
                    slit_id = slit_id - 1
                slitdata = msa_data[slit_id]
                slitdata_model = get_slit_location_model(slitdata)
                msa_transform = slitdata_model | msa_model
                models.append(msa_transform)
                slits.append(slit)
    msa.close()
    s2m = Slit2Msa(slits, models)
    # Identity is for passing the computed wavelength and slit name
    mapping = Mapping((0, 1, 3, 2))
    mapping.inverse = Mapping((0, 1, 3, 2))
    transform = s2m & Identity(1) | mapping
    transform.inputs = ("name", "x_slit", "y_slit", "lam")
    transform.outputs = ("x_msa", "y_msa", "lam", "name")

    return transform


def gwa_to_ifuslit(slits, input_model, disperser, reference_files, slit_y_range):
    """
    Create the transform from ``gwa`` to ``slit_frame``.

    Parameters
    ----------
    slits : list
        A list of slit IDs for all IFU slits 0-29.
    input_model : JwstDataModel
        The input data model.
    disperser : `~jwst.datamodels.DisperserModel`
        A disperser model with the GWA correction applied to it.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'ifufore' reference file.
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.

    Returns
    -------
    model : `~stdatamodels.jwst.transforms.Gwa2Slit` model.
        Transform from ``gwa`` frame to ``slit_frame``.
    """
    ymin, ymax = slit_y_range

    agreq = angle_from_disperser(disperser, input_model)
    lgreq = wavelength_from_disperser(disperser, input_model)

    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    else:
        if velosys is not None:
            velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
            lgreq = lgreq | velocity_corr
            log.info(
                f"Applied Barycentric velocity correction : {velocity_corr[1].amplitude.value}"
            )
    # The wavelength units up to this point are
    # meters as required by the pipeline but the desired output wavelength units is microns.
    # So we are going to Scale the spectral units by 1e6 (meters -> microns)
    is_lamp_exposure = input_model.meta.exposure.type in [
        "NRS_LAMP",
        "NRS_AUTOWAVE",
        "NRS_AUTOFLAT",
    ]
    if input_model.meta.instrument.filter == "OPAQUE" or is_lamp_exposure:
        lgreq = lgreq | Scale(1e6)

    lam_cen = (
        0.5 * (input_model.meta.wcsinfo.waverange_end - input_model.meta.wcsinfo.waverange_start)
        + input_model.meta.wcsinfo.waverange_start
    )
    collimator2gwa = collimator_to_gwa(reference_files, disperser)
    mask = mask_slit(ymin, ymax)

    ifuslicer = IFUSlicerModel(reference_files["ifuslicer"])
    ifupost = IFUPostModel(reference_files["ifupost"])
    slit_models = []
    ifuslicer_model = ifuslicer.model
    for slit in slits:
        slitdata = ifuslicer.data[slit]
        slitdata_model = get_slit_location_model(slitdata)
        ifuslicer_transform = slitdata_model | ifuslicer_model
        ifupost_sl = getattr(ifupost, f"slice_{slit}")
        # construct IFU post transform
        ifupost_transform = _create_ifupost_transform(ifupost_sl)
        msa2gwa = ifuslicer_transform & Const1D(lam_cen) | ifupost_transform | collimator2gwa
        # TODO: Use model sets here
        gwa2slit = gwa_to_ymsa(msa2gwa, lam_cen=lam_cen, slit_y_range=slit_y_range)

        # The comments below list the input coordinates.
        bgwa2msa = (
            # (alpha_out, beta_out, gamma_out), angles at the GWA, coming from the camera
            # (0, - beta_out, alpha_out, beta_out)
            # (0, sy, alpha_out, beta_out)
            # (0, sy, 0, sy, sy, alpha_out, beta_out)
            # ( 0, sy, alpha_in, beta_in, gamma_in, alpha_out, beta_out)
            # (0, sy, alpha_in, beta_in,alpha_out)
            # (0, sy, lambda_computed)
            Mapping((0, 1, 0, 1), n_inputs=3)
            | Const1D(0) * Identity(1) & Const1D(-1) * Identity(1) & Identity(2)
            | Identity(1) & gwa2slit & Identity(2)
            | Mapping((0, 1, 0, 1, 1, 2, 3))
            | Identity(2) & msa2gwa & Identity(2)
            | Mapping((0, 1, 2, 3, 5), n_inputs=7)
            | Identity(2) & lgreq
            | mask
        )

        # transform from ``msa_frame`` to ``gwa`` frame (before the GWA going from detector to sky).
        msa2gwa_out = ifuslicer_transform & Identity(1) | ifupost_transform | collimator2gwa
        msa2bgwa = Mapping((0, 1, 2, 2)) | msa2gwa_out & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
        bgwa2msa.inverse = msa2bgwa
        slit_models.append(bgwa2msa)

    ifuslicer.close()
    ifupost.close()
    return Gwa2Slit(slits, slit_models)


def gwa_to_slit(open_slits, input_model, disperser, reference_files):
    """
    Create the transform from ``gwa`` to ``slit_frame``.

    Parameters
    ----------
    open_slits : list
        A list of slit IDs for all open shutters/slitlets.
    input_model : JwstDataModel
        The input data model.
    disperser : dict
        A corrected disperser ASDF object.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'collimator' and 'msa' reference files.

    Returns
    -------
    model : `~stdatamodels.jwst.transforms.Gwa2Slit` model.
        Transform from ``gwa`` frame to ``slit_frame``.
    """
    agreq = angle_from_disperser(disperser, input_model)
    collimator2gwa = collimator_to_gwa(reference_files, disperser)
    lgreq = wavelength_from_disperser(disperser, input_model)

    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    else:
        if velosys is not None:
            velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
            lgreq = lgreq | velocity_corr
            log.info(
                f"Applied Barycentric velocity correction : {velocity_corr[1].amplitude.value}"
            )

    # The wavelength units up to this point are
    # meters as required by the pipeline but the desired output wavelength units is microns.
    # So we are going to Scale the spectral units by 1e6 (meters -> microns)
    is_lamp_exposure = input_model.meta.exposure.type in [
        "NRS_LAMP",
        "NRS_AUTOWAVE",
        "NRS_AUTOFLAT",
    ]
    if input_model.meta.instrument.filter == "OPAQUE" or is_lamp_exposure:
        lgreq = lgreq | Scale(1e6)

    msa = MSAModel(reference_files["msa"])
    slit_models = []
    slits = []
    for quadrant in range(1, 6):
        slits_in_quadrant = [s for s in open_slits if s.quadrant == quadrant]
        log.info(f"There are {len(slits_in_quadrant)} open slits in quadrant {quadrant}")
        msa_quadrant = getattr(msa, f"Q{quadrant}")

        if any(slits_in_quadrant):
            msa_model = msa_quadrant.model
            msa_data = msa_quadrant.data

            for slit in slits_in_quadrant:
                mask = mask_slit(slit.ymin, slit.ymax)
                slit_id = slit.shutter_id
                # Shutter IDs are numbered starting from 1
                # while FS are numbered starting from 0.
                # "Quadrant 5 is for fixed slits.
                if quadrant != 5:
                    slit_id -= 1
                slitdata = msa_data[slit_id]
                slitdata_model = get_slit_location_model(slitdata)
                msa_transform = slitdata_model | msa_model
                msa2gwa = msa_transform | collimator2gwa
                # TODO: Use model sets here
                gwa2msa = gwa_to_ymsa(msa2gwa, slit=slit, slit_y_range=(slit.ymin, slit.ymax))
                bgwa2msa = (
                    Mapping((0, 1, 0, 1), n_inputs=3)
                    | Const1D(0) * Identity(1) & Const1D(-1) * Identity(1) & Identity(2)
                    | Identity(1) & gwa2msa & Identity(2)
                    | Mapping((0, 1, 0, 1, 2, 3))
                    | Identity(2) & msa2gwa & Identity(2)
                    | Mapping((0, 1, 2, 3, 5), n_inputs=7)
                    | Identity(2) & lgreq
                    | mask
                )
                #   Mapping((0, 1, 2, 5), n_inputs=7) | Identity(2) & lgreq | mask
                # and modify lgreq to accept alpha_in, beta_in, alpha_out
                # msa to before_gwa
                msa2bgwa = msa2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
                bgwa2msa.inverse = msa2bgwa
                slit_models.append(bgwa2msa)
                slits.append(slit)
    msa.close()
    return Gwa2Slit(slits, slit_models)


def angle_from_disperser(disperser, input_model):
    """
    Figure out the angle from the disperser model.

    For gratings this returns a form of the grating equation
    which computes the angle when lambda is known.
    For prism data this returns the Snell model.

    Parameters
    ----------
    disperser : dict
        A corrected disperser ASDF object.
    input_model : JwstDataModel
        The input data model.

    Returns
    -------
    model : `~astropy.modeling.Model`.
        Transform from wavelength to angle.
    """
    sporder = input_model.meta.wcsinfo.spectral_order
    if input_model.meta.instrument.grating.lower() != "prism":
        agreq = AngleFromGratingEquation(disperser.groovedensity, sporder, name="alpha_from_greq")
        return agreq

    system_temperature = input_model.meta.instrument.gwa_tilt
    system_pressure = disperser["pref"]

    snell = Snell(
        disperser["angle"],
        disperser["kcoef"],
        disperser["lcoef"],
        disperser["tcoef"],
        disperser["tref"],
        disperser["pref"],
        system_temperature,
        system_pressure,
        name="snell_law",
    )
    return snell


def wavelength_from_disperser(disperser, input_model):
    """
    Figure out the wavelength from the disperser model.

    For gratings this returns a form of the grating equation
    which computes lambda when all angles are known.

    For prism data this returns a lookup table model
    computing lambda from a known refraction index.

    Parameters
    ----------
    disperser : dict
        A corrected disperser ASDF object.
    input_model : JwstDataModel
        The input data model.

    Returns
    -------
    model : `~astropy.modeling.Model`.
        Transform from angle to wavelength
    """
    sporder = input_model.meta.wcsinfo.spectral_order
    if input_model.meta.instrument.grating.lower() != "prism":
        lgreq = WavelengthFromGratingEquation(
            disperser.groovedensity, sporder, name="lambda_from_gratingeq"
        )
        return lgreq

    lam = np.arange(0.3, 8.005, 0.005) * 1e-6
    system_temperature = input_model.meta.instrument.gwa_tilt
    if system_temperature is None:
        message = "Missing reference temperature (keyword GWA_TILT)."
        log.critical(message)
        raise KeyError(message)
    system_pressure = disperser["pref"]
    tref = disperser["tref"]
    pref = disperser["pref"]
    kcoef = disperser["kcoef"][:]
    lcoef = disperser["lcoef"][:]
    tcoef = disperser["tcoef"][:]
    n = Snell.compute_refraction_index(
        lam, system_temperature, tref, pref, system_pressure, kcoef, lcoef, tcoef
    )
    n = np.flipud(n)
    lam = np.flipud(lam)
    n_from_prism = RefractionIndexFromPrism(disperser["angle"], name="n_prism")

    tab = Tabular1D(points=(n,), lookup_table=lam, bounds_error=False)
    return n_from_prism | tab


def detector_to_gwa(reference_files, detector, disperser):
    """
    Transform from ``sca`` frame to ``gwa`` frame.

    Parameters
    ----------
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'fpa' and 'camera' reference files.
    detector : str
        The detector keyword.
    disperser : dict
        A corrected disperser ASDF object.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from DETECTOR frame to GWA frame.
    """
    with FPAModel(reference_files["fpa"]) as f:
        fpa = getattr(f, detector.lower() + "_model")
    with CameraModel(reference_files["camera"]) as f:
        camera = f.model

    angles = [disperser["theta_x"], disperser["theta_y"], disperser["theta_z"], disperser["tilt_y"]]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name="rotation")
    u2dircos = Unitless2DirCos(name="unitless2directional_cosines")
    # NIRSPEC 1- vs 0- based pixel coordinates issue #1781
    """
    The pipeline works with 0-based pixel coordinates. The Nirspec model,
    stored in reference files, is also 0-based. However, the algorithm specified
    by the IDT team specifies that pixel coordinates are 1-based. This is
    implemented below as a Shift(-1) & Shift(-1) transform. This makes the Nirspec
    instrument WCS pipeline "special" as it requires 1-based inputs.
    As a consequence many steps have to be modified to provide 1-based coordinates
    to the WCS call if the instrument is Nirspec. This is not always easy, especially
    when the step has no knowledge of the instrument.
    This is the reason the algorithm is modified to accept 0-based coordinates.
    This will be discussed in the future with the INS and IDT teams and may be solved
    by changing the algorithm but for now

    model = (models.Shift(-1) & models.Shift(-1) | fpa | camera | u2dircos | rotation)

    is changed to

    model = models.Shift(1) & models.Shift(1) | \
            models.Shift(-1) & models.Shift(-1) | fpa | camera | u2dircos | rotation
    """
    mapping = Mapping((3, 0, 1, 2))
    mapping.inverse = Mapping((1, 2, 3, 0))
    model = (fpa | camera | u2dircos | rotation) & Identity(1) | mapping
    model.inputs = ("x", "y", "name")
    model.outputs = ("name", "angle1", "angle2", "angle3")
    return model


def dms_to_sca(input_model):
    """
    Transform from ``detector`` to ``sca`` coordinates.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.

    Returns
    -------
    model : `~astropy.modeling.core.Model`
        Transform from DMS frame to SCA frame.
    """
    detector = input_model.meta.instrument.detector
    xstart = input_model.meta.subarray.xstart
    ystart = input_model.meta.subarray.ystart
    if xstart is None:
        xstart = 1
    if ystart is None:
        ystart = 1
    # The SCA coordinates are in full frame
    # The inputs are 1-based, remove -1 if/when they are 0-based
    # The outputs must be 1-based because this is what the model expects.
    # If xstart was 0-based and the inputs were 0-based ->
    # Shift(+1)
    subarray2full = models.Shift(xstart - 1) & models.Shift(ystart - 1)
    if detector == "NRS2":
        model = models.Shift(-2047) & models.Shift(-2047) | models.Scale(-1) & models.Scale(-1)
    elif detector == "NRS1":
        model = models.Identity(2)

    # x, y slit_name
    transform = (subarray2full | model) & Identity(1)
    transform.inputs = ("x", "y", "name")
    transform.outputs = ("x", "y", "name")

    return transform


def mask_slit(ymin=-0.55, ymax=0.55):
    """
    Return a model which masks out pixels in a NIRSpec cutout outside the slit.

    Uses ymin, ymax for the slit and the wavelength range to define the location of the slit.

    Parameters
    ----------
    ymin, ymax : float
        The relative min, max boundary of a slit.

    Returns
    -------
    model : `~astropy.modeling.core.Model`
        A model which takes x_slit, y_slit, lam inputs and substitutes the
        values outside the slit with NaN.
    """
    greater_than_ymax = Logical(condition="GT", compareto=ymax, value=np.nan)
    less_than_ymin = Logical(condition="LT", compareto=ymin, value=np.nan)

    model = (
        Mapping((0, 1, 2, 1))
        | Identity(3) & (greater_than_ymax | less_than_ymin | models.Scale(0))
        | Mapping((0, 1, 3, 2, 3))
        | Identity(1)
        & Mapping((0,), n_inputs=2) + Mapping((1,))
        & Mapping((0,), n_inputs=2) + Mapping((1,))
    )
    model.inverse = Identity(3)
    return model


def compute_bounding_box(
    transform, slit_name, wavelength_range, slit_ymin=-0.55, slit_ymax=0.55, refine=True
):
    """
    Compute the bounding box of the projection of a slit/slice on the detector.

    The edges of the slit are used to determine the location
    of the projection of the slit on the detector.
    Because the trace is curved and the wavelength_range may span the
    two detectors, y_min of the projection may be at an arbitrary wavelength.
    The transform is run with a regularly sampled wavelengths to determine y_min.

    Parameters
    ----------
    transform : `astropy.modeling.core.Model`
        The transform from detector to slit.
        `nrs_wcs_set_input` uses "detector to slit", validate_open_slits uses "slit to detector".
    slit_name : int, str, None
        Slit name for which to compute the bounding box.  If None, it is assumed
        the input WCS has already been fixed to a particular slit.
    wavelength_range : tuple
        The wavelength range for the combination of grating and filter.
    slit_ymin : float, optional
        Minimum value for the slit y edge.
    slit_ymax : float, optional
        Maximum value for the slit y edge.
    refine : bool, optional
        If True, the initial bounding box will be refined to limit it to
        ranges with valid wavelengths if possible.

    Returns
    -------
    bbox : tuple
        The bounding box of the projection of the slit on the detector.
    """
    if transform.has_inverse():
        # If the input transform has an inverse, then the transform
        # must be detector to slit and the inverse is slit to detector.
        detector2slit = transform
        slit2detector = detector2slit.inverse
    else:
        # If no inverse is available, the provided transform is slit to detector.
        detector2slit = None
        slit2detector = transform

    lam_min, lam_max = wavelength_range
    step = 1e-10
    nsteps = int((lam_max - lam_min) / step)
    lam_grid = np.linspace(lam_min, lam_max, nsteps)

    def check_range(lower, upper):
        return lower <= upper

    def bbox_from_range(x_range, y_range):
        # The -1 on both is technically because the output of slit2detector is 1-based coordinates.

        # add 10 px margin
        pad_x = (max(0, x_range.min() - 1 - 10) - 0.5, min(2047, x_range.max() - 1 + 10) + 0.5)
        # add 2 px margin
        pad_y = (max(0, y_range.min() - 1 - 2) - 0.5, min(2047, y_range.max() - 1 + 2) + 0.5)

        return pad_x, pad_y

    # Convert FS slit names to numbers
    if isinstance(slit_name, str):
        slit_name = nrs_fs_slit_id(slit_name)

    if slit_name is None:
        x_range_low, y_range_low = slit2detector([0] * nsteps, [slit_ymin] * nsteps, lam_grid)
        x_range_high, y_range_high = slit2detector([0] * nsteps, [slit_ymax] * nsteps, lam_grid)
    else:
        x_range_low, y_range_low, _ = slit2detector(
            slit_name, [0] * nsteps, [slit_ymin] * nsteps, lam_grid
        )
        x_range_high, y_range_high, _ = slit2detector(
            slit_name, [0] * nsteps, [slit_ymax] * nsteps, lam_grid
        )
    x_range = np.hstack((x_range_low, x_range_high))
    y_range = np.hstack((y_range_low, y_range_high))

    # Initial guess for ranges
    bbox = bbox_from_range(x_range, y_range)

    # Run inverse model to narrow range
    if refine and detector2slit is not None and check_range(*bbox[0]) and check_range(*bbox[1]):
        x, y = grid_from_bounding_box(bbox)
        if slit_name is None:
            _, _, lam = detector2slit(x, y)
        else:
            _, _, _, lam = detector2slit(x, y, slit_name)
        valid = np.isfinite(lam)
        if np.any(valid):
            y_range = y[np.isfinite(lam)]
            bbox = bbox_from_range(x_range, y_range)

    return bbox


def generate_compound_bbox(input_model, slits=None, wavelength_range=None, refine=True):
    """
    Generate a compound bounding box for multislit data.

    Parameters
    ----------
    input_model : DataModel
        Input model with WCS set in meta.wcs.
    slits : list or None, optional
        A list of Slit objects.  If not specified, open slits are retrieved
        from the "gwa" to "slit_frame" transform.
    wavelength_range : list or tuple or None, optional
        Wavelength range for the combination of filter and grating.
        If not specified, the range is calculated by calling
        `spectral_order_wrange_from_model`.
    refine : bool, optional
        If True, the initial bounding box will be refined to limit it to
        ranges with valid wavelengths if possible.

    Returns
    -------
    bounding_box : `~astropy.modeling.bounding_box.CompoundBoundingBox`
        Compound bounding box.  Access bounding boxes for individual
        slits with dictionary-like access, keyed on the slit ID.
    """
    if slits is None:
        slits = input_model.meta.wcs.get_transform("gwa", "slit_frame").slits
    wcs_forward_transform = input_model.meta.wcs.forward_transform
    if wavelength_range is None:
        _, wavelength_range = spectral_order_wrange_from_model(input_model)

    bbox_dict = {}

    wcs_forward_transform.inputs = ("x", "y", "name")
    transform = input_model.meta.wcs.get_transform("detector", "slit_frame")
    is_nirspec_ifu = (
        is_nrs_ifu_lamp(input_model) or input_model.meta.exposure.type.lower() == "nrs_ifu"
    )
    if is_nirspec_ifu:
        for slit in slits:
            bb = compute_bounding_box(transform, slit, wavelength_range, refine=refine)
            bbox_dict[slit] = bb
    else:
        for slit in slits:
            bb = compute_bounding_box(
                transform,
                slit.slit_id,
                wavelength_range,
                slit_ymin=slit.ymin,
                slit_ymax=slit.ymax,
                refine=refine,
            )
            bbox_dict[slit.slit_id] = bb

    cbb = mbbox.CompoundBoundingBox.validate(
        wcs_forward_transform, bbox_dict, selector_args=[("name", True)], order="F"
    )

    return cbb


def collimator_to_gwa(reference_files, disperser):
    """
    Transform from collimator to ``gwa`` frame.

    Includes the transforms:
    - through the collimator (going from sky to detector)
    - converting from unitless to directional cosines
    - a 3D rotation before the GWA using th ecorrected disperser angles.

    Parameters
    ----------
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'collimator' reference file.
    disperser : dict
        A corrected disperser ASDF object.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from collimator to ``gwa`` frame.
    """
    with CollimatorModel(reference_files["collimator"]) as f:
        collimator = f.model
    angles = [disperser["theta_x"], disperser["theta_y"], disperser["theta_z"], disperser["tilt_y"]]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name="rotation")
    u2dircos = Unitless2DirCos(name="unitless2directional_cosines")

    return collimator.inverse | u2dircos | rotation


def get_disperser(input_model, disperserfile):
    """
    Return the disperser data model with the GWA correction applied.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    disperserfile : str
        The name of the disperser reference file.

    Returns
    -------
    disperser : dict
        The corrected disperser information.
    """
    disperser = DisperserModel(disperserfile)
    xtilt = input_model.meta.instrument.gwa_xtilt
    ytilt = input_model.meta.instrument.gwa_ytilt
    disperser = correct_tilt(disperser, xtilt, ytilt)
    return disperser


def correct_tilt(disperser, xtilt, ytilt):
    """
    Correct the tilt of the grating by a measured grating tilt angle.

    Parameters
    ----------
    disperser : `~jwst.datamodels.DisperserModel`
        Disperser information.
    xtilt : float
        Value of GWAXTILT keyword - angle in arcsec
    ytilt : float
        Value of GWAYTILT keyword - angle in arcsec

    Returns
    -------
    disp : `~jwst.datamodels.DisperserModel`
        Corrected DisperserModel.

    Notes
    -----
    The GWA_XTILT keyword is used to correct the THETA_Y angle.
    The GWA_YTILT keyword is used to correct the THETA_X angle.
    """

    def _get_correction(gwa_tilt, tilt_angle):
        phi_exposure = gwa_tilt.tilt_model(tilt_angle)
        phi_calibrator = gwa_tilt.tilt_model(gwa_tilt.zeroreadings[0])
        del_theta = 0.5 * (phi_exposure - phi_calibrator) / 3600.0  # in deg
        return del_theta

    disp = disperser.copy()
    disperser.close()
    log.info(f"gwa_ytilt is {ytilt} deg")
    log.info(f"gwa_xtilt is {xtilt} deg")

    if xtilt is not None:
        theta_y_correction = _get_correction(disp.gwa_tiltx, xtilt)
        log.info(f"theta_y correction: {theta_y_correction} deg")
        disp["theta_y"] = disp.theta_y + theta_y_correction
    else:
        log.info("gwa_xtilt not applied")
    if ytilt is not None:
        theta_x_correction = _get_correction(disp.gwa_tilty, ytilt)
        log.info(f"theta_x correction: {theta_x_correction} deg")
        disp.theta_x = disp.theta_x + theta_x_correction
    else:
        log.info("gwa_ytilt not applied")
    return disp


def ifu_msa_to_oteip(reference_files):
    """
    Transform from ``msa_frame`` to ``oteip`` for IFU exposures.

    Parameters
    ----------
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'fore' reference file.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from MSA to OTEIP.
    """
    with FOREModel(reference_files["fore"]) as f:
        fore = f.model

    msa2fore_mapping = Mapping((0, 1, 2, 2), name="msa2fore_mapping")
    msa2fore_mapping.inverse = Mapping((0, 1, 2, 2), name="fore2msa")
    fore_transform = msa2fore_mapping | fore & Identity(1)
    return fore_transform


def msa_to_oteip(reference_files):
    """
    Transform from ``msa_frame`` to ``oteip`` for non IFU exposures.

    Parameters
    ----------
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'fore' reference file.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from MSA to OTEIP.
    """
    with FOREModel(reference_files["fore"]) as f:
        fore = f.model
    msa2fore_mapping = Mapping((0, 1, 2, 2), name="msa2fore_mapping")
    msa2fore_mapping.inverse = Identity(3)
    return msa2fore_mapping | (fore & Identity(1))


def oteip_to_v23(reference_files):
    """
    Transform from ``oteip`` frame to ``v2v3`` frame.

    Parameters
    ----------
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'ote' reference file.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from ``oteip`` to ``v2v3`` frame.
    """
    with OTEModel(reference_files["ote"]) as f:
        ote = f.model

    fore2ote_mapping = Identity(3, name="fore2ote_mapping")
    fore2ote_mapping.inverse = Mapping((0, 1, 2, 2))
    # Create the transform to v2/v3/lambda.  The wavelength units up to this point are
    # meters as required by the pipeline but the desired output wavelength units is microns.
    # So we are going to Scale the spectral units by 1e6 (meters -> microns)
    # The spatial units are currently in deg. Converting to arcsec.
    oteip2v23 = fore2ote_mapping | (ote & Scale(1e6))

    return oteip2v23


def create_frames():
    """
    Create the coordinate frames in the NIRSPEC WCS pipeline.

    Returns
    -------
    det, sca, gwa, slit_frame, slicer_frame, msa_frame, oteip, v2v3, v2v3vacorr, world : tuple
        The coordinate frames. Each is a `~gwcs.coordinate_frames.CoordinateFrame` object.
    """
    det = cf.Frame2D(name="detector", axes_order=(0, 1))
    sca = cf.Frame2D(name="sca", axes_order=(0, 1))
    gwa = cf.Frame2D(
        name="gwa", axes_order=(0, 1), unit=(u.rad, u.rad), axes_names=("alpha_in", "beta_in")
    )
    msa_spatial = cf.Frame2D(
        name="msa_spatial", axes_order=(0, 1), unit=(u.m, u.m), axes_names=("x_msa", "y_msa")
    )
    slit_spatial = cf.Frame2D(
        name="slit_spatial", axes_order=(0, 1), unit=("", ""), axes_names=("x_slit", "y_slit")
    )
    slicer_spatial = cf.Frame2D(
        name="slicer_spatial", axes_order=(0, 1), unit=("", ""), axes_names=("x_slicer", "y_slicer")
    )
    sky = cf.CelestialFrame(name="sky", axes_order=(0, 1), reference_frame=coord.ICRS())
    v2v3_spatial = cf.Frame2D(
        name="v2v3_spatial", axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=("v2", "v3")
    )
    v2v3vacorr_spatial = cf.Frame2D(
        name="v2v3vacorr_spatial",
        axes_order=(0, 1),
        unit=(u.arcsec, u.arcsec),
        axes_names=("v2", "v3"),
    )

    # The oteip_to_v23 incorporates a scale to convert the spectral units from
    # meters to microns.  So the v2v3 output frame will be in u.deg, u.deg, u.micron
    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
    )
    v2v3 = cf.CompositeFrame([v2v3_spatial, spec], name="v2v3")
    v2v3vacorr = cf.CompositeFrame([v2v3vacorr_spatial, spec], name="v2v3vacorr")
    slit_frame = cf.CompositeFrame([slit_spatial, spec], name="slit_frame")
    slicer_frame = cf.CompositeFrame([slicer_spatial, spec], name="slicer")
    msa_frame = cf.CompositeFrame([msa_spatial, spec], name="msa_frame")
    oteip_spatial = cf.Frame2D(
        name="oteip", axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=("X_OTEIP", "Y_OTEIP")
    )
    oteip = cf.CompositeFrame([oteip_spatial, spec], name="oteip")
    world = cf.CompositeFrame([sky, spec], name="world")
    return det, sca, gwa, slit_frame, slicer_frame, msa_frame, oteip, v2v3, v2v3vacorr, world


def create_imaging_frames():
    """
    Create the coordinate frames in the NIRSPEC WCS pipeline.

    Returns
    -------
    det, sca, gwa, msa, oteip, v2v3, v2v3vacorr, world : tuple
        The coordinate frames. Each is a `~gwcs.coordinate_frames.CoordinateFrame` object.
    """
    det = cf.Frame2D(name="detector", axes_order=(0, 1))
    sca = cf.Frame2D(name="sca", axes_order=(0, 1))
    gwa = cf.Frame2D(
        name="gwa", axes_order=(0, 1), unit=(u.rad, u.rad), axes_names=("alpha_in", "beta_in")
    )
    msa = cf.Frame2D(name="msa", axes_order=(0, 1), unit=(u.m, u.m), axes_names=("x_msa", "y_msa"))
    v2v3 = cf.Frame2D(
        name="v2v3", axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=("v2", "v3")
    )
    v2v3vacorr = cf.Frame2D(
        name="v2v3vacorr", axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=("v2", "v3")
    )
    oteip = cf.Frame2D(
        name="oteip", axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=("x_oteip", "y_oteip")
    )
    world = cf.CelestialFrame(name="world", axes_order=(0, 1), reference_frame=coord.ICRS())
    return det, sca, gwa, msa, oteip, v2v3, v2v3vacorr, world


def get_slit_location_model(slitdata):
    """
    Create the transform for the absolute position of a slit on the MSA.

    Parameters
    ----------
    slitdata : ndarray
        An array of shape (5,) with elements:
        slit_id, xcenter, ycenter, xsize, ysize
        This is the slit info in the MSa description file.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        A model which transforms relative position on the slit to
        absolute positions in the quadrant.
        This is later combined with the quadrant model to return
        absolute positions in the MSA.
    """
    num, xcenter, ycenter, xsize, ysize = slitdata
    # fmt: off
    model = models.Scale(xsize) & \
        models.Scale(ysize) | \
        models.Shift(xcenter) & \
        models.Shift(ycenter)
    # fmt: on
    return model


def gwa_to_ymsa(msa2gwa_model, lam_cen=None, slit=None, slit_y_range=None):
    """
    Determine the linear relation d_y(beta_in) for the aperture on the detector.

    Parameters
    ----------
    msa2gwa_model : `astropy.modeling.core.Model`
        The transform from the MSA to the GWA.
    lam_cen : float
        Central wavelength in meters.
    slit : `~stdatamodels.jwst.transforms.models.Slit`
        A Fixed slit or MOS slitlet.
    slit_y_range : tuple or None
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.
        Used for IFU mode only.

    Returns
    -------
    tab : `~astropy.modeling.Tabular1D`
        A 1D lookup table model.
    """
    nstep = 1000
    if slit is not None:
        ymin, ymax = slit.ymin, slit.ymax
    else:
        # The case of IFU data.
        ymin, ymax = slit_y_range
    dy = np.linspace(ymin, ymax, nstep)
    dx = np.zeros(dy.shape)
    if lam_cen is not None:
        # IFU case where IFUPOST has a wavelength dependent distortion
        cosin_grating_k = msa2gwa_model(dx, dy, [lam_cen] * nstep)
    else:
        cosin_grating_k = msa2gwa_model(dx, dy)
    beta_in = cosin_grating_k[1]

    tab = Tabular1D(points=(beta_in,), lookup_table=dy, bounds_error=False, name="tabular")
    return tab


def _nrs_wcs_set_slit_input_legacy(input_model, slit_name):
    """
    Return a WCS object for a specific slit, slice or shutter.

    Does not compute the bounding box.

    Parameters
    ----------
    input_model : JwstDataModel
        A WCS object for the all open slitlets in an observation.
    slit_name : int or str
        Slit.name of an open slit.

    Returns
    -------
    wcsobj : `~gwcs.wcs.WCS`
        WCS object for this slit.
    """
    wcsobj = input_model.meta.wcs

    slit_wcs = copy.deepcopy(wcsobj)
    slit_wcs.set_transform("sca", "gwa", wcsobj.pipeline[1].transform[1:])
    g2s = slit_wcs.pipeline[2].transform
    slit_wcs.set_transform("gwa", "slit_frame", g2s.get_model(slit_name))

    exp_type = input_model.meta.exposure.type
    is_nirspec_ifu = is_nrs_ifu_lamp(input_model) or (exp_type.lower() == "nrs_ifu")
    if is_nirspec_ifu:
        slit_wcs.set_transform(
            "slit_frame", "slicer", wcsobj.pipeline[3].transform.get_model(slit_name) & Identity(1)
        )
    else:
        slit_wcs.set_transform(
            "slit_frame",
            "msa_frame",
            wcsobj.pipeline[3].transform.get_model(slit_name) & Identity(1),
        )
    return slit_wcs


def nrs_wcs_set_input_legacy(
    input_model, slit_name, wavelength_range=None, slit_y_low=None, slit_y_high=None
):
    """
    Return a WCS object for a specific slit, slice or shutter.

    This function is intended to work with old-style NIRSpec WCS implementations,
    to support reading in WCS data from existing datamodels.  For new datamodels,
    produced after v1.18.0, use `nrs_wcs_set_input`.

    Parameters
    ----------
    input_model : JwstDataModel
        A WCS object for the all open slitlets in an observation.
    slit_name : int or str
        Slit.name of an open slit.
    wavelength_range : list
        Wavelength range for the combination of filter and grating.
    slit_y_low, slit_y_high : float
        The lower and upper bounds of the slit. Optional.

    Returns
    -------
    wcsobj : `~gwcs.wcs.WCS`
        WCS object for this slit.
    """
    warnings.warn(
        "The nrs_wcs_set_input_legacy function is intended for use with an "
        "old-style NIRSpec WCS pipeline. "
        "It will be removed in a future build.",
        DeprecationWarning,
        stacklevel=2,
    )

    def _get_y_range(input_model):
        # get the open slits from the model
        # Need them to get the slit ymin,ymax
        g2s = input_model.meta.wcs.get_transform("gwa", "slit_frame")
        open_slits = g2s.slits
        slit = [s for s in open_slits if s.name == slit_name][0]
        return slit.ymin, slit.ymax

    if wavelength_range is None:
        _, wavelength_range = spectral_order_wrange_from_model(input_model)

    slit_wcs = _nrs_wcs_set_slit_input_legacy(input_model, slit_name)
    transform = slit_wcs.get_transform("detector", "slit_frame")
    is_nirspec_ifu = (
        is_nrs_ifu_lamp(input_model) or input_model.meta.exposure.type.lower() == "nrs_ifu"
    )
    if is_nirspec_ifu:
        bb = compute_bounding_box(transform, None, wavelength_range)
    else:
        if slit_y_low is None or slit_y_high is None:
            slit_y_low, slit_y_high = _get_y_range(input_model)
        bb = compute_bounding_box(
            transform, None, wavelength_range, slit_ymin=slit_y_low, slit_ymax=slit_y_high
        )

    slit_wcs.bounding_box = bb
    return slit_wcs


def _fix_slit_name(transform, slit_name):
    """
    Create a new WCS transform by fixing the slit input to a specific value.

    If the inputs for the transform contain "name", then that input is
    fixed to the provided value.  The output value called "name" is also
    dropped from the output for the transform.

    For the inverse transform, if present, the same operations are performed, except
    that the inputs and outputs are not generally explicitly named.  The "name"
    input and output are assumed to be present, to be either the first or last
    input/output, and to match the ordering in the forward transform.  That is,
    if "name" is the first input and the last output in the forward transform,
    it is assumed to be the last input and the first output in the inverse
    transform.

    If "name" is not present in the input transform or if the transform is
    None, then the transform is simply returned unmodified.

    Parameters
    ----------
    transform : `astropy.modeling.core.Model` or None
        The transform to fix.
    slit_name : int or float
        The slit name to fix to.

    Returns
    -------
    new_transform : `astropy.modeling.core.Model`
        A new transform fixed to the input slit.  The "name" input is not
        required for the new transform, and it is not provided on output.
    """
    # Check if anything needs to be done for this transform
    if transform is None:
        return None

    # Directly pick out the slit-specific model to avoid
    # copying all transforms if possible
    if transform.name == "gwa2slit":
        new_transform = transform.get_model(slit_name)
        new_transform.name = transform.name

        # "name" was the first input and first output; it is now dropped.
        new_transform.inputs = transform.inputs[1:]
        new_transform.outputs = transform.outputs[1:]
        return new_transform
    if transform.name == "slit2msa" or transform.name == "slit2slicer":
        new_transform = transform[0].get_model(slit_name) & Identity(1)
        new_transform.name = transform.name

        # "name" was the first input and last output; it is now dropped.
        new_transform.inputs = transform.inputs[1:]
        new_transform.outputs = transform.outputs[:-1]
        return new_transform

    # Fix the "name" input for any other transform if possible
    if "name" not in transform.inputs:
        return transform

    # Name is present: fix it on input
    fixed = {"name": slit_name}

    # Get a mapping to drop the name on output as well
    output_names = []
    output_order = []
    name_idx = None
    for i, output_name in enumerate(transform.outputs):
        if output_name == "name":
            name_idx = i
        else:
            output_names.append(output_name)
            output_order.append(i)
    mapping = Mapping(tuple(output_order), n_inputs=transform.n_outputs)
    mapping.inverse = Identity(len(output_order))

    # New forward transform with fixed input and mapping
    new_transform = fix_inputs(transform, fixed) | mapping
    new_transform.name = transform.name
    new_transform.outputs = output_names

    # Do the same for the inverse transform, except that "name" is usually
    # not explicitly specified
    if transform.has_inverse():
        # get "name" equivalent in inverse input and outputs - it's either first or last
        n_output = transform.inverse.n_outputs
        if name_idx == 0:
            input_name = transform.inverse.inputs[0]
        else:
            input_name = transform.inverse.inputs[-1]
        output_name_idx = transform.inputs.index("name")
        if output_name_idx != 0:
            output_name_idx = n_output - 1

        # Fix for the input
        fix_inverse = {input_name: slit_name}

        # Mapping for the output
        output_order = tuple([i for i in range(n_output) if i != output_name_idx])
        mapping = Mapping(output_order, n_inputs=n_output)
        mapping.inverse = Identity(len(output_order))

        # Make and assign the new inverse transform
        new_inv_transform = fix_inputs(transform.inverse, fix_inverse) | mapping
        new_transform.inverse = new_inv_transform

    return new_transform


def nrs_fs_slit_id(slit_name):
    """
    Get the slit ID corresponding to a fixed slit name.

    Parameters
    ----------
    slit_name : str
        The name of the slit to identify.

    Returns
    -------
    int
        The standard ID for the slit.  Returns -100 if the slit name was not recognized.
    """
    slit_number = -100 + -1 * FIXED_SLIT_NUMS.get(slit_name, 0)
    return slit_number


def nrs_fs_slit_name(slit_id):
    """
    Get the slit name corresponding to a fixed slit ID.

    Parameters
    ----------
    slit_id : int, float, or str
        The slit ID.  Expected values are integers from -100 to -105.

    Returns
    -------
    str
        The name of the slit corresponding to the slit ID.  If the slit ID was not
        recognized, "NONE" is returned.
    """
    try:
        slit_number = -1 * int(slit_id) - 100
    except (ValueError, TypeError):
        return "NONE"
    for key, value in FIXED_SLIT_NUMS.items():
        if value == slit_number:
            return key
    return "NONE"


def nrs_wcs_set_input(input_model, slit_name):
    """
    Make a WCS object for a specific slit or slice.

    Transforms in the new WCS do not require slit IDs on input
    and do not report slit IDs in the output.

    If there is a compound bounding box set for the input model,
    a bounding box will be set for the slit WCS.  Otherwise,
    the bounding box will be None.

    Parameters
    ----------
    input_model : JwstDataModel
        A datamodel that contains a WCS object for the all open slitlets in
        an observation.
    slit_name : int or str
        Slit.name of an open slit.

    Returns
    -------
    wcsobj : `~gwcs.wcs.WCS`
        WCS object for this slit.
    """
    # Get the full WCS object
    full_wcs = input_model.meta.wcs

    # Check for an old-style WCS that needs different handling
    for step in full_wcs.pipeline:
        if isinstance(step.transform, Slit2MsaLegacy):
            # There is a legacy slit2msa transform somewhere in the pipeline.
            # Use the old deepcopy-and-replace method.
            return nrs_wcs_set_input_legacy(input_model, slit_name)

    # Convert FS slit names to numbers
    if str(slit_name).upper() in FIXED_SLIT_NUMS.keys():
        slit_id = nrs_fs_slit_id(str(slit_name).upper())
    else:
        slit_id = slit_name

    # Check slit ID to make sure it's present in the slits
    if "gwa" in full_wcs.available_frames:
        ifu_wcs = False
        slit_ids = full_wcs.get_transform("gwa", "slit_frame").slit_ids
    else:
        # assume it's a coordinate-based IFU WCS
        ifu_wcs = True
        slit_ids = list(range(30))
    if slit_id not in slit_ids:
        raise ValueError(f"Input slit name {slit_name} is not present in the input model.")

    # Fix the slit name input for all transforms
    new_pipeline = []
    for i, step in enumerate(full_wcs.pipeline):
        new_transform = _fix_slit_name(step.transform, slit_id)
        if ifu_wcs and i == 0:
            # No name to fix, so the new_transform is the same as the old
            # transform.  For the first transform in the pipeline, deepcopy
            # it so a slice-specific bounding box can be attached.
            new_transform = copy.deepcopy(new_transform)
        new_pipeline.append((step.frame, new_transform))

    # Make a new WCS with the fixed slit name input
    # Output values from transform will not include the slit name.
    slit_wcs = gwcs.WCS(new_pipeline)

    # Add the bounding box
    if hasattr(full_wcs.pipeline[1].transform, "selector"):
        # For coordinate-based IFU WCS, get the bounding box from the
        # selector that holds the slice map
        det2slicer_selector = full_wcs.pipeline[1].transform.selector
        slit_wcs.bounding_box = det2slicer_selector[slit_id + 1].bounding_box
    elif full_wcs.bounding_box is not None:
        # For all other modern WCS, get the bounding box from the
        # compound bounding box in the full WCS
        slit_wcs.bounding_box = full_wcs.bounding_box[slit_id]

    return slit_wcs


def validate_open_slits(input_model, open_slits, reference_files):
    """
    Remove slits which do not project on the detector from the list of open slits.

    For each slit computes the transform from the slit to the detector and
    determines the bounding box.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    open_slits : list
        List of open slits.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires the 'disperser', 'wavelengthrange', 'msa', 'collimator',
        'fpa', and 'camera' reference files

    Returns
    -------
    open_slits : list
        List of open slits that project onto the detector.
    """

    def _is_valid_slit(domain):
        xlow, xhigh = domain[0]
        ylow, yhigh = domain[1]
        if (
            xlow >= 2048
            or ylow >= 2048
            or xhigh <= 0
            or yhigh <= 0
            or xhigh - xlow < 2
            or yhigh - ylow < 1
        ):
            return False
        else:
            return True

    det2dms = dms_to_sca(input_model).inverse[:-1]
    # read models from reference file
    disperser = DisperserModel(reference_files["disperser"])
    disperser = correct_tilt(
        disperser, input_model.meta.instrument.gwa_xtilt, input_model.meta.instrument.gwa_ytilt
    )

    order, wrange = get_spectral_order_wrange(input_model, reference_files["wavelengthrange"])

    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = order
    agreq = angle_from_disperser(disperser, input_model)
    # GWA to detector
    # model = (fpa | camera | u2dircos | rotation) & Identity(1) | Mapping((3, 0, 1, 2))
    det2gwa = detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)
    gwa2det = det2gwa[:-2].inverse
    # collimator to GWA
    collimator2gwa = collimator_to_gwa(reference_files, disperser)

    # col2det = (Mapping((0, 1, 3)) | collimator2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq |
    #            gwa2det | det2dms)
    col2det = (
        Mapping((0, 1, 3, 2))
        | collimator2gwa & Identity(2)
        | Mapping((3, 0, 1, 2, 4))
        | agreq & Identity(1)
        | gwa2det & Identity(1)
        | det2dms & Identity(1)
    )

    slit2msa = slit_to_msa(open_slits, reference_files["msa"])[0]

    for slit in slit2msa.slits:
        msa2det = slit2msa & Identity(1) | col2det
        bb = compute_bounding_box(msa2det, slit.slit_id, wrange, slit.ymin, slit.ymax)
        valid = _is_valid_slit(bb)
        if not valid:
            log.info(
                f"Removing slit {slit.name} from the list of open slits because the "
                "WCS bounding_box is completely outside the detector."
            )
            idx = np.nonzero([s.name == slit.name for s in open_slits])[0][0]
            open_slits.pop(idx)

    return open_slits


def spectral_order_wrange_from_model(input_model):
    """
    Return the spectral order and wavelength range used in the WCS.

    Parameters
    ----------
    input_model : JwstDataModel
        The data model. Must have been through the assign_wcs step.

    Returns
    -------
    spectral_order : int
        The spectral order.
    wrange : list
        The wavelength range.
    """
    wrange = [input_model.meta.wcsinfo.waverange_start, input_model.meta.wcsinfo.waverange_end]
    spectral_order = input_model.meta.wcsinfo.spectral_order
    return spectral_order, wrange


def nrs_ifu_wcs(input_model):
    """
    Return a list of WCSs for all NIRSPEC IFU slits.

    Parameters
    ----------
    input_model : JwstDataModel
        The data model. Must have passed through the assign_wcs step.

    Returns
    -------
    wcs_list : list
        A list of WCSs for all IFU slits.
    """
    wcs_list = []
    # loop over all IFU slits
    for i in range(30):
        wcs_list.append(nrs_wcs_set_input(input_model, i))
    return wcs_list


def apply_slicemap(input_model, replace_wcs=True):
    """
    Create a pixel-to-slice map and apply it to the WCS.

    By default, the WCS is updated in place in the input model and
    a pixel map is attached in the ``regions`` attribute.  To just
    attach a pixel map, set ``replace_wcs`` to False.

    Parameters
    ----------
    input_model : IFUImageModel
        Input NIRSpec IFU model to be updated.
    replace_wcs : bool
        If False, ``input_model.meta.wcs`` is not modified.
    """
    full_wcs = input_model.meta.wcs

    # Fix the slit name input for initial transforms - the value is not relevant.
    dms2sca = _fix_slit_name(full_wcs.get_transform("detector", "sca"), 0)
    sca2gwa = _fix_slit_name(full_wcs.get_transform("sca", "gwa"), 0)

    # Get the slit-based transforms for slit_frame and slicer, to index into.
    gwa2slit = full_wcs.get_transform("gwa", "slit_frame")
    slit2slicer = full_wcs.get_transform("slit_frame", "slicer")

    # Assume the full detector shape for the regions map and the
    # expected number of IFU slices
    image_shape = (2048, 2048)
    nslice = 30
    regions = np.zeros(image_shape, dtype=np.int32)

    # For each IFU slice, get the bounding box and det2slicer transform
    # and make a pixel to slice map
    transforms = {}
    slice_by_x = {}
    for slit_id in range(nslice):
        bb = full_wcs.bounding_box[slit_id]

        # Construct array indices for pixels in this slice
        x, y = grid_from_bounding_box(bb)
        x_int = x.astype(int)
        y_int = y.astype(int)

        # Get valid wavelengths for the slice
        _, _, lam, _ = full_wcs(x_int, y_int, slit_id)
        valid = np.isfinite(lam)

        # Restrict the slice to valid pixels
        x_int = x_int[valid]
        y_int = y_int[valid]

        # Set the pixel map to the slice value, one-indexed
        # (zero means no transform available)
        regions[y_int, x_int] = slit_id + 1

        # Fix the input for det to slicer transforms by pulling out
        # the specific models needed.
        gwa2slit_by_slice = _fix_slit_name(gwa2slit, slit_id)
        slit2slicer_by_slice = _fix_slit_name(slit2slicer, slit_id)
        new_transform = dms2sca | sca2gwa | gwa2slit_by_slice | slit2slicer_by_slice

        # Bind a bounding box to the transform, including only the valid data
        new_bb = ((x_int.min() - 0.5, x_int.max() + 0.5), (y_int.min() - 0.5, y_int.max() + 0.5))
        bind_bounding_box(new_transform, new_bb, order="F")

        # Keep it for later
        transforms[slit_id + 1] = new_transform

        # Inverse based on slicer x value (detector y in science orientation)
        x_value = np.round(np.nanmean(new_transform(x, y)[0]), 5)
        slice_by_x[x_value] = Mapping((0,), n_inputs=3) | Const1D(slit_id + 1)

    # Attach the slice map to the model
    input_model.regions = regions
    if not replace_wcs:
        # Nothing more to do
        return

    # Slice name mapper for detector pixels
    label_mapper = selector.LabelMapperArray(regions)
    inv_mapper = selector.LabelMapperDict(
        inputs=("x_slice", "y_slice", "lam"),
        mapper=slice_by_x,
        inputs_mapping=Mapping((0,), n_inputs=3),
        atol=1e-4,
    )
    label_mapper.inverse = inv_mapper

    # detector to slicer transform, via a slice label mapper
    det2slicer = selector.RegionsSelector(
        ("x", "y"),
        ("x_slice", "y_slice", "lam"),
        label_mapper=label_mapper,
        selector=transforms,
        name="det2slicer",
    )
    det2slicer.inputs = ("x", "y")
    det2slicer.outputs = ("x_slice", "y_slice", "lam")

    # Make a do-nothing transform to allow a top-level bounding box to be set
    # without a deep copy of the det2slicer transform
    input2det = Identity(2)
    input2det.name = "coord2det"
    input2det.inputs = ("x", "y")
    input2det.outputs = ("x", "y")

    # Fix the slit name input for all further transforms - the value is not relevant.
    slicer_idx = full_wcs.available_frames.index("slicer")
    new_pipeline = [("coordinates", input2det), (full_wcs.pipeline[0].frame, det2slicer)]
    for step in full_wcs.pipeline[slicer_idx:]:
        new_transform = _fix_slit_name(step.transform, 0)
        new_pipeline.append((step.frame, new_transform))

    # Make a new WCS with the slice mapper
    new_wcs = gwcs.WCS(new_pipeline)

    # Attach the new WCS to the model
    input_model.meta.wcs = new_wcs


def _create_ifupost_transform(ifupost_slice):
    """
    Create an IFUPOST transform for a specific slice.

    Parameters
    ----------
    ifupost_slice : `jwst.datamodels.properties.ObjectNode`
        IFUPost transform for a specific slice

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        The transform for this slice.
    """
    linear = ifupost_slice.linear
    polyx = ifupost_slice.xpoly
    polyx_dist = ifupost_slice.xpoly_distortion
    polyy = ifupost_slice.ypoly
    polyy_dist = ifupost_slice.ypoly_distortion

    # the chromatic correction is done here
    # the input is Xslicer, Yslicer, lam
    # The wavelength dependent polynomial is
    # expressed as
    # poly_independent(x, y) + poly_dependent(x, y) * lambda
    model_x = (Mapping((0, 1), n_inputs=3) | polyx) + (
        (Mapping((0, 1), n_inputs=3) | polyx_dist) * (Mapping((2,)) | Identity(1))
    )
    model_y = (Mapping((0, 1), n_inputs=3) | polyy) + (
        (Mapping((0, 1), n_inputs=3) | polyy_dist) * (Mapping((2,)) | Identity(1))
    )

    output2poly_mapping = Identity(2, name="ifupost_outmap")
    output2poly_mapping.inverse = Mapping([0, 1, 2, 0, 1, 2])
    input2poly_mapping = Mapping([0, 1, 2, 0, 1, 2], name="ifupost_inmap")
    input2poly_mapping.inverse = Identity(2)

    model_poly = input2poly_mapping | (model_x & model_y) | output2poly_mapping
    model = linear & Identity(1) | model_poly
    return model


def nrs_lamp(input_model, reference_files, slit_y_range):
    """
    Return the appropriate function for lamp data.

    Parameters
    ----------
    input_model : JwstDataModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Required files depend on the mode.
    slit_y_range : tuple
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    lamp_mode = input_model.meta.instrument.lamp_mode
    if isinstance(lamp_mode, str):
        lamp_mode = lamp_mode.lower()
    else:
        lamp_mode = "none"
    if lamp_mode in ["fixedslit", "brightobj"]:
        return slits_wcs(input_model, reference_files, slit_y_range)
    elif lamp_mode == "ifu":
        return ifu(input_model, reference_files, slit_y_range)
    elif lamp_mode == "msaspec":
        return slits_wcs(input_model, reference_files, slit_y_range)
    else:
        return not_implemented_mode(input_model, reference_files, slit_y_range)


exp_type2transform = {
    "nrs_autoflat": slits_wcs,
    "nrs_autowave": nrs_lamp,
    "nrs_brightobj": slits_wcs,
    "nrs_confirm": imaging,
    "nrs_dark": not_implemented_mode,
    "nrs_fixedslit": slits_wcs,
    "nrs_focus": imaging,
    "nrs_ifu": ifu,
    "nrs_image": imaging,
    "nrs_lamp": nrs_lamp,
    "nrs_mimf": imaging,
    "nrs_msaspec": slits_wcs,
    "nrs_msata": imaging,
    "nrs_taconfirm": imaging,
    "nrs_tacq": imaging,
    "nrs_taslit": imaging,
    "nrs_verify": imaging,
    "nrs_wata": imaging,
}
