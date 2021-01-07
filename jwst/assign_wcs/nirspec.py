"""
Tools to create the WCS pipeline NIRSPEC modes.

Calls create_pipeline() which redirects based on EXP_TYPE.

"""
import logging
import numpy as np

from astropy.modeling import models
from astropy.modeling.models import Mapping, Identity, Const1D, Scale, Tabular1D
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import fits
from gwcs import coordinate_frames as cf

from ..transforms.models import (Rotation3DToGWA, DirCos2Unitless, Slit2Msa,
                                 AngleFromGratingEquation, WavelengthFromGratingEquation,
                                 Gwa2Slit, Unitless2DirCos, Logical, Slit, Snell,
                                 RefractionIndexFromPrism)

from .util import (
    MSAFileError,
    NoDataOnDetectorError,
    not_implemented_mode,
    velocity_correction
)
from . import pointing
from ..datamodels import (CollimatorModel, CameraModel, DisperserModel, FOREModel,
                          IFUFOREModel, MSAModel, OTEModel, IFUPostModel, IFUSlicerModel,
                          WavelengthrangeModel, FPAModel)
from ..lib.exposure_types import is_nrs_ifu_lamp

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging", "ifu", "slits_wcs", "get_open_slits", "nrs_wcs_set_input",
           "nrs_ifu_wcs", "get_spectral_order_wrange"]


def create_pipeline(input_model, reference_files, slit_y_range):
    """
    Create a pipeline list based on EXP_TYPE.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel`, `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.CubeModel`
        The input exposure.
    reference_files : dict
        {reftype: reference_file_name} mapping.
    slit_y_range : list
        The slit Y-range for Nirspec slits, relative to (0, 0) in the center.
    """
    exp_type = input_model.meta.exposure.type.lower()
    if input_model.meta.instrument.grating.lower() == "mirror":
        pipeline = imaging(input_model, reference_files)
    else:
        pipeline = exp_type2transform[exp_type](input_model, reference_files, slit_y_range=slit_y_range)
    if pipeline:
        log.info("Created a NIRSPEC {0} pipeline with references {1}".format(
                exp_type, reference_files))
    return pipeline


def imaging(input_model, reference_files):
    """
    Imaging pipeline.

    It has the following coordinate frames:
    "detector" : the science frame
    "sca" : frame associated with the SCA
    "gwa" " just before the GWA going from detector to sky
    "msa_frame" : at the MSA
    "oteip" : after the FWA
    "v2v3" and "world"

    """
    # Get the corrected disperser model
    disperser = get_disperser(input_model, reference_files['disperser'])

    # DMS to SCA transform
    dms2detector = dms_to_sca(input_model)

    # DETECTOR to GWA transform
    det2gwa = detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)

    gwa_through = Const1D(-1) * Identity(1) & Const1D(-1) * Identity(1) & Identity(1)

    angles = [disperser['theta_x'], disperser['theta_y'],
               disperser['theta_z'], disperser['tilt_y']]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotation').inverse
    dircos2unitless = DirCos2Unitless(name='directional_cosines2unitless')

    col_model = CollimatorModel(reference_files['collimator'])
    col = col_model.model
    col_model.close()

    # Get the default spectral order and wavelength range and record them in the model.
    sporder, wrange = get_spectral_order_wrange(input_model, reference_files['wavelengthrange'])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = sporder

    lam = wrange[0] + (wrange[1] - wrange[0]) * .5

    lam_model = Mapping((0, 1, 1)) | Identity(2) & Const1D(lam)

    gwa2msa = gwa_through | rotation | dircos2unitless | col | lam_model
    gwa2msa.inverse = col.inverse | dircos2unitless.inverse | rotation.inverse | gwa_through

    # Create coordinate frames in the NIRSPEC WCS pipeline
    # "detector", "gwa", "msa", "oteip", "v2v3", "world"
    det, sca, gwa, msa_frame, oteip, v2v3, world = create_imaging_frames()
    if input_model.meta.instrument.filter != 'OPAQUE':
        # MSA to OTEIP transform
        msa2ote = msa_to_oteip(reference_files)
        msa2oteip = msa2ote | Mapping((0, 1), n_inputs=3)
        map1 = Mapping((0, 1, 0, 1))
        minv = msa2ote.inverse
        del minv.inverse
        msa2oteip.inverse = map1 | minv | Mapping((0, 1), n_inputs=3)

        # OTEIP to V2,V3 transform
        with OTEModel(reference_files['ote']) as f:
            oteip2v23 = f.model

        # V2, V3 to world (RA, DEC) transform
        tel2sky = pointing.v23tosky(input_model)

        imaging_pipeline = [(det, dms2detector),
                            (sca, det2gwa),
                            (gwa, gwa2msa),
                            (msa_frame, msa2oteip),
                            (oteip, oteip2v23),
                            (v2v3, tel2sky),
                            (world, None)]
    else:
        # convert to microns if the pipeline ends earlier
        gwa2msa = (gwa2msa | Identity(2) & Scale(10**6)).rename('gwa2msa')
        imaging_pipeline = [(det, dms2detector),
                            (sca, det2gwa),
                            (gwa, gwa2msa),
                            (msa_frame, None)]

    return imaging_pipeline


def ifu(input_model, reference_files, slit_y_range=[-.55, .55]):
    """
    The Nirspec IFU WCS pipeline.

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
    input_model : `~jwst.datamodels.DataModel`
        The input data model.
    reference_files : dict
        The reference files used for this mode.
    slit_y_range : list
        The slit dimensions relative to the center of the slit.
    """
    detector = input_model.meta.instrument.detector
    grating = input_model.meta.instrument.grating
    filter = input_model.meta.instrument.filter
    # Check for ifu reference files
    if reference_files['ifufore'] is None and \
       reference_files['ifuslicer'] is None and \
       reference_files['ifupost'] is None:
        # No ifu reference files, won't be able to create pipeline
        log_message = 'No ifufore, ifuslicer or ifupost reference files'
        log.critical(log_message)
        raise RuntimeError(log_message)
    # Check for data actually being present on NRS2
    log_message = "No IFU slices fall on detector {0}".format(detector)
    if detector == "NRS2" and grating.endswith('M'):
        # Mid-resolution gratings do not project on NRS2.
        log.critical(log_message)
        raise NoDataOnDetectorError(log_message)
    if detector == "NRS2" and grating == "G140H" and filter == "F070LP":
        # This combination of grating and filter does not project on NRS2.
        log.critical(log_message)
        raise NoDataOnDetectorError(log_message)

    slits = np.arange(30)
    # Get the corrected disperser model
    disperser = get_disperser(input_model, reference_files['disperser'])

    # Get the default spectral order and wavelength range and record them in the model.
    sporder, wrange = get_spectral_order_wrange(input_model, reference_files['wavelengthrange'])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = sporder

    # DMS to SCA transform
    dms2detector = dms_to_sca(input_model)
    # DETECTOR to GWA transform
    det2gwa = Identity(2) & detector_to_gwa(reference_files,
                                            input_model.meta.instrument.detector,
                                            disperser)

    # GWA to SLIT
    gwa2slit = gwa_to_ifuslit(slits, input_model, disperser, reference_files, slit_y_range)

    # SLIT to MSA transform
    slit2slicer = ifuslit_to_slicer(slits, reference_files, input_model)

    # SLICER to MSA Entrance
    slicer2msa = slicer_to_msa(reference_files)

    det, sca, gwa, slit_frame, msa_frame, oteip, v2v3, world = create_frames()

    exp_type = input_model.meta.exposure.type.upper()

    is_lamp_exposure = exp_type in ['NRS_LAMP', 'NRS_AUTOWAVE', 'NRS_AUTOFLAT']

    if input_model.meta.instrument.filter == 'OPAQUE' or is_lamp_exposure:
        # If filter is "OPAQUE" or if internal lamp exposure the NIRSPEC WCS pipeline stops at the MSA.
        pipeline = [(det, dms2detector),
                    (sca, det2gwa.rename('detector2gwa')),
                    (gwa, gwa2slit.rename('gwa2slit')),
                    (slit_frame, slit2slicer),
                    ('slicer', slicer2msa),
                    (msa_frame, None)]
    else:
        # MSA to OTEIP transform
        msa2oteip = ifu_msa_to_oteip(reference_files)
        # OTEIP to V2,V3 transform
        # This includes a wavelength unit conversion from meters to microns.
        oteip2v23 = oteip_to_v23(reference_files)

        # V2, V3 to sky
        tel2sky = pointing.v23tosky(input_model) & Identity(1)

        # Create coordinate frames in the NIRSPEC WCS pipeline"
        #
        # The oteip2v2v3 transform converts the wavelength from meters (which is assumed
        # in the whole pipeline) to microns (which is the expected output)
        #
        # "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world"

        pipeline = [(det, dms2detector),
                    (sca, det2gwa.rename('detector2gwa')),
                    (gwa, gwa2slit.rename('gwa2slit')),
                    (slit_frame, slit2slicer),
                    ('slicer', slicer2msa),
                    (msa_frame, msa2oteip.rename('msa2oteip')),
                    (oteip, oteip2v23.rename('oteip2v23')),
                    (v2v3, tel2sky),
                    (world, None)]

    return pipeline


def slits_wcs(input_model, reference_files, slit_y_range):
    """
    The WCS pipeline for MOS and fixed slits.

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
    input_model : `~jwst.datamodels.DataModel`
        The input data model.
    reference_files : dict
        The reference files used for this mode.
    slit_y_range : list
        The slit dimensions relative to the center of the slit.
    """
    open_slits_id = get_open_slits(input_model, reference_files, slit_y_range)
    if not open_slits_id:
        return None
    n_slits = len(open_slits_id)
    log.info("Computing WCS for {0} open slitlets".format(n_slits))

    msa_pipeline = slitlets_wcs(input_model, reference_files, open_slits_id)

    return msa_pipeline


def slitlets_wcs(input_model, reference_files, open_slits_id):
    """
    Create The WCS piepline for MOS and Fixed slits for the
    specific opened shutters/slits. ``slit_y_range`` is taken from
    ``slit.ymin`` and ``slit.ymax``.

    Note: This function is also used by the ``msaflagopen`` step.
    """
    # Get the corrected disperser model
    disperser = get_disperser(input_model, reference_files['disperser'])

    # Get the default spectral order and wavelength range and record them in the model.
    sporder, wrange = get_spectral_order_wrange(input_model, reference_files['wavelengthrange'])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    log.info("SPORDER= {0}, wrange={1}".format(sporder, wrange))
    input_model.meta.wcsinfo.spectral_order = sporder

    # DMS to SCA transform
    dms2detector = dms_to_sca(input_model)
    dms2detector.name = 'dms2sca'
    # DETECTOR to GWA transform
    det2gwa = Identity(2) & detector_to_gwa(reference_files,
                                            input_model.meta.instrument.detector,
                                            disperser)
    det2gwa.name = "det2gwa"

    # GWA to SLIT
    gwa2slit = gwa_to_slit(open_slits_id, input_model, disperser, reference_files)
    gwa2slit.name = "gwa2slit"

    # SLIT to MSA transform
    slit2msa = slit_to_msa(open_slits_id, reference_files['msa'])
    slit2msa.name = "slit2msa"

    # Create coordinate frames in the NIRSPEC WCS pipeline"
    # "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world"
    det, sca, gwa, slit_frame, msa_frame, oteip, v2v3, world = create_frames()

    exp_type = input_model.meta.exposure.type.upper()

    is_lamp_exposure = exp_type in ['NRS_LAMP', 'NRS_AUTOWAVE', 'NRS_AUTOFLAT']

    if input_model.meta.instrument.filter == 'OPAQUE' or is_lamp_exposure:
        # convert to microns if the pipeline ends earlier
        msa_pipeline = [(det, dms2detector),
                        (sca, det2gwa),
                        (gwa, gwa2slit),
                        (slit_frame, slit2msa),
                        (msa_frame, None)]
    else:
        # MSA to OTEIP transform
        msa2oteip = msa_to_oteip(reference_files)
        msa2oteip.name = "msa2oteip"

        # OTEIP to V2,V3 transform
        # This includes a wavelength unit conversion from meters to microns.
        oteip2v23 = oteip_to_v23(reference_files)
        oteip2v23.name = "oteip2v23"

        # V2, V3 to sky
        tel2sky = pointing.v23tosky(input_model) & Identity(1)
        tel2sky.name = "v2v3_to_sky"

        msa_pipeline = [(det, dms2detector),
                        (sca, det2gwa),
                        (gwa, gwa2slit),
                        (slit_frame, slit2msa),
                        (msa_frame, msa2oteip),
                        (oteip, oteip2v23),
                        (v2v3, tel2sky),
                        (world, None)]

    return msa_pipeline


def get_open_slits(input_model, reference_files=None, slit_y_range=[-.55, .55]):
    """Return the opened slits/shutters in a MOS or Fixed Slits exposure.
    """
    exp_type = input_model.meta.exposure.type.lower()
    lamp_mode = input_model.meta.instrument.lamp_mode
    if type(lamp_mode) == str:
        lamp_mode = lamp_mode.lower()
    else:
        lamp_mode = 'none'
    if exp_type in ["nrs_msaspec", "nrs_autoflat"] or ((exp_type in ["nrs_lamp", "nrs_autowave"]) and \
                                                        (lamp_mode == "msaspec")):
        msa_metadata_file, msa_metadata_id, dither_point = get_msa_metadata(
            input_model, reference_files)
        slits = get_open_msa_slits(msa_metadata_file, msa_metadata_id, dither_point, slit_y_range)
    elif exp_type == "nrs_fixedslit":
        slits = get_open_fixed_slits(input_model, slit_y_range)
    elif exp_type == "nrs_brightobj":
        slits = [Slit('S1600A1', 3, 0, 0, 0, slit_y_range[0], slit_y_range[1], 5, 4)]
    elif exp_type in ["nrs_lamp", "nrs_autowave"]:
        if lamp_mode in ['fixedslit', 'brightobj']:
            slits = get_open_fixed_slits(input_model, slit_y_range)
    else:
        raise ValueError("EXP_TYPE {0} is not supported".format(exp_type.upper()))
    if reference_files is not None:
        slits = validate_open_slits(input_model, slits, reference_files)
        log.info("Slits projected on detector {0}: {1}".format(input_model.meta.instrument.detector,
                                                               [sl.name for sl in slits]))
    if not slits:
        log_message = "No open slits fall on detector {0}.".format(input_model.meta.instrument.detector)
        log.critical(log_message)
        raise NoDataOnDetectorError(log_message)
    return slits


def get_open_fixed_slits(input_model, slit_y_range=[-.55, .55]):
    """ Return the opened fixed slits."""
    if input_model.meta.subarray.name is None:
        raise ValueError("Input file is missing SUBARRAY value/keyword.")

    slits = []
    ylow, yhigh = slit_y_range

    s2a1 = Slit('S200A1', 0, 0, 0, 0, ylow, yhigh, 5, 1)
    s2a2 = Slit('S200A2', 1, 0, 0, 0, ylow, yhigh, 5, 2)
    s4a1 = Slit('S400A1', 2, 0, 0, 0, ylow, yhigh, 5, 3)
    s16a1 = Slit('S1600A1', 3, 0, 0, 0, ylow, yhigh, 5, 4)
    s2b1 = Slit('S200B1', 4, 0, 0, 0, ylow, yhigh, 5, 5)

    subarray = input_model.meta.subarray.name.upper()
    if subarray == "SUBS200A1":
        slits.append(s2a1)
    elif subarray == "SUBS200A2":
        slits.append(s2a2)
    elif subarray == "SUBS400A1":
        slits.append(s4a1)
    elif subarray in ("SUB2048", "SUB512", "SUB512S",
                      "SUB1024A", "SUB1024B"):
        slits.append(s16a1)
    elif subarray == "SUBS200B1":
        slits.append(s2b1)
    else:
        slits.extend([s2a1, s2a2, s4a1, s16a1, s2b1])

    return slits


def get_msa_metadata(input_model, reference_files):
    """
    Get the MSA metadata file (MSAMTFL) and the msa metadata ID (MSAMETID).

    """
    try:
        msa_config = reference_files['msametafile']
    except (KeyError, TypeError):
        log.info('MSA metadata file not in reference files dict')
        log.info('Getting MSA metadata file from MSAMETFL keyword')
        msa_config = input_model.meta.instrument.msa_metadata_file
        if msa_config is None:
            message = "msa_metadata_file is None."
            log.critical(message)
            raise MSAFileError(message)
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


def _get_bkg_source_id(bkg_counter, shift_by):
    """
    Compute a ``source_id`` for background slitlets.

    All background slitlets are assigned a source_id of 0.
    A unique ``source_id`` is necessary to keep them separate in exp_to_source.
    A counter is used to assign a unique ``source_id`` that's
    greater than the max ID number of all defined sources.

    Parameters
    ----------
    bkg_counter : int
        The current value of the counter.
    shift_by : int
        The highest of all source_id values.
    """

    return bkg_counter + shift_by


def get_open_msa_slits(msa_file, msa_metadata_id, dither_position,
                       slit_y_range=[-.55, .55]):
    """
    Return the opened MOS slitlets.

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

       column

    Parameters
    ----------
    msa_file : str
        MSA meta data file name, FITS keyword ``MSAMETFL``.
    msa_metadata_id : int
        The MSA meta id for the science file, FITS keyword ``MSAMETID``.
    dither_position : int
        The index in the dither pattern, FITS keyword ``PATT_NUM``.
    slit_y_range : list or tuple of size 2
        The lower and upper limit of the slit.

    Returns
    -------
    slitlets : list
        A list of `~jwst.transforms.models.Slit` objects. Each slitlet is a tuple with
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
        message = "Missing MSA meta (MSAMETFL) file {}".format(msa_file)
        log.error(message)
        raise MSAFileError(message)
    except OSError:
        message = "Unable to read MSA FITS file (MSAMETFL) {0}".format(msa_file)
        log.error(message)
        raise MSAFileError(message)
    except Exception:
        message = "Problem reading MSA metafile (MSAMETFL) {0}".format(msa_file)
        log.error(message)
        raise MSAFileError(message)

    # Get the configuration header from teh _msa.fits file.  The EXTNAME should be 'SHUTTER_INFO'
    msa_conf = msa_file[('SHUTTER_INFO', 1)]
    msa_source = msa_file[("SOURCE_INFO", 1)].data

    # First we are going to filter the msa_file data on the msa_metadata_id
    # and dither_point_index.
    msa_data = [x for x in msa_conf.data if x['msa_metadata_id'] == msa_metadata_id \
                and x['dither_point_index'] == dither_position]

    # Get all source_ids for slitlets with sources.
    # These should not be used when assigning source_id to background slitlets.
    source_ids = set([x[5] for x in msa_conf.data if x['msa_metadata_id'] == msa_metadata_id \
                      and x['dither_point_index'] == dither_position])
    # All BKG shutters in the msa metafile have a source_id value of 0.
    # Remove it from the list of source ids.
    if 0 in source_ids:
        source_ids.remove(0)
    if source_ids:
        max_source_id = max(source_ids) + 1
    else:
        max_source_id = 0

    # define a counter for "all background" slitlets.
    # It will be used to assign a "source_id".
    bkg_counter = 0

    log.debug(f'msa_data with msa_metadata_id = {msa_metadata_id}   {msa_data}')
    log.info(f'Retrieving open MSA slitlets for msa_metadata_id = {msa_metadata_id} '
             f'and dither_index = {dither_position}')

    # Get the unique slitlet_ids
    slitlet_ids_unique = list(set([x['slitlet_id'] for x in msa_data]))

    # SDP may assign a value of "-1" to ``slitlet_id`` - these need to be ignored.
    # JP-436
    if -1 in slitlet_ids_unique:
        slitlet_ids_unique.remove(-1)

    # add a margin to the slit y limits
    margin = 0.05

    # Now lets look at each unique slitlet id
    for slitlet_id in slitlet_ids_unique:
        # Get the rows for the current slitlet_id
        slitlets_sid = [x for x in msa_data if x['slitlet_id'] == slitlet_id]
        open_shutters = [x['shutter_column'] for x in slitlets_sid]

        n_main_shutter = len([s for s in slitlets_sid if s['primary_source'] == 'Y'])
        # In the next part we need to calculate, find, determine 5 things:
        #    quadrant,  xcen, ycen,  ymin, ymax

        # There are no main shutters, all are background
        if n_main_shutter == 0:
            if len(open_shutters) == 1:
                jmin = jmax = j = open_shutters[0]
            else:
                jmin = min([s['shutter_column'] for s in slitlets_sid])
                jmax = max([s['shutter_column'] for s in slitlets_sid])
                j = jmin + (jmax - jmin) // 2 + 1
            ymax = 0.5 + margin + (jmax - j) * 1.15
            ymin = -(-ylow + margin) + (jmin - j) * 1.15
            quadrant = slitlets_sid[0]['shutter_quadrant']
            ycen = j
            xcen = slitlets_sid[0]['shutter_row']  # grab the first as they are all the same
            source_xpos = 0.0
            source_ypos = 0.0
            source_id = _get_bkg_source_id(bkg_counter, max_source_id)
            log.info(f'Slitlet_id {slitlet_id} is background only; assigned source_id = {source_id}')
            bkg_counter += 1
        # There is 1 main shutter, phew, that makes it easier.
        elif n_main_shutter == 1:
            xcen, ycen, quadrant, source_xpos, source_ypos = [
                (s['shutter_row'], s['shutter_column'], s['shutter_quadrant'],
                 s['estimated_source_in_shutter_x'],
                 s['estimated_source_in_shutter_y'])
                for s in slitlets_sid if s['background'] == 'N'][0]

            # y-size
            jmin = min([s['shutter_column'] for s in slitlets_sid])
            jmax = max([s['shutter_column'] for s in slitlets_sid])
            j = ycen
            ymax = yhigh + margin + (jmax - j) * 1.15
            ymin = -(-ylow + margin) + (jmin - j) * 1.15
            source_id = slitlets_sid[0]['source_id']
        # Not allowed....
        else:
            message = ("For slitlet_id = {}, metadata_id = {}, "
                       "and dither_index = {}".format(slitlet_id, msa_metadata_id, dither_position))
            log.warning(message)
            message = ("MSA configuration file has more than 1 shutter with primary source")
            log.warning(message)
            raise MSAFileError(message)

        # subtract 1 because shutter numbers in the MSA reference file are 1-based.
        shutter_id = xcen + (ycen - 1) * 365
        try:
            source_name, source_alias, stellarity, source_ra, source_dec = [
                (s['source_name'], s['alias'], s['stellarity'], s['ra'], s['dec']) \
                for s in msa_source if s['source_id'] == source_id][0]
        except IndexError:
            # all background shutters
            source_name = "background_{}".format(slitlet_id)
            source_alias = "bkg_{}".format(slitlet_id)
            stellarity = 0.0
            source_ra = 0.0
            source_dec = 0.0

        # Create the output list of tuples that contain the required
        # data for further computations
        """
        Convert source positions from PPS to Model coordinate frame.
        The source x,y position in the shutter is given in the msa configuration file,
        columns "estimated_source_in_shutter_x" and "estimated_source_in_shutter_y".
        The source position is in a coordinate system associated with each shutter whose
        origin is the lower left corner of the shutter, positive x is to the right
        and positive y is upwards.
        """
        source_xpos -= 0.5
        source_ypos -= 0.5

        # Create the shutter_state string
        all_shutters = _shutter_id_to_str(open_shutters, ycen)

        slitlets.append(Slit(slitlet_id, shutter_id, dither_position, xcen, ycen, ymin, ymax,
                             quadrant, source_id, all_shutters, source_name, source_alias,
                             stellarity, source_xpos, source_ypos, source_ra, source_dec))
    msa_file.close()
    return slitlets


def _shutter_id_to_str(open_shutters, ycen):
    """
    Return a string representing the open and closed shutters in a slitlet.

    Parameters
    ----------
    open_shutters : list
        List of IDs (shutter_id) of open shutters.
    xcen : int
        X coordinate of main shutter.

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
    all_shutters = all_shutters.astype(np.str)
    all_shutters[cen_ind] = 'x'
    return "".join(all_shutters)


def get_spectral_order_wrange(input_model, wavelengthrange_file):
    """
    Read the spectral order and wavelength range from the reference file.

    Parameters
    ----------
    input_model : `~jwst.datamodels.DataModel`
        The input data model.
    wavelengthrange_file : str
        Reference file of type "wavelengthrange".
    """
    # Nirspec full spectral range
    full_range = [.6e-6, 5.3e-6]

    filter = input_model.meta.instrument.filter
    lamp = input_model.meta.instrument.lamp_state
    grating = input_model.meta.instrument.grating
    exp_type = input_model.meta.exposure.type

    is_lamp_exposure = exp_type in ['NRS_LAMP', 'NRS_AUTOWAVE', 'NRS_AUTOFLAT']

    wave_range_model = WavelengthrangeModel(wavelengthrange_file)
    wrange_selector = wave_range_model.waverange_selector
    if filter == "OPAQUE" or is_lamp_exposure:
        keyword = lamp + '_' + grating
    else:
        keyword = filter + '_' + grating
    try:
        index = wrange_selector.index(keyword)
    except (KeyError, ValueError):
        # Combination of filter_grating is not in wavelengthrange file.
        gratings = [s.split('_')[1] for s in wrange_selector]
        try:
            index = gratings.index(grating)
        except ValueError: # grating not in list
            order = -1
            wrange = full_range
        else:
            order = wave_range_model.order[index]
            wrange = wave_range_model.wavelengthrange[index]
        log.info("Combination {0} missing in wavelengthrange file, setting "
                 "order to {1} and range to {2}.".format(keyword, order, wrange))
    else:
        # Combination of filter_grating is found in wavelengthrange file.
        order = wave_range_model.order[index]
        wrange = wave_range_model.wavelengthrange[index]

    wave_range_model.close()
    return order, wrange


def ifuslit_to_slicer(slits, reference_files, input_model):
    """
    The transform from ``slit_frame`` to ``slicer`` frame.

    Parameters
    ----------
    slits : list
        A list of slit IDs for all slices.
    reference_files : dict
        {reference_type: reference_file_name}
    input_model : `~jwst.datamodels.IFUImageModel`

    Returns
    -------
    model : `~jwst.transforms.Slit2Msa` model.
        Transform from ``slit_frame`` to ``slicer`` frame.
    """
    ifuslicer = IFUSlicerModel(reference_files['ifuslicer'])
    models = []
    ifuslicer_model = ifuslicer.model
    for slit in slits:
        slitdata = ifuslicer.data[slit]
        slitdata_model = (get_slit_location_model(slitdata)).rename('slitdata_model')
        slicer_model = slitdata_model | ifuslicer_model

        msa_transform = slicer_model
        models.append(msa_transform)
    ifuslicer.close()

    return Slit2Msa(slits, models)


def slicer_to_msa(reference_files):
    """
    Trasform from slicer coordinates to MSA entrance.

    Applies the IFUFORE transform.

    """
    with IFUFOREModel(reference_files['ifufore']) as f:
        ifufore = f.model
    slicer2fore_mapping = Mapping((0, 1, 2, 2))
    slicer2fore_mapping.inverse = Identity(3)
    ifufore2fore_mapping = Identity(1)
    ifufore2fore_mapping.inverse = Mapping((0, 1, 2, 2))
    ifu_fore_transform = slicer2fore_mapping | ifufore & Identity(1)
    return ifu_fore_transform


def slit_to_msa(open_slits, msafile):
    """
    The transform from ``slit_frame`` to ``msa_frame``.

    Parameters
    ----------
    open_slits : list
        A list of slit IDs for all open shutters/slitlets.
    msafile : str
        The name of the msa reference file.

    Returns
    -------
    model : `~jwst.transforms.Slit2Msa` model.
        Transform from ``slit_frame`` to ``msa_frame``.
    """
    msa = MSAModel(msafile)
    models = []
    slits = []
    for quadrant in range(1, 6):
        slits_in_quadrant = [s for s in open_slits if s.quadrant == quadrant]
        msa_quadrant = getattr(msa, 'Q{0}'.format(quadrant))
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
    return Slit2Msa(slits, models)


def gwa_to_ifuslit(slits, input_model, disperser, reference_files, slit_y_range):
    """
    The transform from ``gwa`` to ``slit_frame``.

    Parameters
    ----------
    slits : list
        A list of slit IDs for all IFU slits 0-29.
    disperser : `~jwst.datamodels.DisperserModel`
        A disperser model with the GWA correction applied to it.
    filter : str
        The filter used.
    grating : str
        The grating used in the observation.
    reference_files: dict
        Dictionary with reference files returned by CRDS.
    slit_y_range : list or tuple of size 2
        The lower and upper bounds of a slit.

    Returns
    -------
    model : `~jwst.transforms.Gwa2Slit` model.
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
            log.info("Applied Barycentric velocity correction : {}".format(velocity_corr[1].amplitude.value))
    # The wavelength units up to this point are
    # meters as required by the pipeline but the desired output wavelength units is microns.
    # So we are going to Scale the spectral units by 1e6 (meters -> microns)
    is_lamp_exposure = input_model.meta.exposure.type in ['NRS_LAMP', 'NRS_AUTOWAVE', 'NRS_AUTOFLAT']
    if input_model.meta.instrument.filter == 'OPAQUE' or is_lamp_exposure:
        lgreq = lgreq | Scale(1e6)

    lam_cen = 0.5 * (input_model.meta.wcsinfo.waverange_end -
                     input_model.meta.wcsinfo.waverange_start
                     ) + input_model.meta.wcsinfo.waverange_start
    collimator2gwa = collimator_to_gwa(reference_files, disperser)
    mask = mask_slit(ymin, ymax)

    ifuslicer = IFUSlicerModel(reference_files['ifuslicer'])
    ifupost = IFUPostModel(reference_files['ifupost'])
    slit_models = []
    ifuslicer_model = ifuslicer.model
    for slit in slits:
        slitdata = ifuslicer.data[slit]
        slitdata_model = get_slit_location_model(slitdata)
        ifuslicer_transform = (slitdata_model | ifuslicer_model)
        ifupost_sl = getattr(ifupost, "slice_{0}".format(slit))
        # construct IFU post transform
        ifupost_transform = _create_ifupost_transform(ifupost_sl)
        msa2gwa = ifuslicer_transform & Const1D(lam_cen) | ifupost_transform | collimator2gwa
        gwa2slit = gwa_to_ymsa(msa2gwa, lam_cen=lam_cen, slit_y_range=slit_y_range)# TODO: Use model sets here

        # The commnts below list the input coordinates.
        bgwa2msa = (
            # (alpha_out, beta_out, gamma_out), angles at the GWA, coming from the camera
            # (0, - beta_out, alpha_out, beta_out)
            # (0, sy, alpha_out, beta_out)
            # (0, sy, 0, sy, sy, alpha_out, beta_out)
            # ( 0, sy, alpha_in, beta_in, gamma_in, alpha_out, beta_out)
            # (0, sy, alpha_in, beta_in,alpha_out)
            # (0, sy, lambda_computed)
            Mapping((0, 1, 0, 1), n_inputs=3) |
            Const1D(0) * Identity(1) & Const1D(-1) * Identity(1) & Identity(2) | \
            Identity(1) & gwa2slit & Identity(2) | \
            Mapping((0, 1, 0, 1, 1, 2, 3)) | \
            Identity(2) & msa2gwa & Identity(2) | \
            Mapping((0, 1, 2, 3, 5), n_inputs=7) | \
            Identity(2) & lgreq | mask
        )

        # transform from ``msa_frame`` to ``gwa`` frame (before the GWA going from detector to sky).
        msa2gwa_out = ifuslicer_transform & Identity(1) | ifupost_transform | collimator2gwa
        msa2bgwa = Mapping((0, 1, 2, 2)) | msa2gwa_out & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
        bgwa2msa.inverse = msa2bgwa
        slit_models.append(bgwa2msa)

    ifuslicer.close()
    ifupost.close()
    return Gwa2Slit(slits, slit_models)


def gwa_to_slit(open_slits, input_model, disperser,
                reference_files):
    """
    The transform from ``gwa`` to ``slit_frame``.

    Parameters
    ----------
    open_slits : list
        A list of slit IDs for all open shutters/slitlets.
    disperser : dict
        A corrected disperser ASDF object.
    filter : str
        The filter used.
    grating : str
        The grating used in the observation.
    reference_files: dict
        Dictionary with reference files returned by CRDS.

    Returns
    -------
    model : `~jwst.transforms.Gwa2Slit` model.
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
            log.info("Applied Barycentric velocity correction : {}".format(velocity_corr[1].amplitude.value))

    # The wavelength units up to this point are
    # meters as required by the pipeline but the desired output wavelength units is microns.
    # So we are going to Scale the spectral units by 1e6 (meters -> microns)
    is_lamp_exposure = input_model.meta.exposure.type in ['NRS_LAMP', 'NRS_AUTOWAVE', 'NRS_AUTOFLAT']
    if input_model.meta.instrument.filter == 'OPAQUE' or is_lamp_exposure:
        lgreq = lgreq | Scale(1e6)

    msa = MSAModel(reference_files['msa'])
    slit_models = []
    slits = []
    for quadrant in range(1, 6):
        slits_in_quadrant = [s for s in open_slits if s.quadrant == quadrant]
        log.info("There are {0} open slits in quadrant {1}".format(len(slits_in_quadrant), quadrant))
        msa_quadrant = getattr(msa, 'Q{0}'.format(quadrant))

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
                msa_transform = (slitdata_model | msa_model)
                msa2gwa = (msa_transform | collimator2gwa)
                gwa2msa = gwa_to_ymsa(msa2gwa, slit=slit, slit_y_range=(slit.ymin, slit.ymax))# TODO: Use model sets here
                bgwa2msa = Mapping((0, 1, 0, 1), n_inputs=3) | \
                    Const1D(0) * Identity(1) & Const1D(-1) * Identity(1) & Identity(2) | \
                    Identity(1) & gwa2msa & Identity(2) | \
                    Mapping((0, 1, 0, 1, 2, 3)) | Identity(2) & msa2gwa & Identity(2) | \
                    Mapping((0, 1, 2, 3, 5), n_inputs=7) | Identity(2) & lgreq | mask
                    #Mapping((0, 1, 2, 5), n_inputs=7) | Identity(2) & lgreq | mask
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
    For gratings this returns a form of the grating equation
    which computes the angle when lambda is known.

    For prism data this returns the Snell model.
    """
    sporder = input_model.meta.wcsinfo.spectral_order
    if input_model.meta.instrument.grating.lower() != 'prism':
        agreq = AngleFromGratingEquation(disperser.groovedensity,
                                         sporder, name='alpha_from_greq')
        return agreq
    else:
        system_temperature = input_model.meta.instrument.gwa_tilt
        system_pressure = disperser['pref']

        snell = Snell(disperser['angle'], disperser['kcoef'], disperser['lcoef'],
                      disperser['tcoef'], disperser['tref'], disperser['pref'],
                      system_temperature, system_pressure, name="snell_law")
        return snell


def wavelength_from_disperser(disperser, input_model):
    """
    For gratings this returns a form of the grating equation
    which computes lambda when all angles are known.

    For prism data this returns a lookup table model
    computing lambda from a known refraction index.
    """
    sporder = input_model.meta.wcsinfo.spectral_order
    if input_model.meta.instrument.grating.lower() != 'prism':
        lgreq = WavelengthFromGratingEquation(disperser.groovedensity,
                                              sporder, name='lambda_from_gratingeq')
        return lgreq
    else:
        lam = np.arange(0.5, 6.005, 0.005) * 1e-6
        system_temperature = input_model.meta.instrument.gwa_tilt
        if system_temperature is None:
            message = "Missing reference temperature (keyword GWA_TILT)."
            log.critical(message)
            raise KeyError(message)
        system_pressure = disperser['pref']
        tref = disperser['tref']
        pref = disperser['pref']
        kcoef = disperser['kcoef'][:]
        lcoef = disperser['lcoef'][:]
        tcoef = disperser['tcoef'][:]
        n = Snell.compute_refraction_index(lam, system_temperature, tref, pref,
                                           system_pressure, kcoef, lcoef, tcoef
                                           )
        n = np.flipud(n)
        lam = np.flipud(lam)
        n_from_prism = RefractionIndexFromPrism(disperser['angle'], name='n_prism')

        tab = Tabular1D(points=(n,), lookup_table=lam, bounds_error=False)
        return n_from_prism | tab


def detector_to_gwa(reference_files, detector, disperser):
    """
    Transform from ``sca`` frame to ``gwa`` frame.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.
    detector : str
        The detector keyword.
    disperser : dict
        A corrected disperser ASDF object.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from DETECTOR frame to GWA frame.

    """
    with FPAModel(reference_files['fpa']) as f:
        fpa = getattr(f, detector.lower() + '_model')
    with CameraModel(reference_files['camera']) as f:
        camera = f.model

    angles = [disperser['theta_x'], disperser['theta_y'],
               disperser['theta_z'], disperser['tilt_y']]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotation')
    u2dircos = Unitless2DirCos(name='unitless2directional_cosines')
    ## NIRSPEC 1- vs 0- based pixel coordinates issue #1781
    '''
    The pipeline works with 0-based pixel coordinates. The Nirspec model,
    stored in reference files, is also 0-based. However, the algorithm specified
    by the IDT team specifies that pixel coordinates are 1-based. This is
    implemented below as a Shift(-1) & Shift(-1) transform. This makes the Nirspec
    instrument WCS pipeline "special" as it requires 1-based inputs.
    As a consequence many steps have to be modified to provide 1-based coordinates
    to the WCS call if the instrument is Nirspec. This is not always easy, especially
    when the step has no knowledge of the instrument.
    This is the reason the algorithm is modified to acccept 0-based coordinates.
    This will be discussed in the future with the INS and IDT teams and may be solved
    by changing the algorithm but for now

    model = (models.Shift(-1) & models.Shift(-1) | fpa | camera | u2dircos | rotation)

    is changed to

    model = models.Shift(1) & models.Shift(1) | \
            models.Shift(-1) & models.Shift(-1) | fpa | camera | u2dircos | rotation
    '''
    model = fpa | camera | u2dircos | rotation
    return model


def dms_to_sca(input_model):
    """
    Transforms from ``detector`` to ``sca`` coordinates.
    """
    detector = input_model.meta.instrument.detector
    xstart = input_model.meta.subarray.xstart
    ystart = input_model.meta.subarray.ystart
    if xstart is None:
        xstart = 1
    if ystart is None:
        ystart = 1
    # The SCA coordinates are in full frame
    # The inputs are 1-based, remove -1 when'if they are 0-based
    # The outputs must be 1-based because this is what the model expects.
    # If xstart was 0-based and the inputs were 0-based ->
    # Shift(+1)
    subarray2full = models.Shift(xstart - 1) & models.Shift(ystart - 1)
    if detector == 'NRS2':
        model = models.Shift(-2047) & models.Shift(-2047) | models.Scale(-1) & models.Scale(-1)
    elif detector == 'NRS1':
        model = models.Identity(2)
    return subarray2full | model


def mask_slit(ymin=-.55, ymax=.55):
    """
    Returns a model which masks out pixels in a NIRSpec cutout outside the slit.

    Uses ymin, ymax for the slit and the wavelength range to define the location of the slit.

    Parameters
    ----------
    ymin, ymax : float
        ymin and ymax relative boundary of a slit.

    Returns
    -------
    model : `~astropy.modeling.core.Model`
        A model which takes x_slit, y_slit, lam inputs and substitutes the
        values outside the slit with NaN.

    """
    greater_than_ymax = Logical(condition='GT', compareto=ymax, value=np.nan)
    less_than_ymin = Logical(condition='LT', compareto=ymin, value=np.nan)

    model = Mapping((0, 1, 2, 1)) | Identity(3) & (greater_than_ymax | less_than_ymin | models.Scale(0)) | \
          Mapping((0, 1, 3, 2, 3)) | Identity(1) & Mapping((0,), n_inputs=2) + Mapping((1,)) & \
          Mapping((0,), n_inputs=2) + Mapping((1,))
    model.inverse = Identity(3)
    return model


def compute_bounding_box(slit2detector, wavelength_range, slit_ymin=-.55, slit_ymax=.55):
    """
    Compute the bounding box of the projection of a slit/slice on the detector.

    The edges of the slit are used to determine the location
    of the projection of the slit on the detector.
    Because the trace is curved and the wavelength_range may span the
    two detectors, y_min of the projection may be at an arbitrary wavelength.
    The transform is run with a regularly sampled wavelengths to determin y_min.

    Parameters
    ----------
    slit2detector : `astropy.modeling.core.Model`
        The transform from slit to detector.
    wavelength_range : tuple
        The wavelength range for the combination of grating and filter.

    """
    lam_min, lam_max = wavelength_range
    step = 1e-10
    nsteps = int((lam_max - lam_min) / step)
    lam_grid = np.linspace(lam_min, lam_max, nsteps)
    x_range_low, y_range_low = slit2detector([0] * nsteps, [slit_ymin] * nsteps, lam_grid)
    x_range_high, y_range_high = slit2detector([0] * nsteps, [slit_ymax] * nsteps, lam_grid)
    x_range = np.hstack((x_range_low, x_range_high))

    y_range = np.hstack((y_range_low, y_range_high))
    # add 10 px margin
    # The -1 is technically because the output of slit2detector is 1-based coordinates.
    x0 = max(0, x_range.min() - 1 - 10)
    x1 = min(2047, x_range.max() - 1 + 10)
    # add 2 px margin
    y0 = max(0, y_range.min() - 1 - 2)
    y1 = min(2047, y_range.max() - 1 + 2)

    bounding_box = ((x0 - 0.5, x1 + 0.5), (y0 - 0.5, y1 + 0.5))
    return bounding_box


def collimator_to_gwa(reference_files, disperser):
    """
    Transform from collimator to ``gwa`` frame.

    Includes the transforms:
    - through the collimator (going from sky to detector)
    - converting from unitless to directional cosines
    - a 3D rotation before the GWA using th ecorrected disperser angles.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.
    disperser : dict
        A corrected disperser ASDF object.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from collimator to ``gwa`` frame.

    """
    with CollimatorModel(reference_files['collimator']) as f:
        collimator = f.model
    angles = [disperser['theta_x'], disperser['theta_y'],
              disperser['theta_z'], disperser['tilt_y']]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotation')
    u2dircos = Unitless2DirCos(name='unitless2directional_cosines')

    return collimator.inverse | u2dircos | rotation


def get_disperser(input_model, disperserfile):
    """
    Return the disperser data model with the GWA
    correction applied.

    Parameters
    ----------
    input_model : `jwst.datamodels.DataModel`
        The input data model - either an ImageModel or a CubeModel.
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
    xtilt : float
        Value of GWAXTILT keyword - angle in arcsec
    ytilt : float
        Value of GWAYTILT keyword - angle in arcsec
    disperser : `~jwst.datamodels.DisperserModel`
        Disperser information.

    Notes
    -----
    The GWA_XTILT keyword is used to correct the THETA_Y angle.
    The GWA_YTILT keyword is used to correct the THETA_X angle.

    Returns
    -------
    disp : `~jwst.datamodels.DisperserModel`
        Corrected DisperserModel.

    """
    def _get_correction(gwa_tilt, tilt_angle):
        phi_exposure = gwa_tilt.tilt_model(tilt_angle)
        phi_calibrator = gwa_tilt.tilt_model(gwa_tilt.zeroreadings[0])
        del_theta = 0.5 * (phi_exposure - phi_calibrator) / 3600. #in deg
        return del_theta

    disp = disperser.copy()
    disperser.close()
    log.info("gwa_ytilt is {0} deg".format(ytilt))
    log.info("gwa_xtilt is {0} deg".format(xtilt))

    if xtilt is not None:
        theta_y_correction = _get_correction(disp.gwa_tiltx, xtilt)
        log.info('theta_y correction: {0} deg'.format(theta_y_correction))
        disp['theta_y'] = disp.theta_y + theta_y_correction
    else:
        log.info('gwa_xtilt not applied')
    if ytilt is not None:
        theta_x_correction = _get_correction(disp.gwa_tilty, ytilt)
        log.info('theta_x correction: {0} deg'.format(theta_x_correction))
        disp.theta_x = disp.theta_x + theta_x_correction
    else:
        log.info('gwa_ytilt not applied')
    return disp


def ifu_msa_to_oteip(reference_files):
    """
    Transform from ``msa_frame`` to ``oteip`` for IFU exposures.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from MSA to OTEIP.
    """
    with FOREModel(reference_files['fore']) as f:
        fore = f.model

    msa2fore_mapping = Mapping((0, 1, 2, 2), name='msa2fore_mapping')
    msa2fore_mapping.inverse = Mapping((0, 1, 2, 2), name='fore2msa')
    fore_transform = msa2fore_mapping | fore & Identity(1)
    return fore_transform


def msa_to_oteip(reference_files):
    """
    Transform from ``msa_frame`` to ``oteip`` for non IFU exposures.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from MSA to OTEIP.

    """
    with FOREModel(reference_files['fore']) as f:
        fore = f.model
    msa2fore_mapping = Mapping((0, 1, 2, 2), name='msa2fore_mapping')
    msa2fore_mapping.inverse = Identity(3)
    return msa2fore_mapping | (fore & Identity(1))


def oteip_to_v23(reference_files):
    """
    Transform from ``oteip`` frame to ``v2v3`` frame.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from ``oteip`` to ``v2v3`` frame.

    """
    with OTEModel(reference_files['ote']) as f:
        ote = f.model
    fore2ote_mapping = Identity(3, name='fore2ote_mapping')
    fore2ote_mapping.inverse = Mapping((0, 1, 2, 2))
    # Create the transform to v2/v3/lambda.  The wavelength units up to this point are
    # meters as required by the pipeline but the desired output wavelength units is microns.
    # So we are going to Scale the spectral units by 1e6 (meters -> microns)
    # The spatial units are currently in deg. Convertin to arcsec.
    oteip2v23 = fore2ote_mapping | (ote & Scale(1e6))

    return oteip2v23


def create_frames():
    """
    Create the coordinate frames in the NIRSPEC WCS pipeline.

    These are
    "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world".
    """
    det = cf.Frame2D(name='detector', axes_order=(0, 1))
    sca = cf.Frame2D(name='sca', axes_order=(0, 1))
    gwa = cf.Frame2D(name="gwa", axes_order=(0, 1), unit=(u.rad, u.rad),
                      axes_names=('alpha_in', 'beta_in'))
    msa_spatial = cf.Frame2D(name='msa_spatial', axes_order=(0, 1), unit=(u.m, u.m),
                             axes_names=('x_msa', 'y_msa'))
    slit_spatial = cf.Frame2D(name='slit_spatial', axes_order=(0, 1), unit=("", ""),
                             axes_names=('x_slit', 'y_slit'))
    sky = cf.CelestialFrame(name='sky', axes_order=(0, 1), reference_frame=coord.ICRS())
    v2v3_spatial = cf.Frame2D(name='v2v3_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec),
                             axes_names=('V2', 'V3'))

    # The oteip_to_v23 incorporates a scale to convert the spectral units from
    # meters to microns.  So the v2v3 output frame will be in u.deg, u.deg, u.micron
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,),
                            axes_names=('wavelength',))
    v2v3 = cf.CompositeFrame([v2v3_spatial, spec], name='v2v3')
    slit_frame = cf.CompositeFrame([slit_spatial, spec], name='slit_frame')
    msa_frame = cf.CompositeFrame([msa_spatial, spec], name='msa_frame')
    oteip_spatial = cf.Frame2D(name='oteip', axes_order=(0, 1), unit=(u.deg, u.deg),
                               axes_names=('X_OTEIP', 'Y_OTEIP'))
    oteip = cf.CompositeFrame([oteip_spatial, spec], name='oteip')
    world = cf.CompositeFrame([sky, spec], name='world')
    return det, sca, gwa, slit_frame, msa_frame, oteip, v2v3, world


def create_imaging_frames():
    """
    Create the coordinate frames in the NIRSPEC WCS pipeline.
    These are
    "detector", "gwa", "msa_frame", "oteip", "v2v3", "world".
    """
    det = cf.Frame2D(name='detector', axes_order=(0, 1))
    sca = cf.Frame2D(name='sca', axes_order=(0, 1))
    gwa = cf.Frame2D(name="gwa", axes_order=(0, 1), unit=(u.rad, u.rad),
                      axes_names=('alpha_in', 'beta_in'))
    msa = cf.Frame2D(name='msa', axes_order=(0, 1), unit=(u.m, u.m),
                             axes_names=('x_msa', 'y_msa'))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.arcsec, u.arcsec),
                              axes_names=('v2', 'v3'))
    oteip = cf.Frame2D(name='oteip', axes_order=(0, 1), unit=(u.deg, u.deg),
                               axes_names=('x_oteip', 'y_oteip'))
    world = cf.CelestialFrame(name='world', axes_order=(0, 1), reference_frame=coord.ICRS())
    return det, sca, gwa, msa, oteip, v2v3, world


def get_slit_location_model(slitdata):
    """
    The transform for the absolute position of a slit on the MSA.

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
        absolute positions in the quadrant..
        This is later combined with the quadrant model to return
        absolute positions in the MSA.
    """
    num, xcenter, ycenter, xsize, ysize = slitdata
    model = models.Scale(xsize) & models.Scale(ysize) | \
            models.Shift(xcenter) & models.Shift(ycenter)
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
    slit : `~jwst.transforms.models.Slit`
        A Fixed slit or MOS slitlet.
    slit_y_range: list or tuple of size 2
        The lower and upper limit of the slit.
        Used for IFU mode only.
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

    tab = Tabular1D(points=(beta_in,),
                    lookup_table=dy, bounds_error=False, name='tabular')
    return tab


def nrs_wcs_set_input(input_model, slit_name, wavelength_range=None):
    """
    Returns a WCS object for a specific slit, slice or shutter.

    Parameters
    ----------
    input_model : `~jwst.datamodels.DataModel`
        A WCS object for the all open slitlets in an observation.
    slit_name : int or str
        Slit.name of an open slit.
    wavelength_range: list
        Wavelength range for the combination of fliter and grating.

    Returns
    -------
    wcsobj : `~gwcs.wcs.WCS`
        WCS object for this slit.
    """
    import copy # TODO: Add a copy method to gwcs.WCS
    wcsobj = input_model.meta.wcs
    if wavelength_range is None:
        _, wrange = spectral_order_wrange_from_model(input_model)
    else:
        wrange = wavelength_range
    slit_wcs = copy.deepcopy(wcsobj)
    slit_wcs.set_transform('sca', 'gwa', wcsobj.pipeline[1].transform[1:])
    # get the open slits from the model
    # Need them to get the slit ymin,ymax
    g2s = wcsobj.pipeline[2].transform
    open_slits = g2s.slits

    slit_wcs.set_transform('gwa', 'slit_frame', g2s.get_model(slit_name))

    exp_type = input_model.meta.exposure.type
    is_nirspec_ifu = is_nrs_ifu_lamp(input_model) or (exp_type.lower() == 'nrs_ifu')
    if is_nirspec_ifu:
        slit_wcs.set_transform('slit_frame', 'slicer',
                           wcsobj.pipeline[3].transform.get_model(slit_name) & Identity(1))
    else:
        slit_wcs.set_transform('slit_frame', 'msa_frame',
                           wcsobj.pipeline[3].transform.get_model(slit_name) & Identity(1))
    slit2detector = slit_wcs.get_transform('slit_frame', 'detector')

    if is_nirspec_ifu:
        bb = compute_bounding_box(slit2detector, wrange)
    else:
        slit = [s for s in open_slits if s.name == slit_name][0]
        bb = compute_bounding_box(slit2detector, wrange,
                                  slit_ymin=slit.ymin, slit_ymax=slit.ymax)

    slit_wcs.bounding_box = bb
    return slit_wcs


def validate_open_slits(input_model, open_slits, reference_files):
    """
    Remove slits which do not project on the detector from the list of open slits.
    For each slit computes the transform from the slit to the detector and
    determines the bounding box.

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        Input data model

    Returns
    -------
    slit2det : dict
        A dictionary with the slit to detector transform for each slit,
        {slit_id: astropy.modeling.Model}
    """

    def _is_valid_slit(domain):
        xlow, xhigh = domain[0]
        ylow, yhigh = domain[1]
        if (xlow >= 2048 or ylow >= 2048 or
            xhigh <= 0 or yhigh <= 0 or
            xhigh - xlow < 1  or yhigh - ylow < 1):
            return False
        else:
            return True

    det2dms = dms_to_sca(input_model).inverse
    # read models from reference file
    disperser = DisperserModel(reference_files['disperser'])
    disperser = correct_tilt(disperser, input_model.meta.instrument.gwa_xtilt,
                             input_model.meta.instrument.gwa_ytilt)

    order, wrange = get_spectral_order_wrange(input_model,
                                              reference_files['wavelengthrange'])

    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = order
    agreq = angle_from_disperser(disperser, input_model)
    # GWA to detector
    det2gwa = detector_to_gwa(reference_files,
                              input_model.meta.instrument.detector,
                              disperser)
    gwa2det = det2gwa.inverse
    # collimator to GWA
    collimator2gwa = collimator_to_gwa(reference_files, disperser)

    msa = MSAModel(reference_files['msa'])
    col2det = collimator2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq | \
            gwa2det | det2dms
    for quadrant in range(1, 6):
        slits_in_quadrant = [s for s in open_slits if s.quadrant == quadrant]
        if any(slits_in_quadrant):
            msa_quadrant = getattr(msa, "Q{0}".format(quadrant))
            msa_model = msa_quadrant.model
            msa_data = msa_quadrant.data
            for slit in slits_in_quadrant:
                slit_id = slit.shutter_id
                slitdata = msa_data[slit_id]
                slitdata_model = get_slit_location_model(slitdata)
                msa_transform = slitdata_model | msa_model
                msa2det = msa_transform & Identity(1) | col2det
                bb = compute_bounding_box(msa2det, wrange, slit.ymin, slit.ymax)
                valid = _is_valid_slit(bb)
                if not valid:
                    log.info("Removing slit {0} from the list of open slits because the "
                             "WCS bounding_box is completely outside the detector.".format(slit.name))
                    log.debug("Slit bounding_box is {0}".format(bb))
                    idx = np.nonzero([s.name == slit.name for s in open_slits])[0][0]
                    open_slits.pop(idx)

    msa.close()
    return open_slits


def spectral_order_wrange_from_model(input_model):
    """
    Return the spectral order and wavelength range used in the WCS.

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        The data model. Must have been through the assign_wcs step.

    """
    wrange = [input_model.meta.wcsinfo.waverange_start, input_model.meta.wcsinfo.waverange_end]
    spectral_order = input_model.meta.wcsinfo.spectral_order
    return spectral_order, wrange


def nrs_ifu_wcs(input_model):
    """
    Return a list of WCSs for all NIRSPEC IFU slits.

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        The data model. Must have been through the assign_wcs step.
    """
    _, wrange = spectral_order_wrange_from_model(input_model)
    wcs_list = []
    # loop over all IFU slits
    for i in range(30):
        wcs_list.append(nrs_wcs_set_input(input_model, i, wrange))
    return wcs_list


def _create_ifupost_transform(ifupost_slice):
    """
    Create an IFUPOST transform for a specific slice.

    Parameters
    ----------
    ifupost_slice : `jwst.datamodels.properties.ObjectNode`
        IFUPost transform for a specific slice

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
    model_x = ((Mapping((0, 1), n_inputs=3) | polyx) +
                       ((Mapping((0, 1), n_inputs=3) | polyx_dist) *
                        (Mapping((2,)) | Identity(1))))
    model_y = ((Mapping((0, 1), n_inputs=3) | polyy) +
                       ((Mapping((0, 1), n_inputs=3) | polyy_dist) *
                        (Mapping((2,)) | Identity(1))))

    output2poly_mapping = Identity(2, name="{0}_outmap".format('ifupost'))
    output2poly_mapping.inverse = Mapping([0, 1, 2, 0, 1, 2])
    input2poly_mapping = Mapping([0, 1, 2, 0, 1, 2], name="{0}_inmap".format('ifupost'))
    input2poly_mapping.inverse = Identity(2)

    model_poly = input2poly_mapping | (model_x & model_y) | output2poly_mapping
    model = linear & Identity(1) | model_poly
    return model

def nrs_lamp(input_model, reference_files, slit_y_range):
    """Return the appropriate function for lamp data

    Parameters
    ----------
    input_model : `~jwst.datamodels.DataModel`
        The input data model.
    reference_files : dict
        The reference files used for this mode.
    slit_y_range : list
        The slit dimensions relative to the center of the slit.
    """
    lamp_mode = input_model.meta.instrument.lamp_mode
    if type(lamp_mode) == str:
        lamp_mode = lamp_mode.lower()
    else:
        lamp_mode = 'none'
    if lamp_mode in ['fixedslit', 'brightobj']:
        return slits_wcs(input_model, reference_files, slit_y_range)
    elif lamp_mode == 'ifu':
        return ifu(input_model, reference_files, slit_y_range)
    elif lamp_mode == 'msaspec':
        return slits_wcs(input_model, reference_files, slit_y_range)
    else:
        return not_implemented_mode(input_model, reference_files, slit_y_range)

exp_type2transform = {
    'nrs_autoflat':  slits_wcs,
    'nrs_autowave':  nrs_lamp,
    'nrs_brightobj': slits_wcs,
    'nrs_confirm':   imaging,
    'nrs_dark':      not_implemented_mode,
    'nrs_fixedslit': slits_wcs,
    'nrs_focus':     imaging,
    'nrs_ifu':       ifu,
    'nrs_image':     imaging,
    'nrs_lamp':      nrs_lamp,
    'nrs_mimf':      imaging,
    'nrs_msaspec':   slits_wcs,
    'nrs_msata':     imaging,
    'nrs_taconfirm': imaging,
    'nrs_tacq':      imaging,
    'nrs_taslit':    imaging,
    'nrs_verify':    imaging,
    'nrs_wata':      imaging,
}
