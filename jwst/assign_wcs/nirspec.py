"""
Tools to create a WCS pipeline list of steps for NIRSPEC modes.

Call create_pipeline() which redirects based on EXP_TYPE

"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import logging
import numpy as np

from asdf import AsdfFile
from astropy.modeling import models, fitting
from astropy.modeling.models import Mapping, Identity, Const1D, Scale, Shift
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import fits
from gwcs import coordinate_frames as cf

from ..transforms.models import (Rotation3DToGWA, DirCos2Unitless, Slit2Msa,
                                 AngleFromGratingEquation, WavelengthFromGratingEquation,
                                 Gwa2Slit, Unitless2DirCos, Logical, Slit)
from .util import not_implemented_mode
from . import pointing

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_pipeline(input_model, reference_files):
    """
    Create a pipeline list based on EXP_TYPE.

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        Either an ImageModel or a CubeModel
    reference_files : dict
        {reftype: file_name} mapping
        In the pipeline it's returned by CRDS.
    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)
    log.info("Creating a NIRSPEC {0} pipeline with references {1}".format(
        exp_type, reference_files))
    return pipeline


def imaging(input_model, reference_files):
    """
    Imaging pipeline

    frames : detector, gwa, msa, sky
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
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotaton').inverse
    dircos2unitless = DirCos2Unitless(name='directional_cosines2unitless')

    col = AsdfFile.open(reference_files['collimator']).tree['model']

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
        msa2oteip.inverse = Mapping((0, 1, 0, 1)) | msa2ote.inverse | Mapping((0, 1), n_inputs=3)
        # OTEIP to V2,V3 transform
        with AsdfFile.open(reference_files['ote']) as f:
            oteip2v23 = f.tree['model'].copy()

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


def ifu(input_model, reference_files):
    """
    IFU pipeline
    """
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
    det2gwa = Identity(2) & detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)

    # GWA to SLIT
    gwa2slit = gwa_to_ifuslit(slits, disperser, wrange, sporder, reference_files)

    # SLIT to MSA transform
    slit2msa = ifuslit_to_msa(slits, reference_files)

    det, sca, gwa, slit_frame, msa_frame, oteip, v2v3, world = create_frames()
    if input_model.meta.instrument.filter != 'OPAQUE':
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
                    (slit_frame, (Mapping((0, 1, 2, 3)) | slit2msa).rename('slit2msa')),
                    (msa_frame, msa2oteip.rename('msa2oteip')),
                    (oteip, oteip2v23.rename('oteip2v23')),
                    (v2v3, tel2sky),
                    (world, None)]
    else:
        # convert to microns if the pipeline ends earlier
        #slit2msa = (Mapping((0, 1, 2, 3, 4)) | slit2msa | Identity(2) & Scale(10**6)).rename('slit2msa')
        slit2msa = (Mapping((0, 1, 2, 3)) | slit2msa).rename('slit2msa')
        pipeline = [(det, dms2detector),
                    (sca, det2gwa.rename('detector2gwa')),
                    (gwa, gwa2slit.rename('gwa2slit')),
                    (slit_frame, slit2msa),
                    (msa_frame, None)]

    return pipeline


def slits_wcs(input_model, reference_files):
    """
    Create the WCS pipeline for observations using the MSA shutter array or fixed slits.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel`
        The input data model.
    reference_files : dict
        Dictionary with reference files supplied by CRDS.

    """
    open_slits_id = get_open_slits(input_model)
    n_slits = len(open_slits_id)
    log.info("Computing WCS for {0} open slitlets".format(n_slits))
    # Get the corrected disperser model
    disperser = get_disperser(input_model, reference_files['disperser'])

    # Get the default spectral order and wavelength range and record them in the model.
    sporder, wrange = get_spectral_order_wrange(input_model, reference_files['wavelengthrange'])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = sporder

    # DMS to SCA transform
    dms2detector = dms_to_sca(input_model).rename('dms2sca')
    # DETECTOR to GWA transform
    det2gwa = Identity(2) & detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)

    # GWA to SLIT
    gwa2slit = gwa_to_slit(open_slits_id, disperser, wrange, sporder, reference_files)

    # SLIT to MSA transform
    slit2msa = slit_to_msa(open_slits_id, reference_files['msa'])

    # Create coordinate frames in the NIRSPEC WCS pipeline"
    # "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world"
    det, sca, gwa, slit_frame, msa_frame, oteip, v2v3, world = create_frames()
    if input_model.meta.instrument.filter != 'OPAQUE':
        # MSA to OTEIP transform
        msa2oteip = msa_to_oteip(reference_files)

        # OTEIP to V2,V3 transform
        # This includes a wavelength unit conversion from meters to microns.
        oteip2v23 = oteip_to_v23(reference_files)

        # V2, V3 to sky
        tel2sky = pointing.v23tosky(input_model) & Identity(1)

        msa_pipeline = [(det, dms2detector),
                        (sca, det2gwa.rename('det2gwa')),
                        (gwa, gwa2slit.rename('gwa2slit')),
                        (slit_frame, (Mapping((0, 1, 2, 3)) | slit2msa).rename('slit2msa')),
                        (msa_frame, msa2oteip.rename('msa2oteip')),
                        (oteip, oteip2v23.rename('oteip2v23')),
                        (v2v3, tel2sky),
                        (world, None)]
    else:
        # convert to microns if the pipeline ends earlier
        #gwa2slit = (gwa2slit | Identity(2) & Scale(10**6)).rename('gwa2slit')
        gwa2slit = (gwa2slit).rename('gwa2slit')
        msa_pipeline = [(det, dms2detector),
                        (sca, det2gwa),
                        (gwa, gwa2slit),
                        (slit_frame, Mapping((0, 1, 2, 3)) | slit2msa),
                        (msa_frame, None)]

    return msa_pipeline


def get_open_slits(input_model):
    exp_type = input_model.meta.exposure.type.lower()
    if exp_type == "nrs_msaspec":
        msa_config = input_model.meta.instrument.msa_configuration_file
        slits = get_open_msa_slits(msa_config, input_model.meta.instrument.msa_metadata_id)
    elif exp_type == "nrs_fixedslit":
        slits = get_open_fixed_slits(input_model)
    else:
        raise ValueError("EXP_TYPE {0} is not supported".format(exp_type.upper()))
    return slits


def get_open_fixed_slits(input_model):
    slits = []
    slits.append(Slit('S200A1', 0, 0, 0, -.5, .5, 5))
    slits.append(Slit('S200A2', 1, 0, 0, -.5, .5, 5))
    slits.append(Slit('S400A1', 2, 0, 0, -.5, .5, 5))
    slits.append(Slit('S1600A1', 3, 0, 0, -.5, .5, 5))

    if input_model.meta.instrument.detector == 'NRS1':
        if input_model.meta.instrument.filter == 'F070LP' and \
                input_model.meta.instrument.grating == 'G140H':
            slits.append(Slit('S200B1', 4, 0, 0, -.5, .5, 5))
    return slits


def get_open_msa_slits(msa_file, msa_metadata_id):
    """
    Computes (ymin, ymax) of open slitlets.

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
        ('estimated_source_in_shutter_y', '>f4')])

    For example, something like:
        (12, 2, 4, 251, 22, 1, 'Y', 'OPEN', nan, nan),

       column

    Parameters
    ----------
        msa_file : str
            MSA configuration file name, FITS keyword MSACONFL.
        msa_metadata_id : int
            The MSA meta id for the science file, FITS keyword MSAMETID.

    Returns
    -------
    slitlets : list
        A list of slitlets. Each slitlet is a tuple with
        ("name", "shutter_id", "xcen", "ycen", "ymin", "ymax", "quadrant", "source_id", "nshutters")

    """

    slitlets = []

    # If they passed in a string then we shall assume it is the filename
    # of the configuration file.
    with fits.open(msa_file) as msa_file:
        # Get the configuration header from teh _msa.fits file.  The EXTNAME should be 'SHUTTER_INFO'
        msa_conf = msa_file[('SHUTTER_INFO', 1)]
        msa_source = msa_file[("SOURCE_INFO", 1)].data

        # First we are going to filter the msa_file data on the msa_metadata_id
        # as that is all we are interested in for this function.
        msa_data = [x for x in msa_conf.data if x['msa_metadata_id'] == msa_metadata_id]

        log.debug('msa_data with msa_metadata_id = {}   {}'.format(msa_metadata_id, msa_data))
        log.info('Retrieving open slitlets for msa_metadata_id = {}'.format(msa_metadata_id))

        # First thing to do is to get the unique slitlet_ids
        slitlet_ids_unique = list(set([x['slitlet_id'] for x in msa_data]))

        # Now lets look at each unique slitlet id
        for slitlet_id in slitlet_ids_unique:

            # Get the rows for the current slitlet_id
            slitlets_sid = [x for x in msa_data if x['slitlet_id'] == slitlet_id]
            nshutters = len(slitlets_sid)
            # Count the number of backgrounds that have an 'N' (meaning main shutter)
            # This needs to be 0 or 1 and we will have to deal with those differently
            # See: https://github.com/STScI-JWST/jwst/commit/7588668b44b77486cdafb35f7e2eb2dcfa7d1b63#commitcomment-18987564

            n_main_shutter = len([s for s in slitlets_sid if s['background'] == 'N'])

            # In the next part we need to calculate, find, determine 5 things:
            #    quadrant,  xcen, ycen,  ymin, max

            margin = 0.05

            # There are no main shutters, all are background
            if n_main_shutter == 0:
                jmin = min([s['shutter_column'] for s in slitlets_sid])
                jmax = max([s['shutter_column'] for s in slitlets_sid])
                j = (jmax - jmin) // 2 + 1
                ymax = 0.5 + margin + (jmax - j) * 1.15
                ## TODO: check this formula - it is different (assuming it's incorrect in the report).
                ymin = -(0.5 + margin) + (jmin - j) * 1.15
                quadrant = slitlets_sid[0]['shutter_quadrant']
                ycen = j
                xcen = slitlets_sid[0]['shutter_row']  # grab the first as they are all the same
                source_xpos = 0.0
                source_ypos = 0.0
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
                ymax = 0.5 + margin + (jmax - j) * 1.15
                ymin = -(0.5 + margin) + (jmin - j) * 1.15

            # Not allowed....
            else:
                raise ValueError("MSA configuration file has more than 1 shutter with "
                                 "sources for metadata_id = {}".format(msa_metadata_id))

            shutter_id = xcen + (ycen - 1) * 365
            source_id = slitlets_sid[0]['source_id']
            source_name, source_alias, catalog_id, stellarity = [
                (s['source_name'], s['alias'], s['catalog_id'], s['stellarity']) \
                for s in msa_source if s['source_id'] == source_id][0]
            # Create the output list of tuples that contain the required
            # data for further computations
            slitlets.append(Slit(slitlet_id, shutter_id, xcen, ycen, ymin, ymax,
                                 quadrant, source_id, nshutters, source_name, source_alias,
                                 catalog_id, stellarity, source_xpos, source_ypos))

    return slitlets


def get_spectral_order_wrange(input_model, wavelengthrange_file):
    """
    Read the spectral order and wavelength range from the reference file.

    Parameters
    ----------
    filter : str
        The filter used.
    grating : str
        The grating used in the observation.
    wavelength_range_file : str
        Reference file of type "wavelengthrange".
    """
    full_range = [.6e-6, 5.3e-6]

    filter = input_model.meta.instrument.filter
    lamp = input_model.meta.instrument.lamp_state
    grating = input_model.meta.instrument.grating

    wave_range = AsdfFile.open(wavelengthrange_file)
    if filter == "OPAQUE":
        keyword = lamp + '_' + grating
    else:
        keyword = filter + '_' + grating
    try:
        order = wave_range.tree['filter_grating'][keyword]['order']
        wrange = wave_range.tree['filter_grating'][keyword]['range']
    except KeyError:
        order = -1
        wrange = full_range
        log.warning("Combination {0} missing in wavelengthrange file, setting order to -1 and range to {1}.".format(keyword, full_range))
    wave_range.close()
    return order, wrange


def ifuslit_to_msa(slits, reference_files):
    """
    The transform from slit_frame to msa_frame.

    Parameters
    ----------
    slits_id : list
        A list of slit IDs for all open shutters/slitlets.
    msafile : str
        The name of the msa reference file.

    Returns
    -------
    model : `~jwst.transforms.Slit2Msa` model.
        Transform from slit_frame to msa_frame.
    """

    ifuslicer = AsdfFile.open(reference_files['ifuslicer'])
    models = []
    ifuslicer_model = (ifuslicer.tree['model']).rename('ifuslicer_model')
    for slit in slits:
        slitdata = ifuslicer.tree['data'][slit]
        slitdata_model = (get_slit_location_model(slitdata)).rename('slitdata_model')
        msa_transform = slitdata_model | ifuslicer_model
        models.append(msa_transform)
    ifuslicer.close()

    return Slit2Msa(slits, models)


def slit_to_msa(open_slits, msafile):
    """
    The transform from slit_frame to msa_frame.

    Parameters
    ----------
    open_slits : list
        A list of slit IDs for all open shutters/slitlets.
    msafile : str
        The name of the msa reference file.

    Returns
    -------
    model : `~jwst.transforms.Slit2Msa` model.
        Transform from slit_frame to msa_frame.
    """
    msa = AsdfFile.open(msafile)
    models = []
    for quadrant in range(1, 6):
        slits_in_quadrant = [s for s in open_slits if s.quadrant==quadrant]
        if any(slits_in_quadrant):
            msa_data = msa.tree[quadrant]['data']
            msa_model = msa.tree[quadrant]['model']
            for slit in slits_in_quadrant:
                slit_id = slit.shutter_id
                slitdata = msa_data[slit_id]
                slitdata_model = get_slit_location_model(slitdata)
                msa_transform = slitdata_model | msa_model
                models.append(msa_transform)
    msa.close()
    return Slit2Msa(open_slits, models)


def gwa_to_ifuslit(slits, disperser, wrange, order, reference_files):
    """
    GWA to SLIT transform.

    Parameters
    ----------
    slits : list
        A list of slit IDs for all IFU slits 0-29.
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
        Transform from GWA frame to SLIT frame.
   """
    ymin = -.55
    ymax = .55
    agreq = AngleFromGratingEquation(disperser['groove_density'], order, name='alpha_from_greq')
    lgreq = WavelengthFromGratingEquation(disperser['groove_density'], order, name='lambda_from_greq')
    collimator2gwa = collimator_to_gwa(reference_files, disperser)
    mask = mask_slit(ymin, ymax)

    ifuslicer = AsdfFile.open(reference_files['ifuslicer'])
    ifupost = AsdfFile.open(reference_files['ifupost'])
    slit_models = []
    ifuslicer_model = ifuslicer.tree['model']
    for slit in slits:
        slitdata = ifuslicer.tree['data'][slit]
        slitdata_model = get_slit_location_model(slitdata)
        ifuslicer_transform = (slitdata_model | ifuslicer_model)
        ifupost_transform = ifupost.tree[slit]['model']
        msa2gwa = ifuslicer_transform | ifupost_transform | collimator2gwa
        gwa2msa = gwa_to_ymsa(msa2gwa)# TODO: Use model sets here
        bgwa2msa = Mapping((0, 1, 0, 1), n_inputs=3) | \
                 Const1D(0) * Identity(1) & Const1D(-1) * Identity(1) & Identity(2) | \
                 Identity(1) & gwa2msa & Identity(2) | \
                 Mapping((0, 1, 0, 1, 2, 3)) | Identity(2) & msa2gwa & Identity(2) | \
                 Mapping((0, 1, 2, 5), n_inputs=7) | Identity(2) & lgreq | mask

        # msa to before_gwa
        msa2bgwa = msa2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
        bgwa2msa.inverse = msa2bgwa
        slit_models.append(bgwa2msa)

    ifuslicer.close()
    ifupost.close()
    return Gwa2Slit(slits, slit_models)


def gwa_to_slit(open_slits, disperser, wrange, order, reference_files):
    """
    GWA to SLIT transform.

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
        Transform from GWA frame to SLIT frame.
    """
    agreq = AngleFromGratingEquation(disperser['groove_density'], order, name='alpha_from_greq')
    lgreq = WavelengthFromGratingEquation(disperser['groove_density'], order, name='lambda_from_greq')
    collimator2gwa = collimator_to_gwa(reference_files, disperser)

    msa = AsdfFile.open(reference_files['msa'])
    slit_models = []
    for quadrant in range(1, 6):
        slits_in_quadrant = [s for s in open_slits if s.quadrant==quadrant]
        log.info("There are {0} open slits in quadrant {1}".format(len(slits_in_quadrant), quadrant))
        if any(slits_in_quadrant):
            msa_model = msa.tree[quadrant]['model']
            log.info("Getting slits location for quadrant {0}".format(quadrant))
            msa_data = msa.tree[quadrant]['data']
            for slit in slits_in_quadrant:
                mask = mask_slit(slit.ymin, slit.ymax)
                slit_id = slit.shutter_id
                slitdata = msa_data[slit_id]
                slitdata_model = get_slit_location_model(slitdata)
                msa_transform = slitdata_model | msa_model
                msa2gwa = (msa_transform | collimator2gwa)
                gwa2msa = gwa_to_ymsa(msa2gwa)# TODO: Use model sets here
                bgwa2msa = Mapping((0, 1, 0, 1), n_inputs=3) | \
                    Const1D(0) * Identity(1) & Const1D(-1) * Identity(1) & Identity(2) | \
                    Identity(1) & gwa2msa & Identity(2) | \
                    Mapping((0, 1, 0, 1, 2, 3)) | Identity(2) & msa2gwa & Identity(2) | \
                    Mapping((0, 1, 2, 5), n_inputs=7) | Identity(2) & lgreq | mask

                # msa to before_gwa
                msa2bgwa = msa2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
                bgwa2msa.inverse = msa2bgwa
                slit_models.append(bgwa2msa)
    msa.close()
    return Gwa2Slit(open_slits, slit_models)


def detector_to_gwa(reference_files, detector, disperser):
    """
    Transform from DETECTOR frame to GWA frame.

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
    with AsdfFile.open(reference_files['fpa']) as f:
        fpa = f.tree[detector].copy()
    with AsdfFile.open(reference_files['camera']) as f:
        camera = f.tree['model'].copy()

    angles = [disperser['theta_x'], disperser['theta_y'],
               disperser['theta_z'], disperser['tilt_y']]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotaton')
    u2dircos = Unitless2DirCos(name='unitless2directional_cosines')
    model = (models.Shift(-1) & models.Shift(-1) | fpa | camera | u2dircos | rotation)
    return model


def dms_to_sca(input_model):
    """
    Transforms from DMS to SCA coordinates.
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
    # The outputs must be 1-based becaause this is what the model expects.
    # If xstart was 0-based and the inputs were 0-based ->
    # Shift(+1)
    subarray2full = models.Shift(xstart - 1) & models.Shift(ystart - 1)
    if detector == 'NRS2':
        model = models.Shift(-2048) & models.Shift(-2048) | models.Scale(-1) & models.Scale(-1)
    elif detector == 'NRS1':
        model = models.Identity(2)
    return subarray2full | model


def mask_slit(ymin=-.5, ymax=.5):
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


def compute_domain(slit2detector, wavelength_range, slit_ymin=-.5, slit_ymax=.5):
    """
    Compute the projection of a slit/slice on the detector.

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
    lam_min = wavelength_range[0]
    lam_max = wavelength_range[1]

    step = 1e-10
    nsteps = int((lam_max - lam_min) / step)
    lam_grid = np.linspace(lam_min, lam_max, nsteps)
    x_range_low, y_range_low = slit2detector([0] * nsteps, [slit_ymin] * nsteps, lam_grid)
    x_range_high, y_range_high = slit2detector([0] * nsteps, [slit_ymax] * nsteps, lam_grid)
    x_range = np.hstack((x_range_low, x_range_high))
    y_range = np.hstack((y_range_low, y_range_high))
    # add 10 px margin
    # The -1 is technically because the output of slit2detector is 1-based coordinates.
    x0 = int(max(0, x_range.min() -1 -10))
    x1 = int(min(2047, x_range.max() -1 + 10))
    # add 2 px margin
    y0 = int(y_range.min() -1 -2)
    y1 = int(y_range.max() -1 + 2)

    domain = [{'lower': x0, 'upper': x1}, {'lower': y0, 'upper': y1}]
    return domain


def collimator_to_gwa(reference_files, disperser):
    """
    Transform from COLLIMATOR to GWA frame.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.
    disperser : dict
        A corrected disperser ASDF object.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from COLLIMATOR to GWA frame.

    """
    with AsdfFile.open(reference_files['collimator']) as f:
        collimator = f.tree['model'].copy()
    angles = [disperser['theta_x'], disperser['theta_y'],
              disperser['theta_z'], disperser['tilt_y']]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotaton')
    u2dircos = Unitless2DirCos(name='unitless2directional_cosines')

    return collimator.inverse | u2dircos | rotation


def get_disperser(input_model, disperserfile):
    """
    Return the disperser information corrected for the uncertainty in the GWA position.

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
    with AsdfFile.open(disperserfile) as f:
        disperser = f.tree
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
    disperser : dict
        Disperser information.

    """
    def _get_correction(gwa_tilt, tilt_angle):
        phi_exposure = gwa_tilt['tilt_model'](tilt_angle)
        phi_calibrator = gwa_tilt['tilt_model'](gwa_tilt['zeroreadings'][0])
        del_theta = 0.5 * (phi_exposure - phi_calibrator) / 3600. #in deg
        return del_theta

    disp = disperser.copy()
    log.info("gwa_ytilt is {0}".format(ytilt))
    log.info("gwa_xtilt is {0}".format(xtilt))

    if xtilt is not None:
        theta_y_correction = _get_correction(disperser['gwa_tiltx'], xtilt)
        log.info('theta_y correction: {0}'.format(theta_y_correction))
        disp['theta_y'] = disperser['theta_y'] + theta_y_correction
    else:
        log.info('gwa_xtilt not applied')
    if ytilt is not None:
        theta_x_correction = _get_correction(disperser['gwa_tilty'], ytilt)
        log.info('theta_x correction: {0}'.format(theta_x_correction))
        disp['theta_x'] = disperser['theta_x'] + theta_x_correction
    else:
        log.info('gwa_ytilt not applied')
    return disp


def ifu_msa_to_oteip(reference_files):
    """
    Transform from the MSA frame to the OTEIP frame.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from MSA to OTEIP.
    """
    with AsdfFile.open(reference_files['fore']) as f:
        fore = f.tree['model'].copy()
    with AsdfFile.open(reference_files['ifufore']) as f:
        ifufore = f.tree['model'].copy()

    msa2fore_mapping = Mapping((0, 1, 2, 2))
    msa2fore_mapping.inverse = Identity(3)
    ifu_fore_transform = ifufore & Identity(1)
    ifu_fore_transform.inverse = Mapping((0, 1, 2, 2)) | ifufore.inverse & Identity(1)
    fore_transform = msa2fore_mapping | fore & Identity(1)
    return msa2fore_mapping | ifu_fore_transform | fore_transform


def msa_to_oteip(reference_files):
    """
    Transform from the MSA frame to the OTEIP frame.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from MSA to OTEIP.

    """
    with AsdfFile.open(reference_files['fore']) as f:
        fore = f.tree['model'].copy()
    msa2fore_mapping = Mapping((0, 1, 2, 2), name='msa2fore_mapping')
    msa2fore_mapping.inverse = Identity(3)
    return msa2fore_mapping | (fore & Identity(1))


def oteip_to_v23(reference_files):
    """
    Transform from the OTEIP frame to the V2V3 frame.

    Parameters
    ----------
    reference_files: dict
        Dictionary with reference files returned by CRDS.

    Returns
    -------
    model : `~astropy.modeling.core.Model` model.
        Transform from OTEIP to V2V3.

    """
    with AsdfFile.open(reference_files['ote']) as f:
        ote = f.tree['model'].copy()
    fore2ote_mapping = Identity(3, name='fore2ote_mapping')
    fore2ote_mapping.inverse = Mapping((0, 1, 2, 2))
    # Create the transform to v2/v3/lambda.  The wavelength units up to this point are
    # meters as required by the pipeline but the desired output wavelength units is microns.
    # So we are going to Scale the spectral units by 1e6 (meters -> microns)
    # The spatial units are currently in deg. Convertin to arcsec.
    oteip_to_xyan = fore2ote_mapping | (ote & Scale(1e6))
    # Add a shift for the aperture.
    oteip2v23 = oteip_to_xyan | Identity(1) & (Shift(468 / 3600) | Scale(-1)) & Identity(1)
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
    v2v3_spatial = cf.Frame2D(name='v2v3_spatial', axes_order=(0, 1), unit=(u.deg, u.deg),
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
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.deg, u.deg),
                              axes_names=('v2', 'v3'))
    oteip = cf.Frame2D(name='oteip', axes_order=(0, 1), unit=(u.deg, u.deg),
                               axes_names=('x_oteip', 'y_oteip'))
    world = cf.CelestialFrame(name='world', axes_order=(0, 1), reference_frame=coord.ICRS())
    return det, sca, gwa, msa, oteip, v2v3, world


def get_slit_location_model(slitdata):
    """
    Compute the transform for the absolute position of a slit.

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


def gwa_to_ymsa(msa2gwa_model):
    """
    Determine the linear relation d_y(beta_in) for the aperture on the detector.

    Parameters
    ----------
    msa2gwa_model : `astropy.modeling.core.Model`
        The transform from the MSA to the GWA.
    """
    dy = np.linspace(-.55, .55, 1000)
    dx = np.zeros(dy.shape)
    cosin_grating_k = msa2gwa_model(dx, dy)
    fitter = fitting.LinearLSQFitter()
    model = models.Polynomial1D(1)
    poly1d_model = fitter(model, cosin_grating_k[1], dy)
    poly_model = poly1d_model.rename('interpolation')
    return poly_model


def nrs_wcs_set_input(input_model, slit_name, wavelength_range=None):
    """
    Returns a WCS object for this slit.

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
    slit_wcs.set_transform('sca', 'gwa', wcsobj.pipeline[1][1][1:])
    # get the open slits from the model
    # Need them to get the slit ymin,ymax
    g2s = wcsobj.pipeline[2][1]
    open_slits = g2s.slits

    slit_wcs.set_transform('gwa', 'slit_frame', g2s.get_model(slit_name))
    slit_wcs.set_transform('slit_frame', 'msa_frame', wcsobj.pipeline[3][1][1].get_model(slit_name) & Identity(1))
    slit2detector = slit_wcs.get_transform('slit_frame', 'detector')

    if input_model.meta.exposure.type.lower() != 'nrs_ifu':
        slit = [s for s in open_slits if s.name == slit_name][0]
        domain = compute_domain(slit2detector, wrange, slit_ymin=slit.ymin, slit_ymax=slit.ymax)
    else:
        domain = compute_domain(slit2detector, wrange)
    slit_wcs.domain = domain
    return slit_wcs


#def slit_to_detector(input_model, slits_id, lam, reference_files):
    #"""
    #For each slit computes the transform from the MSA to the detector.

    #Used to flag stuck open shutters.
    #Reftypes needed:
    #reftypes = ['fpa', 'camera', 'disperser', 'collimator', 'msa', 'wavelengthrange']

    #Parameters
    #----------
    #input_model : jwst.datamodels.DataModel
        #Input data model
    #slits : list
        #A list of slit IDs. A slit ID is a tuple of (quadrant, slit_nimber)
    #lam : float
        #wavelength, in meters
    #reference_files : dict
        #A dictionary with WCS reference_files, {reftype: refname}

    #Returns
    #-------
    #msa2det : dict
        #A dictionary with the MSA to detector transform for each slit, {slit_id: astropy.modeling.Model}
    #"""
    #slits_msa2det = {}
    ## read models from reference files
    #with AsdfFile.open(reference_files['fpa']) as f:
        #fpa = f.tree[input_model.meta.instrument.detector].copy()
    #with AsdfFile.open(reference_files['camera']) as f:
        #camera = f.tree['model'].copy()
    #with AsdfFile.open(reference_files['disperser']) as f:
        #disperser = f.tree

    #dircos2u = DirCos2Unitless(name='directional2unitless_cosines')
    #disperser = correct_tilt(disperser, input_model.meta.instrument.gwa_xtilt,
                             #input_model.meta.instrument.gwa_ytilt)
    #angles = [disperser['theta_x'], disperser['theta_y'],
               #disperser['theta_z'], disperser['tilt_y']]
    #rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotaton')

    #order, wrange = get_spectral_order_wrange(input_model, reference_files['wavelengthrange'])
    #input_model.meta.wcsinfo.waverange_start = wrange[0]
    #input_model.meta.wcsinfo.waverange_end = wrange[1]
    #input_model.meta.wcsinfo.spectral_order = order

    #agreq = AngleFromGratingEquation(disperser['groove_density'], order, name='alpha_from_greq')
    ## GWA to detector
    #gwa2det = rotation.inverse | dircos2u | camera.inverse | fpa.inverse
    ## collimator to GWA
    #collimator2gwa = collimator_to_gwa(reference_files, disperser)


    #msa = AsdfFile.open(reference_files['msa'])
    #for i in range(1, 5):
        #slit_names = slits_id[slits_id[:, 0] == i]
        #if slit_names.any():
            #msa_model = msa.tree[i]['model']
            #for slit in slit_names:
                #index = slit[1] - 1
                #slitdata = msa.tree[slit[0]]['data'][index]
                #slitdata_model = get_slit_location_model(slitdata)
                ## absolute positions of the slit on the MSA
                #msa_transform = slitdata_model | msa_model
                #s = slitid_to_slit(np.array([slit]))[0]
                #msa2gwa = (msa_transform | collimator2gwa) & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
                #slits_msa2det[s] = msa2gwa | gwa2det
    #msa.close()

    #return slits_msa2det


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
    Return a list of WCS for all NIRSPEC IFU slits.

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


exp_type2transform = {'nrs_tacq': imaging,
                      'nrs_taslit': imaging,
                      'nrs_taconfirm': imaging,
                      'nrs_confirm': imaging,
                      'nrs_fixedslit': slits_wcs,
                      'nrs_ifu': ifu,
                      'nrs_msaspec': slits_wcs,
                      'nrs_image': imaging,
                      'nrs_focus': imaging,
                      'nrs_mimf': imaging,
                      'nrs_bota': imaging,
                      'nrs_autoflat': not_implemented_mode, #TBD
                      'nrs_autowave': not_implemented_mode, #TBD
                      'nrs_lamp': slits_wcs,
                      'nrs_brightobj': slits_wcs
                      }
