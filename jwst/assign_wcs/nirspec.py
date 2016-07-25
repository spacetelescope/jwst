"""
Tools to create a WCS pipeline list of steps for NIRSPEC modes.

Call create_pipeline() which redirects based on EXP_TYPE

"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import logging
import copy
import types
import os.path
import numpy as np

from asdf import AsdfFile
from astropy.modeling import models, fitting
from astropy.modeling.core import Model
from astropy.modeling.models import Mapping, Identity, Const1D
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import fits
from gwcs import coordinate_frames as cf

from ..transforms.models import *
from .util import not_implemented_mode
from . import pointing

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


slit_id2name = {(5, 1): 'S200A1', (5, 2): 'S200A2', (5, 3): 'S400A1',
                (5, 4): 'S1600A1', (5, 5): 'S200B1', (5, 6): 'IFU'}


slit_name2id = {'S200A1': (5, 1), 'S200A2': (5, 2), 'S400A1': (5, 3),
                'S1600A1': (5, 4), 'S200B1': (5, 5), 'IFU': (5, 6)}


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
    det, gwa, msa_frame, oteip, v2v3 = create_imaging_frames()
    if input_model.meta.instrument.filter != 'OPAQUE':
        # MSA to OTEIP transform
        msa2ote = msa_to_oteip(reference_files)
        msa2oteip = msa2ote | Mapping((0, 1), n_inputs=3)
        msa2oteip.inverse = Mapping((0, 1, 0, 1)) | msa2ote.inverse | Mapping((0, 1), n_inputs=3)
        # OTEIP to V2,V3 transform
        with AsdfFile.open(reference_files['ote']) as f:
            oteip2v23 = f.tree['model'].copy()

        # V2, V3 to wprld (RA, DEC) transform
        v23_to_sky = pointing.fitswcs_transform_from_model(input_model)

        imaging_pipeline = [(det, det2gwa),
                            (gwa, gwa2msa),
                            (msa_frame, msa2oteip),
                            (oteip, oteip2v23),
                            (v2v3, v23_to_sky),
                            (world, None)]
    else:
        imaging_pipeline = [(det, det2gwa),
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

    # DETECTOR to GWA transform
    det2gwa = Identity(2) & detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)

    # GWA to SLIT
    gwa2slit = gwa_to_ifuslit(slits, disperser, wrange, sporder, reference_files)

    # SLIT to MSA transform
    slit2msa = ifuslit_to_msa(slits, reference_files)

    det, gwa, slit_frame, msa_frame, oteip, v2v3, world = create_frames()
    if input_model.meta.instrument.filter != 'OPAQUE':
        # MSA to OTEIP transform
        msa2oteip = ifu_msa_to_oteip(reference_files)

        # OTEIP to V2,V3 transform
        oteip2v23 = oteip_to_v23(reference_files)

        # V2, V3 to wprld (RA, DEC ,LAMBDA) transform
        v23_to_sky = pointing.fitswcs_transform_from_model(input_model)
        v23_to_world = v23_to_sky & Identity(1)

        # Create coordinate frames in the NIRSPEC WCS pipeline"
        # "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world"

        pipeline = [(det, det2gwa.rename('detector2gwa')),
                    (gwa, gwa2slit.rename('gwa2slit')),
                    (slit_frame, (Mapping((0, 1, 2, 3, 4)) | slit2msa).rename('slit2msa')),
                    (msa_frame, msa2oteip.rename('msa2oteip')),
                    (oteip, oteip2v23.rename('oteip2v23')),
                    (v2v3, v23_to_world),
                    (world, None)]
    else:

        pipeline = [(det, det2gwa.rename('detector2gwa')),
                    (gwa, gwa2slit.rename('gwa2slit')),
                    (slit_frame, (Mapping((0, 1, 2, 3, 4)) | slit2msa).rename('slit2msa')),
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

    # DETECTOR to GWA transform
    det2gwa = Identity(2) & detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)
    #det2gwa = detector_to_gwa(reference_files, input_model.meta.instrument.detector, disperser)

    # GWA to SLIT
    gwa2slit = gwa_to_slit(open_slits_id, disperser, wrange, sporder, reference_files)

    # SLIT to MSA transform
    slit2msa = slit_to_msa(open_slits_id, reference_files['msa'])

    # Create coordinate frames in the NIRSPEC WCS pipeline"
    # "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world"
    det, gwa, slit_frame, msa_frame, oteip, v2v3, world = create_frames()
    if input_model.meta.instrument.filter != 'OPAQUE':
        # MSA to OTEIP transform
        msa2oteip = msa_to_oteip(reference_files)

        # OTEIP to V2,V3 transform
        oteip2v23 = oteip_to_v23(reference_files)

        # V2, V3 to wprld (RA, DEC ,LAMBDA) transform
        v23_to_sky = pointing.fitswcs_transform_from_model(input_model)
        v23_to_world = v23_to_sky & Identity(1)

        msa_pipeline = [(det, det2gwa),
                        (gwa, gwa2slit),
                        (slit_frame, Mapping((0, 1, 2, 3, 4)) | slit2msa),
                        (msa_frame, msa2oteip),
                        (oteip, oteip2v23),
                        (v2v3, v23_to_world),
                        (world, None)]
    else:
        msa_pipeline = [(det, det2gwa),
                        (gwa, gwa2slit),
                        (slit_frame, Mapping((0, 1, 2, 3, 4)) | slit2msa),
                        (msa_frame, None)]

    return msa_pipeline


def get_open_slits(input_model):
    exp_type = input_model.meta.exposure.type.lower()
    if exp_type == "nrs_msaspec":
        #msa_config = "SPCB-GD-A.msa.fits"
        msa_config = input_model.meta.instrument.msa_configuration_file
        slits = get_open_msa_slits(msa_config)
    elif exp_type == "nrs_fixedslit":
        slits = get_open_fixed_slits(input_model)
    else:
        raise ValueError("EXP_TYPE {0} is not supported".format(exp_type.upper()))
    return slits


def get_open_fixed_slits(input_model):
    if input_model.meta.instrument.detector == 'NRS1':
        if input_model.meta.instrument.filter == 'F070LP' and \
                input_model.meta.instrument.grating == 'G140H':
            slits = np.array([(5, 1), (5, 2), (5, 3), (5, 4), (5, 5)])
        else:
            slits = np.array([(5, 1), (5, 2), (5, 3), (5, 4)])
    else:
        slits = np.array([(5, 1), (5, 2), (5, 3), (5, 4), (5, 5)])
    return slits


def get_open_msa_slits(msa_status):
    """
    Returns slit_id for all open slits.

    Parameters
    ----------
    msa_status : str
        Name of MSA status file.

    Returns
    -------
    msa_slits_id : ndarray
        An array of slit ids.
        A slit_id is a tuple of (quadrant_number, slit_number).

    """
    with fits.open(msa_status) as fmsa:
        msa_slits_id = []
        for i in range(1, 5):
            quadrant = fmsa[i].data
            open_shutters = quadrant[quadrant['status'].nonzero()]
            msa_slits_id.extend([(i, n) for n in open_shutters['NO']])
    return np.array(msa_slits_id)


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
    with AsdfFile.open(reference_files['ifufore']) as f:
        ifufore = f.tree['model']

    ifuslicer = AsdfFile.open(reference_files['ifuslicer'])
    models = {}
    ifuslicer_model = (ifuslicer.tree['model']).rename('ifuslicer_model')
    for slit in slits:
        slitdata = ifuslicer.tree['data'][slit]
        slitdata_model = (get_slit_location_model(slitdata)).rename('slitdata_model')
        msa_transform = slitdata_model | ifuslicer_model
        models[slit] = msa_transform
    ifuslicer.close()

    return Slit2Msa(models)


def slit_to_msa(slits_id, msafile):
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
    msa = AsdfFile.open(msafile)
    models = {}
    for i in range(1, 6):
        slit_names = slits_id[slits_id[:, 0] == i]
        if slit_names.any():
            msa_model = msa.tree[i]['model']
            for slit in slit_names:
                index = slit[1] - 1
                slitdata = msa.tree[slit[0]]['data'][index]
                slitdata_model = get_slit_location_model(slitdata)
                msa_transform = slitdata_model | msa_model
                s = slitid_to_slit(np.array([slit]))[0]
                models[s] = msa_transform
    msa.close()
    return Slit2Msa(models)


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
    agreq = AngleFromGratingEquation(disperser['groove_density'], order, name='alpha_from_greq')
    lgreq = WavelengthFromGratingEquation(disperser['groove_density'], order, name='lambda_from_greq')
    collimator2gwa = collimator_to_gwa(reference_files, disperser)
    lam = (wrange[1] - wrange[0]) / 2 + wrange[0]

    ifuslicer = AsdfFile.open(reference_files['ifuslicer'])
    ifupost = AsdfFile.open(reference_files['ifupost'])
    slit_models = {}
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
                 Mapping((0, 1, 2, 5), n_inputs=7) | Identity(2) & lgreq

        # msa to before_gwa
        #msa2bgwa = Mapping((0, 1, 2, 2)) | msa2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
        msa2bgwa = msa2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
        bgwa2msa.inverse = msa2bgwa
        slit_models[slit] = bgwa2msa

    ifuslicer.close()
    ifupost.close()
    return Gwa2Slit(slit_models)


def gwa_to_slit(slits_id, disperser, wrange, order, reference_files):
    """
    GWA to SLIT transform.

    Parameters
    ----------
    slits_id : list
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
    slit_models = {}
    for i in range(1, 6):
        slit_names = slits_id[slits_id[:, 0] == i]
        if slit_names.any():
            msa_model = msa.tree[i]['model']
            for slit in slit_names:
                index = slit[1] - 1
                slitdata = msa.tree[slit[0]]['data'][index]
                slitdata_model = get_slit_location_model(slitdata)
                msa_transform = slitdata_model | msa_model
                msa2gwa = (msa_transform | collimator2gwa)
                gwa2msa = gwa_to_ymsa(msa2gwa)# TODO: Use model sets here
                bgwa2msa = Mapping((0, 1, 0, 1), n_inputs=3) | \
                    Const1D(0) * Identity(1) & Const1D(-1) * Identity(1) & Identity(2) | \
                    Identity(1) & gwa2msa & Identity(2) | \
                    Mapping((0, 1, 0, 1, 2, 3)) | Identity(2) & msa2gwa & Identity(2) | \
                    Mapping((0, 1, 2, 5), n_inputs=7) | Identity(2) & lgreq

                # msa to before_gwa
                msa2bgwa = msa2gwa & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
                bgwa2msa.inverse = msa2bgwa
                s = slitid_to_slit(np.array([slit]))[0]
                slit_models[s] = bgwa2msa
    msa.close()
    return Gwa2Slit(slit_models)


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
    model = (fpa | camera | u2dircos | rotation)
    return model


def compute_domain(slit2detector, wavelength_range):
    """
    Compute the projection of a slit on the detector.

    Parameters
    ----------
    input_model : `jwst.datamodels.DataModel`
        The input data model - either an ImageModel or a CubeModel.
    open_slits : list
        The open slitlets in the observation.
    reference_files : dict
        {reftype: file_name} mapping
        Reference files for the observation returned by CRDS.

    """
    x = [-.55, .55]
    y = [0, 0]

    lam_min = wavelength_range[0]
    lam_max = wavelength_range[1]

    x_range_low = slit2detector(x, y, lam_min)[0]
    x_range_high = slit2detector(x, y, lam_max)[0]
    x_range_low = np.where(x_range_low > 0, x_range_low, 0)
    x_range_high = np.where(x_range_high < 2048, x_range_high, 2048)
    x0, x1 = x_range_low.min(), x_range_high.max()
    y = [-.55, .55]
    x = [0, 0]
    y_range_min = slit2detector(x, y, lam_min)[1]
    y_range_max = slit2detector(x, y, lam_max)[1]
    y0 = min(y_range_min.min(), y_range_max.min())
    y1 = max(y_range_min.max(), y_range_max.max())

    domain = [{'lower': int(x0), 'upper': int(x1)}, {'lower': int(y0), 'upper': int(y1)}]
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
    # Convert the wavelength to microns
    return fore2ote_mapping | (ote & Identity(1) / Const1D(1e-6))


def create_frames():
    """
    Create the coordinate frames in the NIRSPEC WCS pipeline.

    These are
    "detector", "gwa", "slit_frame", "msa_frame", "oteip", "v2v3", "world".
    """
    det = cf.Frame2D(name='detector', axes_order=(0, 1))
    gwa = cf.Frame2D(name="gwa", axes_order=(0, 1), unit=(u.rad, u.rad),
                      axes_names=('alpha_in', 'beta_in'))
    msa_spatial = cf.Frame2D(name='msa_spatial', axes_order=(0, 1), unit=(u.m, u.m),
                             axes_names=('x_msa', 'y_msa'))
    slit_spatial = cf.Frame2D(name='slit_spatial', axes_order=(0, 1), unit=("", ""),
                             axes_names=('x_slit', 'y_slit'))
    sky = cf.CelestialFrame(name='sky', axes_order=(0, 1), reference_frame=coord.ICRS())
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.m,),
                            axes_names=('wavelength',))
    v2v3_spatial = cf.Frame2D(name='V2V3_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec),
                             axes_names=('V2', 'V3'))
    v2v3 = cf.CompositeFrame([v2v3_spatial, spec], name='v2v3')
    slit_frame = cf.CompositeFrame([slit_spatial, spec], name='slit_frame')
    msa_frame = cf.CompositeFrame([msa_spatial, spec], name='msa_frame')
    oteip_spatial = cf.Frame2D(name='OTEIP_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec),
                               axes_names=('X_OTEIP', 'Y_OTEIP'))
    oteip = cf.CompositeFrame([oteip_spatial, spec], name='oteip')
    world = cf.CompositeFrame([sky, spec], name='world')
    return det, gwa, slit_frame, msa_frame, oteip, v2v3, world


def create_imaging_frames():
    """
    Create the coordinate frames in the NIRSPEC WCS pipeline.
    These are
    "detector", "gwa", "msa_frame", "oteip", "v2v3", "world".
    """
    det = cf.Frame2D(name='detector', axes_order=(0, 1))
    gwa = cf.Frame2D(name="gwa", axes_order=(0, 1), unit=(u.rad, u.rad),
                      axes_names=('alpha_in', 'beta_in'))
    msa = cf.Frame2D(name='msa', axes_order=(0, 1), unit=(u.m, u.m),
                             axes_names=('x_msa', 'y_msa'))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.arcsec, u.arcsec),
                              axes_names=('v2', 'v3'))
    oteip = cf.Frame2D(name='oteip', axes_order=(0, 1), unit=(u.arcsec, u.arcsec),
                               axes_names=('x_oteip', 'y_oteip'))
    return det, gwa, msa, oteip, v2v3


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


def nrs_wcs_set_input(wcsobj, quadrant, slitid, wavelength_range):
    """
    Returns a WCS object for this slit.

    Parameters
    ----------
    wcsobj : `~gwcs.WCS`
        A WCS object for the all open slitlets in an observation.
    quadrant : int
        MSA Quadrant number.
    slit : int
        Slit id within this quadrant.
    wavelength_range: list
        Wavelength range for the combination of fliter and grating.

    Returns
    -------
    wcsobj : `~gwcs.wcs.WCS`
        WCS object for this slit.
    """
    import copy # TODO: Add a copy method to gwcs.WCS
    slit = slitid_to_slit(np.array([(quadrant, slitid)]))[0]
    slit_wcs = copy.deepcopy(wcsobj)
    slit_wcs.set_transform('detector', 'gwa', wcsobj.pipeline[0][1][1:])
    #slit_wcs.set_transform('detector', 'gwa', wcsobj.pipeline[0][1])
    slit_wcs.set_transform('gwa', 'slit_frame', wcsobj.pipeline[1][1].models[slit])
    slit_wcs.set_transform('slit_frame', 'msa_frame', wcsobj.pipeline[2][1][1].models[slit] & Identity(1))

    slit2detector = slit_wcs.get_transform('slit_frame', 'detector')

    domain = compute_domain(slit2detector, wavelength_range)
    slit_wcs.domain = domain
    return slit_wcs


def slit_to_detector(input_model, slits_id, lam, reference_files):
    """
    For each slit computes the transform from the MSA to the detector.

    Used to flag stuck open shutters.
    Reftypes needed:
    reftypes = ['fpa', 'camera', 'disperser', 'collimator', 'msa', 'wavelengthrange']

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        Input data model
    slits : list
        A list of slit IDs. A slit ID is a tuple of (quadrant, slit_nimber)
    lam : float
        wavelength, in meters
    reference_files : dict
        A dictionary with WCS reference_files, {reftype: refname}

    Returns
    -------
    msa2det : dict
        A dictionary with the MSA to detector transform for each slit, {slit_id: astropy.modeling.Model}
    """
    slits_msa2det = {}
    # read models from reference files
    with AsdfFile.open(reference_files['fpa']) as f:
        fpa = f.tree[input_model.meta.instrument.detector].copy()
    with AsdfFile.open(reference_files['camera']) as f:
        camera = f.tree['model'].copy()
    with AsdfFile.open(reference_files['collimator']) as f:
        collimator = f.tree['model'].copy()
    with AsdfFile.open(reference_files['disperser']) as f:
        disperser = f.tree
    with AsdfFile.open(reference_files['wavelengthrange']) as f:
        wave_range = f.tree

    dircos2u = DirCos2Unitless(name='directional2unitless_cosines')
    disperser = correct_tilt(disperser, input_model.meta.instrument.gwa_xtilt,
                             input_model.meta.instrument.gwa_ytilt)
    angles = [disperser['theta_x'], disperser['theta_y'],
               disperser['theta_z'], disperser['tilt_y']]
    rotation = Rotation3DToGWA(angles, axes_order="xyzy", name='rotaton')

    order, wrange = get_spectral_order_wrange(input_model, reference_files['wavelengthrange'])
    input_model.meta.wcsinfo.waverange_start = wrange[0]
    input_model.meta.wcsinfo.waverange_end = wrange[1]
    input_model.meta.wcsinfo.spectral_order = sporder

    agreq = AngleFromGratingEquation(disperser['groove_density'], order, name='alpha_from_greq')
    # GWA to detector
    gwa2det = rotation.inverse | dircos2u | camera.inverse | fpa.inverse
    # collimator to GWA
    collimator2gwa = collimator_to_gwa(reference_files, disperser)


    msa = AsdfFile.open(reference_files['msa'])
    for i in range(1, 5):
        slit_names = slits_id[slits_id[:, 0] == i]
        if slit_names.any():
            msa_model = msa.tree[i]['model']
            for slit in slit_names:
                index = slit[1] - 1
                slitdata = msa.tree[slit[0]]['data'][index]
                slitdata_model = get_slit_location_model(slitdata)
                # absolute positions of the slit on the MSA
                msa_transform = slitdata_model | msa_model
                s = slitid_to_slit(np.array([slit]))[0]
                msa2gwa = (msa_transform | collimator2gwa) & Identity(1) | Mapping((3, 0, 1, 2)) | agreq
                slits_msa2det[s] = msa2gwa | gwa2det
    msa.close()

    return slits_msa2det


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
        wcs_list.append(nrs_wcs_set_input(input_model.meta.wcs, 0, i, wrange))
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
