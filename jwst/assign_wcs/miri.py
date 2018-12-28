import os.path
import logging
import numpy as np
from astropy.modeling import models
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits

import gwcs.coordinate_frames as cf
from gwcs import selector
from gwcs.utils import _toindex
from . import pointing
from ..transforms import models as jwmodels
from .util import (not_implemented_mode, subarray_transform,
                   velocity_correction, bounding_box_from_model, bounding_box_from_subarray)
from ..datamodels import (DistortionModel, FilteroffsetModel,
                          DistortionMRSModel, WavelengthrangeModel,
                          RegionsModel, SpecwcsModel)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging", "lrs", "ifu"]


def create_pipeline(input_model, reference_files):
    """
    Create the WCS pipeline for MIRI modes.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`, `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.CubeModel`
        Data model.
    reference_files : dict
        {reftype: reference file name} mapping.

    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)
    if pipeline:
        log.info("Created a MIRI {0} pipeline with references {1}".format(
            exp_type, reference_files))
    return pipeline


def imaging(input_model, reference_files):
    """
    The MIRI Imaging WCS pipeline.

    It includes three coordinate frames -
    "detector", "v2v3" and "world".

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.
        Uses "distortion" and "filteroffset" reference files.

    """

    # Create the Frames
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')

    # Create the transforms
    distortion = imaging_distortion(input_model, reference_files)
    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        distortion = subarray2full | distortion
        distortion.bounding_box = bounding_box_from_subarray(input_model)
    else:
        # TODO: remove setting the bounding box when it is set in the new ref file.
        try:
            bb = distortion.bounding_box
        except NotImplementedError:
            shape = input_model.data.shape
            # Note: Since bounding_box is attached to the model here it's in reverse order.
            bb = ((-0.5, shape[0] - 0.5), (3.5, shape[1] - 4.5))
            distortion.bounding_box = bb

    tel2sky = pointing.v23tosky(input_model)

    # Create the pipeline
    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)
                ]

    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the "detector" to "v2v3" transform for the MIRI Imager.

    1. Filter dependent shift in (x,y) (!with an oposite
       sign to that delivered by the IT) (uses the "filteroffset" ref file)
    2. Apply MI (uses "distortion" ref file)
    3. Apply Ai and BI matrices (uses "distortion" ref file)
    4. Apply the TI matrix (this gives Xan/Yan coordinates) (uses "distortion" ref file)
    5. Aply the XanYan --> V2V3 transform (uses "distortion" ref file)
    6. Apply V2V3 --> sky transform

    """
    # Read in the distortion.
    with DistortionModel(reference_files['distortion']) as dist:
        distortion = dist.model

    # Check if the transform in the reference file has a ``bounding_box``.
    # If not set a ``bounding_box`` equal to the size of the image.
    try:
        distortion.bounding_box
    except NotImplementedError:
        distortion.bounding_box = bounding_box_from_model(input_model)

    # Add an offset for the filter
    obsfilter = input_model.meta.instrument.filter
    with FilteroffsetModel(reference_files['filteroffset']) as filter_offset:
        filters = filter_offset.filters

    col_offset = None
    row_offset = None
    for f in filters:
        if f.name == obsfilter:
            col_offset = f.column_offset
            row_offset = f.row_offset
            break

    if (col_offset is not None) and (row_offset is not None):
        distortion = models.Shift(col_offset) & models.Shift(row_offset) | distortion

    return distortion


def lrs(input_model, reference_files):
    """
    The LRS-FIXEDSLIT and LRS-SLITLESS WCS pipeline.

    It has two coordinate frames: "detecor" and "world".
    Uses the "specwcs" and "distortion" reference files.

    """
    # Setup the frames.
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    spec = cf.SpectralFrame(name='wavelength', axes_order=(2,), unit=(u.micron,),
                            axes_names=('lambda',))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS(), name='sky')
    v2v3_spatial = cf.Frame2D(name='v2v3_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec))
    v2v3 = cf.CompositeFrame(name="v2v3", frames=[v2v3_spatial, spec])
    world = cf.CompositeFrame(name="world", frames=[sky, spec])

    # Determine the distortion model.
    subarray2full = subarray_transform(input_model)
    with DistortionModel(reference_files['distortion']) as dist:
        distortion = dist.model

    if subarray2full is not None:
        distortion = subarray2full | distortion

    # Incorporate the small rotation
    angle = np.arctan(0.00421924)
    rotation = models.Rotation2D(angle)
    distortion = distortion | rotation

    # Load and process the reference data.
    with fits.open(reference_files['specwcs']) as ref:
        lrsdata = np.array([l for l in ref[1].data])

        # Get the zero point from the reference data.
        # The zero_point is X, Y  (which should be COLUMN, ROW)
        # TODO: Are imx, imy 0- or 1-indexed?  We are treating them here as
        # 0-indexed.  Since they are FITS, they are probably 1-indexed.
        if input_model.meta.exposure.type.lower() == 'mir_lrs-fixedslit':
            zero_point = ref[1].header['imx'], ref[1].header['imy']
        elif input_model.meta.exposure.type.lower() == 'mir_lrs-slitless':
            zero_point = ref[1].header['imxsltl'], ref[1].header['imysltl']
            #zero_point = [35, 442]  # [35, 763] # account for subarray

    # Create the bounding_box
    x0 = lrsdata[:, 3]
    y0 = lrsdata[:, 4]
    x1 = lrsdata[:, 5]

    bb = ((x0.min() - 0.5 + zero_point[0], x1.max() + 0.5 + zero_point[0]),
          (y0.min() - 0.5 + zero_point[1], y0.max() + 0.5 + zero_point[1]))

    # Find the ROW of the zero point which should be the [1] of zero_point
    row_zero_point = zero_point[1]

    # Compute the v2v3 to sky.
    tel2sky = pointing.v23tosky(input_model)

    # Compute the spatial detector to V2V3 transform
    # Take a row centered on zero_point_y and convert it to v2, v3.
    # The forward transform uses constant ``y`` values for each ``x``.
    # The inverse transform uses constant ``v3`` values for each ``v2``.
    x = np.arange(0, input_model.data.shape[1], dtype=np.float)
    y = np.zeros(x.shape, dtype=np.float) + row_zero_point

    zero_point_v2v3 = distortion(*zero_point)

    spatial_forward = models.Identity(1) & models.Const1D(zero_point[1]) | distortion
    spatial_forward.inverse = (models.Identity(1) & models.Const1D(zero_point_v2v3[1]) |
                               distortion.inverse)

    # Create the spectral transforms.
    lrs_wav_model = jwmodels.LRSWavelength(lrsdata, zero_point)

    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    else:
        if velosys is not None:
            velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
            lrs_wav_model = lrs_wav_model | velocity_corr
            log.info("Applied Barycentric velocity correction : {}".format(velocity_corr[1].amplitude.value))

    det_to_v2v3 = models.Mapping((0, 1, 0, 1)) | spatial_forward & lrs_wav_model
    det_to_v2v3.bounding_box = bb[::-1]
    v23_to_world = tel2sky & models.Identity(1)

    # Now the actual pipeline.
    pipeline = [(detector, det_to_v2v3),
                (v2v3, v23_to_world),
                (world, None)
                ]
    return pipeline


def ifu(input_model, reference_files):
    """
    The MIRI MRS WCS pipeline.

    It has the following coordinate frames:
    "detector", "alpha_beta", "v2v3", "world".

    It uses the "distortion", "regions", "specwcs"
    and "wavelengthrange" reference files.
    """
    # Define coordinate frames.
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    alpha_beta = cf.Frame2D(name='alpha_beta_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=('alpha', 'beta'))
    spec_local = cf.SpectralFrame(name='alpha_beta_spectral', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    miri_focal = cf.CompositeFrame([alpha_beta, spec_local], name='alpha_beta')
    v23_spatial = cf.Frame2D(name='V2_V3_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=('v2', 'v3'))
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    v2v3 = cf.CompositeFrame([v23_spatial, spec], name='v2v3')
    icrs = cf.CelestialFrame(name='icrs', reference_frame=coord.ICRS(),
                             axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=('RA', 'DEC'))
    world = cf.CompositeFrame([icrs, spec], name='world')

    # Define the actual transforms
    det2abl = (detector_to_abl(input_model, reference_files)).rename(
        "detector_to_abl")
    abl2v2v3l = (abl_to_v2v3l(input_model, reference_files)).rename("abl_to_v2v3l")

    tel2sky = pointing.v23tosky(input_model) & models.Identity(1)

    # Put the transforms together into a single transform
    shape = input_model.data.shape
    det2abl.bounding_box = ((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5))
    pipeline = [(detector, det2abl),
                (miri_focal, abl2v2v3l),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def detector_to_abl(input_model, reference_files):
    """
    Create the transform from "detector" to "alpha_beta" frame.

    Transform description:
    forward transform
      RegionsSelector
        label_mapper is the regions array
        selector is {slice_number: alpha_model & beta_model & lambda_model}
    backward transform
      RegionsSelector
        label_mapper is LabelMapperDict
           {channel_wave_range (): LabelMapperDict}
                                   {beta: slice_number}
        selector is {slice_number: x_transform & y_transform}
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range

    with DistortionMRSModel(reference_files['distortion']) as dist:
        alpha_model = dist.alpha_model
        beta_model = dist.beta_model
        x_model = dist.x_model
        y_model = dist.y_model
        bzero = dict(zip(dist.bzero.channel_band, dist.bzero.beta_zero))
        bdel = dict(zip(dist.bdel.channel_band, dist.bdel.delta_beta))
        slices = dist.slices

    with SpecwcsModel(reference_files['specwcs']) as f:
        lambda_model = f.model

    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    else:
        if velosys is not None:
            velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
            lambda_model = [m | velocity_corr for m in lambda_model]
            log.info("Applied Barycentric velocity correction : {}".format(velocity_corr[1].amplitude.value))

    with RegionsModel(reference_files['regions']) as f:
        regions = f.regions.copy()

    label_mapper = selector.LabelMapperArray(regions)
    transforms = {}

    for i, sl in enumerate(slices):
        forward = models.Mapping([1, 0, 0, 1, 0]) | \
                alpha_model[i] & beta_model[i] & lambda_model[i]
        inv = models.Mapping([2, 0, 2, 0]) | x_model[i] & y_model[i]
        forward.inverse = inv
        transforms[sl] = forward

    with WavelengthrangeModel(reference_files['wavelengthrange']) as f:
        wr = dict(zip(f.waverange_selector, f.wavelengthrange))

    ch_dict = {}
    for c in channel:
        cb = c + band
        mapper = jwmodels.MIRI_AB2Slice(bzero[cb], bdel[cb], c)
        lm = selector.LabelMapper(inputs=('alpha', 'beta', 'lam'),
                                  mapper=mapper, inputs_mapping=models.Mapping((1,), n_inputs=3))
        ch_dict[tuple(wr[cb])] = lm


    alpha_beta_mapper = selector.LabelMapperRange(('alpha', 'beta', 'lam'), ch_dict,
                                                  models.Mapping((2,)))
    label_mapper.inverse = alpha_beta_mapper


    det2alpha_beta = selector.RegionsSelector(('x', 'y'), ('alpha', 'beta', 'lam'),
                                              label_mapper=label_mapper, selector=transforms)
    return det2alpha_beta


def abl_to_v2v3l(input_model, reference_files):
    """
    Create the transform from "alpha_beta" to "v2v3" frame.

    Transform description:
    forward transform
      RegionsSelector
        label_mapper is LabelMapperDict()
        {channel_wave_range (): channel_number}
        selector is {channel_number: ab2v2 & ab2v3}
    bacward_transform
      RegionsSelector
        label_mapper is LabelMapperDict()
        {channel_wave_range (): channel_number}
        selector is {channel_number: v22ab & v32ab}
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range
    channels = [c + band for c in channel]

    with DistortionMRSModel(reference_files['distortion']) as dist:
        v23 = dict(zip(dist.abv2v3_model.channel_band, dist.abv2v3_model.model))

    with WavelengthrangeModel(reference_files['wavelengthrange']) as f:
        wr = dict(zip(f.waverange_selector, f.wavelengthrange))

    dict_mapper = {}
    sel = {}
    # Since there are two channels in each reference file we need to loop over them
    for c in channels:
        ch = int(c[0])
        dict_mapper[tuple(wr[c])] = models.Mapping((2,), name="mapping_lam") | \
                   models.Const1D(ch, name="channel #")
        ident1 = models.Identity(1, name='identity_lam')
        ident1._inputs = ('lam',)
        chan_v23 = v23[c]
        v23chan_backward = chan_v23.inverse
        del chan_v23.inverse
        v23_spatial = chan_v23
        v23_spatial.inverse = v23chan_backward
        # Tack on passing the third wavelength component
        v23c = v23_spatial & ident1
        sel[ch] = v23c

    wave_range_mapper = selector.LabelMapperRange(('alpha', 'beta', 'lam'), dict_mapper,
                                                  inputs_mapping=models.Mapping([2, ]))
    wave_range_mapper.inverse = wave_range_mapper.copy()
    abl2v2v3l = selector.RegionsSelector(('alpha', 'beta', 'lam'), ('v2', 'v3', 'lam'),
                                      label_mapper=wave_range_mapper,
                                      selector=sel)

    return abl2v2v3l


exp_type2transform = {'mir_image': imaging,
                      'mir_tacq': imaging,
                      'mir_lyot': imaging,
                      'mir_4qpm': imaging,
                      'mir_coroncal': imaging,
                      'mir_lrs-fixedslit': lrs,
                      'mir_lrs-slitless': lrs,
                      'mir_mrs': ifu,
                      'mir_flatmrs': not_implemented_mode,
                      'mir_flatimage': not_implemented_mode,
                      'mir_flat-mrs': not_implemented_mode,
                      'mir_flat-image': not_implemented_mode,
                      'mir_dark': not_implemented_mode,
                      }


def get_wavelength_range(input_model, path=None):
    """
    Return the wavelength range used for computing the WCS.

    Needs access to the reference file used to construct the WCS object.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model after assign_wcs has been run.
    path : str
        Directory where the reference file is. (optional)
    """
    fname = input_model.meta.ref_file.wavelengthrange.name.split('/')[-1]
    if path is None and not os.path.exists(fname):
        raise IOError("Reference file {0} not found. Please specify a path.".format(fname))
    else:
        fname = os.path.join(path, fname)
        f = WavelengthrangeModel(fname)

    wave_range = f.tree['wavelengthrange'].copy()
    wave_channels = f.tree['channels']
    f.close()

    wr = dict(zip(wave_channels, wave_range))
    channel = input_model.meta.instrument.channel
    band = input_model.meta.instrument.band

    return dict([(ch + band, wr[ch + band]) for ch in channel])
