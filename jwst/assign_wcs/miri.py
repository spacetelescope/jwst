import os.path
import logging
import numpy as np
from astropy.modeling import models
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits

from scipy.interpolate import UnivariateSpline
import gwcs.coordinate_frames as cf
from gwcs import selector

from stdatamodels.jwst.datamodels import (DistortionModel, FilteroffsetModel,
                                          DistortionMRSModel, WavelengthrangeModel,
                                          RegionsModel, SpecwcsModel)
from stdatamodels.jwst.transforms.models import (MIRI_AB2Slice, IdealToV2V3)

from . import pointing
from .util import (not_implemented_mode, subarray_transform,
                   velocity_correction, transform_bbox_from_shape,
                   bounding_box_from_subarray)


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging", "lrs", "ifu"]


def create_pipeline(input_model, reference_files):
    """
    Create the WCS pipeline for MIRI modes.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`, `~jwst.datamodels.IFUImageModel`,
                  `~jwst.datamodels.CubeModel`
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
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), axes_names=('v2', 'v3'),
                      unit=(u.arcsec, u.arcsec))
    v2v3vacorr = cf.Frame2D(name='v2v3vacorr', axes_order=(0, 1),
                            axes_names=('v2', 'v3'), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')

    # Create the transforms
    distortion = imaging_distortion(input_model, reference_files)
    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        distortion = subarray2full | distortion
        distortion.bounding_box = bounding_box_from_subarray(input_model)
    else:
        # TODO: remove setting the bounding box if it is set in the new ref file.
        try:
            bb = distortion.bounding_box
        except NotImplementedError:
            shape = input_model.data.shape
            # Note: Since bounding_box is attached to the model here it's in reverse order.
            bb = ((-0.5, shape[0] - 0.5), (3.5, shape[1] - 4.5))
            distortion.bounding_box = bb

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref
    )

    tel2sky = pointing.v23tosky(input_model)

    # Create the pipeline
    pipeline = [(detector, distortion),
                (v2v3, va_corr),
                (v2v3vacorr, tel2sky),
                (world, None)
                ]

    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the "detector" to "v2v3" transform for the MIRI Imager.

    1. Filter dependent shift in (x,y) (!with an opposite
       sign to that delivered by the IT) (uses the "filteroffset" ref file)
    2. Apply MI (uses "distortion" ref file)
    3. Apply Ai and BI matrices (uses "distortion" ref file)
    4. Apply the TI matrix (this gives V2/V3 coordinates) (uses "distortion" ref file)
    5. Apply V2V3 --> sky transform

    """
    # Read in the distortion.
    with DistortionModel(reference_files['distortion']) as dist:
        distortion = dist.model

    # Check if the transform in the reference file has a ``bounding_box``.
    # If not set a ``bounding_box`` equal to the size of the image.
    try:
        bbox = distortion.bounding_box
    except NotImplementedError:
        bbox = distortion.bounding_box = None

    # Add an offset for the filter
    obsfilter = input_model.meta.instrument.filter
    with FilteroffsetModel(reference_files['filteroffset']) as filter_offset:
        filters = filter_offset.filters

    col_offset = None
    row_offset = None
    for f in filters:
        if f.filter == obsfilter:
            col_offset = f.column_offset
            row_offset = f.row_offset
            break

    if col_offset is not None and row_offset is not None:
        distortion = models.Shift(col_offset) & models.Shift(row_offset) | distortion
    if bbox is None:
        distortion.bounding_box = transform_bbox_from_shape(input_model.data.shape)
    else:
        distortion.bounding_box = bbox
    return distortion


def lrs(input_model, reference_files):
    """
    The LRS-FIXEDSLIT and LRS-SLITLESS WCS pipeline.

    Notes
    -----
    It includes three coordinate frames -
    "detector", "v2v3" and "world".
    "v2v3" and "world" each have (spatial, spatial, spectral) components.

    Uses the "specwcs" and "distortion" reference files.

    """
    # Define the various coordinate frames.
    # Original detector frame
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    # Spectral component
    spec = cf.SpectralFrame(name='spec', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    # v2v3 spatial component
    v2v3_spatial = cf.Frame2D(name='v2v3_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec),
                              axes_names=('v2', 'v3'))
    v2v3vacorr_spatial = cf.Frame2D(
        name='v2v3vacorr_spatial',
        axes_order=(0, 1),
        unit=(u.arcsec, u.arcsec),
        axes_names=('v2', 'v3')
    )

    # v2v3 spatial+spectra
    v2v3 = cf.CompositeFrame([v2v3_spatial, spec], name='v2v3')
    v2v3vacorr = cf.CompositeFrame([v2v3vacorr_spatial, spec], name='v2v3vacorr')

    # 'icrs' frame which is the spatial sky component
    icrs = cf.CelestialFrame(name='icrs', reference_frame=coord.ICRS(),
                             axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=('RA', 'DEC'))
    # Final 'world' composite frame with spatial and spectral components
    world = cf.CompositeFrame(name="world", frames=[icrs, spec])

    # Create the transforms
    dettotel = lrs_distortion(input_model, reference_files)
    v2v3tosky = pointing.v23tosky(input_model)
    teltosky = v2v3tosky & models.Identity(1)

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref
    ) & models.Identity(1)

    # Put the transforms together into a single pipeline
    pipeline = [(detector, dettotel),
                (v2v3, va_corr),
                (v2v3vacorr, teltosky),
                (world, None)]

    return pipeline


def lrs_distortion(input_model, reference_files):
    """
    The LRS-FIXEDSLIT and LRS-SLITLESS WCS pipeline.

    Transform from subarray (x, y) to (v2, v3, lambda) using
    the "specwcs" and "distortion" reference files.

    """

    # subarray to full array transform
    subarray2full = subarray_transform(input_model)

    # full array to v2v3 transform for the ordinary imager
    with DistortionModel(reference_files['distortion']) as dist:
        distortion = dist.model

    # Combine models to create subarray to v2v3 distortion
    if subarray2full is not None:
        subarray_dist = subarray2full | distortion
    else:
        subarray_dist = distortion

    ref = fits.open(reference_files['specwcs'])

    with ref:
        lrsdata = np.array([d for d in ref[1].data])
        # Get the zero point from the reference data.
        # The zero_point is X, Y  (which should be COLUMN, ROW)
        # These are 1-indexed in CDP-7 (i.e., SIAF convention) so must be converted to 0-indexed
        if input_model.meta.exposure.type.lower() == 'mir_lrs-fixedslit':
            zero_point = ref[0].header['imx'] - 1, ref[0].header['imy'] - 1
        elif input_model.meta.exposure.type.lower() == 'mir_lrs-slitless':
            zero_point = ref[0].header['imxsltl'] - 1, ref[0].header['imysltl'] - 1
            # Transform to slitless subarray from full array
            zero_point = subarray2full.inverse(zero_point[0], zero_point[1])

    # In the lrsdata reference table, X_center,y_center,wavelength describe the location of the
    # centroid trace along the detector in pixels relative to nominal location.
    # x0,y0(ul) x1,y1 (ur) x2,y2(lr) x3,y3(ll) define corners of the box within which the distortion
    # and wavelength calibration was derived
    xcen = lrsdata[:, 0]
    ycen = lrsdata[:, 1]
    wavetab = lrsdata[:, 2]
    x0 = lrsdata[:, 3]
    y0 = lrsdata[:, 4]
    x1 = lrsdata[:, 5]
    y2 = lrsdata[:, 8]

    # If in fixed slit mode, define the bounding box using the corner locations provided in
    # the CDP reference file.
    if input_model.meta.exposure.type.lower() == 'mir_lrs-fixedslit':

        bb_sub = ((np.floor(x0.min() + zero_point[0]) - 0.5, np.ceil(x1.max() + zero_point[0]) + 0.5),
                  (np.floor(y2.min() + zero_point[1]) - 0.5, np.ceil(y0.max() + zero_point[1]) + 0.5))

    # If in slitless mode, define the bounding box X locations using the subarray x boundaries
    # and the y locations using the corner locations in the CDP reference file.  Make sure to
    # omit the 4 reference pixels on the left edge of slitless subarray.
    if input_model.meta.exposure.type.lower() == 'mir_lrs-slitless':
        bb_sub = ((input_model.meta.subarray.xstart - 1 + 4 - 0.5, input_model.meta.subarray.xsize - 1 + 0.5),
                  (np.floor(y2.min() + zero_point[1]) - 0.5, np.ceil(y0.max() + zero_point[1]) + 0.5))

    # Find the ROW of the zero point
    row_zero_point = zero_point[1]

    # The inputs to the "detector_to_v2v3" transform are
    # - the indices in x spanning the entire image row
    # - y is the y-value of the zero point
    # This is equivalent of making a vector of x, y locations for
    # every pixel in the reference row
    const1d = models.Const1D(row_zero_point)
    const1d.inverse = models.Const1D(row_zero_point)
    det_to_v2v3 = models.Identity(1) & const1d | subarray_dist

    # Now deal with the fact that the spectral trace isn't perfectly up and down along detector.
    # This information is contained in the xcenter/ycenter values in the CDP table, but we'll handle it
    # as a simple rotation using a linear fit to this relation provided by the CDP.

    z = np.polyfit(xcen, ycen, 1)
    slope = 1. / z[0]
    traceangle = np.arctan(slope) * 180. / np.pi  # trace angle in degrees
    rot = models.Rotation2D(traceangle)  # Rotation model

    # Now include this rotation in our overall transform
    # First shift to a frame relative to the trace zeropoint, then apply the rotation
    # to correct for the curved trace.  End in a rotated frame relative to zero at the reference point
    # and where yrot is aligned with the spectral trace)
    xysubtoxyrot = models.Shift(-zero_point[0]) & models.Shift(-zero_point[1]) | rot

    # Next shift back to the subarray frame, and then map to v2v3
    xyrottov2v3 = models.Shift(zero_point[0]) & models.Shift(zero_point[1]) | det_to_v2v3

    # The two models together
    xysubtov2v3 = xysubtoxyrot | xyrottov2v3

    # Work out the spectral component of the transform
    # First compute the reference trace in the rotated-Y frame
    xcenrot, ycenrot = rot(xcen, ycen)
    # The input table of wavelengths isn't perfect, and the delta-wavelength
    # steps show some unphysical behaviour
    # Therefore fit with a spline for the ycenrot->wavelength transform
    # Reverse vectors so that yinv is increasing (needed for spline fitting function)
    yrev = ycenrot[::-1]
    wrev = wavetab[::-1]
    # Spline fit with enforced smoothness
    spl = UnivariateSpline(yrev, wrev, s=0.002)
    # Evaluate the fit at the rotated-y reference points
    wavereference = spl(yrev)
    # wavereference now contains the wavelengths corresponding to regularly-sampled ycenrot, create the model
    wavemodel = models.Tabular1D(lookup_table=wavereference, points=yrev, name='waveref',
                                 bounds_error=False, fill_value=np.nan)

    # Now construct the inverse spectral transform.
    # First we need to create a spline going from wavereference -> ycenrot
    spl2 = UnivariateSpline(wavereference[::-1], ycenrot, s=0.002)
    # Make a uniform grid of wavelength points from min to max, sampled according
    # to the minimum delta in the input table
    dw = np.amin(np.absolute(np.diff(wavereference)))
    wmin = np.amin(wavereference)
    wmax = np.amax(wavereference)
    wgrid = np.arange(wmin, wmax, dw)
    # Evaluate the rotated y locations of the grid
    ygrid = spl2(wgrid)
    # ygrid now contains the rotated y pixel locations corresponding to
    # regularly-sampled wavelengths, create the model
    wavemodel.inverse = models.Tabular1D(lookup_table=ygrid, points=wgrid, name='waverefinv',
                                         bounds_error=False, fill_value=np.nan)

    # Wavelength barycentric correction
    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    else:
        if velosys is not None:
            velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
            wavemodel = wavemodel | velocity_corr
            log.info("Applied Barycentric velocity correction : {}".format(velocity_corr[1].amplitude.value))

    # Construct the full distortion model (xsub,ysub -> v2,v3,wavelength)
    lrs_wav_model = xysubtoxyrot | models.Mapping([1], n_inputs=2) | wavemodel
    dettotel = models.Mapping((0, 1, 0, 1)) | xysubtov2v3 & lrs_wav_model

    # Construct the inverse distortion model (v2,v3,wavelength -> xsub,ysub)
    # Transform to get xrot from v2,v3
    v2v3toxrot = subarray_dist.inverse | xysubtoxyrot | models.Mapping([0], n_inputs=2)
    # wavemodel.inverse gives yrot from wavelength
    # v2,v3,lambda -> xrot,yrot
    xform1 = v2v3toxrot & wavemodel.inverse
    dettotel.inverse = xform1 | xysubtoxyrot.inverse

    # Bounding box is the subarray bounding box, because we're assuming subarray coordinates passed in
    dettotel.bounding_box = bb_sub[::-1]

    return dettotel


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
    alpha_beta = cf.Frame2D(name='alpha_beta_spatial', axes_order=(0, 1),
                            unit=(u.arcsec, u.arcsec), axes_names=('alpha', 'beta'))
    spec_local = cf.SpectralFrame(name='alpha_beta_spectral', axes_order=(2,),
                                  unit=(u.micron,), axes_names=('lambda',))
    miri_focal = cf.CompositeFrame([alpha_beta, spec_local], name='alpha_beta')
    v23_spatial = cf.Frame2D(name='v2v3_spatial', axes_order=(0, 1),
                             unit=(u.arcsec, u.arcsec), axes_names=('v2', 'v3'))
    v2v3vacorr_spatial = cf.Frame2D(name='v2v3vacorr_spatial', axes_order=(0, 1),
                                    unit=(u.arcsec, u.arcsec), axes_names=('v2', 'v3'))

    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,),
                            axes_names=('lambda',))
    v2v3 = cf.CompositeFrame([v23_spatial, spec], name='v2v3')
    v2v3vacorr = cf.CompositeFrame([v2v3vacorr_spatial, spec], name='v2v3vacorr')
    icrs = cf.CelestialFrame(name='icrs', reference_frame=coord.ICRS(),
                             axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=('RA', 'DEC'))
    world = cf.CompositeFrame([icrs, spec], name='world')

    # Define the actual transforms
    det2abl = (detector_to_abl(input_model, reference_files)).rename(
        "detector_to_abl")
    abl2v2v3l = (abl_to_v2v3l(input_model, reference_files)).rename("abl_to_v2v3l")

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref
    ) & models.Identity(1)

    tel2sky = pointing.v23tosky(input_model) & models.Identity(1)

    # Put the transforms together into a single transform
    det2abl.bounding_box = transform_bbox_from_shape(input_model.data.shape)
    pipeline = [(detector, det2abl),
                (miri_focal, abl2v2v3l),
                (v2v3, va_corr),
                (v2v3vacorr, tel2sky),
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
        allregions = f.regions.copy()
        # Use the 80% throughput slice mask
        regions = allregions[7, :, :]

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
        mapper = MIRI_AB2Slice(bzero[cb], bdel[cb], c)
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
    backward_transform
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
        dict_mapper[tuple(wr[c])] = (models.Mapping((2,), name="mapping_lam") |
                                     models.Const1D(ch, name="channel #"))
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
                      'mir_taconfirm': imaging,
                      }


def get_wavelength_range(input_model, path=None):
    """
    Return the wavelength range used for computing the WCS.

    Needs access to the reference file used to construct the WCS object.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImageModel`
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


def store_dithered_position(input_model):
    """Store the location of the dithered pointing
    location in the dither metadata

    Parameters
    ----------
    input_model : `jwst.datamodels.ImageModel`
        Data model containing dither offset information
    """
    # V2_ref and v3_ref should be in arcsec
    idltov23 = IdealToV2V3(
        input_model.meta.wcsinfo.v3yangle,
        input_model.meta.wcsinfo.v2_ref, input_model.meta.wcsinfo.v3_ref,
        input_model.meta.wcsinfo.vparity
    )

    dithered_v2, dithered_v3 = idltov23(input_model.meta.dither.x_offset, input_model.meta.dither.y_offset)

    v23toworld = input_model.meta.wcs.get_transform('v2v3', 'world')
    # v23toworld requires a wavelength along with v2, v3, but value does not affect return
    dithered_ra, dithered_dec, _ = v23toworld(dithered_v2, dithered_v3, 0.0)

    input_model.meta.dither.dithered_ra = dithered_ra
    input_model.meta.dither.dithered_dec = dithered_dec
