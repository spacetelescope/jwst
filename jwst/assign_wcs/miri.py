from pathlib import Path
import logging
import numpy as np
from astropy.modeling import bind_bounding_box
from astropy.modeling import models
from astropy import coordinates as coord
from astropy import units as u
from scipy.interpolate import UnivariateSpline
import gwcs.coordinate_frames as cf
from gwcs import selector

from stdatamodels.jwst.datamodels import (
    DistortionModel,
    FilteroffsetModel,
    DistortionMRSModel,
    WavelengthrangeModel,
    RegionsModel,
    SpecwcsModel,
    MiriLRSSpecwcsModel,
)
from stdatamodels.jwst.transforms.models import MIRI_AB2Slice, IdealToV2V3

from . import pointing
from .util import (
    not_implemented_mode,
    subarray_transform,
    velocity_correction,
    transform_bbox_from_shape,
    bounding_box_from_subarray,
)


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging", "lrs", "ifu"]


def create_pipeline(input_model, reference_files):
    """
    Create the WCS pipeline for MIRI modes.

    Parameters
    ----------
    input_model : ImageModel, IFUImageModel, CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)
    if pipeline:
        log.info(f"Created a MIRI {exp_type} pipeline with references {reference_files}")
    return pipeline


def imaging(input_model, reference_files):
    """
    Create the WCS pipeline for MIRI imaging data.

    It includes three coordinate frames - "detector", "v2v3" and "world".

    Parameters
    ----------
    input_model : ImageModel, IFUImageModel, CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'filteroffset' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    # Create the Frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3", axes_order=(0, 1), axes_names=("v2", "v3"), unit=(u.arcsec, u.arcsec)
    )
    v2v3vacorr = cf.Frame2D(
        name="v2v3vacorr", axes_order=(0, 1), axes_names=("v2", "v3"), unit=(u.arcsec, u.arcsec)
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    # Create the transforms
    distortion = imaging_distortion(input_model, reference_files)
    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        distortion = subarray2full | distortion

        # Bind the bounding box to the distortion model using the bounding box ordering
        # used by GWCS. This makes it clear the bounding box is set correctly to GWCS
        bind_bounding_box(distortion, bounding_box_from_subarray(input_model, order="F"), order="F")
    else:
        # TODO: remove setting the bounding box if it is set in the new ref file.
        try:
            distortion.bounding_box  # noqa: B018
        except NotImplementedError:
            shape = input_model.data.shape
            bind_bounding_box(
                distortion, ((3.5, shape[1] - 4.5), (-0.5, shape[0] - 0.5)), order="F"
            )

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    )

    tel2sky = pointing.v23tosky(input_model)

    # Create the pipeline
    pipeline = [(detector, distortion), (v2v3, va_corr), (v2v3vacorr, tel2sky), (world, None)]

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

    Parameters
    ----------
    input_model : ImageModel, IFUImageModel, or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'filteroffset' reference files.

    Returns
    -------
    distortion : `astropy.modeling.Model`
        The transform from "detector" to "v2v3".
    """
    # Read in the distortion.
    with DistortionModel(reference_files["distortion"]) as dist:
        distortion = dist.model

    # Check if the transform in the reference file has a ``bounding_box``.
    # If not set a ``bounding_box`` equal to the size of the image.
    try:
        bbox = distortion.bounding_box
    except NotImplementedError:
        bbox = distortion.bounding_box = None

    # Add an offset for the filter
    obsfilter = input_model.meta.instrument.filter
    with FilteroffsetModel(reference_files["filteroffset"]) as filter_offset:
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

    # Bind the bounding box to the distortion model using the bounding box ordering used by GWCS.
    # This makes it clear the bounding box is set correctly to GWCS
    bind_bounding_box(
        distortion,
        transform_bbox_from_shape(input_model.data.shape, order="F") if bbox is None else bbox,
        order="F",
    )

    return distortion


def lrs(input_model, reference_files):
    """
    Create the WCS pipeline for LRS-FIXEDSLIT and LRS-SLITLESS data.

    It includes three coordinate frames - "detector", "v2v3" and "world".
    "v2v3" and "world" each have (spatial, spatial, spectral) components.

    Parameters
    ----------
    input_model : ImageModel or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'specwcs' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    # Define the various coordinate frames.
    # Original detector frame
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))

    # Intermediate slit frame
    alpha_beta = cf.Frame2D(
        name="alpha_beta_spatial",
        axes_order=(0, 1),
        unit=(u.arcsec, u.arcsec),
        axes_names=("alpha", "beta"),
    )
    spec_local = cf.SpectralFrame(
        name="alpha_beta_spectral", axes_order=(2,), unit=(u.micron,), axes_names=("lambda",)
    )
    miri_focal = cf.CompositeFrame([alpha_beta, spec_local], name="alpha_beta")

    # Spectral component
    spec = cf.SpectralFrame(name="spec", axes_order=(2,), unit=(u.micron,), axes_names=("lambda",))
    # v2v3 spatial component
    v2v3_spatial = cf.Frame2D(
        name="v2v3_spatial", axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=("v2", "v3")
    )
    v2v3vacorr_spatial = cf.Frame2D(
        name="v2v3vacorr_spatial",
        axes_order=(0, 1),
        unit=(u.arcsec, u.arcsec),
        axes_names=("v2", "v3"),
    )

    # v2v3 spatial+spectra
    v2v3 = cf.CompositeFrame([v2v3_spatial, spec], name="v2v3")
    v2v3vacorr = cf.CompositeFrame([v2v3vacorr_spatial, spec], name="v2v3vacorr")

    # 'icrs' frame which is the spatial sky component
    icrs = cf.CelestialFrame(
        name="icrs",
        reference_frame=coord.ICRS(),
        axes_order=(0, 1),
        unit=(u.deg, u.deg),
        axes_names=("RA", "DEC"),
    )
    # Final 'world' composite frame with spatial and spectral components
    world = cf.CompositeFrame(name="world", frames=[icrs, spec])

    # Create the transforms
    dettoabl = lrs_xytoabl(input_model, reference_files)
    abltov2v3l = lrs_abltov2v3l(input_model, reference_files)
    v2v3tosky = pointing.v23tosky(input_model)
    teltosky = v2v3tosky & models.Identity(1)

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    ) & models.Identity(1)

    # Put the transforms together into a single pipeline
    pipeline = [
        (detector, dettoabl),
        (miri_focal, abltov2v3l),
        (v2v3, va_corr),
        (v2v3vacorr, teltosky),
        (world, None),
    ]

    return pipeline


def lrs_xytoabl(input_model, reference_files):
    """
    Build the first transform for the LRS-FIXEDSLIT and LRS-SLITLESS WCS pipeline.

    Parameters
    ----------
    input_model : ImageModel, IFUImageModel, or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'specwcs' reference files.

    Returns
    -------
    dettoabl : `astropy.modeling.Model`
        The transform from subarray (x, y) to (alpha, beta, lambda).
    """
    # subarray to full array transform
    subarray2full = subarray_transform(input_model)

    # full array to v2v3 transform for the ordinary imager
    with DistortionModel(reference_files["distortion"]) as dist:
        distortion = dist.model

    # Combine models to create subarray to v2v3 distortion
    if subarray2full is not None:
        subarray_dist = subarray2full | distortion
    else:
        subarray_dist = distortion

    refmodel = MiriLRSSpecwcsModel(reference_files["specwcs"])
    if input_model.meta.exposure.type.lower() == "mir_lrs-fixedslit":
        zero_point = refmodel.meta.x_ref - 1, refmodel.meta.y_ref - 1
    elif input_model.meta.exposure.type.lower() == "mir_lrs-slitless":
        zero_point = refmodel.meta.x_ref_slitless - 1, refmodel.meta.y_ref_slitless - 1
        # Transform to slitless subarray from full array
        zero_point = subarray2full.inverse(zero_point[0], zero_point[1])

    # Figure out the typical along-slice pixel scale at the center of the slit
    v2_cen, v3_cen = subarray_dist(zero_point[0], zero_point[1])
    v2_off, v3_off = subarray_dist(zero_point[0] + 1, zero_point[1])
    pscale = np.sqrt(np.power(v2_cen - v2_off, 2) + np.power(v3_cen - v3_off, 2))

    # In the lrsdata reference table, X_center,y_center,wavelength describe the location of the
    # centroid trace along the detector in pixels relative to nominal location.
    # x0,y0(ul) x1,y1 (ur) x2,y2(lr) x3,y3(ll) define corners of the box within which the distortion
    # and wavelength calibration was derived
    xcen = refmodel.wavetable.x_center
    ycen = refmodel.wavetable.y_center
    wavetab = refmodel.wavetable.wavelength
    x0 = refmodel.wavetable.x0
    y0 = refmodel.wavetable.y0
    x1 = refmodel.wavetable.x1
    y2 = refmodel.wavetable.y2
    refmodel.close()
    # If in fixed slit mode, define the bounding box using the corner locations provided in
    # the CDP reference file.
    if input_model.meta.exposure.type.lower() == "mir_lrs-fixedslit":
        bb_sub = (
            (np.floor(x0.min() + zero_point[0]) - 0.5, np.ceil(x1.max() + zero_point[0]) + 0.5),
            (np.floor(y2.min() + zero_point[1]) - 0.5, np.ceil(y0.max() + zero_point[1]) + 0.5),
        )

    # If in slitless mode, define the bounding box X locations using the subarray x boundaries
    # and the y locations using the corner locations in the CDP reference file.  Make sure to
    # omit the 4 reference pixels on the left edge of slitless subarray.
    if input_model.meta.exposure.type.lower() == "mir_lrs-slitless":
        bb_sub = (
            (
                input_model.meta.subarray.xstart - 1 + 4 - 0.5,
                input_model.meta.subarray.xsize - 1 + 0.5,
            ),
            (np.floor(y2.min() + zero_point[1]) - 0.5, np.ceil(y0.max() + zero_point[1]) + 0.5),
        )

    # Now deal with the fact that the spectral trace isn't perfectly up and down along detector.
    # This information is contained in the xcenter/ycenter values in the CDP table,
    # but we'll handle it as a simple x shift using a linear fit
    # to this relation provided by the CDP.
    # First convert the values in CDP table to subarray x/y
    xcen_subarray = xcen + zero_point[0]
    ycen_subarray = ycen + zero_point[1]

    # Fit for X shift as a function of Y
    # Spline fit with enforced smoothness
    spl = UnivariateSpline(ycen_subarray[::-1], xcen_subarray[::-1] - zero_point[0], s=0.002)
    # Evaluate the fit at the y reference points
    xshiftref = spl(ycen_subarray)
    # This function will give slit dX as a function of Y subarray pixel value
    dxmodel = models.Tabular1D(
        lookup_table=xshiftref,
        points=ycen_subarray,
        name="xshiftref",
        bounds_error=False,
        fill_value=np.nan,
    )
    if input_model.meta.exposure.type.lower() == "mir_lrs-fixedslit":
        bb_sub = (bb_sub[0], (dxmodel.points[0].min(), dxmodel.points[0].max()))
    # Fit for the wavelength as a function of Y
    # Reverse the vectors so that yinv is increasing (needed for spline fitting function)
    # Spline fit with enforced smoothness
    spl = UnivariateSpline(ycen_subarray[::-1], wavetab[::-1], s=0.002)
    # Evaluate the fit at the y reference points
    wavereference = spl(ycen_subarray)
    # This model will now give the wavelength corresponding to a given Y subarray pixel value
    wavemodel = models.Tabular1D(
        lookup_table=wavereference,
        points=ycen_subarray,
        name="waveref",
        bounds_error=False,
        fill_value=np.nan,
    )
    # Wavelength barycentric correction
    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    else:
        if velosys is not None:
            velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
            wavemodel = wavemodel | velocity_corr
            log.info(
                f"Applied Barycentric velocity correction : {velocity_corr[1].amplitude.value}"
            )

    # What is the effective slit X as a function of subarray x,y?
    xmodel = models.Mapping([0], n_inputs=2) - (models.Mapping([1], n_inputs=2) | dxmodel)
    # What is the effective Y as a function of subarray x,y?
    ymodel = models.Mapping([1], n_inputs=2)
    # What is the effective XY as a function of subarray x,y?
    xymodel = models.Mapping((0, 1, 0, 1)) | xmodel & ymodel
    # What is the alpha as a function of slit XY?
    alphamodel = (
        models.Mapping([0], n_inputs=2)
        | models.Shift(-zero_point[0])
        | models.Polynomial1D(1, c0=0, c1=pscale)
    )
    # What is the alpha,beta as a function of slit XY? (beta is always zero)
    abmodel = models.Mapping((0, 1, 0)) | alphamodel & models.Const1D(0)

    # Define a shift by the reference point and immediately back again
    # This doesn't do anything effectively,
    # but it stores the reference point for later use in pathloss
    # fmt: off
    reftransform = models.Shift(-zero_point[0]) \
        & models.Shift(-zero_point[1]) \
        | models.Shift(+zero_point[0]) \
        & models.Shift(+zero_point[1])
    # fmt: on
    # Put the transforms together
    xytoab = reftransform | xymodel | abmodel

    # Construct the full distortion model (xsub,ysub -> alpha,beta,wavelength)
    lrs_wav_model = models.Mapping([1], n_inputs=2) | wavemodel
    dettoabl = models.Mapping((0, 1, 0, 1)) | xytoab & lrs_wav_model

    # Construct the inverse distortion model (alpha,beta,wavelength -> xsub,ysub)
    # Go from alpha to slit-X
    slitxmodel = models.Polynomial1D(1, c0=0, c1=1 / pscale) | models.Shift(zero_point[0])
    # Go from lambda to real y
    lam_to_y = wavemodel.inverse
    # Go from slit-x and real y to real-x
    backwards = models.Mapping([0], n_inputs=2) + (models.Mapping([1], n_inputs=2) | dxmodel)

    # Go from alpha,beta,lam to real x
    aa = models.Mapping((0, 2)) | slitxmodel & lam_to_y | backwards
    # Go from alpha,beta,lam to real y
    bb = models.Mapping([2], n_inputs=3) | lam_to_y
    # Go from alpha,beta,lam, to real x,y
    dettoabl.inverse = models.Mapping((0, 1, 2, 0, 1, 2)) | aa & bb

    # Bounding box is the subarray bounding box,
    # because we're assuming subarray coordinates passed in
    bind_bounding_box(dettoabl, bb_sub, order="F")

    return dettoabl


def lrs_abltov2v3l(input_model, reference_files):
    """
    Build the first transform for the LRS-FIXEDSLIT and LRS-SLITLESS WCS pipeline.

    Parameters
    ----------
    input_model : ImageModel, IFUImageModel, or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'specwcs' reference files.

    Returns
    -------
    abl_to_v2v3l : `astropy.modeling.Model`
        The transform from (alpha, beta, lambda) to (v2, v3, lambda).
    """
    # subarray to full array transform
    subarray2full = subarray_transform(input_model)

    # full array to v2v3 transform for the ordinary imager
    with DistortionModel(reference_files["distortion"]) as dist:
        distortion = dist.model

    # Combine models to create subarray to v2v3 distortion
    if subarray2full is not None:
        subarray_dist = subarray2full | distortion
    else:
        subarray_dist = distortion

    refmodel = MiriLRSSpecwcsModel(reference_files["specwcs"])
    if input_model.meta.exposure.type.lower() == "mir_lrs-fixedslit":
        zero_point = refmodel.meta.x_ref - 1, refmodel.meta.y_ref - 1
    elif input_model.meta.exposure.type.lower() == "mir_lrs-slitless":
        zero_point = refmodel.meta.x_ref_slitless - 1, refmodel.meta.y_ref_slitless - 1
        # Transform to slitless subarray from full array
        zero_point = subarray2full.inverse(zero_point[0], zero_point[1])

    refmodel.close()
    # Figure out the typical along-slice pixel scale at the center of the slit
    v2_cen, v3_cen = subarray_dist(zero_point[0], zero_point[1])
    v2_off, v3_off = subarray_dist(zero_point[0] + 1, zero_point[1])
    pscale = np.sqrt(np.power(v2_cen - v2_off, 2) + np.power(v3_cen - v3_off, 2))

    # Go from alpha to slit-X
    slitxmodel = models.Polynomial1D(1, c0=0, c1=1 / pscale) | models.Shift(zero_point[0])
    # Go from beta to slit-Y (row_zero_point plus some offset)
    # Beta should always be zero unless using in a pseudo-ifu mode
    slitymodel = models.Polynomial1D(1, c0=0, c1=1 / pscale) | models.Shift(zero_point[1])
    # Go from alpha-beta to slit xy, and onward to v2v3
    ab_to_v2v3 = slitxmodel & slitymodel | subarray_dist
    # Put it together to pass through wavelength
    abl_to_v2v3l = models.Mapping((0, 1, 2)) | ab_to_v2v3 & models.Identity(1)

    # Define the inverse transform
    # Go from slit X to alpha
    alphamodel = models.Shift(-zero_point[0]) | models.Polynomial1D(1, c0=0, c1=pscale)
    # Go from slit Y to beta
    betamodel = models.Shift(-zero_point[1]) | models.Polynomial1D(1, c0=0, c1=pscale)
    # Go from v2,v3 to slit-x,slit-y
    v2v3_to_xydet = subarray_dist.inverse
    # Go from v2,v3 to alpha, beta
    aa = v2v3_to_xydet | alphamodel & betamodel
    # Go from v2,v3,lambda to alpha,beta,lambda
    abl_to_v2v3l.inverse = models.Mapping((0, 1, 2)) | aa & models.Identity(1)

    return abl_to_v2v3l


def ifu(input_model, reference_files):
    """
    Create the WCS pipeline for MIRI IFU data.

    It has the following coordinate frames: "detector", "alpha_beta", "v2v3", "world".

    Parameters
    ----------
    input_model : ImageModel or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion', 'specwcs', 'regions', and 'wavelengthrange' reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    # Define coordinate frames.
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    alpha_beta = cf.Frame2D(
        name="alpha_beta_spatial",
        axes_order=(0, 1),
        unit=(u.arcsec, u.arcsec),
        axes_names=("alpha", "beta"),
    )
    spec_local = cf.SpectralFrame(
        name="alpha_beta_spectral", axes_order=(2,), unit=(u.micron,), axes_names=("lambda",)
    )
    miri_focal = cf.CompositeFrame([alpha_beta, spec_local], name="alpha_beta")
    v23_spatial = cf.Frame2D(
        name="v2v3_spatial", axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=("v2", "v3")
    )
    v2v3vacorr_spatial = cf.Frame2D(
        name="v2v3vacorr_spatial",
        axes_order=(0, 1),
        unit=(u.arcsec, u.arcsec),
        axes_names=("v2", "v3"),
    )

    spec = cf.SpectralFrame(
        name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("lambda",)
    )
    v2v3 = cf.CompositeFrame([v23_spatial, spec], name="v2v3")
    v2v3vacorr = cf.CompositeFrame([v2v3vacorr_spatial, spec], name="v2v3vacorr")
    icrs = cf.CelestialFrame(
        name="icrs",
        reference_frame=coord.ICRS(),
        axes_order=(0, 1),
        unit=(u.deg, u.deg),
        axes_names=("RA", "DEC"),
    )
    world = cf.CompositeFrame([icrs, spec], name="world")

    # Define the actual transforms
    det2abl = (detector_to_abl(input_model, reference_files)).rename("detector_to_abl")
    abl2v2v3l = (abl_to_v2v3l(input_model, reference_files)).rename("abl_to_v2v3l")

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    ) & models.Identity(1)

    tel2sky = pointing.v23tosky(input_model) & models.Identity(1)

    # Put the transforms together into a single transform
    bind_bounding_box(
        det2abl, transform_bbox_from_shape(input_model.data.shape, order="F"), order="F"
    )
    pipeline = [
        (detector, det2abl),
        (miri_focal, abl2v2v3l),
        (v2v3, va_corr),
        (v2v3vacorr, tel2sky),
        (world, None),
    ]
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

    Parameters
    ----------
    input_model : ImageModel or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion', 'specwcs', 'regions', and 'wavelengthrange' reference files.

    Returns
    -------
    det2alpha_beta : `astropy.modeling.Model`
        The transform from "detector" to "alpha_beta" frame.
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range

    with DistortionMRSModel(reference_files["distortion"]) as dist:
        alpha_model = dist.alpha_model
        beta_model = dist.beta_model
        x_model = dist.x_model
        y_model = dist.y_model
        bzero = dict(zip(dist.bzero.channel_band, dist.bzero.beta_zero, strict=True))
        bdel = dict(zip(dist.bdel.channel_band, dist.bdel.delta_beta, strict=True))
        slices = dist.slices

    with SpecwcsModel(reference_files["specwcs"]) as f:
        lambda_model = f.model

    try:
        velosys = input_model.meta.wcsinfo.velosys
    except AttributeError:
        pass
    else:
        if velosys is not None:
            velocity_corr = velocity_correction(input_model.meta.wcsinfo.velosys)
            lambda_model = [m | velocity_corr for m in lambda_model]
            log.info(
                f"Applied Barycentric velocity correction : {velocity_corr[1].amplitude.value}"
            )

    with RegionsModel(reference_files["regions"]) as f:
        allregions = f.regions.copy()
        # Use the 80% throughput slice mask
        regions = allregions[7, :, :]

    label_mapper = selector.LabelMapperArray(regions)
    transforms = {}

    for i, sl in enumerate(slices):
        forward = models.Mapping([1, 0, 0, 1, 0]) | alpha_model[i] & beta_model[i] & lambda_model[i]
        inv = models.Mapping([2, 0, 2, 0]) | x_model[i] & y_model[i]
        forward.inverse = inv
        transforms[sl] = forward

    with WavelengthrangeModel(reference_files["wavelengthrange"]) as f:
        wr = dict(zip(f.waverange_selector, f.wavelengthrange, strict=True))

    ch_dict = {}
    for c in channel:
        cb = c + band
        mapper = MIRI_AB2Slice(bzero[cb], bdel[cb], c)
        lm = selector.LabelMapper(
            inputs=("alpha", "beta", "lam"),
            mapper=mapper,
            inputs_mapping=models.Mapping((1,), n_inputs=3),
        )
        ch_dict[tuple(wr[cb])] = lm

    alpha_beta_mapper = selector.LabelMapperRange(
        ("alpha", "beta", "lam"), ch_dict, models.Mapping((2,))
    )
    label_mapper.inverse = alpha_beta_mapper

    det2alpha_beta = selector.RegionsSelector(
        ("x", "y"), ("alpha", "beta", "lam"), label_mapper=label_mapper, selector=transforms
    )
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

    Parameters
    ----------
    input_model : ImageModel or CubeModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' and 'wavelengthrange' reference files.

    Returns
    -------
    abl2v2v3l : `astropy.modeling.Model`
        The transform from "alpha_beta" to "v2v3" frame.
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range
    channels = [c + band for c in channel]

    with DistortionMRSModel(reference_files["distortion"]) as dist:
        v23 = dict(zip(dist.abv2v3_model.channel_band, dist.abv2v3_model.model, strict=True))

    with WavelengthrangeModel(reference_files["wavelengthrange"]) as f:
        wr = dict(zip(f.waverange_selector, f.wavelengthrange, strict=True))

    dict_mapper = {}
    sel = {}
    # Since there are two channels in each reference file we need to loop over them
    for c in channels:
        ch = int(c[0])
        dict_mapper[tuple(wr[c])] = models.Mapping((2,), name="mapping_lam") | models.Const1D(
            ch, name="channel #"
        )
        ident1 = models.Identity(1, name="identity_lam")
        ident1._inputs = ("lam",)  # noqa: SLF001
        chan_v23 = v23[c]

        v23chan_backward = chan_v23.inverse
        del chan_v23.inverse
        v23_spatial = chan_v23
        v23_spatial.inverse = v23chan_backward
        # Tack on passing the third wavelength component
        v23c = v23_spatial & ident1
        sel[ch] = v23c

    wave_range_mapper = selector.LabelMapperRange(
        ("alpha", "beta", "lam"),
        dict_mapper,
        inputs_mapping=models.Mapping(
            [
                2,
            ]
        ),
    )
    wave_range_mapper.inverse = wave_range_mapper.copy()
    abl2v2v3l = selector.RegionsSelector(
        ("alpha", "beta", "lam"), ("v2", "v3", "lam"), label_mapper=wave_range_mapper, selector=sel
    )

    return abl2v2v3l


exp_type2transform = {
    "mir_image": imaging,
    "mir_tacq": imaging,
    "mir_lyot": imaging,
    "mir_4qpm": imaging,
    "mir_coroncal": imaging,
    "mir_lrs-fixedslit": lrs,
    "mir_lrs-slitless": lrs,
    "mir_mrs": ifu,
    "mir_flatmrs": not_implemented_mode,
    "mir_flatimage": not_implemented_mode,
    "mir_flat-mrs": not_implemented_mode,
    "mir_flat-image": not_implemented_mode,
    "mir_dark": not_implemented_mode,
    "mir_taconfirm": imaging,
}


def get_wavelength_range(input_model, path=None):
    """
    Return the wavelength range used for computing the WCS.

    Needs access to the reference file used to construct the WCS object.

    Parameters
    ----------
    input_model : ImageModel
        Data model after assign_wcs has been run.
    path : str
        Directory where the reference file is. (optional)

    Returns
    -------
    wave_range : set
        A set of tuples containing the channel and wavelength
        range for each channel used in the WCS.
    """
    fname = Path(input_model.meta.ref_file.wavelengthrange.name.split("/")[-1])
    if path is None and not fname.exists():
        raise OSError(f"Reference file {fname} not found. Please specify a path.")
    else:
        fname = Path(path) / fname
        f = WavelengthrangeModel(fname)

    wave_range = f.tree["wavelengthrange"].copy()
    wave_channels = f.tree["channels"]
    f.close()

    wr = dict(zip(wave_channels, wave_range, strict=True))
    channel = input_model.meta.instrument.channel
    band = input_model.meta.instrument.band

    return {(ch + band, wr[ch + band]) for ch in channel}


def store_dithered_position(input_model):
    """
    Store the location of the dithered pointing location in the dither metadata.

    Parameters
    ----------
    input_model : ImageModel
        Data model containing dither offset information
    """
    # V2_ref and v3_ref should be in arcsec
    idltov23 = IdealToV2V3(
        input_model.meta.wcsinfo.v3yangle,
        input_model.meta.wcsinfo.v2_ref,
        input_model.meta.wcsinfo.v3_ref,
        input_model.meta.wcsinfo.vparity,
    )

    dithered_v2, dithered_v3 = idltov23(
        input_model.meta.dither.x_offset, input_model.meta.dither.y_offset
    )

    v23toworld = input_model.meta.wcs.get_transform("v2v3", "world")
    # v23toworld requires a wavelength along with v2, v3, but value does not affect return
    dithered_ra, dithered_dec, _ = v23toworld(dithered_v2, dithered_v3, 0.0)

    input_model.meta.dither.dithered_ra = dithered_ra
    input_model.meta.dither.dithered_dec = dithered_dec
