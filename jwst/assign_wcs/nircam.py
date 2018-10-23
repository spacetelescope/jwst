import logging

from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import Identity, Const1D, Mapping
import gwcs.coordinate_frames as cf

from . import pointing
from .util import not_implemented_mode, subarray_transform
from ..datamodels import (ImageModel, NIRCAMGrismModel, DistortionModel,
                          CubeModel)
from ..transforms.models import (NIRCAMForwardRowGrismDispersion,
                                 NIRCAMForwardColumnGrismDispersion,
                                 NIRCAMBackwardGrismDispersion)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging", "tsgrism", "wfss"]


def create_pipeline(input_model, reference_files):
    """
    Create the WCS pipeline based on EXP_TYPE.

    Parameters
    ----------
    input_model : `~jwst.datamodel.DataModel`
        Input datamodel for processing
    reference_files : dict {reftype: reference file name}
        The dictionary of reference file names and their associated files.

    Returns
    -------
    pipeline : list
        The pipeline list that is returned is suitable for
        input into  gwcs.wcs.WCS to create a GWCS object.
    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def imaging(input_model, reference_files):
    """
    The NIRCAM imaging WCS pipeline.

    Parameters
    ----------
    input_model : `~jwst.datamodel.DataModel`
        Input datamodel for processing
    reference_files : dict
        The dictionary of reference file names and their associated files
        {reftype: reference file name}.

    Returns
    -------
    pipeline : list
        The pipeline list that is returned is suitable for
        input into  gwcs.wcs.WCS to create a GWCS object.

    Notes
    -----
    It includes three coordinate frames - "detector", "v2v3", and "world",
    and uses the "distortion" reference file.
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')

    subarray2full = subarray_transform(input_model)
    imdistortion = imaging_distortion(input_model, reference_files)
    distortion = subarray2full | imdistortion
    distortion.bounding_box = imdistortion.bounding_box
    del imdistortion.bounding_box
    tel2sky = pointing.v23tosky(input_model)
    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the "detector" to "v2v3" transform for imaging mode.


    Parameters
    ----------
    input_model : `~jwst.datamodel.DataModel`
        Input datamodel for processing
    reference_files : dict
        The dictionary of reference file names and their associated files.

    Returns
    -------
    The transform model

    """
    dist = DistortionModel(reference_files['distortion'])
    transform = dist.model

    try:
        bb = transform.bounding_box
    except NotImplementedError:
        shape = input_model.data.shape
        # Note: Since bounding_box is attached to the model here
        # it's in reverse order.
        """
        A CubeModel is always treated as a stack (in dimension 1)
        of 2D images, as opposed to actual 3D data. In this case
        the bounding box is set to the 2nd and 3rd dimension.
        """
        if isinstance(input_model, CubeModel):
            bb = ((-0.5, shape[1] - 0.5),
                  (-0.5, shape[2] - 0.5))
        elif isinstance(input_model, ImageModel):
            bb = ((-0.5, shape[0] - 0.5),
                  (-0.5, shape[1] - 0.5))
        else:
            raise TypeError("Input is not an ImageModel or CubeModel")

        transform.bounding_box = bb
    dist.close()
    return transform


def tsgrism(input_model, reference_files):
    """Create WCS pipeline for a NIRCAM Time Series Grism observation.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImagingModel`
        The input datamodel, derived from datamodels
    reference_files : dict
        Dictionary of reference file names {reftype: reference file name}.

    Returns
    -------
    pipeline : list
        The pipeline list that is returned is suitable for
        input into  gwcs.wcs.WCS to create a GWCS object.

    Notes
    -----
    The TSGRISM mode should function effectively like the grism mode
    except that subarrays will be allowed. Since the transform models
    depend on the original full frame coordinates of the observation,
    the regular grism transforms will need to be shifted to the full
    frame coordinates around the trace transform.

    TSGRISM is only slated to work with GRISMR and Mod A
    """

    # The input is the grism image
    if not isinstance(input_model, CubeModel):
        raise TypeError('The input data model must be a CubeModel.')

    # make sure this is a grism image
    if "NRC_TSGRISM" != input_model.meta.exposure.type:
        raise ValueError('The input exposure is not a NIRCAM time series grism')

    if input_model.meta.instrument.module != "A":
        raise ValueError('NRC_TSGRISM mode only supports module A')

    if input_model.meta.instrument.pupil != "GRISMR":
        raise ValueError('NRC_TSGRIM mode only supports GRISMR')

    gdetector = cf.Frame2D(name='grism_detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    detector = cf.Frame2D(name='full_detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.deg, u.deg))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')

    # translate the x,y detector-in to x,y detector out coordinates
    # Get the disperser parameters which are defined as a model for each
    # spectral order
    with NIRCAMGrismModel(reference_files['specwcs']) as f:
        displ = f.displ
        dispx = f.dispx
        dispy = f.dispy
        invdispx = f.invdispx
        invdispl = f.invdispl
        orders = f.orders

    # now create the appropriate model for the grismr
    det2det = NIRCAMForwardRowGrismDispersion(orders,
                                              lmodels=displ,
                                              xmodels=invdispx,
                                              ymodels=dispy)

    det2det.inverse = NIRCAMBackwardGrismDispersion(orders,
                                                    lmodels=invdispl,
                                                    xmodels=dispx,
                                                    ymodels=dispy)

    # input into the forward transform is x,y,x0,y0,order
    # where x,y is the pixel location in the grism image
    # and x0,y0 is the source location in the "direct" image
    # For this mode, the source is always at crpix1,crpis2
    # discussion with nadia that wcsinfo might not be available
    # here but crpix info could be in wcs.source_location or similar
    # TSGRISM mode places the sources at crpix, and all subarrays
    # begin at 0,0, so no need to translate the crpix to full frame
    # because they already are in full frame coordinates.
    xc, yc = (input_model.meta.wcsinfo.crpix1, input_model.meta.wcsinfo.crpix2)
    xcenter = Const1D(xc)
    xcenter.inverse = Const1D(xc)
    ycenter = Const1D(yc)
    ycenter.inverse = Const1D(yc)

    setra = Const1D(input_model.meta.wcsinfo.crval1)
    setra.inverse = Const1D(input_model.meta.wcsinfo.crval1)
    setdec = Const1D(input_model.meta.wcsinfo.crval2)
    setdec.inverse = Const1D(input_model.meta.wcsinfo.crval2)

    # x, y, order in goes to transform to full array location and order
    # get the shift to full frame coordinates
    sub2full = subarray_transform(input_model) & Identity(1)
    sub2direct = (sub2full | Mapping((0, 1, 0, 1, 2)) |
                  (Identity(2) & xcenter & ycenter & Identity(1)) |
                  det2det)

    # take us from full frame detector to v2v3
    distortion = imaging_distortion(input_model, reference_files) & Identity(2)

    # v2v3 to the sky
    # remap the tel2sky inverse as well since we can feed it the values of
    # crval1, crval2 which correspond to crpix1, crpix2. This leaves
    # us with a calling structure:
    #  (x, y, order) <-> (wavelength, order)
    tel2sky = pointing.v23tosky(input_model) & Identity(2)
    t2skyinverse = tel2sky.inverse
    newinverse = Mapping((0, 1, 0, 1)) | setra & setdec & Identity(2) | t2skyinverse
    tel2sky.inverse = newinverse

    pipeline = [(gdetector, sub2direct),
                (detector, distortion),
                (v2v3, tel2sky),
                (world, None)]

    return pipeline


def wfss(input_model, reference_files):
    """
    Create the WCS pipeline for a NIRCAM grism observation.

    Parameters
    ----------
    input_model: `~jwst.datamodels.ImagingModel`
        The input datamodel, derived from datamodels
    reference_files: dict
        Dictionary {reftype: reference file name}.

    Returns
    -------
    pipeline : list
        The pipeline list that is returned is suitable for
        input into  gwcs.wcs.WCS to create a GWCS object.

    Notes
    -----
    The tree in the grism reference file has a section for each order.
    Not sure if there will be a separate passband reference file needed for
    the wavelength scaling or wedge offsets.

    The direct image the catalog has been created from was processed through
    resample, but the dispersed images have not been resampled. This is OK if
    the trace and dispersion solutions are defined with respect to the
    distortion-corrected image. The catalog from the combined direct image
    has object locations in in detector space and the RA DEC of the object on
    the sky.

    The WCS information for the grism image  plus the observed filter will be
    used to translate these to pixel locations for each of the objects.
    The grism images will then use their grism trace information to translate
    to detector space. The translation is assumed to be one-to-one for purposes
    of identifying the center of the object trace.

    The extent of the trace for each object can then be calculated based on
    the grism in use (row or column). Where the left/bottom of the trace starts
    at t = 0 and the right/top of the trace ends at t = 1, as long as they
    have been defined as such by the team.

    The extraction box is calculated to be the minimum bounding box of the
    object extent in the segmentation map associated with the direct image.
    The values of the min and max corners, taken from the computed minimum
    bounding box are saved in the photometry catalog in units of RA, DEC
    so they can be translated to pixels by the dispersed image's imaging-wcs.
    """

    # The input is the grism image
    if not isinstance(input_model, ImageModel):
        raise TypeError('The input data model must be an ImageModel.')

    # make sure this is a grism image
    if "NRC_WFSS" not in input_model.meta.exposure.type:
            raise ValueError('The input exposure is not a NIRCAM grism')

    # Create the empty detector as a 2D coordinate frame in pixel units
    gdetector = cf.Frame2D(name='grism_detector',
                           axes_order=(0, 1),
                           unit=(u.pix, u.pix))

    # translate the x,y detector-in to x,y detector out coordinates
    # Get the disperser parameters which are defined as a model for each
    # spectral order
    with NIRCAMGrismModel(reference_files['specwcs']) as f:
        displ = f.displ
        dispx = f.dispx
        dispy = f.dispy
        invdispx = f.invdispx
        invdispy = f.invdispy
        invdispl = f.invdispl
        orders = f.orders

    # now create the appropriate model for the grism[R/C]
    if "GRISMR" in input_model.meta.instrument.pupil:
        det2det = NIRCAMForwardRowGrismDispersion(orders,
                                                  lmodels=displ,
                                                  xmodels=invdispx,
                                                  ymodels=dispy)

    elif "GRISMC" in input_model.meta.instrument.pupil:
        det2det = NIRCAMForwardColumnGrismDispersion(orders,
                                                     lmodels=displ,
                                                     xmodels=dispx,
                                                     ymodels=invdispy)

    det2det.inverse = NIRCAMBackwardGrismDispersion(orders,
                                                    lmodels=invdispl,
                                                    xmodels=dispx,
                                                    ymodels=dispy)

    # create the pipeline to construct a WCS object for the whole image
    # which can translate ra,dec to image frame reference pixels
    # it also needs to be part of the grism image wcs pipeline to
    # go from detector to world coordinates. However, the grism image
    # will be effectively translating pixel->world coordinates in a
    # manner that gives you the originating 'imaging' pixels ra and dec,
    # not the ra/dec on the sky from the pointing wcs of the grism image.
    image_pipeline = imaging(input_model, reference_files)

    # input is (x,y,x0,y0,order) -> x0, y0, wave, order
    # x0, y0 is in the image-detector reference frame already
    # and are fed to the wcs to calculate the ra,dec, pix offsets
    # and order are used to calculate the wavelength of the pixel
    grism_pipeline = [(gdetector, det2det)]

    # pass the x0,y0, wave, order, through the pipeline
    imagepipe = []
    world = image_pipeline.pop()
    for cframe, trans in image_pipeline:
        trans = trans & (Identity(2))
        imagepipe.append((cframe, trans))
    imagepipe.append((world))
    grism_pipeline.extend(imagepipe)
    return grism_pipeline


exp_type2transform = {'nrc_image': imaging,
                      'nrc_wfss': wfss,
                      'nrc_tacq': imaging,
                      'nrc_taconfirm': imaging,
                      'nrc_coron': imaging,
                      'nrc_focus': imaging,
                      'nrc_tsimage': imaging,
                      'nrc_tsgrism': tsgrism,
                      'nrc_led': not_implemented_mode,
                      'nrc_dark': not_implemented_mode,
                      'nrc_flat': not_implemented_mode,
                      }
