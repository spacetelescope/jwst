"""
Utility function for assign_wcs.

"""
import warnings
import logging
import functools
import numpy as np

from astropy.utils.misc import isiterable
from astropy.io import fits
from astropy.modeling import models as astmodels
from astropy.table import QTable
from astropy.constants import c

from gwcs import WCS
from gwcs.wcstools import wcs_from_fiducial, grid_from_bounding_box
from gwcs import utils as gwutils

from . import pointing
from ..lib.catalog_utils import SkyObject
from ..transforms.models import GrismObject
from ..datamodels import WavelengthrangeModel, DataModel, CubeModel, IFUCubeModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["reproject", "wcs_from_footprints", "velocity_correction"]


class MissingMSAFileError(Exception):

    def __init__(self, message):
        super(MissingMSAFileError, self).__init__(message)


class NoDataOnDetectorError(Exception):
    """WCS solution indicates no data on detector

    When WCS solutions are available, the solutions indicate that no data
    will be present, raise this exception.

    Specific example is for NIRSpec and the NRS2 detector. For various
    configurations of the MSA, it is possible that no dispersed spectra will
    appear on NRS2. This is not a failure of calibration, but needs to be
    called out in order for the calling architecture to be aware of this.

    """

    def __init__(self, message=None):
        if message is None:
            message = 'WCS solution indicate that no science is in the data.'
        super(NoDataOnDetectorError, self).__init__(message)


def _domain_to_bounding_box(domain):
    # TODO: remove this when domain is completely removed
    bb = tuple([(item['lower'], item['upper']) for item in domain])
    if len(bb) == 1:
        bb = bb[0]
    return bb


def reproject(wcs1, wcs2, origin=0):
    """
    Given two WCSs return a function which takes pixel coordinates in
    the first WCS and computes their location in the second one.

    It performs the forward transformation of ``wcs1`` followed by the
    inverse of ``wcs2``.

    Parameters
    ----------
    wcs1, wcs2 : `~gwcs.wcs.WCS`
        WCS objects.

    Returns
    -------
    _reproject : func
        Function to compute the transformations.  It takes x, y
        positions in ``wcs1`` and returns x, y positions in ``wcs2``.
    """

    def _reproject(x, y):
        sky = wcs1.forward_transform(x, y)
        return wcs2.backward_transform(*sky)
    return _reproject


def wcs_from_footprints(dmodels, refmodel=None, transform=None, bounding_box=None, domain=None):
    """
    Create a WCS from a list of input data models.

    A fiducial point in the output coordinate frame is created from  the
    footprints of all WCS objects. For a spatial frame this is the center
    of the union of the footprints. For a spectral frame the fiducial is in
    the beginning of the footprint range.
    If ``refmodel`` is None, the first WCS object in the list is considered
    a reference. The output coordinate frame and projection (for celestial frames)
    is taken from ``refmodel``.
    If ``transform`` is not suplied, a compound transform is created using
    CDELTs and PC.
    If ``bounding_box`` is not supplied, the bounding_box of the new WCS is computed
    from bounding_box of all input WCSs.

    Parameters
    ----------
    dmodels : list of `~jwst.datamodels.DataModel`
        A list of data models.
    refmodel : `~jwst.datamodels.DataModel`, optional
        This model's WCS is used as a reference.
        WCS. The output coordinate frame, the projection and a
        scaling and rotation transform is created from it. If not supplied
        the first model in the list is used as ``refmodel``.
    transform : `~astropy.modeling.core.Model`, optional
        A transform, passed to :meth:`~gwcs.wcstools.wcs_from_fiducial`
        If not supplied Scaling | Rotation is computed from ``refmodel``.
    bounding_box : tuple, optional
        Bounding_box of the new WCS.
        If not supplied it is computed from the bounding_box of all inputs.
    """
    if domain is not None:
        warnings.warning("'domain' was deprecated in 0.8 and will be removed from next"
                         "version. Use 'bounding_box' instead.")
        bb = _domain_to_bounding_box(domain)
    else:
        bb = bounding_box
    wcslist = [im.meta.wcs for im in dmodels]
    if not isiterable(wcslist):
        raise ValueError("Expected 'wcslist' to be an iterable of WCS objects.")
    if not all([isinstance(w, WCS) for w in wcslist]):
        raise TypeError("All items in wcslist are to be instances of gwcs.WCS.")
    if refmodel is None:
        refmodel = dmodels[0]
    else:
        if not isinstance(refmodel, DataModel):
            raise TypeError("Expected refmodel to be an instance of DataModel.")

    fiducial = compute_fiducial(wcslist, bb)

    prj = astmodels.Pix2Sky_TAN()

    if transform is None:
        transform = []
        wcsinfo = pointing.wcsinfo_from_model(refmodel)
        sky_axes, spec, other = gwutils.get_axes(wcsinfo)
        rotation = astmodels.AffineTransformation2D(wcsinfo['PC'])
        transform.append(rotation)
        if sky_axes:
            cdelt1, cdelt2 = wcsinfo['CDELT'][sky_axes]
            scale = np.sqrt(np.abs(cdelt1 * cdelt2))
            scales = astmodels.Scale(scale) & astmodels.Scale(scale)
            transform.append(scales)

        if transform:
            transform = functools.reduce(lambda x, y: x | y, transform)

    out_frame = refmodel.meta.wcs.output_frame
    wnew = wcs_from_fiducial(fiducial, coordinate_frame=out_frame,
                             projection=prj, transform=transform)

    footprints = [w.footprint().T for w in wcslist]
    domain_bounds = np.hstack([wnew.backward_transform(*f) for f in footprints])
    for axs in domain_bounds:
        axs -= axs.min()
    bounding_box = []
    for axis in out_frame.axes_order:
        axis_min, axis_max = domain_bounds[axis].min(), domain_bounds[axis].max()
        bounding_box.append((axis_min, axis_max))
    bounding_box = tuple(bounding_box)
    ax1, ax2 = np.array(bounding_box)[sky_axes]
    offset1 = (ax1[1] - ax1[0]) / 2
    offset2 = (ax2[1] - ax2[0]) / 2
    offsets = astmodels.Shift(-offset1) & astmodels.Shift(-offset2)

    wnew.insert_transform('detector', offsets, after=True)
    wnew.bounding_box = bounding_box
    return wnew


def compute_fiducial(wcslist, bounding_box=None, domain=None):
    """
    For a celestial footprint this is the center.
    For a spectral footprint, it is the beginning of the range.

    This function assumes all WCSs have the same output coordinate frame.
    """
    if domain is not None:
        warnings.warning("'domain' was deprecated in 0.8 and will be removed from next"
                         "version. Use 'bounding_box' instead.")
    axes_types = wcslist[0].output_frame.axes_type
    spatial_axes = np.array(axes_types) == 'SPATIAL'
    spectral_axes = np.array(axes_types) == 'SPECTRAL'
    footprints = np.hstack([w.footprint(bounding_box=bounding_box).T for w in wcslist])
    spatial_footprint = footprints[spatial_axes]
    spectral_footprint = footprints[spectral_axes]

    fiducial = np.empty(len(axes_types))
    if (spatial_footprint).any():
        lon, lat = spatial_footprint
        lon, lat = np.deg2rad(lon), np.deg2rad(lat)
        x_mean = np.mean(np.cos(lat) * np.cos(lon))
        y_mean = np.mean(np.cos(lat) * np.sin(lon))
        z_mean = np.mean(np.sin(lat))
        lon_fiducial = np.rad2deg(np.arctan2(y_mean, x_mean)) % 360.0
        lat_fiducial = np.rad2deg(np.arctan2(z_mean, np.sqrt(x_mean ** 2 +
            y_mean ** 2)))
        fiducial[spatial_axes] = lon_fiducial, lat_fiducial
    if (spectral_footprint).any():
        fiducial[spectral_axes] = spectral_footprint.min()
    return fiducial


def is_fits(input):
    """
    Returns
    --------
    isFits: tuple
        An ``(isfits, fitstype)`` tuple.  The values of ``isfits`` and
        ``fitstype`` are specified as:

         - ``isfits``: True|False
         - ``fitstype``: if True, one of 'waiver', 'mef', 'simple'; if False, None

    Notes
    -----
    Input images which do not have a valid FITS filename will automatically
    result in a return of (False, None).

    In the case that the input has a valid FITS filename but runs into some
    error upon opening, this routine will raise that exception for the calling
    routine/user to handle.
    """

    isfits = False
    fitstype = None
    names = ['fits', 'fit', 'FITS', 'FIT']
    #determine if input is a fits file based on extension
    # Only check type of FITS file if filename ends in valid FITS string
    f = None
    fileclose = False
    if isinstance(input, fits.HDUList):
        isfits = True
        f = input
    else:
        isfits = True in [input.endswith(l) for l in names]

    # if input is a fits file determine what kind of fits it is
    # waiver fits len(shape) == 3
    if isfits:
        if not f:
            try:
                f = fits.open(input, mode='readonly')
                fileclose = True
            except Exception:
                if f is not None:
                    f.close()
                raise
        data0 = f[0].data
        if data0 is not None:
            try:
                if isinstance(f[1], fits.TableHDU):
                    fitstype = 'waiver'
            except IndexError:
                fitstype = 'simple'

        else:
            fitstype = 'mef'
        if fileclose:
            f.close()

    return isfits, fitstype


def subarray_transform(input_model):
    """
    Inputs are in full frame coordinates.
    If a subarray observation - shift the inputs.

    """
    xstart = input_model.meta.subarray.xstart
    ystart = input_model.meta.subarray.ystart
    if xstart is None:
        xstart = 1
    if ystart is None:
        ystart = 1
    subarray2full = astmodels.Shift(xstart - 1) & astmodels.Shift(ystart - 1)
    return subarray2full


def not_implemented_mode(input_model, ref):
    """
    Return ``None`` if assign_wcs has not been implemented for a mode.
    """
    exp_type = input_model.meta.exposure.type
    message = "WCS for EXP_TYPE of {0} is not implemented.".format(exp_type)
    log.critical(message)
    return None


def get_object_info(catalog_name=None):
    """Return a list of SkyObjects from the direct image

    the source_catalog step are read into a list of  SkyObjects
    which can be referenced by catalog id.

    Parameters
    ----------
    catalog_name : str
        The name of the photutils catalog

    Returns
    -------
    objects : list[jwst.transforms.models.SkyObject]
        A list of SkyObject tuples

    Notes
    -----

    """
    if catalog_name is None:
        raise TypeError("Expected name of the catalog file")
    objects = []
    catalog = QTable.read(catalog_name, format='ascii.ecsv')

    # validate that the expected columns are there
    # id is just a bad name for a param, but it's used in the catalog
    required_fields = list(SkyObject()._fields)
    if "sid" in required_fields:
        required_fields[required_fields.index("sid")] = "id"

    try:
        if not set(required_fields).issubset(set(catalog.colnames)):
            difference = set(catalog.colnames).difference(required_fields)
            raise KeyError("Missing required columns in source catalog ({}): {}"
                           .format(catalog_name, difference))
    except AttributeError as e:
        print("Problem validating object catalog columns {0:s}: {1}"
              .format(catalog_name, e))

    # The columns are named sky_bbox_ll, sky_bbox_ul, sky_bbox_lr, and sky_bbox_ur, each of
    # which is a SkyCoord (i.e. RA & Dec & frame) at one corner of the minimal bounding box.
    # There will also be a sky_bbox property as a 4-tuple of SkyCoord, but that is not
    # serializable (hence, the four separate columns). This is not yet merged in photutils
    # -- I discovered some bugs with SkyCoord and serialization, some of which have just been
    # fixed in astropy.

    for row in catalog:
        objects.append(SkyObject(sid=row['id'],
                                 xcentroid=row['xcentroid'],
                                 ycentroid=row['ycentroid'],
                                 sky_centroid=row['sky_centroid'],
                                 abmag=row['abmag'],
                                 abmag_error=row['abmag_error'],
                                 sky_bbox_ll=row['sky_bbox_ll'],
                                 sky_bbox_lr=row['sky_bbox_lr'],
                                 sky_bbox_ul=row['sky_bbox_ul'],
                                 sky_bbox_ur=row['sky_bbox_ur'],
                                 )
                       )
    return objects


def create_grism_bbox(input_model, reference_files,
                      mmag_extract=99.0):
    """Create bounding boxes for each object in the catalog

    The sky coordinates in the catalog image are first related
    to the grism image. They need to go through the WCS object
    in order to find the "direct image" pixel location, which is
    also in a detector pixel coordinate frame. This "direct image"
    location can then be sent through the trace polynomials to find
    the spectral location on the grism image for that wavelength and order.


    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model which holds the grism image
    reference_files : dict
        Dictionary of reference files
    mmag_extract : float
        The faintest magnitude to extract from the catalog

    Returns
    -------
    A list of GrismObject(s) for every source in the catalog
    Each grism object contains information about it's
    spectral extent

    Notes
    -----
    The wavelengthrange reference file is used to govern the extent of the
    bounding box for each object. The name of the catalog has been stored
    in the input models meta information under the source_catalog key.

    It's left to the calling routine to cut the bounding boxes at the extent
    of the detector (for example, extract 2d would only extract the on-detector
    portion of the bounding box)

    Bounding box dispersion direction is dependent on the filter and module for
    NIRCAM and changes for GRISMR, but is consistent for GRISMC, see
    https://jwst-docs.stsci.edu/display/JTI/NIRCam+Wide+Field+Slitless+Spectroscopy

    NIRISS only has one detector, but GRISMC disperses along rows and GRISMR disperses
    along columns.

    """
    # figure out the dispersion direction, shouldn't this be in the polynomials already?
    disperse_row_right = True  # disperse to increasing x
    disperse_column = False  # column always disperses to increasing y

    instr_name = input_model.meta.instrument.name
    if instr_name == "NIRCAM":
        module, grism, filter_name = (input_model.meta.instrument.module,
                                      input_model.meta.instrument.pupil,
                                      input_model.meta.instrument.filter)
        if ((module == "B") and (grism == "GRISMR")):
            disperse_row_right = False
        elif (grism == "GRISMC"):
            disperse_column = True
            disperse_row_right = False

    elif instr_name == "NIRISS":
        grism, filter_name = (input_model.meta.instrument.filter,
                              input_model.meta.instrument.pupil)

        if "R" == grism[-1]:
            disperse_column = True
            disperse_row_right = False
    else:
        raise ValueError("Input model is from unexpected instrument")

    # get the array extent to exclude boxes not contained on the detector
    xsize = input_model.meta.subarray.xsize
    ysize = input_model.meta.subarray.ysize

    skyobject_list = get_object_info(input_model.meta.source_catalog.filename)

    # get the imaging transform to record the center of the object in the image
    # here, image is in the imaging reference frame, before going through the
    # dispersion coefficients
    sky_to_detector = input_model.meta.wcs.get_transform('world', 'detector')
    sky_to_grism = input_model.meta.wcs.get_transform('world', 'grism_detector')

    # Get the disperser parameters which have the wave limits
    with WavelengthrangeModel(reference_files['wavelengthrange']) as f:
        waverange_selector = f.waverange_selector
        orders = f.order

    # All objects in the catalog will use the same filter for translation
    # that filter is the one that was used in front of the grism
    fselect = waverange_selector.index(filter_name)

    grism_objects = []  # the return list of GrismObjects
    for obj in skyobject_list:
        if obj.abmag < mmag_extract:
            # could add logic to ignore object if too far off image,

            # save the image frame center of the object
            # takes in ra, dec, wavelength, order but wave and order
            # don't get used until the detector->grism_detector transform
            xcenter, ycenter, _, _ = sky_to_detector(obj.sky_centroid.icrs.ra.value,
                                                     obj.sky_centroid.icrs.dec.value,
                                                     1, 1)

            order_bounding = {}
            waverange = {}
            for oidx, order in enumerate(orders):
                # The orders of the bounding box in the non-dispersed image
                # drive the extraction extent. The location of the min and
                # max wavelengths for each order are used to get the location
                # of the +/- sides of the bounding box in the grism image
                lmin, lmax = waverange[oidx][fselect]

                # we need to be specific with dispersion direction here?
                # I think this should be taken care of in the trace polys
                wave_min = lmax
                wave_max = lmin
                if (disperse_row_right or disperse_column):
                    wave_min = lmin
                    wave_max = lmax


                xmin, ymin, _, _, _ = sky_to_grism(obj.sky_bbox_ll.ra.value, obj.sky_bbox_ll.dec.value, lmin, order)
                xmax, ymax, _, _, _ = sky_to_grism(obj.sky_bbox_ur.ra.value, obj.sky_bbox_ur.dec.value, lmax, order)

                # convert to integer pixels, making use of python3 round to integer, 2.7 rounds to float
                # if disperse_column:
                #     cdisp = abs(round(bxmax)-round(bxmin)) // 2
                #     xmin, ymin, xmax, ymax = map(round,[xmin, ymin-cdisp, xmax, ymax+cdisp])
                # else:
                #     cdisp = abs(round(bymax)-round(bymin)) // 2
                #     xmin, ymin, xmax, ymax = map(round,[xmin-cdisp, ymin, xmax+cdisp, ymax])


                # don't add objects and orders which are entirely off the detector
                # this could also live in extract_2d
                # partial_order marks partial off-detector objects which are near enough to cause
                # spectra to be observed on the detector. This is usefull because the catalog often is
                # created from a resampled direct image that is bigger than the detector FOV for a single
                # grism exposure.
                exclude = False
                partial_order = False

                if ((ymin < 0) or (ymax > ysize)):
                    partial_order = True
                if ((ymin < 0) and (ymax < 0)):
                    exclude = True
                if (ymin > ysize):
                    exclude = True

                if ((xmin < 0) or (xmax > xsize)):
                    partial_order = True
                if ((xmin < 0) and (xmax < 0)):
                    exclude = True
                if (xmin > xsize):
                    exclude = True

                if partial_order:
                    log.info("Partial order on detector for obj: {} order: {}".format(obj.sid, order))
                if exclude:
                    log.info("Excluding off-image object: {}, order {}".format(obj.sid, order))
                else:
                    order_bounding[order] = ((round(ymin), round(ymax)), (round(xmin), round(xmax)))
                    waverange[order] = ((lmin, lmax))
            # add lmin and lmax used for the orders here?
            # input_model.meta.wcsinfo.waverange_start keys covers the
            # full range of all the orders

            if len(order_bounding) > 0:
                grism_objects.append(GrismObject(sid=obj.sid,
                                                 order_bounding=order_bounding,
                                                 sky_centroid=obj.sky_centroid,
                                                 partial_order=partial_order,
                                                 waverange=waverange,
                                                 sky_bbox_ll=obj.sky_bbox_ll,
                                                 sky_bbox_lr=obj.sky_bbox_lr,
                                                 sky_bbox_ul=obj.sky_bbox_ul,
                                                 sky_bbox_ur=obj.sky_bbox_ur,
                                                 xcentroid=xcenter,
                                                 ycentroid=ycenter))

    return grism_objects


def get_num_msa_open_shutters(shutter_state):
    """
    Return the number of open shutters in a slitlet.

    Parameters
    ----------
    shutter_state : str
        ``Slit.shutter_state`` attribute - a combination of
        ``1`` - open shutter, ``0`` - closed shutter, ``x`` - main shutter.
    """
    num = shutter_state.count('1')
    if 'x' in shutter_state:
        num += 1
    return num


def bounding_box_from_shape(model):
    """Create a bounding box from the shape of the data.

    Note: The bounding box of a ``CubeModel`` is the bounding_box of one
    of the stacked images.
    """
    if isinstance(model, (CubeModel, IFUCubeModel)):
        shape = model.data[0].shape
    else:
        shape = model.data.shape

    bbox = ((-0.5, shape[1] - 0.5),
            (-0.5, shape[0] - 0.5))
    return bbox


def update_s_region_imaging(model):
    """
    Update the ``S_REGION`` keyword using ``WCS.footprint``.
    """

    bbox = model.meta.wcs.bounding_box

    if bbox is None:
        bbox = bounding_box_from_shape(model)

    # footprint is an array of shape (2, 4) as we
    # are interested only in the footprint on the sky
    footprint = model.meta.wcs.footprint(bbox, center=True, axis_type="spatial").T
    # take only imaging footprint
    footprint = footprint[:2, :]

    # Make sure RA values are all positive
    negative_ind = footprint[0] < 0
    if negative_ind.any():
        footprint[0][negative_ind] = 360 + footprint[0][negative_ind]

    footprint = footprint.T
    update_s_region_keyword(model, footprint)


def update_s_region_spectral(model):
    swcs = model.meta.wcs

    bbox = swcs.bounding_box
    if bbox is None:
        bbox = bounding_box_from_shape(model)

    x, y = grid_from_bounding_box(bbox)
    ra, dec, lam = swcs(x, y)
    footprint = np.array([[np.nanmin(ra), np.nanmin(dec)],
                 [np.nanmax(ra), np.nanmin(dec)],
                 [np.nanmax(ra), np.nanmax(dec)],
                 [np.nanmin(ra), np.nanmax(dec)]])
    update_s_region_keyword(model, footprint)


def update_s_region_keyword(model, footprint):
    """ Update the S_REGION keyword.
    """
    s_region = (
        "POLYGON ICRS "
        " {0:.9f} {1:.9f}"
        " {2:.9f} {3:.9f}"
        " {4:.9f} {5:.9f}"
        " {6:.9f} {7:.9f}".format(*footprint.flatten()))
    if "nan" in s_region:
        # do not update s_region if there are NaNs.
        log.info("There are NaNs in s_region, S_REGION not updated.")
    else:
        model.meta.wcsinfo.s_region = s_region
        log.info("Update S_REGION to {}".format(model.meta.wcsinfo.s_region))


def velocity_correction(velosys):
    """
    Compute wavelength correction to Barycentric reference frame.

    Parameters
    ----------
    velosys : float
        Radial velocity wrt Barycenter [m / s].
    """
    correction = (1 / (1 + velosys / c.value))
    model =  astmodels.Identity(1) * astmodels.Const1D(correction, name="velocity_correction")
    model.inverse = astmodels.Identity(1) / astmodels.Const1D(correction, name="inv_vel_correciton")

    return model
