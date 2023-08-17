"""
Utility function for assign_wcs.

"""
import logging
import functools
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.utils.misc import isiterable
from astropy.io import fits
from astropy.modeling import models as astmodels
from astropy.table import QTable
from astropy.constants import c
from typing import Union, List

from gwcs import WCS
from gwcs.wcstools import wcs_from_fiducial, grid_from_bounding_box
from gwcs import utils as gwutils
from stpipe.exceptions import StpipeExitException

from stdatamodels.jwst.datamodels import JwstDataModel
from stdatamodels.jwst.datamodels import WavelengthrangeModel
from stdatamodels.jwst.transforms.models import GrismObject

from . import pointing
from ..lib.catalog_utils import SkyObject


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


_MAX_SIP_DEGREE = 6


__all__ = ["reproject", "wcs_from_footprints", "velocity_correction",
           "MSAFileError", "NoDataOnDetectorError", "compute_scale",
           "calc_rotation_matrix", "wrap_ra", "update_fits_wcsinfo"]


class MSAFileError(Exception):

    def __init__(self, message):
        super(MSAFileError, self).__init__(message)


class NoDataOnDetectorError(StpipeExitException):
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
        # The first argument instructs stpipe CLI tools to exit with status
        # 64 when this exception is raised.
        super().__init__(64, message)


def _domain_to_bounding_box(domain):
    # TODO: remove this when domain is completely removed
    bb = tuple([(item['lower'], item['upper']) for item in domain])
    if len(bb) == 1:
        bb = bb[0]
    return bb


def reproject(wcs1, wcs2):
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


def compute_scale(wcs: WCS, fiducial: Union[tuple, np.ndarray],
                  disp_axis: int = None, pscale_ratio: float = None) -> float:
    """Compute scaling transform.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        Reference WCS object from which to compute a scaling factor.

    fiducial : tuple
        Input fiducial of (RA, DEC) or (RA, DEC, Wavelength) used in calculating reference points.

    disp_axis : int
        Dispersion axis integer. Assumes the same convention as `wcsinfo.dispersion_direction`

    pscale_ratio : int
        Ratio of input to output pixel scale

    Returns
    -------
    scale : float
        Scaling factor for x and y or cross-dispersion direction.

    """
    spectral = 'SPECTRAL' in wcs.output_frame.axes_type

    if spectral and disp_axis is None:
        raise ValueError('If input WCS is spectral, a disp_axis must be given')

    crpix = np.array(wcs.invert(*fiducial))

    delta = np.zeros_like(crpix)
    spatial_idx = np.where(np.array(wcs.output_frame.axes_type) == 'SPATIAL')[0]
    delta[spatial_idx[0]] = 1

    crpix_with_offsets = np.vstack((crpix, crpix + delta, crpix + np.roll(delta, 1))).T
    crval_with_offsets = wcs(*crpix_with_offsets, with_bounding_box=False)

    coords = SkyCoord(ra=crval_with_offsets[spatial_idx[0]], dec=crval_with_offsets[spatial_idx[1]], unit="deg")
    xscale = np.abs(coords[0].separation(coords[1]).value)
    yscale = np.abs(coords[0].separation(coords[2]).value)

    if pscale_ratio is not None:
        xscale *= pscale_ratio
        yscale *= pscale_ratio

    if spectral:
        # Assuming scale doesn't change with wavelength
        # Assuming disp_axis is consistent with DataModel.meta.wcsinfo.dispersion.direction
        return yscale if disp_axis == 1 else xscale

    return np.sqrt(xscale * yscale)


def calc_rotation_matrix(roll_ref: float, v3i_yang: float, vparity: int = 1) -> List[float]:
    """Calculate the rotation matrix.

    Parameters
    ----------
    roll_ref : float
        Telescope roll angle of V3 North over East at the ref. point in radians

    v3i_yang : float
        The angle between ideal Y-axis and V3 in radians.

    vparity : int
        The x-axis parity, usually taken from the JWST SIAF parameter VIdlParity.
        Value should be "1" or "-1".

    Returns
    -------
    matrix: [pc1_1, pc1_2, pc2_1, pc2_2]
        The rotation matrix

    Notes
    -----
    The rotation is

       ----------------
       | pc1_1  pc2_1 |
       | pc1_2  pc2_2 |
       ----------------

    """
    if vparity not in (1, -1):
        raise ValueError(f'vparity should be 1 or -1. Input was: {vparity}')

    rel_angle = roll_ref - (vparity * v3i_yang)

    pc1_1 = vparity * np.cos(rel_angle)
    pc1_2 = np.sin(rel_angle)
    pc2_1 = vparity * -np.sin(rel_angle)
    pc2_2 = np.cos(rel_angle)

    return [pc1_1, pc1_2, pc2_1, pc2_2]


def wcs_from_footprints(dmodels, refmodel=None, transform=None, bounding_box=None,
                        pscale_ratio=None, pscale=None, rotation=None,
                        shape=None, crpix=None, crval=None):
    """
    Create a WCS from a list of input data models.

    A fiducial point in the output coordinate frame is created from  the
    footprints of all WCS objects. For a spatial frame this is the center
    of the union of the footprints. For a spectral frame the fiducial is in
    the beginning of the footprint range.
    If ``refmodel`` is None, the first WCS object in the list is considered
    a reference. The output coordinate frame and projection (for celestial frames)
    is taken from ``refmodel``.
    If ``transform`` is not supplied, a compound transform is created using
    CDELTs and PC.
    If ``bounding_box`` is not supplied, the bounding_box of the new WCS is computed
    from bounding_box of all input WCSs.

    Parameters
    ----------
    dmodels : list of `~jwst.datamodels.JwstDataModel`
        A list of data models.
    refmodel : `~jwst.datamodels.JwstDataModel`, optional
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
    pscale_ratio : float, None, optional
        Ratio of input to output pixel scale. Ignored when either
        ``transform`` or ``pscale`` are provided.
    pscale : float, None, optional
        Absolute pixel scale in degrees. When provided, overrides
        ``pscale_ratio``. Ignored when ``transform`` is provided.
    rotation : float, None, optional
        Position angle of output image’s Y-axis relative to North.
        A value of 0.0 would orient the final output image to be North up.
        The default of `None` specifies that the images will not be rotated,
        but will instead be resampled in the default orientation for the camera
        with the x and y axes of the resampled image corresponding
        approximately to the detector axes. Ignored when ``transform`` is
        provided.
    shape : tuple of int, None, optional
        Shape of the image (data array) using ``numpy.ndarray`` convention
        (``ny`` first and ``nx`` second). This value will be assigned to
        ``pixel_shape`` and ``array_shape`` properties of the returned
        WCS object.
    crpix : tuple of float, None, optional
        Position of the reference pixel in the image array.  If ``crpix`` is not
        specified, it will be set to the center of the bounding box of the
        returned WCS object.
    crval : tuple of float, None, optional
        Right ascension and declination of the reference pixel. Automatically
        computed if not provided.

    """
    bb = bounding_box
    wcslist = [im.meta.wcs for im in dmodels]

    if not isiterable(wcslist):
        raise ValueError("Expected 'wcslist' to be an iterable of WCS objects.")

    if not all([isinstance(w, WCS) for w in wcslist]):
        raise TypeError("All items in wcslist are to be instances of gwcs.WCS.")

    if refmodel is None:
        refmodel = dmodels[0]
    else:
        if not isinstance(refmodel, JwstDataModel):
            raise TypeError("Expected refmodel to be an instance of DataModel.")

    fiducial = compute_fiducial(wcslist, bb)
    if crval is not None:
        # overwrite spatial axes with user-provided CRVAL:
        i = 0
        for k, axt in enumerate(wcslist[0].output_frame.axes_type):
            if axt == 'SPATIAL':
                fiducial[k] = crval[i]
                i += 1

    ref_fiducial = np.array([refmodel.meta.wcsinfo.ra_ref, refmodel.meta.wcsinfo.dec_ref])

    prj = astmodels.Pix2Sky_TAN()

    if transform is None:
        transform = []
        wcsinfo = pointing.wcsinfo_from_model(refmodel)
        sky_axes, spec, other = gwutils.get_axes(wcsinfo)

        # Need to put the rotation matrix (List[float, float, float, float])
        # returned from calc_rotation_matrix into the correct shape for
        # constructing the transformation
        v3yangle = np.deg2rad(refmodel.meta.wcsinfo.v3yangle)
        vparity = refmodel.meta.wcsinfo.vparity
        if rotation is None:
            roll_ref = np.deg2rad(refmodel.meta.wcsinfo.roll_ref)
        else:
            roll_ref = np.deg2rad(rotation) + (vparity * v3yangle)

        pc = np.reshape(
            calc_rotation_matrix(roll_ref, v3yangle, vparity=vparity),
            (2, 2)
        )

        rotation = astmodels.AffineTransformation2D(pc, name='pc_rotation_matrix')
        transform.append(rotation)

        if sky_axes:
            if not pscale:
                pscale = compute_scale(refmodel.meta.wcs, ref_fiducial,
                                       pscale_ratio=pscale_ratio)
            transform.append(astmodels.Scale(pscale, name='cdelt1') & astmodels.Scale(pscale, name='cdelt2'))

        if transform:
            transform = functools.reduce(lambda x, y: x | y, transform)

    out_frame = refmodel.meta.wcs.output_frame
    input_frame = refmodel.meta.wcs.input_frame
    wnew = wcs_from_fiducial(fiducial, coordinate_frame=out_frame, projection=prj,
                             transform=transform, input_frame=input_frame)

    footprints = [w.footprint().T for w in wcslist]
    domain_bounds = np.hstack([wnew.backward_transform(*f) for f in footprints])
    axis_min_values = np.min(domain_bounds, axis=1)
    domain_bounds = (domain_bounds.T - axis_min_values).T

    output_bounding_box = []
    for axis in out_frame.axes_order:
        axis_min, axis_max = domain_bounds[axis].min(), domain_bounds[axis].max()
        output_bounding_box.append((axis_min, axis_max))

    output_bounding_box = tuple(output_bounding_box)
    if crpix is None:
        offset1, offset2 = wnew.backward_transform(*fiducial)
        offset1 -= axis_min_values[0]
        offset2 -= axis_min_values[1]
    else:
        offset1, offset2 = crpix
    offsets = astmodels.Shift(-offset1, name='crpix1') & astmodels.Shift(-offset2, name='crpix2')

    wnew.insert_transform('detector', offsets, after=True)
    wnew.bounding_box = output_bounding_box

    if shape is None:
        shape = [int(axs[1] - axs[0] + 0.5) for axs in output_bounding_box[::-1]]

    wnew.pixel_shape = shape[::-1]
    wnew.array_shape = shape

    return wnew


def compute_fiducial(wcslist, bounding_box=None):
    """
    For a celestial footprint this is the center.
    For a spectral footprint, it is the beginning of the range.

    This function assumes all WCSs have the same output coordinate frame.
    """

    axes_types = wcslist[0].output_frame.axes_type
    spatial_axes = np.array(axes_types) == 'SPATIAL'
    spectral_axes = np.array(axes_types) == 'SPECTRAL'
    footprints = np.hstack([w.footprint(bounding_box=bounding_box).T for w in wcslist])
    spatial_footprint = footprints[spatial_axes]
    spectral_footprint = footprints[spectral_axes]

    fiducial = np.empty(len(axes_types))
    if spatial_footprint.any():
        lon, lat = spatial_footprint
        lon, lat = np.deg2rad(lon), np.deg2rad(lat)
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)

        x_mid = (np.max(x) + np.min(x)) / 2.
        y_mid = (np.max(y) + np.min(y)) / 2.
        z_mid = (np.max(z) + np.min(z)) / 2.
        lon_fiducial = np.rad2deg(np.arctan2(y_mid, x_mid)) % 360.0
        lat_fiducial = np.rad2deg(np.arctan2(z_mid, np.sqrt(x_mid ** 2 + y_mid ** 2)))
        fiducial[spatial_axes] = lon_fiducial, lat_fiducial
    if spectral_footprint.any():
        fiducial[spectral_axes] = spectral_footprint.min()
    return fiducial


def is_fits(input_img):
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
    # determine if input is a fits file based on extension
    # Only check type of FITS file if filename ends in valid FITS string
    f = None
    fileclose = False
    if isinstance(input_img, fits.HDUList):
        isfits = True
        f = input_img
    else:
        isfits = True in [input_img.endswith(suffix) for suffix in names]

    # if input is a fits file determine what kind of fits it is
    # waiver fits len(shape) == 3
    if isfits:
        if not f:
            try:
                f = fits.open(input_img, mode='readonly')
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
    Return an offset model if the observation uses a subarray.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        Data model.

    Returns
    -------
    subarray2full : `~astropy.modeling.core.Model` or ``None``
        Returns a (combination of ) ``Shift`` models if a subarray is used.
        Returns ``None`` if a full frame observation.
    """
    tr_xstart = astmodels.Identity(1)
    tr_ystart = astmodels.Identity(1)

    # These quantities are 1-based
    xstart = input_model.meta.subarray.xstart
    ystart = input_model.meta.subarray.ystart

    if xstart is not None and xstart != 1:
        tr_xstart = astmodels.Shift(xstart - 1)

    if ystart is not None and ystart != 1:
        tr_ystart = astmodels.Shift(ystart - 1)

    if (isinstance(tr_xstart, astmodels.Identity) and
            isinstance(tr_ystart, astmodels.Identity)):
        # the case of a full frame observation
        return None
    else:
        subarray2full = tr_xstart & tr_ystart
        return subarray2full


def not_implemented_mode(input_model, ref, slit_y_range=None):
    """
    Return ``None`` if assign_wcs has not been implemented for a mode.
    """
    exp_type = input_model.meta.exposure.type
    message = "WCS for EXP_TYPE of {0} is not implemented.".format(exp_type)
    log.critical(message)
    return None


def get_object_info(catalog_name=None):
    """Return a list of SkyObjects from the direct image

    The source_catalog step catalog items are read into a list
    of  SkyObjects which can be referenced by catalog id. Only
    the columns needed by the WFSS code are saved.

    Parameters
    ----------
    catalog_name : str, astropy.table.table.Qtable
        The name of the photutils catalog or its quantities table

    Returns
    -------
    objects : list[jwst.transforms.models.SkyObject]
        A list of SkyObject tuples

    Notes
    -----

    """
    if isinstance(catalog_name, str):
        if len(catalog_name) == 0:
            err_text = "Empty catalog filename"
            log.error(err_text)
            raise ValueError(err_text)
        try:
            catalog = QTable.read(catalog_name, format='ascii.ecsv')
        except FileNotFoundError as e:
            log.error("Could not find catalog file: {0}".format(e))
            raise FileNotFoundError("Could not find catalog: {0}".format(e))
    elif isinstance(catalog_name, QTable):
        catalog = catalog_name
    else:
        err_text = "Need to input string name of catalog or astropy.table.table.QTable instance"
        log.error(err_text)
        raise TypeError(err_text)

    objects = []

    # validate that the expected columns are there
    required_fields = set(SkyObject()._fields)

    try:
        if not set(required_fields).issubset(set(catalog.colnames)):
            difference = set(required_fields).difference(set(catalog.colnames))
            err_text = "Missing required columns in source catalog: {0}".format(difference)
            log.error(err_text)
            raise KeyError(err_text)
    except AttributeError as e:
        err_text = "Problem validating object catalog columns: {0}".format(e)
        log.error(err_text)
        raise AttributeError

    # The columns are named sky_bbox_ll, sky_bbox_ul, sky_bbox_lr,
    # and sky_bbox_ur, each of which is a SkyCoord (i.e. RA & Dec & frame) at
    # one corner of the minimal bounding box. There will also be a sky_bbox
    # property as a 4-tuple of SkyCoord, but that is not serializable
    # (hence, the four separate columns).

    for row in catalog:
        objects.append(SkyObject(label=row['label'],
                                 xcentroid=row['xcentroid'],
                                 ycentroid=row['ycentroid'],
                                 sky_centroid=row['sky_centroid'],
                                 isophotal_abmag=row['isophotal_abmag'],
                                 isophotal_abmag_err=row['isophotal_abmag_err'],
                                 sky_bbox_ll=row['sky_bbox_ll'],
                                 sky_bbox_lr=row['sky_bbox_lr'],
                                 sky_bbox_ul=row['sky_bbox_ul'],
                                 sky_bbox_ur=row['sky_bbox_ur'],
                                 is_extended=row['is_extended']
                                 )
                       )
    return objects


def create_grism_bbox(input_model,
                      reference_files=None,
                      mmag_extract=None,
                      extract_orders=None,
                      wfss_extract_half_height=None,
                      wavelength_range=None,
                      nbright=None):
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
    reference_files : dict, optional
        Dictionary of reference file names.
        If ``None``, ``wavelength_range`` must be supplied to specify
        the orders and corresponding wavelength ranges to be used in extraction.
    mmag_extract : float, optional
        The faintest magnitude to extract from the catalog.
    extract_orders : list, optional
        The list of orders to extract, if specified this will
        override the orders listed in the wavelengthrange reference file.
        If ``None``, the default one in the wavelengthrange reference file is used.
    wfss_extract_half_height : int, optional
        Cross-dispersion extraction half height in pixels, WFSS mode.
        Overwrites the computed extraction height in ``GrismObject.order_bounding.``
        If ``None``, it's computed from the segmentation map,
        using the min and max wavelength for each of the orders that
        are available.
    wavelength_range : dict, optional
        Pairs of {spectral_order: (wave_min, wave_max)} for each order.
        If ``None``, the default one in the wavelengthrange reference file is used.
    nbright : int, optional
        The number of brightest objects to extract from the catalog.

    Returns
    -------
    grism_objects : list
        A list of GrismObject(s) for every source in the catalog.
        Each grism object contains information about its
        spectral extent.

    Notes
    -----
    The wavelengthrange reference file is used to govern
    the extent of the bounding box for each object. The name of the
    catalog has been stored in the input models meta information under
    the source_catalog key.

    It's left to the calling routine to cut the bounding boxes at the
    extent of the detector (for example, extract 2d would only extract
    the on-detector portion of the bounding box)

    Bounding box dispersion direction is dependent on the filter and
    module for NIRCAM and changes for GRISMR, but is consistent for GRISMC,
    see https://jwst-docs.stsci.edu/display/JTI/NIRCam+Wide+Field+Slitless+Spectroscopy

    NIRISS has one detector.  GRISMC disperses along rows and
    GRISMR disperses along columns.

    If ``wfss_extract_half_height`` is specified it is used to compute the extent in
    the cross-dispersion direction, which becomes ``2 * wfss_extract_half_height + 1``.
    ``wfss_extract_half_height`` can only be applied to point source objects.

    """
    instr_name = input_model.meta.instrument.name
    if instr_name == "NIRCAM":
        filter_name = input_model.meta.instrument.filter
    elif instr_name == "NIRISS":
        filter_name = input_model.meta.instrument.pupil
    else:
        raise ValueError("create_grism_object works with NIRCAM and NIRISS WFSS exposures only.")

    if reference_files is None:
        # Get the list of extract_orders and lmin, lmax from wavelength_range.
        if wavelength_range is None:
            message = "If reference files are not supplied, ``wavelength_range`` must be provided."
            raise TypeError(message)
    else:
        # Get the list of extract_orders and lmin, lmax from the ``wavelengthrange`` reference file.
        with WavelengthrangeModel(reference_files['wavelengthrange']) as f:
            if 'WFSS' not in f.meta.exposure.type:
                err_text = "Wavelengthrange reference file not for WFSS"
                log.error(err_text)
                raise ValueError(err_text)
            ref_extract_orders = f.extract_orders
            if extract_orders is None:
                # ref_extract_orders = extract_orders
                extract_orders = [x[1] for x in ref_extract_orders if x[0] == filter_name].pop()

            wavelength_range = f.get_wfss_wavelength_range(filter_name, extract_orders)

    if mmag_extract is None:
        mmag_extract = 999.  # extract all objects, regardless of magnitude
    else:
        log.info("Extracting objects < abmag = {0}".format(mmag_extract))
    if not isinstance(mmag_extract, (int, float)):
        raise TypeError(f"Expected mmag_extract to be a number, got {mmag_extract}")

    # extract the catalog objects
    if input_model.meta.source_catalog is None:
        err_text = "No source catalog listed in datamodel."
        log.error(err_text)
        raise ValueError(err_text)

    log.info(f"Getting objects from {input_model.meta.source_catalog}")

    return _create_grism_bbox(input_model, mmag_extract, wfss_extract_half_height, wavelength_range,
                              nbright)


def _create_grism_bbox(input_model, mmag_extract=None, wfss_extract_half_height=None,
                       wavelength_range=None, nbright=None):

    log.debug(f'Extracting with wavelength_range {wavelength_range}')

    # this contains the pure information from the catalog with no translations
    skyobject_list = get_object_info(input_model.meta.source_catalog)
    # get the imaging transform to record the center of the object in the image
    # here, image is in the imaging reference frame, before going through the
    # dispersion coefficients

    sky_to_detector = input_model.meta.wcs.get_transform('world', 'detector')
    sky_to_grism = input_model.meta.wcs.backward_transform

    grism_objects = []  # the return list of GrismObjects
    for obj in skyobject_list:
        if obj.isophotal_abmag is not None:
            if obj.isophotal_abmag < mmag_extract:
                # could add logic to ignore object if too far off image,

                # save the image frame center of the object
                # takes in ra, dec, wavelength, order but wave and order
                # don't get used until the detector->grism_detector transform
                xcenter, ycenter, _, _ = sky_to_detector(obj.sky_centroid.icrs.ra.value,
                                                         obj.sky_centroid.icrs.dec.value,
                                                         1, 1)

                order_bounding = {}
                waverange = {}
                partial_order = {}
                for order in wavelength_range:
                    # range_select = [(x[2], x[3]) for x in wavelengthrange if (x[0] == order and x[1] == filter_name)]
                    # The orders of the bounding box in the non-dispersed image
                    # drive the extraction extent. The location of the min and
                    # max wavelengths for each order are used to get the
                    # location of the +/- sides of the bounding box in the
                    # grism image
                    lmin, lmax = wavelength_range[order]
                    ra = np.array([obj.sky_bbox_ll.ra.value, obj.sky_bbox_lr.ra.value,
                                   obj.sky_bbox_ul.ra.value, obj.sky_bbox_ur.ra.value])
                    dec = np.array([obj.sky_bbox_ll.dec.value, obj.sky_bbox_lr.dec.value,
                                    obj.sky_bbox_ul.dec.value, obj.sky_bbox_ur.dec.value])
                    x1, y1, _, _, _ = sky_to_grism(ra, dec, [lmin] * 4, [order] * 4)
                    x2, y2, _, _, _ = sky_to_grism(ra, dec, [lmax] * 4, [order] * 4)

                    xstack = np.hstack([x1, x2])
                    ystack = np.hstack([y1, y2])

                    # Subarrays are only allowed in nircam tsgrism mode. The polynomial transforms
                    # only work with the full frame coordinates. The code here is called during extract_2d,
                    # and is creating bounding boxes which should be in the full frame coordinates, it just
                    # uses the input catalog and the magnitude to limit the objects that need bounding boxes.

                    # Tsgrism is always supposed to have the source object at the same pixel, and that is
                    # hardcoded into the transforms. At least a while ago, the 2d extraction for tsgrism mode
                    # didn't call this bounding box code. So I think it's safe to leave the subarray
                    # subtraction out, i.e. do not subtract x/ystart.

                    xmin = np.nanmin(xstack)
                    xmax = np.nanmax(xstack)
                    ymin = np.nanmin(ystack)
                    ymax = np.nanmax(ystack)

                    if wfss_extract_half_height is not None and not obj.is_extended:
                        if input_model.meta.wcsinfo.dispersion_direction == 2:
                            ra_center, dec_center = obj.sky_centroid.ra.value, obj.sky_centroid.dec.value
                            center, _, _, _, _ = sky_to_grism(ra_center, dec_center, (lmin + lmax) / 2, order)
                            xmin = center - wfss_extract_half_height
                            xmax = center + wfss_extract_half_height
                        elif input_model.meta.wcsinfo.dispersion_direction == 1:
                            ra_center, dec_center = obj.sky_centroid.ra.value, obj.sky_centroid.dec.value
                            _, center, _, _, _ = sky_to_grism(ra_center, dec_center, (lmin + lmax) / 2, order)
                            ymin = center - wfss_extract_half_height
                            ymax = center + wfss_extract_half_height
                        else:
                            raise ValueError("Cannot determine dispersion direction.")

                    # Convert floating-point corner values to whole pixel indexes
                    xmin = gwutils._toindex(xmin)
                    xmax = gwutils._toindex(xmax)
                    ymin = gwutils._toindex(ymin)
                    ymax = gwutils._toindex(ymax)

                    # Don't add objects and orders that are entirely off the detector.
                    # "partial_order" marks objects that are near enough to the detector
                    # edge to have some spectrum on the detector.
                    # This is useful because the catalog often is created from a resampled direct
                    # image that is bigger than the detector FOV for a single grism exposure.
                    exclude = False
                    ispartial = False

                    # Here we check to ensure that the extraction region `pts`
                    # has at least two pixels of width in the dispersion
                    # direction, and one in the cross-dispersed direction when
                    # placed into the subarray extent.
                    pts = np.array([[ymin, xmin], [ymax, xmax]])
                    subarr_extent = np.array([[0, 0],
                                             [input_model.meta.subarray.ysize - 1,
                                              input_model.meta.subarray.xsize - 1]])

                    if input_model.meta.wcsinfo.dispersion_direction == 1:
                        # X-axis is dispersion direction
                        disp_col = 1
                        xdisp_col = 0
                    else:
                        # Y-axis is dispersion direction
                        disp_col = 0
                        xdisp_col = 1

                    dispaxis_check = (pts[1, disp_col] - subarr_extent[0, disp_col] > 0) and \
                                     (subarr_extent[1, disp_col] - pts[0, disp_col] > 0)
                    xdispaxis_check = (pts[1, xdisp_col] - subarr_extent[0, xdisp_col] >= 0) and \
                                      (subarr_extent[1, xdisp_col] - pts[0, xdisp_col] >= 0)

                    contained = dispaxis_check and xdispaxis_check

                    inidx = np.all(np.logical_and(subarr_extent[0] <= pts, pts <= subarr_extent[1]), axis=1)

                    if not contained:
                        exclude = True
                        log.info("Excluding off-image object: {}, order {}".format(obj.label, order))
                    elif contained >= 1:
                        outbox = pts[np.logical_not(inidx)]
                        if len(outbox) > 0:
                            ispartial = True
                            log.info("Partial order on detector for obj: {} order: {}".format(obj.label, order))

                    if not exclude:
                        order_bounding[order] = ((ymin, ymax), (xmin, xmax))
                        waverange[order] = ((lmin, lmax))
                        partial_order[order] = ispartial

                if len(order_bounding) > 0:
                    grism_objects.append(GrismObject(sid=obj.label,
                                                     order_bounding=order_bounding,
                                                     sky_centroid=obj.sky_centroid,
                                                     partial_order=partial_order,
                                                     waverange=waverange,
                                                     sky_bbox_ll=obj.sky_bbox_ll,
                                                     sky_bbox_lr=obj.sky_bbox_lr,
                                                     sky_bbox_ul=obj.sky_bbox_ul,
                                                     sky_bbox_ur=obj.sky_bbox_ur,
                                                     xcentroid=xcenter,
                                                     ycentroid=ycenter,
                                                     is_extended=obj.is_extended,
                                                     isophotal_abmag=obj.isophotal_abmag))

    # At this point we have a list of grism objects limited to
    # isophotal_abmag < mmag_extract. We now need to further restrict
    # the list to the N brightest objects, as given by nbright.
    if nbright is None:
        # Include all objects, regardless of brightness
        final_objects = grism_objects
    else:
        # grism_objects is a list of objects, so it's not easy or practical
        # to sort it directly. So create a list of the isophotal_abmags, which
        # we'll then use to find the N brightest objects.
        indxs = np.argsort([obj.isophotal_abmag for obj in grism_objects])

        # Create a final grism object list containing only the N brightest objects
        final_objects = []
        final_objects = [grism_objects[i] for i in indxs[:nbright]]
        del grism_objects

    log.info(f"Total of {len(final_objects)} grism objects defined")
    if len(final_objects) == 0:
        log.warning("No grism objects saved; check catalog or step params")

    return final_objects


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


def transform_bbox_from_shape(shape):
    """Create a bounding box from the shape of the data.

    This is appropriate to attached to a transform.

    Parameters
    ----------
    shape : tuple
        The shape attribute from a `numpy.ndarray` array

    Returns
    -------
    bbox : tuple
        Bounding box in y, x order.
    """
    bbox = ((-0.5, shape[-2] - 0.5),
            (-0.5, shape[-1] - 0.5))
    return bbox


def wcs_bbox_from_shape(shape):
    """Create a bounding box from the shape of the data.

    This is appropriate to attach to a wcs object
    Parameters
    ----------
    shape : tuple
        The shape attribute from a `numpy.ndarray` array

    Returns
    -------
    bbox : tuple
        Bounding box in x, y order.
    """
    bbox = ((-0.5, shape[-1] - 0.5),
            (-0.5, shape[-2] - 0.5))
    return bbox


def bounding_box_from_subarray(input_model):
    """Create a bounding box from the subarray size.

    Note: The bounding_box assumes full frame coordinates.
    It is set to ((ystart, ystart + xsize), (xstart, xstart + xsize)).
    It is in 0-based coordinates.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        The data model.

    Returns
    -------
    bbox : tuple
        Bounding box in y, x order.
    """
    bb_xstart = -0.5
    bb_xend = -0.5
    bb_ystart = -0.5
    bb_yend = -0.5

    if input_model.meta.subarray.xsize is not None:
        bb_xend = input_model.meta.subarray.xsize - 0.5
    if input_model.meta.subarray.ysize is not None:
        bb_yend = input_model.meta.subarray.ysize - 0.5

    bbox = ((bb_ystart, bb_yend), (bb_xstart, bb_xend))
    return bbox


def update_s_region_imaging(model):
    """
    Update the ``S_REGION`` keyword using ``WCS.footprint``.
    """

    bbox = model.meta.wcs.bounding_box

    if bbox is None:
        bbox = wcs_bbox_from_shape(model.data.shape)
        model.meta.wcs.bounding_box = bbox

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


def compute_footprint_spectral(model):
    """
    Determine spatial footprint for spectral observations using the instrument model.

    Parameters
    ----------
    model : `~jwst.datamodels.IFUImageModel`
        The output of assign_wcs.
    """
    swcs = model.meta.wcs
    bbox = swcs.bounding_box
    if bbox is None:
        bbox = wcs_bbox_from_shape(model.data.shape)

    x, y = grid_from_bounding_box(bbox)
    ra, dec, lam = swcs(x, y)

    # the wrapped ra values are forced to be on one side of ra-border
    # the wrapped ra are used to determine the correct  min and max ra
    ra = wrap_ra(ra)
    min_ra = np.nanmin(ra)
    max_ra = np.nanmax(ra)

    # for the footprint we want the ra values to fall between 0 to 360
    if min_ra < 0:
        min_ra = min_ra + 360.0
    if max_ra >= 360.0:
        max_ra = max_ra - 360.0
    footprint = np.array([[min_ra, np.nanmin(dec)],
                          [max_ra, np.nanmin(dec)],
                          [max_ra, np.nanmax(dec)],
                          [min_ra, np.nanmax(dec)]])
    lam_min = np.nanmin(lam)
    lam_max = np.nanmax(lam)
    return footprint, (lam_min, lam_max)


def update_s_region_spectral(model):
    """ Update the S_REGION keyword.
    """
    footprint, spectral_region = compute_footprint_spectral(model)
    update_s_region_keyword(model, footprint)
    model.meta.wcsinfo.spectral_region = spectral_region


def compute_footprint_nrs_slit(slit):
    """ Compute the footprint of a Nirspec slit using the instrument model.

    Parameters
    ----------
    slit : `~jwst.datamodels.SlitModel`
    """
    slit2world = slit.meta.wcs.get_transform("slit_frame", "world")
    # Define the corners of a virtual slit. The center of the slit is (0, 0).
    virtual_corners_x = [-.5, -.5, .5, .5]
    virtual_corners_y = [slit.slit_ymin, slit.slit_ymax, slit.slit_ymax, slit.slit_ymin]
    # Use a default wavelength or 2 microns as input to the transform.
    input_lam = [2e-6] * 4
    ra, dec, lam = slit2world(virtual_corners_x,
                              virtual_corners_y,
                              input_lam)
    footprint = np.array([ra, dec]).T
    lam_min = np.nanmin(lam)
    lam_max = np.nanmax(lam)
    return footprint, (lam_min, lam_max)


def update_s_region_nrs_slit(slit):
    footprint, spectral_region = compute_footprint_nrs_slit(slit)
    update_s_region_keyword(slit, footprint)
    slit.meta.wcsinfo.spectral_region = spectral_region


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


def compute_footprint_nrs_ifu(dmodel, mod):
    """
    Determine NIRSPEC IFU footprint using the instrument model.

    For efficiency this function uses the transforms directly,
    instead of the WCS object. The common transforms in the WCS
    model chain are referenced and reused; only the slice specific
    transforms are computed.

    If the transforms change this function should be revised.

    Parameters
    ----------
    output_model : `~jwst.datamodels.IFUImageModel`
        The output of assign_wcs.
    mod : module
        The imported ``nirspec`` module.

    Returns
    -------
    footprint : ndarray
        The spatial footprint
    spectral_region : tuple
        The wavelength range for the observation.
    """
    ra_total = []
    dec_total = []
    lam_total = []
    _, wrange = mod.spectral_order_wrange_from_model(dmodel)
    pipe = dmodel.meta.wcs.pipeline

    # Get the GWA to slit_frame transform
    g2s = pipe[2].transform

    # Construct a list of the transforms between coordinate frames.
    # Set a place holder ``Identity`` transform at index 2 and 3.
    # Update them with slice specific transforms.
    transforms = [pipe[0].transform]
    transforms.append(pipe[1].transform[1:])
    transforms.append(astmodels.Identity(1))
    transforms.append(astmodels.Identity(1))
    transforms.extend([step.transform for step in pipe[4:-1]])

    for sl in range(30):
        transforms[2] = g2s.get_model(sl)
        # Create the full transform from ``slit_frame`` to ``detector``.
        # It is used to compute the bounding box.
        m = functools.reduce(lambda x, y: x | y, [tr.inverse for tr in transforms[:3][::-1]])
        bbox = mod.compute_bounding_box(m, wrange)
        # Add the remaining transforms - from ``sli_frame`` to ``world``
        transforms[3] = pipe[3].transform.get_model(sl) & astmodels.Identity(1)
        mforw = functools.reduce(lambda x, y: x | y, transforms)
        x1, y1 = grid_from_bounding_box(bbox)
        ra, dec, lam = mforw(x1, y1)
        ra_total.extend(np.ravel(ra))
        dec_total.extend(np.ravel(dec))
        lam_total.extend(np.ravel(lam))
    # the wrapped ra values are forced to be on one side of ra-border
    # the wrapped ra are used to determine the correct  min and max ra
    ra_total = wrap_ra(ra_total)
    ra_max = np.nanmax(ra_total)
    ra_min = np.nanmin(ra_total)
    # for the footprint we want ra to be between 0 to 360
    if ra_min < 0:
        ra_min = ra_min + 360.0
    if ra_max >= 360.0:
        ra_max = ra_max - 360.0

    dec_max = np.nanmax(dec_total)
    dec_min = np.nanmin(dec_total)
    lam_max = np.nanmax(lam_total)
    lam_min = np.nanmin(lam_total)
    footprint = np.array([ra_min, dec_min, ra_max, dec_min, ra_max, dec_max, ra_min, dec_max])
    return footprint, (lam_min, lam_max)


def update_s_region_nrs_ifu(output_model, mod):
    """
    Update S_REGION for NRS_IFU observations using calculated footprint.

    Parameters
    ----------
    output_model : `~jwst.datamodels.IFUImageModel`
        The output of assign_wcs.
    mod : module
        The imported ``nirspec`` module.
    """
    footprint, spectral_region = compute_footprint_nrs_ifu(output_model, mod)
    update_s_region_keyword(output_model, footprint)
    output_model.meta.wcsinfo.spectral_region = spectral_region


def update_s_region_mrs(output_model):
    """
    Update S_REGION for MIRI_MRS observations using the WCS transforms.

    Parameters
    ----------
    output_model : `~jwst.datamodels.IFUImageModel`
        The output of assign_wcs.
    """
    footprint, spectral_region = compute_footprint_spectral(output_model)
    update_s_region_keyword(output_model, footprint)
    output_model.meta.wcsinfo.spectral_region = spectral_region


def velocity_correction(velosys):
    """
    Compute wavelength correction to Barycentric reference frame.

    Parameters
    ----------
    velosys : float
        Radial velocity wrt Barycenter [m / s].
    """
    correction = (1 / (1 + velosys / c.value))
    model = astmodels.Identity(1) * astmodels.Const1D(correction, name="velocity_correction")
    model.inverse = astmodels.Identity(1) / astmodels.Const1D(correction, name="inv_vel_correction")

    return model


def wrap_ra(ravalues):
    """Test for 0/360 wrapping in ra values.

    If exists it makes it difficult to determine
    ra range of a region on the sky. This problem is solved by putting them all
    on "one side" of 0/360 border

    Parameters
    ----------
    ravalues : numpy.ndarray
        input RA values

    Returns
    ------
    a numpy array of ra values all on "same side" of 0/360 border
    """

    ravalues_array = np.array(ravalues)
    index_good = np.where(np.isfinite(ravalues_array))
    ravalues_wrap = ravalues_array[index_good].copy()
    median_ra = np.nanmedian(ravalues_wrap)

    # using median to test if there is any wrapping going on
    wrap_index = np.where(np.fabs(ravalues_wrap - median_ra) > 180.0)
    nwrap = wrap_index[0].size

    # get all the ra on the same "side" of 0/360
    if nwrap != 0 and median_ra < 180:
        ravalues_wrap[wrap_index] = ravalues_wrap[wrap_index] - 360.0

    if nwrap != 0 and median_ra > 180:
        ravalues_wrap[wrap_index] = ravalues_wrap[wrap_index] + 360.0

    # if the input ravaules are a list - return a list
    if isinstance(ravalues, list):
        ravalues = ravalues_wrap.tolist()

    return ravalues_wrap


def in_ifu_slice(slice_wcs, ra, dec, lam):
    """
    Given RA, DEC and LAM return the x, y positions within a slice.

    Parameters
    ----------
    slice_wcs : `~gwcs.WCS`
        Slice WCS object.
    ra, dec, lam : float, ndarray
        Physical Coordinates.

    Returns
    -------
    x, y : float, ndarray
        x, y locations within the slice.
    """
    slicer2world = slice_wcs.get_transform('slicer', 'world')
    slx, sly, sllam = slicer2world.inverse(ra, dec, lam)

    # Compute the slice X coordinate using the center of the slit.
    SLX, _, _ = slice_wcs.get_transform('slit_frame', 'slicer')(0, 0, 2e-6)
    onslice_ind = np.isclose(slx, SLX, atol=5e-4)

    return onslice_ind


def update_fits_wcsinfo(datamodel, max_pix_error=0.01, degree=None, npoints=32,
                        crpix=None, projection='TAN', imwcs=None, **kwargs):
    """
    Update ``datamodel.meta.wcsinfo`` based on a FITS WCS + SIP approximation
    of a GWCS object. By default, this function will approximate
    the datamodel's GWCS object stored in ``datamodel.meta.wcs`` but it can
    also approximate a user-supplied GWCS object when provided via
    the ``imwcs`` parameter.

    The default mode in using this attempts to achieve roughly 0.01 pixel
    accuracy over the entire image.

    This function uses the :py:meth:`~gwcs.wcs.WCS.to_fits_sip` to
    create FITS WCS representations of GWCS objects. Only most important
    :py:meth:`~gwcs.wcs.WCS.to_fits_sip` parameters are exposed here. Other
    arguments to :py:meth:`~gwcs.wcs.WCS.to_fits_sip` can be passed via
    ``kwargs`` - see "Other Parameters" section below.
    Please refer to the documentation of :py:meth:`~gwcs.wcs.WCS.to_fits_sip`
    for more details.

    .. warning::
        This function modifies input data model's ``datamodel.meta.wcsinfo``
        members.


    Parameters
    ----------

    datamodel : `ImageModel`
        The input data model for imaging or WFSS mode whose ``meta.wcsinfo``
        field should be updated from GWCS. By default, ``datamodel.meta.wcs``
        is used to compute FITS WCS + SIP approximation. When ``imwcs`` is
        not `None` then computed FITS WCS will be an approximation of the WCS
        provided through the ``imwcs`` parameter.

    max_pix_error : float, optional
        Maximum allowed error over the domain of the pixel array. This
        error is the equivalent pixel error that corresponds to the maximum
        error in the output coordinate resulting from the fit based on
        a nominal plate scale.

    degree : int, iterable, None, optional
        Degree of the SIP polynomial. Default value `None` indicates that
        all allowed degree values (``[1...6]``) will be considered and
        the lowest degree that meets accuracy requerements set by
        ``max_pix_error`` will be returned. Alternatively, ``degree`` can be
        an iterable containing allowed values for the SIP polynomial degree.
        This option is similar to default `None` but it allows caller to
        restrict the range of allowed SIP degrees used for fitting.
        Finally, ``degree`` can be an integer indicating the exact SIP degree
        to be fit to the WCS transformation. In this case
        ``max_pixel_error`` is ignored.

    npoints : int, optional
        The number of points in each dimension to sample the bounding box
        for use in the SIP fit. Minimum number of points is 3.

    crpix : list of float, None, optional
        Coordinates (1-based) of the reference point for the new FITS WCS.
        When not provided, i.e., when set to `None` (default) the reference
        pixel already specified in ``wcsinfo`` will be re-used. If
        ``wcsinfo`` does not contain ``crpix`` information, then the
        reference pixel will be chosen near the center of the bounding box
        for axes corresponding to the celestial frame.

    projection : str, `~astropy.modeling.projections.Pix2SkyProjection`, optional
        Projection to be used for the created FITS WCS. It can be specified
        as a string of three characters specifying a FITS projection code
        from Table 13 in
        `Representations of World Coordinates in FITS \
        <https://doi.org/10.1051/0004-6361:20021326>`_
        (Paper I), Greisen, E. W., and Calabretta, M. R., A & A, 395,
        1061-1075, 2002. Alternatively, it can be an instance of one of the
        `astropy's Pix2Sky_* <https://docs.astropy.org/en/stable/modeling/\
        reference_api.html#module-astropy.modeling.projections>`_
        projection models inherited from
        :py:class:`~astropy.modeling.projections.Pix2SkyProjection`.

    imwcs : `gwcs.WCS`, None, optional
        Imaging GWCS object for WFSS mode whose FITS WCS approximation should
        be computed and stored in the ``datamodel.meta.wcsinfo`` field.
        When ``imwcs`` is `None` then WCS from ``datamodel.meta.wcs``
        will be used.

        .. warning::

            Used with WFSS modes only. For other modes, supplying a different
            WCS from ``datamodel.meta.wcs`` will result in the GWCS and
            FITS WCS descriptions to diverge.

    Other Parameters
    ----------------

    max_inv_pix_error : float, None, optional
        Maximum allowed inverse error over the domain of the pixel array
        in pixel units. With the default value of `None` no inverse
        is generated.

    inv_degree : int, iterable, None, optional
        Degree of the SIP polynomial. Default value `None` indicates that
        all allowed degree values (``[1...6]``) will be considered and
        the lowest degree that meets accuracy requerements set by
        ``max_pix_error`` will be returned. Alternatively, ``degree`` can be
        an iterable containing allowed values for the SIP polynomial degree.
        This option is similar to default `None` but it allows caller to
        restrict the range of allowed SIP degrees used for fitting.
        Finally, ``degree`` can be an integer indicating the exact SIP degree
        to be fit to the WCS transformation. In this case
        ``max_inv_pixel_error`` is ignored.

    bounding_box : tuple, None, optional
        A pair of tuples, each consisting of two numbers
        Represents the range of pixel values in both dimensions
        ((xmin, xmax), (ymin, ymax))

    verbose : bool, optional
        Print progress of fits.


    Returns
    -------
    FITS header with all SIP WCS keywords


    Raises
    ------
    ValueError
        If the WCS is not at least 2D, an exception will be raised. If the
        specified accuracy (both forward and inverse, both rms and maximum)
        is not achieved an exception will be raised.


    Notes
    -----

    Use of this requires a judicious choice of required accuracies.
    Attempts to use higher degrees (~7 or higher) will typically fail due
    to floating point problems that arise with high powers.

    For more details, see :py:meth:`~gwcs.wcs.WCS.to_fits_sip`.

    """
    if crpix is None:
        crpix = [datamodel.meta.wcsinfo.crpix1, datamodel.meta.wcsinfo.crpix2]
    if None in crpix:
        crpix = None

    # For WFSS modes the imaging WCS is passed as an argument.
    # For imaging modes it is retrieved from the datamodel.
    if imwcs is None:
        imwcs = datamodel.meta.wcs

    # make a copy of kwargs:
    kwargs = {k: v for k, v in kwargs.items()}

    # override default values for "other parameters":
    max_inv_pix_error = kwargs.pop('max_inv_pix_error', None)
    inv_degree = kwargs.pop('inv_degree', None)
    if inv_degree is None:
        inv_degree = range(1, _MAX_SIP_DEGREE)

    # limit default 'degree' range to _MAX_SIP_DEGREE:
    if degree is None:
        degree = range(1, _MAX_SIP_DEGREE)

    hdr = imwcs.to_fits_sip(
        max_pix_error=max_pix_error,
        degree=degree,
        max_inv_pix_error=max_inv_pix_error,
        inv_degree=inv_degree,
        npoints=npoints,
        crpix=crpix,
        **kwargs
    )

    # update meta.wcsinfo with FITS keywords except for naxis*
    del hdr['naxis*']

    # maintain convention of lowercase keys
    hdr_dict = {k.lower(): v for k, v in hdr.items()}

    # delete naxis, cdelt, pc from wcsinfo
    rm_keys = ['naxis', 'cdelt1', 'cdelt2',
               'pc1_1', 'pc1_2', 'pc2_1', 'pc2_2',
               'a_order', 'b_order', 'ap_order', 'bp_order']

    rm_keys.extend(f"{s}_{i}_{j}" for i in range(10) for j in range(10)
                   for s in ['a', 'b', 'ap', 'bp'])

    for key in rm_keys:
        if key in datamodel.meta.wcsinfo.instance:
            del datamodel.meta.wcsinfo.instance[key]

    # update meta.wcs_info with fit keywords
    datamodel.meta.wcsinfo.instance.update(hdr_dict)

    return hdr


def wfss_imaging_wcs(wfss_model, imaging, bbox=None, **kwargs):
    """ Add a FITS WCS approximation for imaging mode to WFSS headers.

    Parameters
    ----------

    wfss_model : `~ImageModel`
        Input WFSS model (NRC or NIS).
    imaging : func, callable
        The ``imaging`` function in the ``niriss`` or ``nircam`` modules.
    bbox : tuple or None
        The bounding box over which to approximate the distortion solution.
        Typically this is based on the shape of the direct image.

    """
    xstart = wfss_model.meta.subarray.xstart
    ystart = wfss_model.meta.subarray.ystart
    reference_files = get_wcs_reference_files(wfss_model)
    image_pipeline = imaging(wfss_model, reference_files)
    imwcs = WCS(image_pipeline)
    if bbox is not None:
        imwcs.bounding_box = bbox
    elif xstart is not None and ystart is not None and (xstart != 1 or ystart != 1):
        imwcs.bounding_box = bounding_box_from_subarray(wfss_model)
    else:
        imwcs.bounding_box = wcs_bbox_from_shape(wfss_model.data.shape)

    _ = update_fits_wcsinfo(wfss_model, projection='TAN', imwcs=imwcs, bounding_box=None, **kwargs)


def get_wcs_reference_files(datamodel):
    """Retrieve names of WCS reference files for NIS_WFSS and NRC_WFSS modes.

    Parameters
    ----------

    datamodel : `~ImageModel`
        Input WFSS file (NRC or NIS).

    """
    from jwst.assign_wcs import AssignWcsStep
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        val = step.get_reference_file(datamodel, reftype)
        if val.strip() == 'N/A':
            refs[reftype] = None
        else:
            refs[reftype] = val
    return refs
