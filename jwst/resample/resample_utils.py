from copy import deepcopy
import logging
import math

import asdf
import numpy as np
from astropy import units as u

from stdatamodels.jwst.datamodels.dqflags import pixel
from astropy.coordinates import SkyCoord

from stcal.alignment.util import (
    compute_scale,
    wcs_from_sregions,
)
from gwcs import wcstools

from stcal.alignment.util import compute_s_region_keyword
from stcal.resample import UnsupportedWCSError
from stcal.resample.utils import compute_mean_pixel_area
from stcal.resample.utils import build_mask as _stcal_build_mask


__all__ = ["build_mask", "resampled_wcs_from_models"]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def resampled_wcs_from_models(
    input_models,
    pixel_scale_ratio=1.0,
    pixel_scale=None,
    output_shape=None,
    rotation=None,
    crpix=None,
    crval=None,
):
    """
    Compute the WCS of the resampled image from input models and specified WCS parameters.

    Parameters
    ----------
    input_models : `~jwst.datamodel.ModelLibrary`
        Each datamodel must have a ``model.meta.wcs`` set to a ~gwcs.WCS object.
    pixel_scale_ratio : float, optional
        Desired pixel scale ratio defined as the ratio of the desired output
        pixel scale to the first input model's pixel scale computed from this
        model's WCS at the fiducial point (taken as the ``ref_ra`` and
        ``ref_dec`` from the ``wcsinfo`` meta attribute of the first input
        image). Ignored when ``pixel_scale`` is specified.
    pixel_scale : float, None, optional
        Desired pixel scale (in degrees) of the output WCS. When provided,
        overrides ``pixel_scale_ratio``.
    output_shape : tuple of two integers (int, int), None, optional
        Shape of the image (data array) using ``np.ndarray`` convention
        (``ny`` first and ``nx`` second). This value will be assigned to
        ``pixel_shape`` and ``array_shape`` properties of the returned
        WCS object.
    rotation : float, None, optional
        Position angle of output image's Y-axis relative to North.
        A value of 0.0 would orient the final output image to be North up.
        The default of `None` specifies that the images will not be rotated,
        but will instead be resampled in the default orientation for the
        camera with the x and y axes of the resampled image corresponding
        approximately to the detector axes. Ignored when ``transform`` is
        provided.
    crpix : tuple of float, None, optional
        Position of the reference pixel in the resampled image array.
        If ``crpix`` is not specified, it will be set to the center of the
        bounding box of the returned WCS object.
    crval : tuple of float, None, optional
        Right ascension and declination of the reference pixel.
        Automatically computed if not provided.

    Returns
    -------
    wcs : ~gwcs.wcs.WCS
        The WCS object corresponding to the combined input footprints.
    pscale_in : float
        Computed pixel scale (in degrees) of the first input image.
    pscale_out : float
        Computed pixel scale (in degrees) of the output image.
    pixel_scale_ratio : float
        Pixel scale ratio (output to input).
    """
    # build a list of WCS of all input models:
    sregion_list = []
    ref_wcs = None

    with input_models:
        for model in input_models:
            if ref_wcs is None:
                ref_wcsinfo = model.meta.wcsinfo.instance
                ref_wcs = deepcopy(model.meta.wcs)
                shape = model.data.shape

            sregion_list.append(model.meta.wcsinfo.s_region)
            input_models.shelve(model, modify=False)

    if not sregion_list:
        raise ValueError("No input models.")

    naxes = ref_wcs.output_frame.naxes
    if naxes != 2:
        raise UnsupportedWCSError(
            f"Output WCS needs 2 coordinate axes but the supplied WCS has {naxes} axes."
        )

    if pixel_scale is None:
        # TODO: at some point we should switch to compute_mean_pixel_area
        #       instead of compute_scale.
        pscale_in0 = compute_scale(
            ref_wcs, fiducial=np.array([ref_wcsinfo["ra_ref"], ref_wcsinfo["dec_ref"]])
        )
        pixel_scale = pscale_in0 * pixel_scale_ratio
        log.info(f"Pixel scale ratio (pscale_out/pscale_in): {pixel_scale_ratio}")
        log.info(f"Computed output pixel scale: {3600 * pixel_scale} arcsec.")
    else:
        pscale_in0 = np.rad2deg(math.sqrt(compute_mean_pixel_area(ref_wcs, shape=shape)))

        pixel_scale_ratio = pixel_scale / pscale_in0
        log.info(f"Output pixel scale: {3600 * pixel_scale} arcsec.")
        log.info(f"Computed pixel scale ratio (pscale_out/pscale_in): {pixel_scale_ratio}.")

    wcs = wcs_from_sregions(
        sregion_list,
        ref_wcs=ref_wcs,
        ref_wcsinfo=ref_wcsinfo,
        pscale_ratio=pixel_scale_ratio,
        pscale=pixel_scale,
        rotation=rotation,
        shape=output_shape,
        crpix=crpix,
        crval=crval,
    )

    return wcs, pscale_in0, pixel_scale, pixel_scale_ratio


def build_mask(dqarr, bitvalue):
    """
    Build a bit mask from an input DQ array and a bitvalue flag.

    Parameters
    ----------
    dqarr : numpy.ndarray
        Data quality array.
    bitvalue : int
        Bit value to be used for flagging good pixels.

    Returns
    -------
    numpy.ndarray
        Bit mask, where 1 is good and 0 is bad.
    """
    return _stcal_build_mask(dqarr=dqarr, good_bits=bitvalue, flag_name_map=pixel)


def is_sky_like(frame):
    """
    Check that a frame is a sky-like frame by looking at its output units.

    If output units are either ``deg`` or ``arcsec`` the frame is considered
    a sky-like frame (as opposite to, e.g., a Cartesian frame.)

    Parameters
    ----------
    frame : gwcs.WCS
        WCS object to check.

    Returns
    -------
    bool
        ``True`` if the frame is sky-like, ``False`` otherwise.
    """
    return u.Unit("deg") in frame.unit or u.Unit("arcsec") in frame.unit


def load_custom_wcs(asdf_wcs_file, output_shape=None):
    """
    Load a custom output WCS from an ASDF file.

    Parameters
    ----------
    asdf_wcs_file : str
        Path to an ASDF file containing a GWCS structure. The WCS object
        must be under the ``"wcs"`` key. Additional keys recognized by
        :py:func:`load_custom_wcs` are: ``"pixel_area"``, ``"pixel_scale"``,
        ``"pixel_shape"``, and ``"array_shape"``. The latter two are used only
        when the WCS object does not have the corresponding attributes set.
        Pixel scale and pixel area should be provided in units of ``arcsec``
        and ``arcsec**2``.
    output_shape : tuple of int, optional
        Array shape for the output data.  If not provided,
        the custom WCS must specify one of (in order of priority):
        ``array_shape``, ``pixel_shape``, or ``bounding_box``.

    Returns
    -------
    wcs_dict : dict
        A dictionary with three key/value pairs:
        ``"wcs"``, ``"pixel_scale"``, and ``"pixel_area"``.
        ``"pixel_scale"``, and ``"pixel_area"`` may be `None` if not stored in
        the ASDF file. If ``"pixel_area"`` is provided but ``"pixel_scale"``
        is not then pixel scale will be computed from pixel area assuming
        square pixels: ``pixel_scale = sqrt(pixel_area)``.
    """
    if not asdf_wcs_file:
        return None

    with asdf.open(asdf_wcs_file) as af:
        wcs = deepcopy(af.tree["wcs"])
        user_pixel_area = af.tree.get("pixel_area", None)
        user_pixel_scale = af.tree.get("pixel_scale", None)
        if user_pixel_scale is None and user_pixel_area is not None:
            user_pixel_scale = np.sqrt(user_pixel_area)

        user_pixel_shape = af.tree.get("pixel_shape", None)
        user_array_shape = af.tree.get(
            "array_shape", None if user_pixel_shape is None else user_pixel_shape[::-1]
        )

    if output_shape is not None:
        wcs.array_shape = output_shape[::-1]
    elif wcs.array_shape is None:
        if user_array_shape is not None:
            wcs.array_shape = user_array_shape
        elif getattr(wcs, "bounding_box", None) is not None:
            wcs.array_shape = tuple(
                int(axs[1] + 0.5) for axs in wcs.bounding_box.bounding_box(order="C")
            )
        else:
            raise ValueError(
                "Step argument 'output_shape' is required when custom WCS "
                "does not have 'array_shape', 'pixel_shape', or "
                "'bounding_box' attributes set."
            )

    wcs_dict = {
        "wcs": wcs,
        "pixel_area": user_pixel_area,
        "pixel_scale": user_pixel_scale,
    }
    return wcs_dict


def find_miri_lrs_sregion(sregion_model1, wcs):
    """
    Find s region for MIRI LRS resampled data.

    Parameters
    ----------
    sregion_model1 : str
        The s_regions of the first input model
    wcs : gwcs.WCS
        Spatial/spectral WCS.

    Returns
    -------
    sregion : str
        The s_region for the resample data.
    """
    # use the first sregion to set the width of the slit
    spatial_box = sregion_model1
    s = spatial_box.split(" ")
    a1 = float(s[3])
    b1 = float(s[4])
    a2 = float(s[5])
    b2 = float(s[6])
    a3 = float(s[7])
    b3 = float(s[8])
    a4 = float(s[9])
    b4 = float(s[10])

    # convert each corner to SkyCoord
    coord1 = SkyCoord(a1, b1, unit="deg")
    coord2 = SkyCoord(a2, b2, unit="deg")
    coord3 = SkyCoord(a3, b3, unit="deg")
    coord4 = SkyCoord(a4, b4, unit="deg")

    # Find the distance between the corners
    # corners are counterclockwise from 1,2,3,4
    sep1 = coord1.separation(coord2)
    sep2 = coord2.separation(coord3)
    sep3 = coord3.separation(coord4)
    sep4 = coord4.separation(coord1)

    # use the separation values so we can find the min value later
    sep = [sep1.value, sep2.value, sep3.value, sep4.value]

    # the minimum separation is the slit width
    min_sep = np.min(sep)
    min_sep = min_sep * u.deg  # set the units to degrees

    log.info(f"Estimated MIRI LRS slit width: {min_sep * 3600} arcsec.")
    # now use the combined WCS to map all pixels to the slit center
    bbox = wcs.bounding_box
    grid = wcstools.grid_from_bounding_box(bbox)
    ra, dec, _ = np.array(wcs(*grid))
    ra = ra.flatten()
    dec = dec.flatten()
    # ra and dec are the values along the output resampled slit center
    # using the first point and last point find the position angle
    star1 = SkyCoord(ra[0] * u.deg, dec[0] * u.deg, frame="icrs")
    star2 = SkyCoord(ra[-1] * u.deg, dec[-1] * u.deg, frame="icrs")
    position_angle = star1.position_angle(star2).to(u.deg)

    # 90 degrees to the position angle of the slit will define s_region
    pos_angle = position_angle - 90.0 * u.deg

    star_c1 = star1.directional_offset_by(pos_angle, min_sep / 2)
    star_c2 = star1.directional_offset_by(pos_angle, -min_sep / 2)
    star_c3 = star2.directional_offset_by(pos_angle, min_sep / 2)
    star_c4 = star2.directional_offset_by(pos_angle, -min_sep / 2)

    # set  these values to footprint
    # ra,dec corners are in counter-clockwise direction
    footprint = [
        star_c1.ra.value,
        star_c1.dec.value,
        star_c3.ra.value,
        star_c3.dec.value,
        star_c4.ra.value,
        star_c4.dec.value,
        star_c2.ra.value,
        star_c2.dec.value,
    ]
    footprint = np.array(footprint)
    s_region = compute_s_region_keyword(footprint)
    return s_region
