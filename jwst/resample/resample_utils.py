from copy import deepcopy
import logging
import math
import warnings

import asdf
import numpy as np
from astropy import units as u
from drizzle.utils import decode_context as _drizzle_decode_context

from stdatamodels.jwst.datamodels.dqflags import pixel

from stcal.alignment.util import (
    compute_scale,
    wcs_from_sregions,
)
from stcal.resample import UnsupportedWCSError
from stcal.resample.utils import compute_mean_pixel_area
from stcal.resample.utils import (
    build_driz_weight as _stcal_build_driz_weight,
    build_mask as _stcal_build_mask,
)


__all__ = [
    "build_mask",
    "decode_context",
    "make_output_wcs",
    "resampled_wcs_from_models"
]

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
    Computes the WCS of the resampled image from input models and
    specified WCS parameters.

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
            "Output WCS needs 2 coordinate axes but the "
            f"supplied WCS has {naxes} axes."
        )

    if pixel_scale is None:
        # TODO: at some point we should switch to compute_mean_pixel_area
        #       instead of compute_scale.
        pscale_in0 = compute_scale(
            ref_wcs,
            fiducial=np.array([ref_wcsinfo["ra_ref"], ref_wcsinfo["dec_ref"]])
        )
        pixel_scale = pscale_in0 * pixel_scale_ratio
        log.info(
            f"Pixel scale ratio (pscale_out/pscale_in): {pixel_scale_ratio}"
        )
        log.info(f"Computed output pixel scale: {3600 * pixel_scale} arcsec.")
    else:
        pscale_in0 = np.rad2deg(
            math.sqrt(compute_mean_pixel_area(ref_wcs, shape=shape))
        )

        pixel_scale_ratio = pixel_scale / pscale_in0
        log.info(f"Output pixel scale: {3600 * pixel_scale} arcsec.")
        log.info(
            "Computed pixel scale ratio (pscale_out/pscale_in): "
            f"{pixel_scale_ratio}."
        )

    wcs = wcs_from_sregions(
        sregion_list,
        ref_wcs=ref_wcs,
        ref_wcsinfo=ref_wcsinfo,
        pscale_ratio=pixel_scale_ratio,
        pscale=pixel_scale,
        rotation=rotation,
        shape=output_shape,
        crpix=crpix,
        crval=crval
    )

    return wcs, pscale_in0, pixel_scale, pixel_scale_ratio


def make_output_wcs(input_models, ref_wcs=None,
                    pscale_ratio=None, pscale=None, rotation=None, shape=None,
                    crpix=None, crval=None):
    """
    .. deprecated:: 1.17.2
        :py:func:`make_output_wcs` has been deprecated and will be removed
        in a future release. Use :py:func:`resampled_wcs_from_models` instead.


    Generate output WCS here based on footprints of all input WCS objects.

    Parameters
    ----------
    input_models : `~jwst.datamodel.ModelLibrary`
        The datamodels to combine into a single output WCS. Each datamodel must
        have a ``meta.wcs.s_region`` attribute.

    ref_wcs : gwcs.WCS, None, optional
        Custom WCS to use as the output WCS. If not provided,
        the reference WCS will be taken as the WCS of the first input model, with
        its bounding box adjusted to encompass all input frames.

    pscale_ratio : float, None, optional
        Ratio of input to output pixel scale. Ignored when ``pscale``
        is provided.

    pscale : float, None, optional
        Absolute pixel scale in degrees. When provided, overrides
        ``pscale_ratio``.

    rotation : float, None, optional
        Position angle of output image Y-axis relative to North.
        A value of 0.0 would orient the final output image to be North up.
        The default of `None` specifies that the images will not be rotated,
        but will instead be resampled in the default orientation for the camera
        with the x and y axes of the resampled image corresponding
        approximately to the detector axes.

    shape : tuple of int, None, optional
        Shape of the image (data array) using ``numpy.ndarray`` convention
        (``ny`` first and ``nx`` second). This value will be assigned to
        ``pixel_shape`` and ``array_shape`` properties of the returned
        WCS object.

    crpix : tuple of float, None, optional
        Position of the reference pixel in the image array. If ``crpix`` is not
        specified, it will be set to the center of the bounding box of the
        returned WCS object.

    crval : tuple of float, None, optional
        Right ascension and declination of the reference pixel. Automatically
        computed if not provided.

    Returns
    -------
    output_wcs : object
        WCS object, with defined domain, covering entire set of input frames
    """
    warnings.warn(
        "'make_output_wcs()' has been deprecated since 1.17.2 and "
        "will be removed in a future release. "
        "Use 'resampled_wcs_from_models()' instead.",
        DeprecationWarning,
        stacklevel=2
    )

    if ref_wcs is None:
        sregion_list = []
        with input_models:
            for i, model in enumerate(input_models):
                sregion_list.append(model.meta.wcsinfo.s_region)
                if i == 0:
                    example_model = model
                    ref_wcs = example_model.meta.wcs
                    ref_wcsinfo = example_model.meta.wcsinfo.instance
                input_models.shelve(model)
        naxes = ref_wcs.output_frame.naxes

        if naxes != 2:
            msg = ("Output WCS needs 2 spatial axes "
                   f"but the supplied WCS has {naxes} axes.")
            raise RuntimeError(msg)

        output_wcs = wcs_from_sregions(
            sregion_list,
            ref_wcs=ref_wcs,
            ref_wcsinfo=ref_wcsinfo,
            pscale_ratio=pscale_ratio,
            pscale=pscale,
            rotation=rotation,
            shape=shape,
            crpix=crpix,
            crval=crval
        )
        del example_model

    else:
        naxes = ref_wcs.output_frame.naxes
        if naxes != 2:
            msg = ("Output WCS needs 2 spatial axes "
                   f"but the supplied WCS has {naxes} axes.")
            raise RuntimeError(msg)
        output_wcs = deepcopy(ref_wcs)
        if shape is not None:
            output_wcs.array_shape = shape

    # Check that the output data shape has no zero-length dimensions
    if not np.prod(output_wcs.array_shape):
        msg = f"Invalid output frame shape: {tuple(output_wcs.array_shape)}"
        raise ValueError(msg)

    return output_wcs


def build_driz_weight(model, weight_type=None, good_bits=None):
    """
    .. deprecated:: 1.17.2
        :py:func:`build_driz_weight` has been deprecated and will be removed
        in a future release. Use :py:func:`stcal.utils.build_driz_weight`
        instead.


    Create a weight map for use by drizzle
    """
    warnings.warn(
        "'build_driz_weight()' has been deprecated since 1.17.2 and "
        "will be removed in a future release. "
        "Use 'stcal.utils.build_driz_weight()' instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return _stcal_build_driz_weight(
        model=model,
        weight_type=weight_type,
        good_bits=good_bits
    )


def build_mask(dqarr, bitvalue):
    """Build a bit mask from an input DQ array and a bitvalue flag

    In the returned bit mask, 1 is good, 0 is bad
    """
    return _stcal_build_mask(
        dqarr=dqarr,
        good_bits=bitvalue,
        flag_name_map=pixel
    )


def is_sky_like(frame):
    """ Checks that a frame is a sky-like frame by looking at its output units.
    If output units are either ``deg`` or ``arcsec`` the frame is considered
    a sky-like frame (as opposite to, e.g., a Cartesian frame.)
    """
    return u.Unit("deg") in frame.unit or u.Unit("arcsec") in frame.unit


def decode_context(context, x, y):
    """
    .. deprecated:: 1.17.2
        :py:func:`decode_context` has been deprecated and will be removed in a
        future release. Use :py:func:`drizzle.utils.decode_context` instead.


    Get 0-based indices of input images that contributed to (resampled)
    output pixel with coordinates ``x`` and ``y``.

    Parameters
    ----------
    context: numpy.ndarray
        A 3D `~numpy.ndarray` of integral data type.

    x: int, list of integers, numpy.ndarray of integers
        X-coordinate of pixels to decode (3rd index into the ``context`` array)

    y: int, list of integers, numpy.ndarray of integers
        Y-coordinate of pixels to decode (2nd index into the ``context`` array)

    Returns
    -------
    A list of `numpy.ndarray` objects each containing indices of input images
    that have contributed to an output pixel with coordinates ``x`` and ``y``.
    The length of returned list is equal to the number of input coordinate
    arrays ``x`` and ``y``.

    Examples
    --------
    An example context array for an output image of array shape ``(5, 6)``
    obtained by resampling 80 input images.

    >>> import numpy as np
    >>> from jwst.resample.resample_utils import decode_context
    >>> con = np.array(
    ...     [[[0, 0, 0, 0, 0, 0],
    ...       [0, 0, 0, 36196864, 0, 0],
    ...       [0, 0, 0, 0, 0, 0],
    ...       [0, 0, 0, 0, 0, 0],
    ...       [0, 0, 537920000, 0, 0, 0]],
    ...      [[0, 0, 0, 0, 0, 0,],
    ...       [0, 0, 0, 67125536, 0, 0],
    ...       [0, 0, 0, 0, 0, 0],
    ...       [0, 0, 0, 0, 0, 0],
    ...       [0, 0, 163856, 0, 0, 0]],
    ...      [[0, 0, 0, 0, 0, 0],
    ...       [0, 0, 0, 8203, 0, 0],
    ...       [0, 0, 0, 0, 0, 0],
    ...       [0, 0, 0, 0, 0, 0],
    ...       [0, 0, 32865, 0, 0, 0]]],
    ...     dtype=np.int32
    ... )
    >>> decode_context(con, [3, 2], [1, 4])
    [array([ 9, 12, 14, 19, 21, 25, 37, 40, 46, 58, 64, 65, 67, 77]),
     array([ 9, 20, 29, 36, 47, 49, 64, 69, 70, 79])]

    """
    warnings.warn(
        "'decode_context()' has been deprecated since 1.17.2 and "
        "will be removed in a future release. "
        "Use 'drizzle.utils.decode_context()' instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return _drizzle_decode_context(context, x, y)


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
            "array_shape",
            None if user_pixel_shape is None else user_pixel_shape[::-1]
        )

    if output_shape is not None:
        wcs.array_shape = output_shape[::-1]
    elif wcs.array_shape is None:
        if user_array_shape is not None:
            wcs.array_shape = user_array_shape
        elif getattr(wcs, "bounding_box", None) is not None:
            wcs.array_shape = tuple(
                int(axs[1] + 0.5)
                for axs in wcs.bounding_box.bounding_box(order="C")
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
