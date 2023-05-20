from copy import deepcopy

from astropy.modeling.rotations import RotationSequence3D
from gwcs.wcs import WCS
import numpy as np
from tweakwcs.correctors import JWSTWCSCorrector
from tweakwcs.linearfit import build_fit_matrix

from stdatamodels.jwst.datamodels import ImageModel
from ..assign_wcs.util import update_fits_wcsinfo
from ..assign_wcs.pointing import _v23tosky


_RAD2ARCSEC = 3600.0 * np.rad2deg(1.0)


def _wcsinfo_from_wcs_transform(wcs):
    frames = wcs.available_frames
    if 'v2v3' not in frames or 'world' not in frames or frames[-1] != 'world':
        raise ValueError(
            "Unsupported WCS structure."
        )

    # Initially get v2_ref, v3_ref, and roll_ref from
    # the v2v3 to world transform. Also get ra_ref, dec_ref
    t = wcs.get_transform(frames[-2], 'world')
    for m in t:
        if isinstance(m, RotationSequence3D) and m.parameters.size == 5:
            v2_ref, nv3_ref, roll_ref, dec_ref, nra_ref = m.angles.value
            break
    else:
        raise ValueError(
            "Unsupported WCS structure."
        )

    # overwrite v2_ref, v3_ref, and roll_ref with
    # values from the tangent plane when available:
    if 'v2v3corr' in frames:
        # get v2_ref, v3_ref, and roll_ref from
        # the v2v3 to v2v3corr transform:
        frm1 = 'v2v3vacorr' if 'v2v3vacorr' in frames else 'v2v3'
        tpcorr = wcs.get_transform(frm1, 'v2v3corr')
        v2_ref, nv3_ref, roll_ref = tpcorr['det_to_optic_axis'].angles.value

    wcsinfo = {
        'v2_ref': 3600 * v2_ref,
        'v3_ref': -3600 * nv3_ref,
        'roll_ref': roll_ref,
        'ra_ref': -nra_ref,
        'dec_ref': dec_ref
    }

    return wcsinfo


def adjust_wcs(wcs, delta_ra=0.0, delta_dec=0.0, delta_roll=0.0,
               scale_factor=1.0):
    """
    Apply corrections to an imaging WCS of 'cal' data models.

    .. warning::
        This function is not designed to handle neither FITS WCS nor
        GWCS of resampled images. It is designed specifically for GWCS of
        calibrated imaging data models that can be used as input to Stage 3
        of the JWST pipeline (with suffixes '_cal', '_tweakreg', '_skymatch').

    Parameters
    ----------
    wcs : `gwcs.WCS`
        WCS object to be adjusted. Must be an imaging JWST WCS of a calibrated
        data model.

    delta_ra : float, optional
        Additional rotation (in degrees) to be applied along the longitude
        direction.

    delta_dec : float, optional
        Additional rotation (in degrees) to be applied along the latitude
        direction.

    delta_roll : float, optional
        Additional rotation (in degrees) to be applied to the telescope roll
        angle (rotation about V1 axis).

    scale_factor : float, optional
        A multiplicative scale factor to be applied to the current scale
        (if any) in the WCS. If input ``wcs`` does not have a scale factor
        applied, it is assumed to be 1. The scale factor is applied in
        a tangent plane perpendicular to the V1 axis of the telescope.

    Returns
    -------

    wcs : `gwcs.WCS`
        Adjusted WCS object.
    """
    # find the last frame in the pipeline that starts with 'v2v3':
    pipeline = deepcopy(wcs.pipeline)
    for step in pipeline[::-1]:
        if (
            step.frame.name and
            step.frame.name.startswith('v2v3') and
            step.transform.name == 'v23tosky'
        ):
            s = step
            break
    else:
        raise ValueError("Unknown WCS structure")

    v2, v3, roll, dec, ra = s.transform.parameters[-5:]
    v3 = -v3
    ra = -ra + delta_ra
    dec += delta_dec
    roll += delta_roll

    s.transform = _v23tosky(v2, v3, roll, ra, dec)

    wcs = WCS(pipeline)
    if scale_factor != 1.0:
        # apply scale factor in the tangent plane:
        corr = JWSTWCSCorrector(
            wcs,
            {
                'v2_ref': 3600.0 * v2,
                'v3_ref': 3600.0 * v3,
                'roll_ref': 0.0
            }
        )
        corr.set_correction(matrix=build_fit_matrix(0.0, scale_factor))
        wcs = corr.wcs

    return wcs


def transfer_wcs_correction(to_image, from_image, matrix=None, shift=None):
    """
    Applies the same *total* WCS correction that was applied to the WCS in the
    ``from_image`` data model to the WCS of the ``to_image`` data model.
    In some ways this function is analogous function to the ``tweakback``
    function for HST available in the
    `drizzlepac package <https://github.com/spacetelescope/drizzlepac>`_.

    One fundamental difference between this function and ``tweakback``
    is that JWST data models do not keep a history
    of data's WCS via alternative WCS as it is done in HST data and so it is
    impossible to select and apply only one particular WCS correction if
    there were multiple corrections previously applied to a WCS. The
    tangent-plane correction in JWST WCS is cumulative/total correction.
    If you would like to apply a specific/custom correction, you can do that
    via ``matrix`` and ``shift`` arguments which is defined in the
    reference tangent plane provided by the ``from_image``'s WCS.

    When providing your own corrections via ``matrix`` and ``shift`` arguments,
    this function is similar to the :py:func:`adjust_wcs` function
    but provides an alternative way of specifying corrections via affine
    transformations in a reference tangent plane.

    .. warning::
        Upon return, if the ``to_image`` argument is an `ImageModel` it will be
        modified with an updated ``ImageModel.meta.wcs`` WCS model.
        If ``to_image`` argument is a file name of an ``ImageModel``, that
        model will be read in, its WCS will be updated, and the updated model
        will be written out to the same file. BACKUP the file in ``to_image``
        argument before calling this function.

    Parameters
    ----------

    to_image : str, ImageModel
        Image model to which the correction should be applied/transferred to.

        .. warning::
            If it is a string file name then, upon return, this file
            will be **overwritten** with a data model with an updated WCS.

    from_image : str, ImageModel, gwcs.wcs.WCS
        A data model whose WCS was previously corrected.
        This data model plays two roles: 1) it is the reference WCS which
        provides a tangent plane in which corrections have been defined, and
        2) it provides WCS corrections to be applied to ``to_image``'s
        WCS.

        If the WCS of the ``from_image`` data model does not contain
        corections, then *both* ``matrix`` *and* ``shift`` arguments *must be
        supplied*.

    matrix : 2D list, 2D numpy.ndarray, None, optional
        A 2D matrix part of an affine transformation defined in the tangent
        plane derived from the ``from_image``'s WCS.

        .. note::
            When provided, ``shift`` argument *must also be provided* in
            which case ``matrix`` and ``shift`` arguments override the
            correction (if present) from the  ``from_file``'s WCS.

    shift : list, numpy.ndarray, None, optional
        A list of length 2 representing the translational part of an affine
        transformation (in arcsec) defined in the tangent plane derived
        from the ``from_image``'s WCS.

        .. note::
            When provided, ``matrix`` argument *must also be provided* in
            which case ``matrix`` and ``shift`` arguments override the
            correction (if present) from the  ``from_file``'s WCS.

    Returns
    -------

    Upon return, if the ``to_image`` argument is an `ImageModel` it will be
    modified with an updated ``ImageModel.meta.wcs`` WCS model.
    If ``to_image`` argument is a file name of an ``ImageModel``, that
    model will be read in, its WCS will be updated, and the updated model
    will be written to the same file. BACKUP the file in ``to_image``
    argument before calling this function.

    """
    if isinstance(to_image, str):
        to_file = to_image
        to_image = ImageModel(to_file)
    else:
        to_file = None

    if isinstance(from_image, str):
        from_image = ImageModel(from_image)
        refwcs = from_image.meta.wcs
        ref_wcsinfo = from_image.meta.wcsinfo.instance

    elif isinstance(from_image, WCS):
        refwcs = from_image
        ref_wcsinfo = _wcsinfo_from_wcs_transform(refwcs)

    else:  # an ImageModel object
        refwcs = from_image.meta.wcs
        ref_wcsinfo = from_image.meta.wcsinfo.instance

    if ((matrix is None and shift is not None) or
            (matrix is not None and shift is None)):
        raise ValueError(
            "Both 'matrix' and 'shift' must be either None or not None."
        )

    elif matrix is None and shift is None:
        # use the correction that was applied to the "reference" WCS - the
        # WCS from the "from_image":
        if 'v2v3corr' not in refwcs.available_frames:
            raise ValueError(
                "The WCS of the 'from_image' data model has not been "
                "previously corrected. A corrected WCS is needed in order to "
                "\"transfer\" it to the 'to_file' when both 'matrix' and "
                "'shift' are None."
            )

        t = refwcs.get_transform('v2v3', 'v2v3corr')
        affine = t['tp_affine']
        matrix = affine.matrix.value
        shift = _RAD2ARCSEC * affine.translation.value  # convert to arcsec

    from_corr = JWSTWCSCorrector(
        refwcs,
        ref_wcsinfo
    )

    to_corr = JWSTWCSCorrector(
        to_image.meta.wcs,
        to_image.meta.wcsinfo.instance
    )

    to_corr.set_correction(matrix=matrix, shift=shift, ref_tpwcs=from_corr)
    to_image.meta.wcs = to_corr.wcs

    update_fits_wcsinfo(to_image)

    if to_file:
        to_image.write(to_file)
