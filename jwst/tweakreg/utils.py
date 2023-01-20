from copy import deepcopy

from gwcs.wcs import WCS
from tweakwcs.correctors import JWSTWCSCorrector
from tweakwcs.linearfit import build_fit_matrix

from ..assign_wcs.pointing import _v23tosky


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
