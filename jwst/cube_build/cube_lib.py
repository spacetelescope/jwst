import logging
import math
from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy.stats import circmean
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import pointing

log = logging.getLogger(__name__)


@dataclass
class CubeWCS:
    """A class for the WCS information of 3D spectral cubes."""

    naxis1: int
    naxis2: int
    naxis3: int
    cdelt1: float
    cdelt2: float
    cdelt3: float
    crpix1: float
    crpix2: float
    crpix3: float
    crval1: float
    crval2: float
    crval3: float
    cdelt3_normal: np.ndarray
    rot_angle: float
    xcoord: np.ndarray
    ycoord: np.ndarray
    zcoord: np.ndarray
    a_min: float
    a_max: float
    b_min: float
    b_max: float
    lambda_min: float
    lambda_max: float
    linear_wave: bool
    wavelength_table: np.ndarray | None = None


def radec2std(crval1, crval2, ra, dec, rot_angle=None):
    """
    Compute the tangent projection coordinates ``(xi,eta)`` from ``ra,dec`` using ``crval1`` and ``crval2``.

    Parameters
    ----------
    crval1 : float
      RA value of tangent point
    crval2 : float
      Dec value of tangent point
    ra : ndarray or float
      A list (or single value) of RA points to convert
    dec : ndarray  or float
      A list (or single value) of Dec points to convert
    rot_angle : float or None, optional
      Rotation angle given in degrees

    Returns
    -------
    xi : float
        X-axis tangent plane coordinate of ``ra,dec``
    eta : float
        Y-axis tangent plane coordinate of ``ra,dec``
    """  # noqa: E501
    if np.isscalar(ra):
        ra = np.asarray([ra])
        dec = np.asarray([dec])

    rad2arcsec = (180.0 * 3600.0) / math.pi
    deg2rad = math.pi / 180.0

    ra0 = crval1 * deg2rad
    dec0 = crval2 * deg2rad
    radiff = ra * deg2rad - ra0
    decr = dec * deg2rad

    h = np.sin(decr) * math.sin(dec0) + np.cos(decr) * math.cos(dec0) * np.cos(radiff)

    xi = np.cos(decr) * np.sin(radiff) / h
    eta = (np.sin(decr) * math.cos(dec0) - np.cos(decr) * math.sin(dec0) * np.cos(radiff)) / h

    xi = -xi
    xi = xi * rad2arcsec
    eta = eta * rad2arcsec

    # xi is made negative so it increases in the opposite direction
    # of ra to match the images the Parity of the ifu_cube.

    if rot_angle is not None:
        temp1 = xi * np.cos(-rot_angle * deg2rad) - eta * np.sin(-rot_angle * deg2rad)
        temp2 = xi * np.sin(-rot_angle * deg2rad) + eta * np.cos(-rot_angle * deg2rad)
        xi = temp1
        eta = temp2

    return xi, eta


def std2radec(crval1, crval2, xi, eta):
    """
    Compute ``ra,dec`` from the tangent plane rectangular coordinates.

    Compute the ``ra,dec`` values of tangent plane rectangular coordinates using
    ``crval1, crval2`` (the tangent point). This routine takes the rectangular
    plane and projects it onto the spherical plane using ``crval1, crval2`` as
    the tangent plane.

    Parameters
    ----------
    crval1 : float
      RA value of tangent point
    crval2 : float
      Dec value of tangent point
    xi : float
      Xi rectangular coordinate of tangent plane projected ``ra,dec``
    eta : float
      Eta rectangular coordinate of tangent plane projected ``ra,dec``

    Returns
    -------
    ra : float
      List (or single value) of RA computed values
    dec : float
      List (or single value) of Dec computed values
    """
    if np.isscalar(xi):
        eta = np.asarray([eta])
        xi = np.asarray([xi])

    rad2arcsec = (180.0 * 3600.0) / math.pi
    deg2rad = math.pi / 180.0

    ra0 = crval1 * deg2rad
    dec0 = crval2 * deg2rad

    # tangent projection
    xi = xi / rad2arcsec
    eta = eta / rad2arcsec
    xi = -xi  # xi is made negative so it increases in the opposite direction
    # of ra to match the images the Parity of the ifu_cube is
    # for ra is PC1_1 = -1

    ra0 = crval1 * deg2rad
    dec0 = crval2 * deg2rad

    beta = math.cos(dec0) - eta * math.sin(dec0)
    angle = xi / beta

    ra = np.arctan(angle) + ra0
    gamma = np.sqrt(xi * xi + beta * beta)
    angle = eta * math.cos(dec0) + math.sin(dec0)
    angle = angle / gamma
    dec = np.arctan(angle)
    ra = ra / deg2rad
    dec = dec / deg2rad

    mask = ra < 0
    ra[mask] += 360.0
    mask = ra > 360.0
    ra[mask] -= 360.0
    return ra, dec


def check_cube_size(naxis1, naxis2, naxis3, limit):
    """Check the size of cube and throw an exception if too large."""
    total = naxis1 * naxis2 * naxis3
    if total > limit:
        msg = f"naxis1 * naxis2 * naxis3 = {total} exceeds limit of {limit}"
        log.error(msg)
        raise ValueError(msg)


def set_wcs(
    corner_ra,
    corner_dec,
    lambda_min,
    lambda_max,
    ra_center,
    dec_center,
    rot_angle,
    nspax_x,
    nspax_y,
    cdelt1,
    cdelt2,
    cdelt3,
    linear_wave,
    wavelength_table,
):
    """
    Set up the 3D spatial/spectral World Coordinate System (WCS).

    This method calculates the spatial field of view (RA/Dec mapped to tangent plane
    coordinates xi and eta) and sets the corresponding FITS WCS standard attributes
    (CRVAL, CRPIX, NAXIS) along with the spaxel center arrays for spatial (x, y)
    and spectral (z) axes.

    Parameters
    ----------
    corner_ra : array-like of float
      Right Ascension values at the footprint corners, in degrees.
    corner_dec : array-like of float
      Declination values at the footprint corners, in degrees.
    lambda_min : float
      Minimum wavelength threshold for the spectral axis.
    lambda_max : float
      Maximum wavelength threshold for the spectral axis.
    rot_angle : float
      Position or rotation angle for the sky-to-tangent plane transformation, in degrees.

    Returns
    -------
    crval1 : float
        Reference point RA coordinate.
    crval2 : float
        Reference point Dec coordinate.
    crval3 : float
        Reference point wavelength coordinate.
    crpix1 : float
        Reference point pixel index along x-axis (1-indexed).
    crpix2 : float
        Reference point pixel index along y-axis (1-indexed).
    crpix3 : float
        Reference point pixel index along wavelength axis (1-indexed).
    naxis1 : int
        Cube pixel dimension along x.
    naxis2 : int
        Cube pixel dimension along y.
    naxis3 : int
        Cube pixel dimension along wavelength.
    a_min : float
        Minimum spatial bound in standard coordinates.
    a_max : float
        Maximum spatial bound in standard coordinates.
    b_min : float
        Minimum spatial bound in standard coordinates.
    b_max : float
        Maximum spatial bound in standard coordinates.
    xcoord : ndarray
        1D array of spaxel center positions along x.
    ycoord : ndarray
        1D array of spaxel center positions along y.
    zcoord : ndarray
        1D array of spaxel center positions along wavelength.
    lambda_min : float
        Adjusted minimum bound for linear wavelength axis.
    lambda_max : float
        Adjusted maximum bound for linear wavelength axis.
    cdelt3_normal : ndarray
        Array of wavelength bin step sizes used for flux normalization.
    """
    ra_min = np.min(corner_ra)
    ra_max = np.max(corner_ra)

    dec_min = np.min(corner_dec)
    dec_max = np.max(corner_dec)
    dec_ave = (dec_min + dec_max) / 2.0

    # we can not average ra values because of the convergence
    # of hour angles.
    ravalues = np.zeros(2)
    ravalues[0] = ra_min
    ravalues[1] = ra_max

    # astropy circmean assumes angles are in radians,
    # we have angles in degrees
    ra_ave = circmean(ravalues * u.deg).value % 360

    if ra_center is not None:
        crval1 = ra_center
    else:
        crval1 = ra_ave

    if dec_center is not None:
        crval2 = dec_center
    else:
        crval2 = dec_ave

    # find the 4 corners
    xi_corner = []
    eta_corner = []
    num = len(corner_ra)

    for i in range(num):
        xi, eta = radec2std(crval1, crval2, corner_ra[i], corner_dec[i], rot_angle)
        xi_corner.append(xi)
        eta_corner.append(eta)
    xi_min = min(xi_corner)
    xi_max = max(xi_corner)
    eta_min = min(eta_corner)
    eta_max = max(eta_corner)

    # find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
    # to find location of center abs of min values is how many pixels
    # we want a symmetric cube centered on xi,eta = 0
    xilimit = max(np.abs(xi_min), np.abs(xi_max))
    etalimit = max(np.abs(eta_min), np.abs(eta_max))

    na = math.ceil((xilimit / cdelt1).item()) + 1
    nb = math.ceil((etalimit / cdelt2).item()) + 1

    # if the user set the nspax_x or nspax_y then redefine na, nb
    # it is assumed that both values are ODD numbers
    # We want the central pixel to be the tangent point with na/nb pixels on either
    # side of central pixel.
    if nspax_x is not None:
        na = int(nspax_x / 2)
    if nspax_y is not None:
        nb = int(nspax_y / 2)

    xi_min = 0.0 - (na * cdelt1) - (cdelt1 / 2.0)
    xi_max = (na * cdelt1) + (cdelt1 / 2.0)

    eta_min = 0.0 - (nb * cdelt2) - (cdelt2 / 2.0)
    eta_max = (nb * cdelt2) + (cdelt2 / 2.0)

    crpix1 = float(na) + 1.0
    crpix2 = float(nb) + 1.0

    naxis1 = na * 2 + 1
    naxis2 = nb * 2 + 1

    a_min = xi_min
    a_max = xi_max
    b_min = eta_min
    b_max = eta_max
    # center of spaxels
    xcoord = np.zeros(naxis1)
    xstart = xi_min + cdelt1 / 2.0
    xcoord = np.arange(start=xstart, stop=xstart + naxis1 * cdelt1, step=cdelt1)
    ycoord = np.zeros(naxis2)
    ystart = eta_min + cdelt2 / 2.0
    ycoord = np.arange(start=ystart, stop=ystart + naxis2 * cdelt2, step=cdelt2)
    # depending on the naxis and cdelt values the x,ycoord can have 1 more element than naxis.
    # Clean up arrays dropping extra values at the end.
    xcoord = xcoord[0:naxis1]
    ycoord = ycoord[0:naxis2]

    # xv, yv = np.meshgrid(xcoord, ycoord)
    # xcenters = xv.flatten()
    # ycenters = yv.flatten()

    # set up the lambda (z) coordinate of the cube
    cdelt3_normal = None
    if linear_wave:
        range_lambda = lambda_max - lambda_min
        naxis3 = int(math.ceil(range_lambda / cdelt3))

        # adjust max based on integer value of naxis3
        lambda_max = lambda_min + (naxis3) * cdelt3
        zcoord = np.zeros(naxis3)
        # CRPIX3 for FITS is 1 (center of first pixel)
        # CRVAL3 then is lambda_min + cdelt3/ 2.0, which is also zcoord[0]
        # Note that these are all values at the center of a spaxel
        crval3 = lambda_min + cdelt3 / 2.0
        crpix3 = 1.0
        zstart = lambda_min + cdelt3 / 2.0
        zcoord = np.arange(start=zstart, stop=lambda_max, step=cdelt3)
        zcoord = zcoord[0:naxis3]
    else:
        naxis3 = len(wavelength_table)
        zcoord = np.asarray(wavelength_table)
        crval3 = wavelength_table[0]
        crpix3 = 1.0
    # set up the cdelt3_normal normalizing array used
    cdelt3_normal = np.zeros(naxis3)
    for j in range(naxis3 - 1):
        cdelt3_normal[j] = zcoord[j + 1] - zcoord[j]

    cdelt3_normal[naxis3 - 1] = cdelt3_normal[naxis3 - 2]

    return CubeWCS(
        naxis1=naxis1,
        naxis2=naxis2,
        naxis3=naxis3,
        cdelt1=cdelt1,
        cdelt2=cdelt2,
        cdelt3=cdelt3,
        crpix1=crpix1,
        crpix2=crpix2,
        crpix3=crpix3,
        crval1=crval1,
        crval2=crval2,
        crval3=crval3,
        cdelt3_normal=cdelt3_normal,
        rot_angle=rot_angle,
        xcoord=xcoord,
        ycoord=ycoord,
        zcoord=zcoord,
        a_min=a_min,
        a_max=a_max,
        b_min=b_min,
        b_max=b_max,
        lambda_min=lambda_min,
        lambda_max=lambda_max,
        linear_wave=linear_wave,
        wavelength_table=wavelength_table,
    )


def print_cube_geometry(coord_system, wcs):
    """Print to log the general properties of the size of the IFU cube."""
    log.info("Cube Geometry:")
    if coord_system == "internal_cal":
        log.info(
            "axis#  Naxis  CRPIX    CRVAL      CDELT(arcsec)  "
            " Min & Max (along slice, across slice)"
        )
    else:
        log.info("axis#  Naxis  CRPIX    CRVAL       CDELT(arcsec)  Min & Max (xi, eta arcsec)")
    log.info(
        "Axis 1 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f",
        wcs.naxis1,
        wcs.crpix1,
        wcs.crval1,
        wcs.cdelt1,
        wcs.a_min,
        wcs.a_max,
    )
    log.info(
        "Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f",
        wcs.naxis2,
        wcs.crpix2,
        wcs.crval2,
        wcs.cdelt2,
        wcs.b_min,
        wcs.b_max,
    )
    if wcs.linear_wave:
        log.info("axis#  Naxis  CRPIX    CRVAL      CDELT(microns)  Min & Max (microns)")
        log.info(
            "Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f",
            wcs.naxis3,
            wcs.crpix3,
            wcs.crval3,
            wcs.cdelt3,
            wcs.lambda_min,
            wcs.lambda_max,
        )

    if not wcs.linear_wave:
        log.info("Non-linear wavelength dimension; CDELT3 variable")
        log.info("axis#  Naxis  CRPIX    CRVAL     Min & Max (microns)")
        wave_array = wcs.wavelength_table  # if wcs.wavelength_table is not None else wcs.zcoord
        log.info(
            "Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f",
            wcs.naxis3,
            wcs.crpix3,
            wcs.crval3,
            wave_array[0],
            wave_array[wcs.naxis3 - 1],
        )

    if wcs.rot_angle is not None:
        log.info("Rotation angle between RA-Dec and Slicer-Plane %12.8f", wcs.rot_angle)


def find_spaxel_flux(spaxel_iflux, spaxel_flux, spaxel_weight, spaxel_var):
    """
    Calculate the final normalized flux and variance for each spaxel.

    Normalizes the weighted flux and variance arrays after mapping overlap
    between the detector array and the sky. Spaxels with no contributing points
    (where spaxel_iflux == 0) are ignored to prevent division-by-zero errors.

    Parameters
    ----------
    spaxel_iflux : numpy.ndarray
        Array containing the counts/number of detector elements contributing to each spaxel.
    spaxel_flux : numpy.ndarray
         Array of weighted cumulative flux values. Modified in-place.
    spaxel_weight : numpy.ndarray
         Array of total weights accumulated for each spaxel.Modified in-place.
    spaxel_var : numpy.ndarray
         Array of cumulative variance values. Modified in-place.

    Returns
    -------
    spaxel_flux : numpy.ndarray
         Array of weighted cumulative flux values. Modified in-place.
    spaxel_var : numpy.ndarray
         Array of cumulative variance values. Modified in-place.
    """
    # after overlapping is performed between the detector mapped array and sky find the
    # final spaxel values for flux and variance.

    # Don't apply any normalization if no points contributed to a spaxel
    # (i.e., don't divide by zero)

    good = spaxel_iflux > 0

    # Normalize the weighted sum of pixel fluxes by the sum of the weights
    spaxel_flux[good] = spaxel_flux[good] / spaxel_weight[good]
    # Normalize the variance by  the weights
    spaxel_var[good] = spaxel_var[good] / (spaxel_weight[good] * spaxel_weight[good])

    return (spaxel_flux, spaxel_var)


def set_final_dq_flags(spaxel_dq, spaxel_flux, spaxel_var, spaxel_weight):
    """
    Set up the final DQ flags.

    These flags include:
    * Good data (0)
    * NON_SCIENCE
    * DO_NOT_USE.

    Spaxels with zero total weight received no valid detector data — typically
    because all contributing pixels were saturated. For these spaxels,
    DO_NOT_USE and NON_SCIENCE flags are added to the DQ array, and the
    corresponding flux and variance values are set to NaN.

    Parameters
    ----------
    spaxel_dq : numpy.ndarray
        1D array of DQ flags for each spaxel. Modified in-place.
    spaxel_flux : numpy.ndarray
        1D array of flux values for each spaxel. Modified in-place.
    spaxel_var : numpy.ndarray
        1D array of variance values for each spaxel. Modified in-place.
    spaxel_weight : numpy.ndarray
        1D array of accumulated weights for each spaxel. Used read-only.

    Returns
    -------
    spaxel_flux : numpy.ndarray
        Flux array with zero-weight spaxels set to NaN.
    spaxel_var : numpy.ndarray
        Variance array with zero-weight spaxels set to NaN.
    spaxel_dq : numpy.ndarray
        DQ array with DO_NOT_USE and NON_SCIENCE bits set for zero-weight spaxels.
    """
    # convert all remaining spaxel_weight = 0 to DO_NOT_USE
    weight_is_zero = spaxel_weight == 0

    # If the data was only from saturating data, then the weight = 0.
    # Flux is nan, we need to ADD DO_NOT_USE

    spaxel_dq[weight_is_zero] |= (
        datamodels.dqflags.pixel["DO_NOT_USE"] | datamodels.dqflags.pixel["NON_SCIENCE"]
    )

    spaxel_flux[weight_is_zero] = np.nan
    spaxel_var[weight_is_zero] = np.nan
    return (spaxel_flux, spaxel_var, spaxel_dq)


def setup_final_model(
    model_ref,
    cube_wcs,
    spaxel_flux,
    spaxel_iflux,
    spaxel_var,
    spaxel_dq,
    instrument,
    rot_angle,
    coord_system,
    interpolation,
    weighting,
    weight_power,
    rois,
    roiw,
    roiw_table,
    output_name,
):
    """
    Set up the final meta WCS info of IFU cube along with other FITS keywords.

    Parameters
    ----------
    model_ref : `~stdatamodels.jwst.datamodels.IFUImageModel`
        The first IFU image model to use to fill in basic header values.

    Returns
    -------
    result : `~stdatamodels.jwst.datamodels.IFUCubeModel`
        IFU cube datamodel with data arrays filled in.
    """
    status = 0
    # loop over the wavelength planes to confirm each plane has some data
    # for initial or final planes that do not have any data - eliminated them
    # from the IFUcube
    # Rearrange values from 1d vectors into 3d cubes

    flux = spaxel_flux.reshape((cube_wcs.naxis3, cube_wcs.naxis2, cube_wcs.naxis1))
    wmap = spaxel_iflux.reshape((cube_wcs.naxis3, cube_wcs.naxis2, cube_wcs.naxis1))

    var = spaxel_var.reshape((cube_wcs.naxis3, cube_wcs.naxis2, cube_wcs.naxis1))
    dq = spaxel_dq.reshape((cube_wcs.naxis3, cube_wcs.naxis2, cube_wcs.naxis1))

    # fluxes can be Nan if all the input data was from saturated data.
    nan_flux = np.where(np.isnan(flux))
    dq[nan_flux] |= datamodels.dqflags.pixel["DO_NOT_USE"]

    var = np.sqrt(var)
    if cube_wcs.linear_wave:
        pixels = np.arange(cube_wcs.naxis3)

        # Calculate wavelengths
        # We add 1 to 'pixels' to convert 0-based Python indexing to 1-based FITS indexing
        wavelength_table = cube_wcs.crval3 + (pixels + 1 - cube_wcs.crpix3) * cube_wcs.cdelt3
        wave = np.asarray(wavelength_table, dtype=np.float32)
        num = len(wave)
        alldata = np.array([(wave[None].T,)], dtype=[("wavelength", "<f4", (num, 1))])
        # always write the wavetable
        cube_model = datamodels.IFUCubeModel(
            data=flux, dq=dq, err=var, weightmap=wmap, wavetable=alldata
        )
    else:
        wave = np.asarray(cube_wcs.wavelength_table, dtype=np.float32)
        num = len(wave)
        alldata = np.array([(wave[None].T,)], dtype=[("wavelength", "<f4", (num, 1))])

        cube_model = datamodels.IFUCubeModel(
            data=flux, dq=dq, err=var, weightmap=wmap, wavetable=alldata
        )

    cube_model.update(model_ref)
    if output_name is not None:
        cube_model.meta.filename = output_name

    cube_model.meta.wcsinfo.crval1 = cube_wcs.crval1
    cube_model.meta.wcsinfo.crval2 = cube_wcs.crval2
    cube_model.meta.wcsinfo.crpix1 = cube_wcs.crpix1
    cube_model.meta.wcsinfo.crpix2 = cube_wcs.crpix2

    cube_model.meta.wcsinfo.cdelt1 = cube_wcs.cdelt1 / 3600.0  # convert to degrees
    cube_model.meta.wcsinfo.cdelt2 = cube_wcs.cdelt2 / 3600.0
    # Now that we've got a pixel scale, set photometric area keywords
    # cube_model.meta.photometry.pixelarea_arcsecsq = cube_wcs.cdelt1 * cube_wcs.cdelt2
    # cube_model.meta.photometry.pixelarea_steradians = (
    #    cube_model.meta.photometry.pixelarea_arcsecsq * 2.3504e-11
    # )
    if cube_wcs.linear_wave:
        cube_model.meta.wcsinfo.crval3 = cube_wcs.crval3
        cube_model.meta.wcsinfo.cdelt3 = cube_wcs.cdelt3
        cube_model.meta.wcsinfo.ctype3 = "WAVE"
        cube_model.meta.wcsinfo.crpix3 = cube_wcs.crpix3
        cube_model.meta.ifu.roi_spatial = float(rois)
        cube_model.meta.ifu.roi_wave = float(roiw)
        # even though we are writing a WAVE-TAB we
        # do not want to set these parameters:
        #   ctype3="WAVE-TAB",
        #   ps3_0 = 'WCS-TABLE'
        #   ps3_1 = 'wavelength'
        # because some viewers (e.g. DS9) report an incorrect wavelength range
    else:
        cube_model.meta.wcsinfo.ctype3 = "WAVE-TAB"
        cube_model.meta.wcsinfo.ps3_0 = "WCS-TABLE"
        cube_model.meta.wcsinfo.ps3_1 = "wavelength"
        cube_model.meta.wcsinfo.crval3 = 1.0
        cube_model.meta.wcsinfo.crpix3 = 1.0
        cube_model.meta.wcsinfo.cdelt3 = None
        cube_model.meta.ifu.roi_wave = np.mean(roiw_table)
        cube_model.wavedim = f"(1,{num:d})"

    cube_model.meta.wcsinfo.ctype1 = "RA---TAN"
    cube_model.meta.wcsinfo.ctype2 = "DEC--TAN"
    cube_model.meta.wcsinfo.cunit1 = "deg"
    cube_model.meta.wcsinfo.cunit2 = "deg"

    cube_model.meta.wcsinfo.cunit3 = "um"
    cube_model.meta.wcsinfo.wcsaxes = 3
    cube_model.meta.wcsinfo.pc1_1 = -1
    cube_model.meta.wcsinfo.pc1_2 = 0
    cube_model.meta.wcsinfo.pc1_3 = 0

    cube_model.meta.wcsinfo.pc2_1 = 0
    cube_model.meta.wcsinfo.pc2_2 = 1
    cube_model.meta.wcsinfo.pc2_3 = 0

    cube_model.meta.wcsinfo.pc3_1 = 0
    cube_model.meta.wcsinfo.pc3_2 = 0
    cube_model.meta.wcsinfo.pc3_3 = 1

    if rot_angle is None:
        rot_angle = 0.0
    cube_model.meta.wcsinfo.pc1_1 = -np.cos(rot_angle * np.pi / 180.0)
    cube_model.meta.wcsinfo.pc1_2 = np.sin(rot_angle * np.pi / 180.0)
    cube_model.meta.wcsinfo.pc2_1 = np.sin(rot_angle * np.pi / 180.0)
    cube_model.meta.wcsinfo.pc2_2 = np.cos(rot_angle * np.pi / 180.0)

    cube_model.meta.ifu.flux_extension = "SCI"
    cube_model.meta.ifu.error_extension = "ERR"
    cube_model.meta.ifu.error_type = "ERR"
    cube_model.meta.ifu.dq_extension = "DQ"
    cube_model.meta.ifu.weighting = str(weighting)
    cube_model.meta.ifu.weight_power = weight_power

    if interpolation == "drizzle":
        # stick in values of 0, otherwise it is NaN and
        # fits file can not be written because these
        # values are defined in ifucube.schema.yaml
        cube_model.meta.ifu.weight_power = 0
        cube_model.meta.ifu.roi_wave = 0
        cube_model.meta.ifu.roi_spatial = 0
        cube_model.meta.ifu.weighting = str(interpolation)

    if coord_system == "internal_cal":
        # stick in values of 0, otherwise it is NaN and
        # fits file can not be written because these
        # values are defined in ifucube.schema.yaml
        cube_model.meta.ifu.weight_power = 0
        cube_model.meta.ifu.roi_wave = 0
        cube_model.meta.ifu.roi_spatial = 0

        # un-correct cdelt for degree conversion
        cube_model.meta.wcsinfo.cdelt1 = cube_wcs.cdelt1
        cube_model.meta.wcsinfo.cdelt2 = cube_wcs.cdelt2

        # correct "RA" axis orientation
        cube_model.meta.wcsinfo.pc1_1 *= -1.0

        if instrument == "MIRI":
            cube_model.meta.wcsinfo.cunit1 = "arcsec"
            cube_model.meta.wcsinfo.cunit2 = "arcsec"
            cube_model.meta.wcsinfo.ctype1 = "MRSALPHA"
            cube_model.meta.wcsinfo.ctype2 = "MRSBETA"

        if instrument == "NIRSPEC":
            cube_model.meta.wcsinfo.cunit1 = "meter"
            cube_model.meta.wcsinfo.cunit2 = "meter"
            cube_model.meta.wcsinfo.ctype1 = "NRSSLICERX"
            cube_model.meta.wcsinfo.ctype2 = "NRSSLICERY"

    # set WCS information
    wcsobj = pointing.create_fitswcs(cube_model)
    cube_model.meta.wcs = wcsobj
    cube_model.meta.wcs.bounding_box = (
        (0, cube_wcs.naxis1 - 1),
        (0, cube_wcs.naxis2 - 1),
        (0, cube_wcs.naxis3 - 1),
    )

    cube_model.meta.cal_step.cube_build = "COMPLETE"
    # problem with cube_build - contains only 0 data
    if status == 1:
        cube_model.meta.cal_step.cube_build = "SKIPPED"

    result = (cube_model, status)
    return result
