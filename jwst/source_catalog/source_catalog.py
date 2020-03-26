import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
import astropy.units as u
import photutils

from ..datamodels import ImageModel


def make_source_catalog(model, kernel_fwhm, kernel_xsize, kernel_ysize,
                        snr_threshold, npixels, deblend_nlevels=32,
                        deblend_contrast=0.001, deblend_mode='exponential',
                        connectivity=8, deblend=False):
    """
    Create a final catalog of source photometry and morphologies.

    Parameters
    ----------
    model : `ImageModel`
        The input `ImageModel` of a single drizzled image.  The
        input image is assumed to be background subtracted.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel
        used to filter the image before thresholding.  Filtering the
        image will smooth the noise and maximize detectability of
        objects with a shape similar to the kernel.

    kernel_xsize : odd int
        The size in the x dimension (columns) of the kernel array.

    kernel_ysize : odd int
        The size in the y dimension (row) of the kernel array.

    snr_threshold : float
        The signal-to-noise ratio per pixel above the ``background`` for
        which to consider a pixel as possibly being part of a source.

    npixels : int
        The number of connected pixels, each greater than the threshold
        that an object must have to be detected.  ``npixels`` must be a
        positive integer.

    deblend_nlevels : int, optional
        The number of multi-thresholding levels to use for deblending
        sources.  Each source will be re-thresholded at
        ``deblend_nlevels``, spaced exponentially or linearly (see the
        ``deblend_mode`` keyword), between its minimum and maximum
        values within the source segment.

    deblend_contrast : float, optional
        The fraction of the total (blended) source flux that a local
        peak must have to be considered as a separate object.
        ``deblend_contrast`` must be between 0 and 1, inclusive.  If
        ``deblend_contrast = 0`` then every local peak will be made a
        separate object (maximum deblending).  If ``deblend_contrast =
        1`` then no deblending will occur.  The default is 0.001, which
        will deblend sources with a magnitude differences of about 7.5.

    deblend_mode : {'exponential', 'linear'}, optional
        The mode used in defining the spacing between the
        multi-thresholding levels (see the ``deblend_nlevels`` keyword)
        when deblending sources.

    connectivity : {4, 8}, optional
        The type of pixel connectivity used in determining how pixels
        are grouped into a detected source.  The options are 4 or 8
        (default).  4-connected pixels touch along their edges.
        8-connected pixels touch along their edges or corners.  For
        reference, SExtractor uses 8-connected pixels.

    deblend : bool, optional
        Whether to deblend overlapping sources.  Source deblending
        requires scikit-image.

    Returns
    -------
    catalog : `~astropy.Table` or `None`
        An astropy Table containing the source photometry and
        morphologies.  If no sources are detected then `None` is
        returned.
    """

    if not isinstance(model, ImageModel):
        raise ValueError('The input model must be a ImageModel.')

    # Use this when model.wht contains an IVM map
    # Calculate "background-only" error assuming the weight image is an
    # inverse-variance map (IVM).  The weight image is clipped because it
    # may contain zeros.
    # bkg_error = np.sqrt(1.0 / np.clip(model.wht, 1.0e-20, 1.0e20))
    # threshold = snr_threshold * bkg_error

    # Estimate the 1-sigma noise in the image empirically because model.wht
    # does not yet contain an IVM map
    mask = (model.wht == 0)
    data_mean, data_median, data_std = sigma_clipped_stats(
        model.data, mask=mask, sigma=3.0, maxiters=10)
    threshold = data_median + (data_std * snr_threshold)

    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=kernel_xsize, y_size=kernel_ysize)
    kernel.normalize()

    segm = photutils.detect_sources(model.data, threshold, npixels=npixels,
                                    filter_kernel=kernel,
                                    connectivity=connectivity)

    # segm=None for photutils >= 0.7; segm.nlabels == 0 for photutils < 0.7
    if segm is None or segm.nlabels == 0:
        return None

    # source deblending requires scikit-image
    if deblend:
        segm = photutils.deblend_sources(model.data, segm, npixels=npixels,
                                         filter_kernel=kernel,
                                         nlevels=deblend_nlevels,
                                         contrast=deblend_contrast,
                                         mode=deblend_mode,
                                         connectivity=connectivity,
                                         relabel=True)

    # Calculate total error, including source Poisson noise.
    # This calculation assumes that the data and bkg_error images are in
    # units of electron/s.  Poisson noise is not included for pixels
    # where data < 0.
    exptime = model.meta.resample.product_exposure_time    # total exptime
    # total_error = np.sqrt(bkg_error**2 +
    #                       np.maximum(model.data / exptime, 0))
    total_error = np.sqrt(data_std**2 + np.maximum(model.data / exptime, 0))

    wcs = model.get_fits_wcs()
    source_props = photutils.source_properties(
        model.data, segm, error=total_error, filter_kernel=kernel, wcs=wcs)

    columns = ['id', 'xcentroid', 'ycentroid', 'sky_centroid', 'area',
               'source_sum', 'source_sum_err', 'semimajor_axis_sigma',
               'semiminor_axis_sigma', 'orientation',
               'sky_bbox_ll', 'sky_bbox_ul', 'sky_bbox_lr', 'sky_bbox_ur']
    catalog = source_props.to_table(columns=columns)

    # convert orientation to degrees
    orient_deg = catalog['orientation'].to(u.deg)
    catalog.replace_column('orientation', orient_deg)

    # define orientation position angle
    rot = _get_rotation(wcs)
    catalog['orientation_sky'] = ((270. - rot +
                                   catalog['orientation'].value) * u.deg)

    # define flux in microJanskys
    nsources = len(catalog)
    pixelarea = model.meta.photometry.pixelarea_arcsecsq
    if pixelarea is None:
        micro_Jy = np.full(nsources, np.nan)
    else:
        micro_Jy = (catalog['source_sum'] *
                    model.meta.photometry.conversion_microjanskys *
                    model.meta.photometry.pixelarea_arcsecsq)

    # define AB mag
    abmag = np.full(nsources, np.nan)
    mask = np.isfinite(micro_Jy)
    abmag[mask] = -2.5 * np.log10(micro_Jy[mask]) + 23.9
    catalog['abmag'] = abmag

    # define AB mag error
    # assuming SNR >> 1 (otherwise abmag_error is asymmetric)
    abmag_error = (2.5 * np.log10(np.e) * catalog['source_sum_err'] /
                   catalog['source_sum'])
    abmag_error[~mask] = np.nan
    catalog['abmag_error'] = abmag_error

    return catalog


def _get_rotation(wcs, skew_tolerance=0.01):
    """
    Get the rotation of an image from its WCS.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS` instance
        The image WCS.

    skew_tolerance : float, optional
        The absolute allowable difference between the x and y axis
        rotations.

    Returns
    -------
    rot : float
        The counterclockwise rotation angle of the image.
    """

    if hasattr(wcs.wcs, 'pc'):
        rot_matrix = wcs.wcs.pc
    elif hasattr(wcs.wcs, 'cd'):
        rot_matrix = wcs.wcs.cd
    else:
        raise ValueError("wcs must have a PC or CD matrix.")

    if rot_matrix[1, 0] == 0 and rot_matrix[0, 1] == 0:    # no rotation
        return 0.

    if np.linalg.det(rot_matrix) < 0:
        sgn = -1
    else:
        sgn = 1

    xrot = np.arctan2(sgn * rot_matrix[1, 0],
                      sgn * rot_matrix[0, 0]) * 180. / np.pi
    yrot = np.arctan2(-rot_matrix[0, 1], rot_matrix[1, 1]) * 180. / np.pi

    if np.abs(xrot - yrot) > skew_tolerance:
        raise ValueError('x and y axes are skewed: x_rot={0:.2f} deg, '
                         'y_rot={1:.2f} deg'.format(xrot, yrot))

    if wcs.wcs.lonpole != 180:
        xrot += (180.0 - wcs.wcs.lonpole)

    return xrot
