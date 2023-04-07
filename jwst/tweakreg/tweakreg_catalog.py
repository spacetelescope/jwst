import logging

from astropy.table import Table
import numpy as np
from photutils.detection import DAOStarFinder

from stdatamodels.jwst.datamodels import dqflags, ImageModel

from ..source_catalog.detection import JWSTBackground

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def make_tweakreg_catalog(model, kernel_fwhm, snr_threshold, sharplo=0.2,
                          sharphi=1.0, roundlo=-1.0, roundhi=1.0,
                          brightest=None, peakmax=None, bkg_boxsize=400):
    """
    Create a catalog of point-line sources to be used for image
    alignment in tweakreg.

    Parameters
    ----------
    model : `ImageModel`
        The input `ImageModel` of a single image.  The input image is
        assumed to be background subtracted.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel
        used to filter the image before thresholding.  Filtering the
        image will smooth the noise and maximize detectability of
        objects with a shape similar to the kernel.

    snr_threshold : float
        The signal-to-noise ratio per pixel above the ``background`` for
        which to consider a pixel as possibly being part of a source.

    sharplo : float, optional
        The lower bound on sharpness for object detection.

    sharphi : float, optional
        The upper bound on sharpness for object detection.

    roundlo : float, optional
        The lower bound on roundness for object detection.

    roundhi : float, optional
        The upper bound on roundness for object detection.

    brightest : int, None, optional
        Number of brightest objects to keep after sorting the full object list.
        If ``brightest`` is set to `None`, all objects will be selected.

    peakmax : float, None, optional
        Maximum peak pixel value in an object. Only objects whose peak pixel
        values are *strictly smaller* than ``peakmax`` will be selected.
        This may be used to exclude saturated sources. By default, when
        ``peakmax`` is set to `None`, all objects will be selected.

        .. warning::
            `DAOStarFinder` automatically excludes objects whose peak
            pixel values are negative. Therefore, setting ``peakmax`` to a
            non-positive value would result in exclusion of all objects.

    bkg_boxsize : float, optional
        The background mesh box size in pixels.

    Returns
    -------
    catalog : `~astropy.Table`
        An astropy Table containing the source catalog.
    """
    if not isinstance(model, ImageModel):
        raise TypeError('The input model must be an ImageModel.')

    # Mask the non-imaging area (e.g. MIRI)
    coverage_mask = ((dqflags.pixel['NON_SCIENCE'] +
                      dqflags.pixel['DO_NOT_USE']) &
                     model.dq).astype(bool)

    columns = ['id', 'xcentroid', 'ycentroid', 'flux']
    try:
        bkg = JWSTBackground(model.data, box_size=bkg_boxsize,
                             coverage_mask=coverage_mask)
        threshold_img = bkg.background + (snr_threshold * bkg.background_rms)
    except ValueError as e:
        log.warning(f"Error determining sky background: {e.args[0]}")
        # return an empty catalog
        catalog = Table(names=columns, dtype=(int, float, float, float))
        return catalog

    threshold = np.median(threshold_img)  # DAOStarFinder requires float

    daofind = DAOStarFinder(fwhm=kernel_fwhm, threshold=threshold,
                            sharplo=sharplo, sharphi=sharphi, roundlo=roundlo,
                            roundhi=roundhi, brightest=brightest,
                            peakmax=peakmax)
    sources = daofind(model.data, mask=coverage_mask)

    if sources:
        catalog = sources[columns]
    else:
        # return an empty table
        catalog = Table(names=columns, dtype=(int, float, float, float))

    return catalog
