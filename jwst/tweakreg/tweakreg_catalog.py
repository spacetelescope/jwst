from astropy.table import Table
import numpy as np
from photutils import detect_threshold, DAOStarFinder

from ..datamodels import ImageModel


def make_tweakreg_catalog(model, kernel_fwhm, snr_threshold, sharplo=0.2,
                          sharphi=1.0, roundlo=-1.0, roundhi=1.0):
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
        The lower bound on roundess for object detection.

    roundhi : float, optional
        The upper bound on roundess for object detection.

    Returns
    -------
    catalog : `~astropy.Table`
        An astropy Table containing the source catalog.
    """

    if not isinstance(model, ImageModel):
        raise ValueError('The input model must be a ImageModel.')

    threshold_img = detect_threshold(model.data, snr=snr_threshold)
    # TODO:  use threshold image based on error array
    threshold = threshold_img[0, 0]     # constant image

    daofind = DAOStarFinder(fwhm=kernel_fwhm, threshold=threshold,
                            sharplo=sharplo, sharphi=sharphi, roundlo=roundlo,
                            roundhi=roundhi)
    sources = daofind(model.data)

    columns = ['id', 'xcentroid', 'ycentroid', 'flux']
    if len(sources) > 0:
        catalog = sources[columns]
    else:
        catalog = Table(names=columns, dtype=(np.int_, np.float_, np.float_,
                                              np.float_))

    return catalog
