import logging
import inspect
import warnings

from astropy.table import Table
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
from photutils.detection import DAOStarFinder, IRAFStarFinder
from photutils.segmentation import SourceFinder, SourceCatalog

from stdatamodels.jwst.datamodels import dqflags, ImageModel

from ..source_catalog.detection import JWSTBackground

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def _SourceFinderWrapper(data, threshold, mask=None, **kwargs):
    """
    Wrapper function for photutils.source_finder.SourceFinder to make input
    and output consistent with DAOStarFinder and IRAFStarFinder.
    see `photutils segmentation tutorial <https://photutils.readthedocs.io/en/stable/segmentation.html>`_.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.
    threshold : float
        The absolute image value above which to select sources.
    mask : array_like (bool), optional
        The image mask
    **kwargs : Additional keyword arguments are passed to `photutils.segmentation.SourceFinder`
        and/or `photutils.segmentation.SourceCatalog`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.
    """

    default_kwargs = {'npixels': 10,
                      'progress_bar': False,
                      }
    kwargs = {**default_kwargs, **kwargs}

    # handle passing kwargs into SourceFinder and SourceCatalog
    # note that this suppresses TypeError: unexpected keyword arguments
    # so user must be careful to know which kwargs are passed in here
    finder_args = list(inspect.signature(SourceFinder).parameters)
    catalog_args = list(inspect.signature(SourceCatalog).parameters)
    finder_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in finder_args}
    catalog_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in catalog_args}
    if ('kron_params' in catalog_dict.keys()) and (catalog_dict['kron_params'] is None):
        catalog_dict['kron_params'] = (2.5, 1.4, 0.0)  # necessary because cannot specify default in Step spec string

    finder = SourceFinder(**finder_dict)
    segment_map = finder(data, threshold, mask=mask)
    sources = SourceCatalog(data, segment_map, mask=mask, **catalog_dict).to_table()
    sources.rename_column('label', 'id')
    sources.rename_column('segment_flux', 'flux')

    return sources


def _IRAFStarFinderWrapper(data, threshold, mask=None, **kwargs):
    """
    Wrapper function for `photutils.detection.IRAFStarFinder` to make inputs
    and outputs consistent across the three allowed detection methods

    Parameters
    ----------
    data : array_like
        The 2D array of the image.
    threshold : float
        The absolute image value above which to select sources.
    mask : array_like (bool), optional
        The image mask
    **kwargs : Additional keyword arguments are passed to
        `photutils.detection.IRAFStarFinder`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.
    """

    # defaults are not necessary to repeat here when running full pipeline step
    # but direct call to make_tweakreg_catalog will fail without 'fwhm' specified
    default_kwargs = {'fwhm': 2.5,}
    kwargs = {**default_kwargs, **kwargs}

    # note that this suppresses TypeError: unexpected keyword arguments
    # so user must be careful to know which kwargs are passed in here
    finder_args = list(inspect.signature(IRAFStarFinder).parameters)
    finder_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in finder_args}
    fwhm = finder_dict.pop('fwhm')

    starfind = IRAFStarFinder(threshold, fwhm, **finder_dict)
    sources = starfind(data, mask=mask)

    return sources


def _DaoStarFinderWrapper(data, threshold, mask=None, **kwargs):
    """
    Wrapper function for `photutils.detection.DAOStarFinder` to make inputs
    and outputs consistent across the three allowed detection methods

    Parameters
    ----------
    data : array_like
        The 2D array of the image.
    threshold : float
        The absolute image value above which to select sources.
    mask : array_like (bool), optional
        The image mask
    **kwargs : Additional keyword arguments are passed to
        a`photutils.detection.DAOStarFinder`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.
    """

    # defaults are not necessary to repeat here when running full pipeline step
    # but direct call to make_tweakreg_catalog will fail without 'fwhm' specified
    default_kwargs = {'fwhm': 2.5,}
    kwargs = {**default_kwargs, **kwargs}

    # for consistency with IRAFStarFinder, allow minsep_fwhm to be passed in
    # and convert to pixels in the same way that IRAFStarFinder does
    # see IRAFStarFinder readthedocs page and also
    # https://github.com/astropy/photutils/issues/1561
    if "minsep_fwhm" in kwargs:
        min_sep_pix = max(2, int(kwargs["minsep_fwhm"] * kwargs["fwhm"] + 0.5))
        kwargs["min_separation"] = min_sep_pix

    # note that this suppresses TypeError: unexpected keyword arguments
    # so user must be careful to know which kwargs are passed in here
    finder_args = list(inspect.signature(DAOStarFinder).parameters)
    finder_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in finder_args}
    fwhm = finder_dict.pop('fwhm')
    starfind = DAOStarFinder(threshold, fwhm, **finder_dict)
    sources = starfind(data, mask=mask)

    return sources


def make_tweakreg_catalog(model, snr_threshold, bkg_boxsize=400, starfinder='iraf', starfinder_kwargs={}):
    """
    Create a catalog of point-line sources to be used for image
    alignment in tweakreg.

    Parameters
    ----------
    model : `ImageModel`
        The input `ImageModel` of a single image.  The input image is
        assumed to be background subtracted.

    snr_threshold : float
        The signal-to-noise ratio per pixel above the ``background`` for
        which to consider a pixel as possibly being part of a source.

    bkg_boxsize : float, optional
        The background mesh box size in pixels.

    starfinder : str, optional
        The `photutils` star finder to use.  Options are 'dao', 'iraf', or 'segmentation'.

            - 'dao': `photutils.detection.DAOStarFinder`
            - 'iraf': `photutils.detection.IRAFStarFinder`
            - 'segmentation': `photutils.segmentation.SourceFinder` and `photutils.segmentation.SourceCatalog`

    starfinder_kwargs : dict, optional
        additional keyword arguments to be passed to the star finder.
        for 'segmentation', these can be kwargs to `photutils.segmentation.SourceFinder`
        and/or `photutils.segmentation.SourceCatalog`.
        for 'dao' or 'iraf', these are kwargs to `photutils.detection.DAOStarFinder`
        or `photutils.detection.IRAFStarFinder`, respectively.
        Defaults are as stated in the docstrings of those functions unless noted here:

            - 'dao': fwhm=2.5
            - 'iraf': fwhm=2.5
            - 'segmentation': npixels=10, progress_bar=False

    Returns
    -------
    catalog : `~astropy.Table`
        An astropy Table containing the source catalog.
    """
    if not isinstance(model, ImageModel):
        raise TypeError('The input model must be an ImageModel.')

    if starfinder.lower() in ['dao', 'daostarfinder']:
        StarFinder = _DaoStarFinderWrapper
    elif starfinder.lower() in ['iraf', 'irafstarfinder']:
        StarFinder = _IRAFStarFinderWrapper
    elif starfinder.lower() in ['segmentation', 'sourcefinder']:
        StarFinder = _SourceFinderWrapper
    else:
        raise ValueError(f"Unknown starfinder type: {starfinder}")

    # Mask the non-imaging area, e.g. reference pixels and MIRI non-science area
    coverage_mask = ((dqflags.pixel['NON_SCIENCE'] +
                      dqflags.pixel['DO_NOT_USE']) &
                     model.dq).astype(bool)

    columns = ['id', 'xcentroid', 'ycentroid', 'flux']
    try:
        bkg = JWSTBackground(model.data, box_size=bkg_boxsize,
                             coverage_mask=coverage_mask)
        with warnings.catch_warnings():
            # suppress warning about NaNs being automatically masked - this is desired
            warnings.simplefilter('ignore', AstropyUserWarning) 
            threshold_img = bkg.background + (snr_threshold * bkg.background_rms)
    except ValueError as e:
        log.warning(f"Error determining sky background: {e.args[0]}")
        # return an empty catalog
        catalog = Table(names=columns, dtype=(int, float, float, float))
        return catalog

    threshold = np.median(threshold_img)  # DAOStarFinder requires float
    sources = StarFinder(model.data, threshold, mask=coverage_mask, **starfinder_kwargs)

    if sources:
        catalog = sources[columns]
    else:
        # return an empty table
        catalog = Table(names=columns, dtype=(int, float, float, float))

    return catalog
