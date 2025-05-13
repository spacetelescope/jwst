import logging
import inspect
import warnings

from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils import lazyproperty
from astropy.stats import SigmaClip

import numpy as np
from photutils.detection import DAOStarFinder, IRAFStarFinder
from photutils.segmentation import SourceFinder, SourceCatalog
from photutils.background import Background2D, MedianBackground

from stdatamodels.jwst.datamodels import dqflags, ImageModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["make_tweakreg_catalog"]


class NoCatalogError(Exception):
    """Exception raised when background not determined or no sources are found in the image."""

    pass


class JWSTBackground:
    """Class to estimate a 2D background and background RMS noise in an image."""

    def __init__(self, data, box_size=100, coverage_mask=None):
        """
        Initialize the class.

        Parameters
        ----------
        data : `~numpy.ndarray`
            The input 2D image for which to estimate the background.

        box_size : int or array_like (int)
            The box size along each axis.  If ``box_size`` is a scalar then
            a square box of size ``box_size`` will be used.  If ``box_size``
            has two elements, they should be in ``(ny, nx)`` order.

        coverage_mask : array_like (bool), optional
            A boolean mask, with the same shape as ``data``, where a `True`
            value indicates the corresponding element of ``data`` is masked.
            Masked data are excluded from calculations. ``coverage_mask``
            should be `True` where there is no coverage (i.e., no data) for
            a given pixel (e.g., blank areas in a mosaic image). It should
            not be used for bad pixels.
        """
        self.data = data
        self.box_size = np.asarray(box_size).astype(int)  # must be integer
        self.coverage_mask = coverage_mask

    @lazyproperty
    def _background2d(self):
        """
        Estimate the 2D background and background RMS noise in an image.

        Returns
        -------
        background : `photutils.background.Background2D`
            A Background2D object containing the 2D background and
            background RMS noise estimates.
        """
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        filter_size = (3, 3)

        # All data have NaNs.  Suppress warnings about them.
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=AstropyUserWarning)
            try:
                bkg = Background2D(
                    self.data,
                    self.box_size,
                    filter_size=filter_size,
                    coverage_mask=self.coverage_mask,
                    sigma_clip=sigma_clip,
                    bkg_estimator=bkg_estimator,
                )
            except ValueError:
                # use the entire unmasked array
                bkg = Background2D(
                    self.data,
                    self.data.shape,
                    filter_size=filter_size,
                    coverage_mask=self.coverage_mask,
                    sigma_clip=sigma_clip,
                    bkg_estimator=bkg_estimator,
                    exclude_percentile=100.0,
                )
                log.info(
                    "Background could not be estimated in meshes. "
                    "Using the entire unmasked array for background "
                    f"estimation: bkg_boxsize={self.data.shape}."
                )

        return bkg

    @lazyproperty
    def background(self):
        """
        Compute the 2-D background if it has not been computed yet, then return it.

        Returns
        -------
        background : `~numpy.ndarray`
            The 2D background image.
        """
        return self._background2d.background

    @lazyproperty
    def background_rms(self):
        """
        Compute the 2-D background RMS image if it has not been computed yet, then return it.

        Returns
        -------
        background_rms : `~numpy.ndarray`
            The 2D background RMS image.
        """
        return self._background2d.background_rms


def _sourcefinder_wrapper(data, threshold, mask=None, **kwargs):
    """
    Make input and output of SourceFinder consistent with IRAFStarFinder and DAOStarFinder.

    Wrapper function for photutils.source_finder.SourceFinder to make input
    and output consistent with DAOStarFinder and IRAFStarFinder.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.
    threshold : float
        The absolute image value above which to select sources.
    mask : array_like (bool), optional
        The image mask
    **kwargs : dict
        Additional keyword arguments passed to `photutils.segmentation.SourceFinder`
        and/or `photutils.segmentation.SourceCatalog`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.

    References
    ----------
    `photutils segmentation tutorial <https://photutils.readthedocs.io/en/stable/segmentation.html>`_.
    """
    default_kwargs = {
        "npixels": 10,
        "progress_bar": False,
    }
    kwargs = {**default_kwargs, **kwargs}

    # handle passing kwargs into SourceFinder and SourceCatalog
    # note that this suppresses TypeError: unexpected keyword arguments
    # so user must be careful to know which kwargs are passed in here
    finder_args = list(inspect.signature(SourceFinder).parameters)
    catalog_args = list(inspect.signature(SourceCatalog).parameters)
    finder_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in finder_args}
    catalog_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in catalog_args}
    if ("kron_params" in catalog_dict.keys()) and (catalog_dict["kron_params"] is None):
        # necessary because cannot specify default in Step spec string
        catalog_dict["kron_params"] = (2.5, 1.4, 0.0)

    finder = SourceFinder(**finder_dict)
    segment_map = finder(data, threshold, mask=mask)
    sources = SourceCatalog(data, segment_map, mask=mask, **catalog_dict).to_table()
    sources.rename_column("label", "id")
    sources.rename_column("segment_flux", "flux")

    return sources


def _iraf_starfinder_wrapper(data, threshold, mask=None, **kwargs):
    """
    Make input and output of IRAFStarFinder consistent with SourceFinder and DAOStarFinder.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.
    threshold : float
        The absolute image value above which to select sources.
    mask : array_like (bool), optional
        The image mask
    **kwargs : dict
        Additional keyword arguments passed to `photutils.detection.IRAFStarFinder`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.
    """
    # defaults are not necessary to repeat here when running full pipeline step
    # but direct call to make_tweakreg_catalog will fail without 'fwhm' specified
    default_kwargs = {
        "fwhm": 2.5,
    }
    kwargs = {**default_kwargs, **kwargs}

    # note that this suppresses TypeError: unexpected keyword arguments
    # so user must be careful to know which kwargs are passed in here
    finder_args = list(inspect.signature(IRAFStarFinder).parameters)
    finder_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in finder_args}
    fwhm = finder_dict.pop("fwhm")

    starfind = IRAFStarFinder(threshold, fwhm, **finder_dict)
    sources = starfind(data, mask=mask)

    return sources


def _dao_starfinder_wrapper(data, threshold, mask=None, **kwargs):
    """
    Make input and output of DAOStarFinder consistent with SourceFinder and IRAFStarFinder.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.
    threshold : float
        The absolute image value above which to select sources.
    mask : array_like (bool), optional
        The image mask
    **kwargs : dict
        Additional keyword arguments passed to `photutils.detection.DAOStarFinder`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.
    """
    # defaults are not necessary to repeat here when running full pipeline step
    # but direct call to make_tweakreg_catalog will fail without 'fwhm' specified
    default_kwargs = {
        "fwhm": 2.5,
    }
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
    fwhm = finder_dict.pop("fwhm")
    starfind = DAOStarFinder(threshold, fwhm, **finder_dict)
    sources = starfind(data, mask=mask)

    return sources


def make_tweakreg_catalog(
    model,
    snr_threshold,
    bkg_boxsize=400,
    coverage_mask=None,
    starfinder_name="iraf",
    starfinder_kwargs=None,
):
    """
    Create a catalog of point-line sources to be used for image alignment in tweakreg.

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
    coverage_mask : array_like (bool), optional
        A boolean mask with the same shape as ``model.data``, where a `True`
        value indicates the corresponding element of ``model.data`` is masked.
        Masked pixels will not be included in any source.
    starfinder_name : str, optional
        The `photutils` star finder to use.  Options are 'dao', 'iraf', or 'segmentation'.
        - 'dao': `photutils.detection.DAOStarFinder`
        - 'iraf': `photutils.detection.IRAFStarFinder`
        - 'segmentation': `photutils.segmentation.SourceFinder`
    starfinder_kwargs : dict, optional
        Additional keyword arguments to be passed to the star finder.
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
        raise TypeError("The input model must be an ImageModel.")
    if starfinder_kwargs is None:
        starfinder_kwargs = {}

    if starfinder_name.lower() in ["dao", "daostarfinder"]:
        starfinder = _dao_starfinder_wrapper
    elif starfinder_name.lower() in ["iraf", "irafstarfinder"]:
        starfinder = _iraf_starfinder_wrapper
    elif starfinder_name.lower() in ["segmentation", "sourcefinder"]:
        starfinder = _sourcefinder_wrapper
    else:
        raise ValueError(f"Unknown starfinder type: {starfinder_name}")

    # Mask the non-imaging area, e.g. reference pixels and MIRI non-science area
    if coverage_mask is None:
        coverage_mask = (
            (dqflags.pixel["NON_SCIENCE"] + dqflags.pixel["DO_NOT_USE"]) & model.dq
        ).astype(bool)

    try:
        bkg = JWSTBackground(model.data, box_size=bkg_boxsize, coverage_mask=coverage_mask)
        with warnings.catch_warnings():
            # suppress warning about NaNs being automatically masked - this is desired
            warnings.simplefilter("ignore", AstropyUserWarning)
            threshold_img = bkg.background + (snr_threshold * bkg.background_rms)
    except ValueError as e:
        raise NoCatalogError(f"Error determining sky background: {e.args[0]}") from None

    threshold = np.median(threshold_img)  # DAOStarFinder requires float
    sources = starfinder(model.data, threshold, mask=coverage_mask, **starfinder_kwargs)
    if not sources:
        raise NoCatalogError("No sources found in the image.")

    return sources
