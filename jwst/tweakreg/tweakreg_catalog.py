import inspect
import logging
import warnings

import astropy.units as u
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import SigmaClip, gaussian_fwhm_to_sigma
from astropy.table import Table
from astropy.utils import lazyproperty
from astropy.utils.exceptions import AstropyUserWarning
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder, IRAFStarFinder
from photutils.segmentation import SourceCatalog, SourceFinder
from photutils.segmentation.catalog import DEFAULT_COLUMNS
from photutils.utils import NoDetectionsWarning
from stdatamodels.jwst.datamodels import ImageModel, dqflags

log = logging.getLogger(__name__)


__all__ = ["make_tweakreg_catalog"]


SOURCECAT_COLUMNS = DEFAULT_COLUMNS + [
    "ellipticity",
    "sky_bbox_ll",
    "sky_bbox_ul",
    "sky_bbox_lr",
    "sky_bbox_ur",
]


class JWSTBackground:
    """Class to estimate a 2D background and background RMS noise in an image."""

    def __init__(self, data, box_size=100, coverage_mask=None):
        """
        Initialize the class.

        Parameters
        ----------
        data : ndarray
            The input 2D image for which to estimate the background.

        box_size : int or array-like (int)
            The box size along each axis.  If ``box_size`` is a scalar then
            a square box of size ``box_size`` will be used.  If ``box_size``
            has two elements, they should be in ``(ny, nx)`` order.

        coverage_mask : array-like (bool), optional
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
        background : ndarray
            The 2D background image.
        """
        return self._background2d.background

    @lazyproperty
    def background_rms(self):
        """
        Compute the 2-D background RMS image if it has not been computed yet, then return it.

        Returns
        -------
        background_rms : ndarray
            The 2D background RMS image.
        """
        return self._background2d.background_rms


def _convolve_data(data, kernel_fwhm, mask=None):
    """
    Convolve the data with a Gaussian2D kernel.

    Parameters
    ----------
    data : ndarray
        The 2D array to convolve.
    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.
    mask : array-like, bool, optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.

    Returns
    -------
    convolved_data : ndarray
        The convolved 2D array.
    """
    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma)

    # All data have NaNs.  Suppress warnings about them.
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=AstropyUserWarning)
        return convolve(data, kernel, mask=mask)


def _rename_columns(sources):
    """
    Rename catalog columns and add astropy units to be consistent between the three star finders.

    Table is modified in place.

    Parameters
    ----------
    sources : `~astropy.table.QTable`
        The source catalog table for which to rename columns.
    """
    rename_map = {
        "pa": "orientation",  # for iraf and dao star finder compatibility with source_catalog step
        "npix": "area",  # for iraf and dao star finder compatibility with source_catalog step
        "segment_flux": "flux",  # for sourcefinder compatibility with tweakreg step
        "label": "id",  # for sourcefinder compatibility with tweakreg step
    }
    for old_col, new_col in rename_map.items():
        if old_col in sources.colnames:
            sources.rename_column(old_col, new_col)
    units_map = {"orientation": u.deg}
    for col, unit in units_map.items():
        if col in sources.colnames:
            sources[col] = u.Quantity(sources[col], unit=unit)


def _sourcefinder_wrapper(data, threshold_img, kernel_fwhm, mask=None, **kwargs):
    """
    Make input and output of SourceFinder consistent with IRAFStarFinder and DAOStarFinder.

    Wrapper function for photutils.source_finder.SourceFinder to make input
    and output consistent with DAOStarFinder and IRAFStarFinder.

    Parameters
    ----------
    data : array-like
        The 2D array of the image.
    threshold_img : ndarray
        The per-pixel absolute image value above which to select sources.
    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.
    mask : array-like (bool), optional
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
    :ref:`photutils segmentation tutorial <photutils:image_segmentation>`.
    """
    default_kwargs = {
        "npixels": 10,
        "progress_bar": False,
    }
    kwargs = {**default_kwargs, **kwargs}

    # convolve the data with a Gaussian kernel
    if kernel_fwhm > 0:
        conv_data = _convolve_data(data, kernel_fwhm, mask=mask)
    else:
        conv_data = data

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
    segment_map = finder(conv_data, threshold_img, mask=mask)
    if segment_map is None:
        return None, None
    sources = SourceCatalog(
        data,
        segment_map,
        mask=mask,
        convolved_data=conv_data,
        **catalog_dict,
    ).to_table(columns=SOURCECAT_COLUMNS)

    return sources, segment_map


def _iraf_starfinder_wrapper(data, threshold_img, kernel_fwhm, mask=None, **kwargs):
    """
    Make input and output of IRAFStarFinder consistent with SourceFinder and DAOStarFinder.

    Parameters
    ----------
    data : array-like
        The 2D array of the image.
    threshold_img : ndarray
        The per-pixel absolute image value above which to select sources.
    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the Gaussian kernel
    mask : array-like (bool), optional
        The image mask
    **kwargs : dict
        Additional keyword arguments passed to `photutils.detection.IRAFStarFinder`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.
    segmentation_image : ndarray or None
        The segmentation image, or None if not applicable.
    """
    # note that this suppresses TypeError: unexpected keyword arguments
    # so user must be careful to know which kwargs are passed in here
    finder_args = list(inspect.signature(IRAFStarFinder).parameters)
    finder_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in finder_args}

    threshold = np.median(threshold_img)  # only float is supported, not per-pixel value
    starfind = IRAFStarFinder(threshold, kernel_fwhm, **finder_dict)
    sources = starfind(data, mask=mask)
    return sources, None


def _dao_starfinder_wrapper(data, threshold_img, kernel_fwhm, mask=None, **kwargs):
    """
    Make input and output of DAOStarFinder consistent with SourceFinder and IRAFStarFinder.

    Parameters
    ----------
    data : array-like
        The 2D array of the image.
    threshold_img : ndarray
        The per-pixel absolute image value above which to select sources.
    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the Gaussian kernel
    mask : array-like (bool), optional
        The image mask
    **kwargs : dict
        Additional keyword arguments passed to `photutils.detection.DAOStarFinder`.

    Returns
    -------
    sources : `~astropy.table.QTable`
        A table containing the found sources.
    segmentation_image : ndarray or None
        The segmentation image, or None if not applicable.
    """
    # for consistency with IRAFStarFinder, allow minsep_fwhm to be passed in
    # and convert to pixels in the same way that IRAFStarFinder does
    # see IRAFStarFinder readthedocs page and also
    # https://github.com/astropy/photutils/issues/1561
    if "minsep_fwhm" in kwargs:
        min_sep_pix = max(2, int(kwargs["minsep_fwhm"] * kernel_fwhm + 0.5))
        kwargs["min_separation"] = min_sep_pix

    # note that this suppresses TypeError: unexpected keyword arguments
    # so user must be careful to know which kwargs are passed in here
    finder_args = list(inspect.signature(DAOStarFinder).parameters)
    finder_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in finder_args}

    threshold = np.median(threshold_img)  # only float is supported, not per-pixel value
    starfind = DAOStarFinder(threshold, kernel_fwhm, **finder_dict)
    sources = starfind(data, mask=mask)
    return sources, None


def make_tweakreg_catalog(
    model,
    snr_threshold,
    kernel_fwhm,
    bkg_boxsize=400,
    coverage_mask=None,
    starfinder_name="iraf",
    starfinder_kwargs=None,
):
    """
    Create a catalog of point-line sources to be used for image alignment in tweakreg.

    Parameters
    ----------
    model : `~stdatamodels.jwst.datamodels.ImageModel`
        The input `~stdatamodels.jwst.datamodels.ImageModel` of a single image.  The input image is
        assumed to be background subtracted.
    snr_threshold : float
        The signal-to-noise ratio per pixel above the ``background`` for
        which to consider a pixel as possibly being part of a source.
    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the Gaussian kernel
        used to convolve the image.
    bkg_boxsize : float, optional
        The background mesh box size in pixels.
    coverage_mask : array-like (bool), optional
        A boolean mask with the same shape as ``model.data``, where a `True`
        value indicates the corresponding element of ``model.data`` is masked.
        Masked pixels will not be included in any source.
    starfinder_name : str, optional
        The ``photutils`` star finder to use.  Options are 'dao', 'iraf', or 'segmentation':

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
    catalog : `~astropy.table.Table`
        An astropy Table containing the source catalog.
    segmentation_image : ndarray or None
        The segmentation image, or None if not applicable.
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

    # Compute the background and threshold image
    try:
        bkg = JWSTBackground(model.data, box_size=bkg_boxsize, coverage_mask=coverage_mask)
        with warnings.catch_warnings():
            # suppress warning about NaNs being automatically masked - this is desired
            warnings.simplefilter("ignore", AstropyUserWarning)
            threshold_img = snr_threshold * bkg.background_rms
        data = model.data - bkg.background
    except ValueError as e:
        log.warning(f"Error determining sky background: {e.args[0]}")
        sources = _empty_table()
        _rename_columns(sources)
        return sources, None

    # Run the star finder
    with warnings.catch_warnings():  # handle lack of detections later
        warnings.filterwarnings(
            "ignore", category=NoDetectionsWarning, message="No sources were found"
        )
        sources, segmentation_image = starfinder(
            data,
            threshold_img,
            kernel_fwhm,
            mask=coverage_mask,
            **starfinder_kwargs,
        )
    if not sources:
        log.warning("No sources found in the image.")
        sources = _empty_table()

    _rename_columns(sources)
    return sources, segmentation_image


def _empty_table():
    """
    Return an empty table with the correct column names and dtypes.

    Returns
    -------
    `~astropy.table.Table`
        An empty table with the correct column names and dtypes.
    """
    default_names = ["id", "xcentroid", "ycentroid", "flux"]
    default_dtypes = (int, float, float, float)
    return Table(names=default_names, dtype=default_dtypes).copy()
