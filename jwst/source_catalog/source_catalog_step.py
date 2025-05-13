"""Module for the source catalog step."""

import warnings
from pathlib import Path

from crds.core.exceptions import CrdsLookupError
import numpy as np

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.utils.exceptions import AstropyUserWarning

from stdatamodels.jwst import datamodels

from .reference_data import ReferenceData
from .source_catalog import JWSTSourceCatalog
from ..stpipe import Step
from jwst.tweakreg.tweakreg_catalog import make_tweakreg_catalog, NoCatalogError


__all__ = ["SourceCatalogStep"]


class SourceCatalogStep(Step):
    """Create a final catalog of source photometry and morphologies."""

    class_alias = "source_catalog"

    spec = """
        bkg_boxsize = integer(default=1000)   # background mesh box size in pixels
        kernel_fwhm = float(default=2.0)      # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=3.0)    # SNR threshold above the bkg
        npixels = integer(default=25)         # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        aperture_ee1 = integer(default=30)    # aperture encircled energy 1
        aperture_ee2 = integer(default=50)    # aperture encircled energy 2
        aperture_ee3 = integer(default=70)    # aperture encircled energy 3
        ci1_star_threshold = float(default=2.0)  # CI 1 star threshold
        ci2_star_threshold = float(default=1.8)  # CI 2 star threshold
        suffix = string(default='cat')        # Default suffix for output files
    """  # noqa: E501

    reference_file_types = ["apcorr", "abvegaoffset"]

    def _get_reffile_paths(self, model):
        filepaths = []
        for reffile_type in self.reference_file_types:
            try:
                filepath = self.get_reference_file(model, reffile_type)
                self.log.info(f"Using {reffile_type.upper()} reference file: {filepath}")
            except CrdsLookupError as err:
                msg = f"{err} Source catalog will not be created."
                self.log.warning(msg)
                return None

            filepaths.append(filepath)
        return filepaths

    def process(self, input_model):
        """
        Create the catalog from the input datamodel.

        Parameters
        ----------
        input_model : str or `ImageModel`
            A FITS filename or an `ImageModel` of a drizzled image.

        Returns
        -------
        catalog : `astropy.table.Table` or None
            The source catalog, or None if no sources were found.
        """
        with datamodels.open(input_model) as model:
            reffile_paths = self._get_reffile_paths(model)
            aperture_ee = (self.aperture_ee1, self.aperture_ee2, self.aperture_ee3)

            try:
                refdata = ReferenceData(model, reffile_paths, aperture_ee)
                aperture_params = refdata.aperture_params
                abvega_offset = refdata.abvega_offset
            except RuntimeError as err:
                msg = f"{err} Source catalog will not be created."
                self.log.warning(msg)
                return None

            coverage_mask = np.isnan(model.err) | (model.wht == 0)
            convolved_data = _convolve_data(model.data, self.kernel_fwhm, mask=coverage_mask)

            starfinder_kwargs = {
                "npixels": self.npixels,
                "deblend": self.deblend,
                "connectivity": 8,
                "nlevels": 32,
                "contrast": 0.001,
                "mode": "exponential",
                "relabel": True,
                "convolved_data": convolved_data,
                "error": model.err,
                "wcs": model.meta.wcs,
            }
            try:
                catalog = make_tweakreg_catalog(
                    model,
                    self.snr_threshold,
                    bkg_boxsize=self.bkg_boxsize,
                    coverage_mask=coverage_mask,
                    starfinder_name="segmentation",
                    starfinder_kwargs=starfinder_kwargs,
                )
            except NoCatalogError as err:
                msg = f"{err} Source catalog will not be created."
                self.log.warning(msg)
                return None

            JWSTSourceCatalog.convert_to_jy(model)
            ci_star_thresholds = (self.ci1_star_threshold, self.ci2_star_threshold)
            catobj = JWSTSourceCatalog(
                model,
                catalog,
                self.kernel_fwhm,
                aperture_params,
                abvega_offset,
                ci_star_thresholds,
            )
            catalog = catobj.catalog

            if self.save_results:
                cat_filepath = self.make_output_path(ext=".ecsv")
                catalog.write(cat_filepath, format="ascii.ecsv", overwrite=True)
                model.meta.source_catalog = Path(cat_filepath).name
                self.log.info(f"Wrote source catalog: {cat_filepath}")

                # segm_model = datamodels.SegmentationMapModel(segment_img.data)
                # segm_model.update(model, only="PRIMARY")
                # segm_model.meta.wcs = model.meta.wcs
                # segm_model.meta.wcsinfo = model.meta.wcsinfo
                # self.save_model(segm_model, suffix="segm")
                # model.meta.segmentation_map = segm_model.meta.filename
                # self.log.info(f"Wrote segmentation map: {segm_model.meta.filename}")

        return catalog


def _convolve_data(data, kernel_fwhm, mask=None):
    """
    Convolve the data with a Gaussian2D kernel.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The 2D array to convolve.
    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.
    mask : array_like, bool, optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.

    Returns
    -------
    convolved_data : `~numpy.ndarray`
        The convolved 2D array.
    """
    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma)
    kernel.normalize(mode="integral")

    # All data have NaNs.  Suppress warnings about them.
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=AstropyUserWarning)
        return convolve(data, kernel, mask=mask, normalize_kernel=True)
