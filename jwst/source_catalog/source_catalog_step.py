"""Module for the source catalog step."""

import logging
from pathlib import Path

import numpy as np
from crds.core.exceptions import CrdsLookupError
from stdatamodels.jwst import datamodels

from jwst.source_catalog.reference_data import ReferenceData
from jwst.source_catalog.source_catalog import JWSTSourceCatalog
from jwst.stpipe import Step
from jwst.tweakreg.tweakreg_catalog import make_tweakreg_catalog

__all__ = ["SourceCatalogStep"]

log = logging.getLogger(__name__)


class SourceCatalogStep(Step):
    """Create a final catalog of source photometry and morphologies."""

    class_alias = "source_catalog"

    spec = """

        aperture_ee1 = integer(default=30)    # aperture encircled energy 1
        aperture_ee2 = integer(default=50)    # aperture encircled energy 2
        aperture_ee3 = integer(default=70)    # aperture encircled energy 3
        ci1_star_threshold = float(default=2.0)  # CI 1 star threshold
        ci2_star_threshold = float(default=1.8)  # CI 2 star threshold
        suffix = string(default='cat')        # Default suffix for output files
        starfinder = option('dao', 'iraf', 'segmentation', default='segmentation') # Star finder to use.

        # general starfinder options
        snr_threshold = float(default=3.0) # SNR threshold above the bkg for star finder
        bkg_boxsize = integer(default=1000) # The background mesh box size in pixels.
        kernel_fwhm = float(default=2.0) # Gaussian kernel FWHM in pixels

        # kwargs for DAOStarFinder and IRAFStarFinder, only used if starfinder is 'dao' or 'iraf'
        minsep_fwhm = float(default=0.0) # Minimum separation between detected objects in FWHM
        sigma_radius = float(default=1.5) # Truncation radius of the Gaussian kernel, units of sigma
        sharplo = float(default=0.5) # The lower bound on sharpness for object detection.
        sharphi = float(default=2.0) # The upper bound on sharpness for object detection.
        roundlo = float(default=0.0) # The lower bound on roundness for object detection.
        roundhi = float(default=0.2) # The upper bound on roundness for object detection.
        brightest = integer(default=200) # Keep top ``brightest`` objects
        peakmax = float(default=None) # Filter out objects with pixel values >= ``peakmax``

        # kwargs for SourceCatalog and SourceFinder, only used if starfinder is 'segmentation'
        npixels = integer(default=25) # Minimum number of connected pixels
        connectivity = option(4, 8, default=8) # The connectivity defining the neighborhood of a pixel
        nlevels = integer(default=32) # Number of multi-thresholding levels for deblending
        contrast = float(default=0.001) # Fraction of total source flux an object must have to be deblended
        multithresh_mode = option('exponential', 'linear', 'sinh', default='exponential') # Multi-thresholding mode
        localbkg_width = integer(default=0) # Width of rectangular annulus used to compute local background around each source
        apermask_method = option('correct', 'mask', 'none', default='correct') # How to handle neighboring sources
        kron_params = float_list(min=2, max=3, default=None) # Parameters defining Kron aperture
        deblend = boolean(default=False) # deblend sources?
    """  # noqa: E501

    reference_file_types = ["apcorr", "abvegaoffset"]

    def _get_reffile_paths(self, model):
        filepaths = []
        for reffile_type in self.reference_file_types:
            try:
                filepath = self.get_reference_file(model, reffile_type)
                log.info(f"Using {reffile_type.upper()} reference file: {filepath}")
            except CrdsLookupError as err:
                msg = f"{err} Source catalog will not be created."
                log.warning(msg)
                return None

            filepaths.append(filepath)
        return filepaths

    def process(self, input_model):
        """
        Create the catalog from the input datamodel.

        Parameters
        ----------
        input_model : str or `~stdatamodels.jwst.datamodels.ImageModel`
            A FITS filename or an `~stdatamodels.jwst.datamodels.ImageModel` of a drizzled image.

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
                log.warning(msg)
                return None

            coverage_mask = np.isnan(model.err) | (model.wht == 0)

            # convert to Jy before calling make_tweakreg_catalog so the outputs end up in Jy
            JWSTSourceCatalog.convert_mjysr_to_jy(model)

            starfinder_kwargs = {
                "sigma_radius": self.sigma_radius,
                "minsep_fwhm": self.minsep_fwhm,
                "sharplo": self.sharplo,
                "sharphi": self.sharphi,
                "roundlo": self.roundlo,
                "roundhi": self.roundhi,
                "peakmax": self.peakmax,
                "brightest": self.brightest,
                "npixels": self.npixels,
                "connectivity": int(self.connectivity),  # option returns a string, so cast to int
                "nlevels": self.nlevels,
                "contrast": self.contrast,
                "mode": self.multithresh_mode,
                "localbkg_width": self.localbkg_width,
                "apermask_method": self.apermask_method,
                "kron_params": self.kron_params,
                "deblend": self.deblend,
                "error": model.err,
                "wcs": model.meta.wcs,
                "relabel": True,
            }
            catalog, segment_img = make_tweakreg_catalog(
                model,
                self.snr_threshold,
                self.kernel_fwhm,
                bkg_boxsize=self.bkg_boxsize,
                coverage_mask=coverage_mask,
                starfinder_name=self.starfinder,
                starfinder_kwargs=starfinder_kwargs,
            )
            if len(catalog) == 0:
                log.warning("No sources found in the image. Catalog will be empty.")
                return None

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
                log.info(f"Wrote source catalog: {cat_filepath}")

                if segment_img is not None:
                    segm_model = datamodels.SegmentationMapModel(segment_img.data)
                    segm_model.update(model, only="PRIMARY")
                    segm_model.meta.wcs = model.meta.wcs
                    segm_model.meta.wcsinfo = model.meta.wcsinfo
                    self.save_model(segm_model, suffix="segm")
                    model.meta.segmentation_map = segm_model.meta.filename
                    log.info(f"Wrote segmentation map: {segm_model.meta.filename}")

        return catalog
