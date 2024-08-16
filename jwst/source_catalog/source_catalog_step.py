"""
Module for the source catalog step.
"""

import os

from crds.core.exceptions import CrdsLookupError
import numpy as np

from stdatamodels.jwst import datamodels

from .detection import convolve_data, JWSTBackground, JWSTSourceFinder
from .reference_data import ReferenceData
from .source_catalog import JWSTSourceCatalog
from ..stpipe import Step

__all__ = ["SourceCatalogStep"]


class SourceCatalogStep(Step):
    """
    Create a final catalog of source photometry and morphologies.

    Parameters
    -----------
    input : str or `ImageModel`
        A FITS filename or an `ImageModel` of a drizzled image.
    """

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
    """

    reference_file_types = ['apcorr', 'abvegaoffset']

    def _get_reffile_paths(self, model):
        filepaths = []
        for reffile_type in self.reference_file_types:
            try:
                filepath = self.get_reference_file(model, reffile_type)
                self.log.info(f'Using {reffile_type.upper()} reference file: '
                              f'{filepath}')
            except CrdsLookupError as err:
                msg = f'{err} Source catalog will not be created.'
                self.log.warning(msg)
                return None

            filepaths.append(filepath)
        return filepaths

    def process(self, input_model):
        with datamodels.open(input_model) as model:
            reffile_paths = self._get_reffile_paths(model)
            aperture_ee = (self.aperture_ee1, self.aperture_ee2,
                           self.aperture_ee3)

            try:
                refdata = ReferenceData(model, reffile_paths,
                                        aperture_ee)
                aperture_params = refdata.aperture_params
                abvega_offset = refdata.abvega_offset
            except RuntimeError as err:
                msg = f'{err} Source catalog will not be created.'
                self.log.warning(msg)
                return None

            coverage_mask = np.isnan(model.err) | (model.wht == 0)
            bkg = JWSTBackground(model.data, box_size=self.bkg_boxsize,
                                 coverage_mask=coverage_mask)
            model.data -= bkg.background

            threshold = self.snr_threshold * bkg.background_rms
            finder = JWSTSourceFinder(threshold, self.npixels,
                                      deblend=self.deblend)

            convolved_data = convolve_data(model.data, self.kernel_fwhm,
                                           mask=coverage_mask)
            segment_img = finder(convolved_data, mask=coverage_mask)
            if segment_img is None:
                return None

            ci_star_thresholds = (self.ci1_star_threshold,
                                  self.ci2_star_threshold)
            catobj = JWSTSourceCatalog(model, segment_img, convolved_data,
                                       self.kernel_fwhm, aperture_params,
                                       abvega_offset, ci_star_thresholds)
            catalog = catobj.catalog

            # add back background to data so input model is unchanged
            model.data += bkg.background

            if self.save_results:
                cat_filepath = self.make_output_path(ext='.ecsv')
                catalog.write(cat_filepath, format='ascii.ecsv',
                              overwrite=True)
                model.meta.source_catalog = os.path.basename(cat_filepath)
                self.log.info(f'Wrote source catalog: {cat_filepath}')

                segm_model = datamodels.SegmentationMapModel(segment_img.data)
                segm_model.update(model, only="PRIMARY")
                segm_model.meta.wcs = model.meta.wcs
                segm_model.meta.wcsinfo = model.meta.wcsinfo
                self.save_model(segm_model, suffix='segm')
                model.meta.segmentation_map = segm_model.meta.filename
                self.log.info('Wrote segmentation map: '
                              f'{segm_model.meta.filename}')

        return catalog
