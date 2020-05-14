"""
Module for the source catalog step.
"""

import os
import warnings

from crds.core.exceptions import CrdsLookupError
from photutils.utils.exceptions import NoDetectionsWarning

from .source_catalog import (ReferenceData, Background, make_kernel,
                             make_segment_img, calc_total_error,
                             SourceCatalog)
from .. import datamodels
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

    spec = """
        bkg_boxsize = float(default=100)      # background mesh box size in pixels
        kernel_fwhm = float(default=2.0)      # Gaussian kernel FWHM in pixels
        kernel_xsize = float(default=None)    # Kernel x size in pixels
        kernel_ysize = float(default=None)    # Kernel y size in pixels
        snr_threshold = float(default=3.0)    # SNR threshold above the bkg
        npixels = float(default=5.0)          # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        aperture_ee1 = float(default=30)      # aperture encircled energy 1
        aperture_ee2 = float(default=50)      # aperture encircled energy 2
        aperture_ee3 = float(default=70)      # aperture encircled energy 3
        suffix = string(default='cat')        # Default suffix for output files
    """

    reference_file_types = ['apcorr', 'abvegaoffset']

    def process(self, input_model):
        if self.kernel_xsize is not None or self.kernel_ysize is not None:
            warnings.simplefilter('default')
            warnings.warn('kernel_xsize and kernel_ysize are deprecated and '
                          'no longer used', DeprecationWarning)

        with datamodels.open(input_model) as model:
            try:
                apcorr_fn = self.get_reference_file(input_model, 'apcorr')
            except CrdsLookupError:
                apcorr_fn = None
            self.log.info(f'Using APCORR reference file {apcorr_fn}')

            try:
                abvegaoffset_fn = self.get_reference_file(input_model,
                                                          'abvegaoffset')

            except CrdsLookupError:
                abvegaoffset_fn = None
            self.log.info('Using ABVEGAOFFSET reference file '
                          f'{abvegaoffset_fn}')

            aperture_ee = (self.aperture_ee1, self.aperture_ee2,
                           self.aperture_ee3)
            try:
                refdata = ReferenceData(model, aperture_ee=aperture_ee,
                                        apcorr_filename=apcorr_fn,
                                        abvegaoffset_filename=abvegaoffset_fn)
                aperture_params = refdata.aperture_params
                abvega_offset = refdata.abvega_offset
            except RuntimeError as err:
                msg = f'{err} Source catalog will not be created.'
                self.log.warn(msg)
                return

            coverage_mask = (model.wht == 0)
            if coverage_mask.all():
                self.log.warn('There are no pixels with non-zero weight. '
                              'Source catalog will not be created.')
                return

            bkg = Background(model.data, box_size=self.bkg_boxsize,
                             mask=coverage_mask)
            model.data -= bkg.background

            threshold = self.snr_threshold * bkg.background_rms
            kernel = make_kernel(self.kernel_fwhm)
            with warnings.catch_warnings():
                # suppress NoDetectionsWarning from photutils
                warnings.filterwarnings('ignore',
                                        category=NoDetectionsWarning)
                segment_img = make_segment_img(model.data, threshold,
                                               npixels=self.npixels,
                                               kernel=kernel,
                                               mask=coverage_mask,
                                               deblend=self.deblend)
            if segment_img is None:
                self.log.warn('No sources were found. Source catalog will '
                              'not be created.')
                return
            self.log.info(f'Detected {segment_img.nlabels} sources')

            # TODO: update when model contains errors
            total_error = calc_total_error(model)

            catobj = SourceCatalog(model, segment_img, error=total_error,
                                   kernel=kernel,
                                   kernel_fwhm=self.kernel_fwhm,
                                   aperture_params=aperture_params,
                                   abvega_offset=abvega_offset)
            catalog = catobj.catalog

            if self.save_results:
                cat_filepath = self.make_output_path(ext='.ecsv')
                catalog.write(cat_filepath, format='ascii.ecsv',
                              overwrite=True)
                model.meta.source_catalog = os.path.basename(cat_filepath)
                self.log.info(f'Wrote source catalog: {cat_filepath}')

        return catalog
