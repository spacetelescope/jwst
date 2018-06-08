#!/usr/bin/env python

from ..stpipe import Step
from ..datamodels import DrizProductModel
from . import source_catalog

__all__ = ["SourceCatalogStep"]


class SourceCatalogStep(Step):
    """
    Create a final catalog of source photometry and morphologies.

    Parameters
    -----------
    input : str or `DrizProductModel`
        A FITS filename or a `DrizProductModel` of a single drizzled
        image.  The input image is assumed to be background subtracted.
    """

    spec = """
        kernel_fwhm = float(default=2.0)      # Gaussian kernel FWHM in pixels
        kernel_xsize = float(default=5)       # Kernel x size in pixels
        kernel_ysize = float(default=5)       # Kernel y size in pixels
        snr_threshold = float(default=3.0)    # SNR threshold above the bkg
        npixels = float(default=5.0)          # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        output_ext = string(default='.ecsv')  # Default type of output
        suffix = string(default='cat')        # Default suffix for output files
        save_results = boolean(default=True)  # Save the output catalog
    """

    def process(self, input):
        kernel_fwhm = self.kernel_fwhm
        kernel_xsize = self.kernel_xsize
        kernel_ysize = self.kernel_ysize
        snr_threshold = self.snr_threshold
        npixels = self.npixels
        deblend = self.deblend

        with DrizProductModel(input) as model:
            catalog = source_catalog.make_source_catalog(
                model, kernel_fwhm, kernel_xsize, kernel_ysize, snr_threshold,
                npixels, deblend=deblend)

            if len(catalog) == 0:
                self.log.info('No sources were found.  Source catalog will '
                              'not be written.')
                return

            self.log.info('Detected {0} sources'.format(len(catalog)))

            if self.save_results:
                cat_filepath = self.make_output_path(basepath=model.meta.filename)
                catalog.write(
                    cat_filepath, format='ascii.ecsv', overwrite=True
                )
                self.log.info('Wrote source catalog: {0}'
                              .format(cat_filepath))
                model.meta.source_catalog.filename = cat_filepath

        # nothing is returned because this is the last step
        return
