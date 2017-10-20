#!/usr/bin/env python

from ..stpipe import Step, cmdline
from ..datamodels import DrizProductModel
from ..lib.catalog_utils import replace_suffix_ext
from . import source_catalog


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
        kernel_fwhm = float(default=2.0)    # Gaussian kernel FWHM in pixels
        kernel_xsize = float(default=5)     # Kernel x size in pixels
        kernel_ysize = float(default=5)     # Kernel y size in pixels
        snr_threshold = float(default=3.0)  # SNR threshold above the bkg
        npixels = float(default=5.0)        # min number of pixels in source
        deblend = boolean(default=False)    # deblend sources?
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

            old_suffixes = ['i2d']
            output_dir = self.search_attr('output_dir')
            cat_filepath = replace_suffix_ext(model.meta.filename,
                                              old_suffixes, 'cat',
                                              output_ext='ecsv',
                                              output_dir=output_dir)
            catalog.write(cat_filepath, format='ascii.ecsv', overwrite=True)
            self.log.info('Wrote source catalog: {0}'
                          .format(cat_filepath))
            model.meta.source_catalog.filename = cat_filepath

        # nothing is returned because this is the last step
        return


if __name__ == '__main__':
    cmdline.step_script(source_catalog_step)
