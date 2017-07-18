#!/usr/bin/env python

from ..stpipe import Step, cmdline
from ..datamodels import DrizProductModel
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
            self.log.info('Detected {0} sources'.format(len(catalog)))

            catalog_filename = model.meta.filename.replace('_i2d.fits',
                                                           '_cat.ecsv')
            catalog.write(catalog_filename, format='ascii.ecsv')
            self.log.info('Wrote source catalog: {0}'
                          .format(catalog_filename))
            model.meta.source_catalog.filename = catalog_filename

        # nothing is returned because this is the last step
        return


if __name__ == '__main__':
    cmdline.step_script(source_catalog_step)
