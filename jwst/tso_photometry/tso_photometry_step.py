#!/usr/bin/env python

from ..stpipe import Step, cmdline
from ..datamodels import DrizProductModel
from . import source_catalog


class TSOPhotometryStep(Step):
    """
    Perform circular aperture photometry on Time-Series (Imaging)
    Observations (TSO).

    Parameters
    -----------
    input : str or `DrizProductModel`
        A FITS filename or an `DrizProductModel` of a single drizzled
        image.  The input image is assumed to be background subtracted.
    """

    spec = """
        catalog_format = string(default='ecsv')   # Catalog output file format
        kernel_fwhm = float(default=2.0)    # Gaussian kernel FWHM in pixels
        kernel_xsize = float(default=5)     # Kernel x size in pixels
        kernel_ysize = float(default=5)     # Kernel y size in pixels
        snr_threshold = float(default=3.0)  # SNR threshold above the bkg
        npixels = float(default=5.0)        # min number of pixels in source
        deblend = boolean(default=False)    # deblend sources?
    """

    def process(self, input):

        catalog_format = self.catalog_format
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

            catalog_filename = model.meta.filename.replace(
                '.fits', '_cat.{0}'.format(catalog_format))
            if catalog_format == 'ecsv':
                fmt = 'ascii.ecsv'
            elif catalog_format == 'fits':
                # NOTE: The catalog must not contain any 'None' values.
                #       FITS will also not clobber existing files.
                fmt = 'fits'
            else:
                raise ValueError('catalog_format must be "ecsv" or "fits".')
            catalog.write(catalog_filename, format=fmt)
            self.log.info('Wrote source catalog: {0}'.
                          format(catalog_filename))
            model.meta.source_catalog.filename = catalog_filename

        # because the is the last CALIMAGE3 step, nothing is returned
        return


if __name__ == '__main__':
    cmdline.step_script(source_catalog_step)
