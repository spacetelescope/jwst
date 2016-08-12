#!/usr/bin/env python

from ..stpipe import Step, cmdline
from ..datamodels import ModelContainer
from .tweakreg_catalog import make_tweakreg_catalog


class TweakregCatalogStep(Step):
    """
    Create a catalog of compact source positions for tweakreg.

    Parameters
    -----------
    input : str or `ModelContainer`
        An association table filename or a `ModelContainer`
        instance.
    """

    spec = """
        catalog_format = string(default='ecsv')   # Catalog output file format
        kernel_fwhm = float(default=2.5)    # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=5.0)  # SNR threshold above the bkg
    """

    def process(self, input):

        catalog_format = self.catalog_format
        kernel_fwhm = self.kernel_fwhm
        snr_threshold = self.snr_threshold

        model = ModelContainer(input)

        for image_model in model:
            catalog = make_tweakreg_catalog(image_model, kernel_fwhm,
                                            snr_threshold)
            filename = image_model.meta.filename
            self.log.info('Detected {0} sources in {1}.'.
                          format(len(catalog), filename))

            catalog_filename = filename.replace('.fits', '_cat.{0}'.
                                                format(catalog_format))
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
            image_model.meta.tweakreg_catalog.filename = catalog_filename
            image_model.catalog = catalog

        model.close()
        return model


if __name__ == '__main__':
    cmdline.step_script(tweakreg_catalog_step)
