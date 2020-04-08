#!/usr/bin/env python

import os

from .source_catalog import SourceCatalog
from .. import datamodels
from ..stpipe import Step

__all__ = ["SourceCatalogStep"]


class SourceCatalogStep(Step):
    """
    Create a final catalog of source photometry and morphologies.

    Parameters
    -----------
    input : str or `ImageModel`
        A FITS filename or a `ImageModel` of a single drizzled
        image.  The input image is assumed to be background subtracted.
    """

    spec = """
        bkg_boxsize = float(default=100)      # background mesh box size in pixels
        kernel_fwhm = float(default=2.0)      # Gaussian kernel FWHM in pixels
        kernel_xsize = float(default=5)       # Kernel x size in pixels
        kernel_ysize = float(default=5)       # Kernel y size in pixels
        snr_threshold = float(default=3.0)    # SNR threshold above the bkg
        npixels = float(default=5.0)          # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        aperture_ee1 = float(default=30)      # aperture encircled energy 1
        aperture_ee2 = float(default=50)      # aperture encircled energy 2
        aperture_ee3 = float(default=70)      # aperture encircled energy 3
        output_ext = string(default='.ecsv')  # Default file extension
        suffix = string(default='cat')        # Default suffix for output files
    """

    reference_file_types = ['apcorr', 'abvega_offset']

    def process(self, input_model):
        with datamodels.open(input_model) as model:
            apcorr_filename = self.get_reference_file(input_model, 'apcorr')
            self.log.info(f'Using apcorr reference file {apcorr_filename}')

            # TODO: fix after reference file is in CRDS
            #abvega_offset_filename = self.get_reference_file(
            #    input_model, 'abvega_offset')
            abvega_offset_filename = None
            self.log.info('Using abvega_offset reference file ' +
                          f'{abvega_offset_filename}')

            if self.aperture_ee2 < self.aperture_ee1:
                raise ValueError('aperture_ee2 value must be less than '
                                 'aperture_ee1')
            if self.aperture_ee3 < self.aperture_ee2:
                raise ValueError('aperture_ee3 value must be less than '
                                 'aperture_ee2')
            aperture_ee = (int(self.aperture_ee1), int(self.aperture_ee2),
                           int(self.aperture_ee3))

            catobj = SourceCatalog(model, bkg_boxsize=self.bkg_boxsize,
                                   kernel_fwhm=self.kernel_fwhm,
                                   snr_threshold=self.snr_threshold,
                                   npixels=self.npixels, deblend=self.deblend,
                                   aperture_ee=aperture_ee,
                                   apcorr_filename=apcorr_filename,
                                   abvega_offset_filename=abvega_offset_filename)
            catalog = catobj.run()

            if self.save_results:
                cat_filepath = self.make_output_path()
                catalog.write(cat_filepath, format='ascii.ecsv',
                              overwrite=True)
                model.meta.source_catalog = os.path.basename(cat_filepath)
                self.log.info(f'Wrote source catalog: {cat_filepath}')

        return catalog
