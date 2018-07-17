#!/usr/bin/env python

from ..stpipe import Step, cmdline
from ..datamodels import CubeModel
from ..lib.catalog_utils import replace_suffix_ext
from .tso_photometry import tso_aperture_photometry

__all__ = ["TSOPhotometryStep"]


class TSOPhotometryStep(Step):
    """
    Perform circular aperture photometry on imaging Time Series
    Observations (TSO).

    Parameters
    -----------
    input : str or `CubeModel`
        A filename for either a FITS image or and association table or a
        `CubeModel`.
    """

    spec = """
        save_catalog = boolean(default=False)  # save exposure-level catalog
    """

    def process(self, input):
        with CubeModel(input) as model:
            xcenter = model.meta.wcsinfo.crpix1 - 1    # 1-based origin
            ycenter = model.meta.wcsinfo.crpix2 - 1    # 1-based origin

            # Define the circular aperture and circular annulus for
            # photometry.  All radii are in pixel units.

            # TODO:  move these instrument/pupil-dependent aperture
            # definitions to a reference file
            if model.meta.instrument.name == 'NIRCAM':
                if model.meta.instrument.pupil == 'WLP8':
                    # weak-lens exposure
                    radius = 50
                    radius_inner = 60
                    radius_outer = 70
                else:
                    # focused exposure
                    radius = 3
                    radius_inner = 4
                    radius_outer = 5
            elif model.meta.instrument.name == 'MIRI':
                    # focused exposure
                    # (there is no weak-lens option for MIRI)
                    radius = 3.5    # 0.385 arcsec
                    radius_inner = 9.09    # 1.0 arcsec
                    radius_outer = 10.91    # 1.2 arcsec

            catalog = tso_aperture_photometry(model, xcenter, ycenter,
                                              radius, radius_inner,
                                              radius_outer)

            if self.save_catalog:
                old_suffixes = ['calints', 'crfints']
                output_dir = self.search_attr('output_dir')
                cat_filepath = replace_suffix_ext(model.meta.filename,
                                                  old_suffixes, 'phot',
                                                  output_ext='ecsv',
                                                  output_dir=output_dir)
                catalog.write(cat_filepath, format='ascii.ecsv',
                              overwrite=True)
                self.log.info('Wrote TSO photometry catalog: {0}'.
                              format(cat_filepath))

        return catalog
