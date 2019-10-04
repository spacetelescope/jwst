#!/usr/bin/env python

from ..stpipe import Step
from ..datamodels import CubeModel
from ..datamodels import TsoPhotModel
from ..lib.catalog_utils import replace_suffix_ext
from .tso_photometry import tso_aperture_photometry

__all__ = ['TSOPhotometryStep']


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

    reference_file_types = ['tsophot']

    def process(self, input):
        with CubeModel(input) as model:

            xcenter = model.meta.wcsinfo.crpix1 - 1    # 1-based origin
            ycenter = model.meta.wcsinfo.crpix2 - 1    # 1-based origin

            tsophot_filename = self.get_reference_file(model, 'tsophot')
            self.log.debug('Reference file name = {}'.format(tsophot_filename))
            if tsophot_filename == 'N/A':
                self.log.warning('No TSOPHOT reference file found;')
                self.log.warning('the tso_photometry step will be skipped.')
                return None

            pupil_name = 'ANY'
            if model.meta.instrument.pupil is not None:
                pupil_name = model.meta.instrument.pupil

            (radius, radius_inner, radius_outer) = get_ref_data(
                        tsophot_filename, pupil=pupil_name)
            self.log.debug('Using reference file {}'.format(tsophot_filename))
            self.log.debug('radius = {}'.format(radius))
            self.log.debug('radius_inner = {}'.format(radius_inner))
            self.log.debug('radius_outer = {}'.format(radius_outer))

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


def get_ref_data(reffile, pupil='ANY'):

    ref_model = TsoPhotModel(reffile)
    radii = ref_model.radii
    value = None
    val_any_pupil = None
    for item in radii:
        if item.pupil == pupil.upper():
            value = (item.radius,
                     item.radius_inner, item.radius_outer)
            break
        elif item.pupil == 'ANY' and val_any_pupil is None:
            # Save this value as a fallback, in case we don't find a match
            # to an actual pupil name.
            val_any_pupil = (item.radius,
                             item.radius_inner, item.radius_outer)

    if value is not None:
        (radius, radius_inner, radius_outer) = value
    elif val_any_pupil is not None:
        (radius, radius_inner, radius_outer) = val_any_pupil
    else:
        (radius, radius_inner, radius_outer) = (0., 0., 0.)

    ref_model.close()

    return (radius, radius_inner, radius_outer)
