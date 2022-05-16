"""
Module for parsing APCORR and ABVEGAOFFSET reference file data.
"""

import logging

from astropy.utils import lazyproperty
import numpy as np

from .. import datamodels
from ..datamodels import ImageModel, ABVegaOffsetModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ReferenceData:
    """
    Class for APCORR and ABVEGAOFFSET reference file data needed by
    `SourceCatalogStep`.

    Parameters
    ----------
    model : `ImageModel`
        An `ImageModel` of drizzled image.

    reffile_paths : list of 2 str
        The full path filename of the APCORR and ABVEGAOFFSET reference
        files.

    aperture_ee : tuple of 3 int
        The aperture encircled energies to be used for aperture
        photometry.  The values must be 3 strictly-increasing integers.
        Valid values are defined in the APCORR reference files (20, 30,
        40, 50, 60, 70, or 80).

    Attributes
    ----------
    aperture_params : `dict`
        A dictionary containing the aperture parameters (radii, aperture
        corrections, and background annulus inner and outer radii).

    abvega_offset : float
        Offset to convert from AB to Vega magnitudes.  The value
        represents m_AB - m_Vega.
    """

    def __init__(self, model, reffile_paths, aperture_ee):
        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')
        self.model = model

        self.aperture_ee = self._validate_aperture_ee(aperture_ee)
        self.apcorr_filename = reffile_paths[0]
        self.abvegaoffset_filename = reffile_paths[1]

        self.instrument = self.model.meta.instrument.name
        self.detector = self.model.meta.instrument.detector
        self.filtername = self.model.meta.instrument.filter
        self.pupil = model.meta.instrument.pupil
        self.subarray = self.model.meta.subarray.name

        log.info(f'Instrument: {self.instrument}')
        if self.detector is not None:
            log.info(f'Detector: {self.detector}')
        if self.filtername is not None:
            log.info(f'Filter: {self.filtername}')
        if self.pupil is not None:
            log.info(f'Pupil: {self.pupil}')
        if self.subarray is not None:
            log.info(f'Subarray: {self.subarray}')

    @staticmethod
    def _validate_aperture_ee(aperture_ee):
        """
        Validate the input ``aperture_ee``.
        """
        aperture_ee = np.array(aperture_ee).astype(int)
        if not np.all(aperture_ee[1:] > aperture_ee[:-1]):
            raise ValueError('aperture_ee values must be strictly '
                             'increasing')
        if len(aperture_ee) != 3:
            raise ValueError('aperture_ee must contain only 3 values')
        if np.any(np.logical_or(aperture_ee <= 0, aperture_ee >= 100)):
            raise ValueError('aperture_ee values must be between 0 and 100')
        return aperture_ee

    @lazyproperty
    def _aperture_ee_table(self):
        """
        Get the encircled energy table for the given instrument
        configuration.
        """
        if self.instrument in ('NIRCAM', 'NIRISS'):
            selector = {'filter': self.filtername, 'pupil': self.pupil}
        elif self.instrument == 'MIRI':
            selector = {'filter': self.filtername, 'subarray': self.subarray}
        elif self.instrument == 'FGS':
            selector = None
        else:
            raise RuntimeError(f'{self.instrument} is not a valid instrument')

        apcorr_model = datamodels.open(self.apcorr_filename)
        apcorr = apcorr_model.apcorr_table
        if selector is None:  # FGS
            ee_table = apcorr
        else:
            mask_idx = [apcorr[key] == value
                        for key, value in selector.items()]
            ee_table = apcorr[np.logical_and.reduce(mask_idx)]

        if len(ee_table) == 0:
            raise RuntimeError('APCORR reference file data is missing for '
                               f'{selector}.')

        return ee_table

    def _get_ee_table_row(self, aperture_ee):
        """
        Get the encircled energy row for the input ``aperture_ee``.
        """
        ee_percent = np.round(self._aperture_ee_table['eefraction'] * 100)
        row_mask = (ee_percent == aperture_ee)
        ee_row = self._aperture_ee_table[row_mask]
        if len(ee_row) == 0:
            raise RuntimeError('Aperture encircled energy value of '
                               f'{aperture_ee} appears to be invalid. No '
                               'matching row was found in the APCORR '
                               'reference file {self.apcorr_filename}')
        if len(ee_row) > 1:
            raise RuntimeError('More than one matching row was found in '
                               'the APCORR reference file '
                               f'{self.apcorr_filename}')
        return ee_row

    @lazyproperty
    def aperture_params(self):
        """
        A dictionary containing the aperture parameters (radii, aperture
        corrections, and background annulus inner and outer radii).
        """
        if self.apcorr_filename is None:
            log.warning('APCorrModel reference file was not input. Using '
                        'fallback aperture sizes without any aperture '
                        'corrections.')

            params = {'aperture_radii': np.array((1.0, 2.0, 3.0)),
                      'aperture_corrections': np.array((1.0, 1.0, 1.0)),
                      'aperture_ee': np.array((1, 2, 3)),
                      'bkg_aperture_inner_radius': 5.0,
                      'bkg_aperture_outer_radius': 10.0}
            return params

        params = {}
        radii = []
        apcorrs = []
        skyins = []
        skyouts = []
        for aper_ee in self.aperture_ee:
            row = self._get_ee_table_row(aper_ee)
            radii.append(row['radius'][0])
            apcorrs.append(row['apcorr'][0])
            skyins.append(row['skyin'][0])
            skyouts.append(row['skyout'][0])

        if self.model.meta.resample.pixel_scale_ratio is not None:
            # pixel_scale_ratio is the ratio of the resampled to the native
            # pixel scale (values < 1 have smaller resampled pixels)
            pixel_scale_ratio = self.model.meta.resample.pixel_scale_ratio
        else:
            log.warning('model.meta.resample.pixel_scale_ratio was not '
                        'found. Assuming the native detector pixel scale '
                        '(i.e., pixel_scale_ratio = 1)')
            pixel_scale_ratio = 1.0

        params['aperture_ee'] = self.aperture_ee
        params['aperture_radii'] = np.array(radii) / pixel_scale_ratio
        params['aperture_corrections'] = np.array(apcorrs)

        skyins = np.unique(skyins)
        skyouts = np.unique(skyouts)
        if len(skyins) != 1 or len(skyouts) != 1:
            raise RuntimeError('Expected to find only one value for skyin '
                               'and skyout in the APCORR reference file for '
                               'a given selector.')
        params['bkg_aperture_inner_radius'] = skyins[0] / pixel_scale_ratio
        params['bkg_aperture_outer_radius'] = skyouts[0] / pixel_scale_ratio

        return params

    @lazyproperty
    def abvega_offset(self):
        """
        Offset to convert from AB to Vega magnitudes.

        The value represents m_AB - m_Vega.
        """
        if self.abvegaoffset_filename is None:
            log.warning('ABVEGAOFFSET reference file was not input. '
                        'Catalog Vega magnitudes are not correct.')
            return 0.0

        if self.instrument in ('NIRCAM', 'NIRISS'):
            selector = {'filter': self.filtername, 'pupil': self.pupil}
        elif self.instrument == 'MIRI':
            selector = {'filter': self.filtername}
        elif self.instrument == 'FGS':
            selector = {'detector': self.detector}
        else:
            raise RuntimeError(f'{self.instrument} is not a valid instrument')

        abvegaoffset_model = ABVegaOffsetModel(self.abvegaoffset_filename)
        offsets_table = abvegaoffset_model.abvega_offset

        try:
            mask_idx = [offsets_table[key] == value
                        for key, value in selector.items()]
        except KeyError as badkey:
            raise KeyError(f'{badkey} not found in ABVEGAOFFSET reference '
                           f'file {self.abvegaoffset_filename}')

        row = offsets_table[np.logical_and.reduce(mask_idx)]

        if len(row) == 0:
            raise RuntimeError('Did not find matching row in ABVEGAOFFSET '
                               f'reference file {self.abvegaoffset_filename}')
        if len(row) > 1:
            raise RuntimeError('Found more than one matching row in '
                               'ABVEGAOFFSET reference file '
                               f'{self.abvegaoffset_filename}')

        abvega_offset = row['abvega_offset'][0]
        log.info(f'AB to Vega magnitude offset {abvega_offset:.5f}')
        abvegaoffset_model.close()
        return abvega_offset
