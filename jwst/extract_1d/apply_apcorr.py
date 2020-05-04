from typing import Tuple
from scipy.interpolate import interp2d, interp1d
from astropy.io import fits

from ..assign_wcs.util import compute_scale
from ..datamodels import DataModel


class ApCorr:
    """Perform aperture correction on input extraction data.

    """
    def __init__(self, input_model: DataModel, appcorr_table: fits.FITS_rec, slitname: str = None,
                 location: Tuple[float, float] = None):
        self.correction = None

        self.model = input_model
        self._reference_table = appcorr_table
        self.slitname = slitname
        self.location = location

        self.match_keys = self._get_match_keys()
        self.match_pars = self._get_match_pars()

        if self.slitname and self.slitname != self.model.meta.exposure.type:
            self.match_pars['slit'] = self.slitname

        self.reference = self._reduce_reftable()

        if 'size' in self.reference and self.reference.columns['size'].unit.startswith('arcsec'):
            if self.location is not None:
                self.reference['size'][0] /= compute_scale(self.model.wcs, location)
            else:
                raise ValueError(
                    f'If the size column for the input APCORR reference file are in units with arcseconds, a location '
                    f'(RA, DEC) must be provided in order to convert to pixels.'
                )

        self.apcorr_func = self._approx_apcorr_fn()

    def _get_match_keys(self):
        """Get column keys needed for reducing the reference table."""
        match_pars = {
            'MIRI': {
                'LRS': {'subarray': ['name']},
                'MRS': {'instrument': []},
            },
            'NIRSPEC': {
                'IFU': {'instrument': ['filter', 'grating']},
                'MSASPEC': {'instrument': ['filter', 'grating']},
                'FIXEDSLIT': {'instrument': ['filter', 'grating', 'slit']},
            },
            'NIRCAM': {
                'WFSS': {'insturment': ['filter', 'pupil']}
            },
            'NIRISS': {
                'WFSS': {'instrument': ['filter', 'pupil']}
            }
        }

        instrument = self.model.meta.instrument.name.upper()
        exptype = self.model.meta.exposure.type.upper()

        relevant_pars = match_pars[instrument]

        for key in relevant_pars.keys():
            if key in exptype:
                return relevant_pars[key]

    def _get_match_pars(self):
        match_pars = {}

        for node, keys in self.match_keys.items():
            meta_node = self.model.meta[node]

            for key in keys:
                match_pars[key if key != 'name' else node] = getattr(meta_node, key)

        return match_pars

    def _reduce_reftable(self):
        """Reduce full reference table to a matched row."""
        table = self._reference_table.copy()

        for key, value in self.match_pars.items():
            table = table[table[key] == value]

        return table

    def _approx_apcorr_fn(self):
        """Generate an approximate function for interpolating apcorr values to input wavelength and size."""
        wavelength = self.reference['wavelength'][0][:self.reference['nelem_wl'][0]]

        if self.model.meta.exposure.type == 'MIR_MRS':
            size = self.reference['radius'][0][:self.reference['nelem_wl'][0]]
            apcorr = self.reference['apcorr'][0][:self.reference['nelem_wl'][0]]

            return interp1d(size, wavelength, apcorr)

        size = self.reference['size'][0][:self.reference['nelem_size'][0]]
        apcorr = self.reference['apcorr'][0][:self.reference['nelem_size'][0], :self.reference['nelem_wl'][0]]

        return interp2d(size, wavelength, apcorr)

    def apply_apcorr(self, spec_table: fits.FITS_rec):
        """Apply interpolated aperture correction value to source-related extraction results.

        Parameters
        ----------
        spec_table : `~fits.FITS_rec`
            Table of aperture corrections values from apcorr reference file.

        """
        cols_to_correct = ('flux', 'surf_bright', 'error', 'sb_error')

        for row in spec_table:
            correction = self.apcorr_func(row['npixels'], row['wavelength'])

            for col in cols_to_correct:
                row[col] *= correction
