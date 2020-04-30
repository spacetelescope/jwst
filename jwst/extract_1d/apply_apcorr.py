from typing import Tuple
from scipy.interpolate import interp2d
from astropy.io import fits

from ..assign_wcs.util import compute_scale
from ..datamodels import DataModel, ReferenceFileModel


class ApCorr:
    """Perform aperture correction on input extraction data.

    """
    def __init__(self, input_model: DataModel, appcorr_table: fits.FITS_rec, slitname: str = None,
                 location: Tuple[float, float] = None):
        self.reference = None
        self.correction = None

        self.model = input_model
        self._reference_table = appcorr_table
        self.slitname = slitname
        self.location = location
        self.match_keys = self._get_match_keys()

        self.match_pars = {
            key: getattr(input_model.meta, key) for key in self.match_keys if hasattr(input_model.meta, key)
        }

        if self.slitname:
            self.match_pars['slit'] = self.slitname

        self._reduce_reftable()

        if self.reference['size'].unit.startswith('arcsec'):
            if self.location is not None:
                self.reference['size'] /= compute_scale(self.model.wcs, location)
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
                'LRS': ['subarray'],
                'MRS': ['radius', 'axis_ratio', 'axis_pa']
            },
            'NIRSPEC': {
                'IFU/MOS': ['filter', 'grating', 'pixphase'],
                'FS': ['filter', 'grating', 'slit', 'pixphase'],
            },
            'NIRCAM': {
                'WFSS': ['filter', 'pupil']
            },
            'NIRISS': {
                'WFSS': ['filter', 'pupil']
            }
        }

        instrument = self.model.meta.instrument.name.upper()
        mode = self.model.meta.mode.name.upper()

        relevant_pars = match_pars[instrument]

        if mode in ('IFU', 'MOS'):
            return relevant_pars['IFU/MOS']

        return relevant_pars[mode]

    def _reduce_reftable(self):
        """Reduce full reference table to a matched row."""
        table = self._reference_table.copy()

        for key, value in self.match_pars:
            table = table[table['key'] == value]

        self.reference = table

    def _approx_apcorr_fn(self):
        """Generate an approximate function for interpolating apcorr values to input wavelength and size."""
        wavelength = self.reference['wavelength'][0][:self.reference['nelem_wl']]
        size = self.reference['size'][0][:self.reference['nelem_size']]

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
            correction = self.apcorr_func(row['wavelength'], row['npixels'])

            for col in cols_to_correct:
                row[col] *= correction


# def inperpolate_apcorr(wavelength, apcorr_wl, size, size_units, **compute_scale_kwargs):
#
#     if size_units.startswith('arcsec'):
#         # This assumes that the image scale is the same in both the dispersion
#         # and cross-dispersion directions, which is probably not valid.
#         # This is the pixel scale in arcseconds per pixel.
#         # pixel_scale = compute_scale(**compute_scale_kwargs)
#         # size /= pixel_scale                     # convert to pixels
#
#     wavelength[np.isnan(wavelength)] = -1.
#     extr_size = np.interp(wavelength, apcorr_wl, size,
#                           left=np.NaN, right=np.NaN)
#     no_cal = np.isnan(extr_size)
