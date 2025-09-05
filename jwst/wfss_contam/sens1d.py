import logging

import numpy as np

from jwst.photom.photom import find_row

log = logging.getLogger(__name__)

__all__ = ["get_photom_data", "create_1d_sens"]


def get_photom_data(phot_model, filter_name, pupil, order):
    """
    Retrieve wavelength and response data from photom ref file.

    Wavelengths from the reference file are expected to be in units of microns.

    Parameters
    ----------
    phot_model : `jwst.datamodels.NrcWfssPhotomModel` or `jwst.datamodels.NisWfssPhotomModel`
        Photom ref file data model
    filter_name : str
        Filter value
    pupil : str
        Pupil value
    order : int
        Spectral order number

    Returns
    -------
    ref_waves : float array
        Wavelengths from the ref file.
    relresps : float array
        Wavelength-dependent response (flux calibration) values from the ref file,
        same shape as `ref_waves`.
    """
    # Get the appropriate row of data from the reference table
    phot_table = phot_model.phot_table
    fields_to_match = {"filter": filter_name, "pupil": pupil, "order": order}
    row = find_row(phot_table, fields_to_match)
    tabdata = phot_table[row]

    # Scalar conversion factor
    scalar_conversion = tabdata["photmjsr"]  # unit is MJy / sr

    # Get the length of the relative response arrays in this row
    nelem = tabdata["nelem"]

    # Load the wavelength and relative response arrays
    ref_waves = tabdata["wavelength"][:nelem]
    relresps = scalar_conversion * tabdata["relresponse"][:nelem]

    # Make sure waves and relresps are in increasing wavelength order
    if not np.all(np.diff(ref_waves) > 0):
        index = np.argsort(ref_waves)
        ref_waves = ref_waves[index].copy()
        relresps = relresps[index].copy()

    return ref_waves, relresps


def create_1d_sens(data_waves, ref_waves, relresps):
    """
    Find photometric conversion values based on per-pixel wavelength-dependent response.

    Parameters
    ----------
    data_waves : float array
        Input data wavelength values
    ref_waves : float array
        1D wavelength vector on which relative response values are sampled
    relresps : float array
        1D photometric response values, as a function of wavelength

    Returns
    -------
    sens_1d : float array
        1D array of computed photometric conversion values
    no_cal : int array
        1D mask indicating where no conversion is available
    """
    # Interpolate the photometric response values onto the
    # 1D wavelength grid of the data
    sens_1d = np.interp(data_waves, ref_waves, relresps, left=np.nan, right=np.nan)

    sens_1d[np.isnan(sens_1d)] = 0
    no_cal = sens_1d == 0

    return sens_1d, no_cal
