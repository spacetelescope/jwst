"""Get integrated flux as a function of time for a multi-integration spectroscopic observation."""

import logging

import numpy as np
from astropy.table import Table
from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from jwst.white_light.white_light import white_light

log = logging.getLogger(__name__)

__all__ = ["WhiteLightStep"]


class WhiteLightStep(Step):
    """Sum the spectroscopic flux over all wavelengths in each integration."""

    class_alias = "white_light"

    spec = """
    min_wavelength     = float(default=None)      # Default wavelength minimum for integration
    max_wavelength     = float(default=None)      # Default wavelength maximum for integration
    output_ext         = string(default='.ecsv')  # Output file type
    suffix             = string(default='whtlt')  # Default suffix for output files
    """  # noqa: E501

    reference_file_types = ["wavelengthrange"]

    def process(self, input_data):
        """
        Sum the flux over all wavelengths in each integration.

        Produce an integrated (“white”) flux as a function of time for the target. This
        is to be applied to the _x1dints product in a spectroscopic
        Time-Series Observation (TSO).

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.MultiSpecModel`
            Either the path to the file or the science data model for the sum.

        Returns
        -------
        result : `~astropy.table.table.QTable` or None
            Table containing the integrated flux as a function of time.
            If input is invalid, None is returned.
        """
        # This step will not modify the input, so always skip making a copy
        input_model = self.prepare_output(input_data, make_copy=False)

        # First check for valid input
        if not input_model.hasattr("spec") or len(input_model.spec) == 0:
            log.error("No valid input spectra found.")
            return None

        # Load the wavelength range reference file
        waverange_table = self._get_reference_wavelength_range(input_model)

        # Call the white light curve generation routine
        result = white_light(
            input_model,
            waverange_table=waverange_table,
            min_wave=self.min_wavelength,
            max_wave=self.max_wavelength,
        )

        # Write the output catalog
        if self.save_results:
            output_path = self.make_output_path()
            result.write(output_path, format="ascii.ecsv", overwrite=True)

        return result

    def _get_reference_wavelength_range(self, input_model):
        """
        Get the wavelength range reference file and convert it to an astropy table.

        Parameters
        ----------
        input_model : datamodels.Model
            The input data model from which to extract the wavelength range.

        Returns
        -------
        astropy.table.Table or None
            If the input model is of type NIS_SOSS, returns a table with the wavelength
            range information. Otherwise, returns None.
        """
        if input_model.meta.exposure.type != "NIS_SOSS":
            return None
        wavelengthrange_file = self.get_reference_file(input_model, "wavelengthrange")
        if wavelengthrange_file == "N/A":
            log.warning(
                "No wavelength range reference file found. "
                "The entire wavelength range in the input spectral tables will be used "
                "for all spectral orders."
            )
            return None
        log.info(f"Using wavelength range reference file: {wavelengthrange_file}")
        with datamodels.WavelengthrangeModel(wavelengthrange_file) as f:
            return Table(
                np.array(f.wavelengthrange.instance),
                names=["order", "filter", "min_wave", "max_wave"],
                dtype=[int, str, float, float],
            )
