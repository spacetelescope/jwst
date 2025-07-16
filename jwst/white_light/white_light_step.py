"""Get integrated flux as a function of time for a multi-integration spectroscopic observation."""

from astropy.table import Table
import numpy as np
from stdatamodels.jwst import datamodels

from jwst.stpipe import Step

from .white_light import white_light

__all__ = ["WhiteLightStep"]


class WhiteLightStep(Step):
    """
    Sum the spectroscopic flux over all wavelengths in each integration.

    Produce an integrated (“white”) flux as a function of time for the target. This
    is to be applied to the _x1dints product in a spectroscopic
    Time-Series Observation (TSO).
    """

    class_alias = "white_light"

    spec = """
    min_wavelength     = float(default=None)      # Default wavelength minimum for integration
    max_wavelength     = float(default=None)      # Default wavelength maximum for integration
    output_ext         = string(default='.ecsv')  # Output file type
    suffix             = string(default='whtlt')  # Default suffix for output files
    """  # noqa: E501

    def process(self, step_input):
        """
        Sum the flux over all wavelengths in each integration as a function of time for the target.

        Parameters
        ----------
        step_input : str or MultiSpecModel
            Either the path to the file or the science data model for the sum.

        Returns
        -------
        result : astropy.table.table.QTable
            Table containing the integrated flux as a function of time.
        """
        # load the wavelength range reference file

        # Load the input
        with datamodels.open(step_input) as input_model:
            wr = self._get_reference_wavelength_range(input_model)
            # Call the white light curve generation routine
            result = white_light(
                input_model, wr=wr, min_wave=self.min_wavelength, max_wave=self.max_wavelength
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
        with datamodels.WavelengthrangeModel(wavelengthrange_file) as f:
            return Table(
                np.array(f.wavelengthrange.instance),
                names=["order", "filter", "min_wave", "max_wave"],
                dtype=[int, str, float, float],
            )
