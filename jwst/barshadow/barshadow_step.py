"""Barshadow correction step."""

import logging

from stdatamodels.jwst import datamodels

from jwst.barshadow import bar_shadow
from jwst.stpipe import Step

__all__ = ["BarShadowStep"]

log = logging.getLogger(__name__)


class BarShadowStep(Step):
    """
    Correct NIRSpec MOS data for barshadow effects.

    Bar shadow correction depends on the position of a pixel along the slit
    and the wavelength. It is only applied to uniform sources and only for
    NRS MSA data.
    """

    class_alias = "barshadow"

    spec = """
        inverse = boolean(default=False)    # Invert the operation
        source_type = string(default=None)  # Process as specified source type.
    """  # noqa: E501

    reference_file_types = ["barshadow"]

    def process(self, input_data):
        """
        Perform the barshadow correction step.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.MultiSlitModel`
            Input JWST file name or datamodel object.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.MultiSlitModel`
            JWST datamodel object with corrected science data and barshadow
            extension(s) added.
        """
        # Open the input data model
        output_model = self.prepare_output(input_data)

        # Check for MOS data
        exp_type = output_model.meta.exposure.type
        if exp_type != "NRS_MSASPEC":
            output_model.meta.cal_step.barshadow = "SKIPPED"
            return output_model

        if self.use_correction_pars:
            correction_pars = self.correction_pars
            barshadow_model = None
        else:
            correction_pars = None

            # Get the name of the bar shadow reference file to use
            barshadow_name = self.get_reference_file(output_model, "barshadow")
            log.info(f"Using BARSHADOW reference file {barshadow_name}")

            # Check for a valid reference file
            if barshadow_name == "N/A":
                log.warning("No BARSHADOW reference file found")
                log.warning("Bar shadow step will be skipped")
                output_model.meta.cal_step.barshadow = "SKIPPED"
                return output_model

            # Open the barshadow ref file data model
            barshadow_model = datamodels.BarshadowModel(barshadow_name)

        # Do the bar shadow correction
        output_model, self.correction_pars = bar_shadow.do_correction(
            output_model,
            barshadow_model,
            inverse=self.inverse,
            source_type=self.source_type,
            correction_pars=correction_pars,
        )

        if barshadow_model:
            barshadow_model.close()
        output_model.meta.cal_step.barshadow = "COMPLETE"

        return output_model
