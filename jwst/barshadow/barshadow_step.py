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
        input_data : DataModel
            Input JWST datamodel object.

        Returns
        -------
        result : DataModel
            JWST datamodel object with corrected science data and barshadow
            extension(s) added.
        """
        # Open the input data model
        with datamodels.open(input_data) as input_model:
            exp_type = input_model.meta.exposure.type
            if exp_type == "NRS_MSASPEC":
                if self.use_correction_pars:
                    correction_pars = self.correction_pars
                    barshadow_model = None

                else:
                    correction_pars = None

                    # Get the name of the bar shadow reference file to use
                    self.barshadow_name = self.get_reference_file(input_model, "barshadow")
                    log.info(f"Using BARSHADOW reference file {self.barshadow_name}")

                    # Check for a valid reference file
                    if self.barshadow_name == "N/A":
                        log.warning("No BARSHADOW reference file found")
                        log.warning("Bar shadow step will be skipped")
                        result = input_model.copy()
                        result.meta.cal_step.barshadow = "SKIPPED"
                        return result

                    # Open the barshadow ref file data model
                    barshadow_model = datamodels.BarshadowModel(self.barshadow_name)

                # Do the bar shadow correction
                result, self.correction_pars = bar_shadow.do_correction(
                    input_model,
                    barshadow_model,
                    inverse=self.inverse,
                    source_type=self.source_type,
                    correction_pars=correction_pars,
                )

                if barshadow_model:
                    barshadow_model.close()
                result.meta.cal_step.barshadow = "COMPLETE"
            else:
                input_model.meta.cal_step.barshadow = "SKIPPED"
                result = input_model
        return result
