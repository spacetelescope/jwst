import logging

from stdatamodels.jwst import datamodels

from jwst.pathloss import pathloss
from jwst.stpipe import Step

__all__ = ["PathLossStep"]

log = logging.getLogger(__name__)


class PathLossStep(Step):
    """
    Apply the path loss correction to a science exposure.

    Pathloss depends on the centering of the source in the aperture if the
    source is a point source.
    """

    class_alias = "pathloss"

    spec = """
        inverse = boolean(default=False)    # Invert the operation
        source_type = string(default=None)  # Process as specified source type
        user_slit_loc = float(default=None)   # User-provided correction to MIRI LRS source location
    """  # noqa: E501

    reference_file_types = ["pathloss"]

    def process(self, input_data):
        """
        Execute the pathloss step.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.JwstDataModel`
            The input file name or datamodel to be corrected.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
            The updated datamodel..
        """
        # Open the input data model
        output_model = self.prepare_output(input_data)

        if self.use_correction_pars:
            correction_pars = self.correction_pars
            pathloss_model = None
        else:
            correction_pars = None

            # Get the name of the pathloss reference file to use
            pathloss_name = self.get_reference_file(output_model, "pathloss")
            log.info(f"Using PATHLOSS reference file {pathloss_name}")

            # Check for a valid reference file
            if pathloss_name == "N/A":
                log.warning("No PATHLOSS reference file found")
                log.warning("Pathloss step will be skipped")
                output_model.meta.cal_step.pathloss = "SKIPPED"
                return output_model

            # Open the pathloss ref file data model
            if output_model.meta.exposure.type.upper() in ["MIR_LRS-FIXEDSLIT", "MIR_WFSS"]:
                pathloss_model = datamodels.MirLrsPathlossModel(pathloss_name)
            else:
                pathloss_model = datamodels.PathlossModel(pathloss_name)

        # Do the pathloss correction
        output_model, self.correction_pars = pathloss.do_correction(
            output_model,
            pathloss_model,
            inverse=self.inverse,
            source_type=self.source_type,
            correction_pars=correction_pars,
            user_slit_loc=self.user_slit_loc,
        )

        if pathloss_model:
            pathloss_model.close()

        return output_model
