import logging

from stdatamodels.jwst import datamodels

from jwst.residual_fringe import residual_fringe
from jwst.stpipe import Step

__all__ = ["ResidualFringeStep"]

log = logging.getLogger(__name__)


class ResidualFringeStep(Step):
    """
    Apply residual fringe correction to a MIRI MRS image.

    Requires frequency parameters provided in the FRINGEFREQ reference file.
    """

    class_alias = "residual_fringe"

    spec = """
        skip = boolean(default=True)
        save_intermediate_results  = boolean(default = False)
        search_output_file = boolean(default = False)
        ignore_region_min = list(default = None)
        ignore_region_max = list(default = None)
        suffix = string(default = 'residual_fringe')
    """  # noqa: E501

    reference_file_types = ["fringefreq", "regions"]

    def process(self, input_data):
        """
        Perform the residual fringe correction.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.IFUImageModel`
            Input data to correct.  Must be a MIRI MRS IFU image.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.IFUImageModel`
            The corrected datamodel.
        """
        # Sets the transmission level to use in the regions file.
        # 80% is standard.
        transmission_level = 80

        # set up the dictionary to ignore wavelength regions in the residual fringe correction
        ignore_regions = {"num": 0, "min": [], "max": []}
        if self.ignore_region_min is not None:
            for region in self.ignore_region_min:
                ignore_regions["min"].append(float(region))
        min_num = len(ignore_regions["min"])

        if self.ignore_region_max is not None:
            for region in self.ignore_region_max:
                ignore_regions["max"].append(float(region))
        max_num = len(ignore_regions["max"])

        if max_num != min_num:
            log.error("Number of minimum and maximum wavelengths to ignore are not the same")
            raise ValueError("Number of ignore_region_min does not match ignore_region_max")

        ignore_regions["num"] = min_num
        if min_num > 0:
            log.info(f"Ignoring {min_num} wavelength regions")

        output_model = self.prepare_output(input_data)

        if isinstance(output_model, datamodels.IFUImageModel):
            exptype = output_model.meta.exposure.type
        else:
            raise TypeError(f"Failed to process input type: {type(output_model)}")

        # Set up residual fringe correction parameters
        pars = {
            "transmission_level": transmission_level,
            "save_intermediate_results": self.save_intermediate_results,
            "make_output_path": self.make_output_path,
        }

        if exptype != "MIR_MRS":
            log.warning("Residual fringe correction is only for MIRI MRS data")
            log.warning(f"Input is: {exptype}")
            output_model.meta.cal_step.residual_fringe = "SKIPPED"
            return output_model

        # 1. set up the reference files
        # 2. correct the  model
        # 3. return from step

        residual_fringe_filename = self.get_reference_file(output_model, "fringefreq")
        log.info(f"Using FRINGEFREQ reference file: {residual_fringe_filename}")

        # set up regions reference file
        regions_filename = self.get_reference_file(output_model, "regions")
        log.info(f"Using MRS regions reference file: {regions_filename}")

        # Check for a valid reference files. If they are not found skip step
        if residual_fringe_filename == "N/A" or regions_filename == "N/A":
            if residual_fringe_filename == "N/A":
                log.warning("No FRINGEFREQ reference file found")
                log.warning("Residual Fringe step will be skipped")

            if regions_filename == "N/A":
                log.warning("No MRS regions reference file found")
                log.warning("Residual Fringe step will be skipped")

            output_model.meta.cal_step.residual_fringe = "SKIPPED"
            return output_model

        # Do the correction
        rfc = residual_fringe.ResidualFringeCorrection(
            output_model,
            residual_fringe_filename,
            regions_filename,
            ignore_regions,
            **pars,
        )
        output_model = rfc.do_correction()
        output_model.meta.cal_step.residual_fringe = "COMPLETE"
        return output_model
