from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from jwst.residual_fringe import residual_fringe

__all__ = ["ResidualFringeStep"]


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
        input_data : str or IFUImageModel
            Input data to correct.  Must be a MIRI MRS IFU image.

        Returns
        -------
        IFUImageModel
            The corrected datamodel.
        """
        self.transmission_level = 80  # sets the transmission level to use in the regions file
        # 80% is what other steps use.

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
            self.log.error("Number of minimum and maximum wavelengths to ignore are not the same")
            raise ValueError("Number of ignore_region_min does not match ignore_region_max")

        ignore_regions["num"] = min_num

        if min_num > 0:
            self.log.info(f"Ignoring {min_num} wavelength regions")

        self.ignore_regions = ignore_regions

        input_data = datamodels.open(input_data)

        if isinstance(input_data, datamodels.IFUImageModel):
            exptype = input_data.meta.exposure.type
        else:
            raise TypeError(f"Failed to process input type: {type(input_data)}")

        # Set up residual fringe correction parameters
        pars = {
            "transmission_level": self.transmission_level,
            "save_intermediate_results": self.save_intermediate_results,
            "make_output_path": self.make_output_path,
        }

        if exptype != "MIR_MRS":
            self.log.warning("Residual fringe correction is only for MIRI MRS data")
            self.log.warning(f"Input is: {exptype}")
            input_data.meta.cal_step.residual_fringe = "SKIPPED"
            return input_data

        # 1. set up the reference files
        # 2. correct the  model
        # 3. return from step

        self.residual_fringe_filename = self.get_reference_file(input_data, "fringefreq")
        self.log.info(f"Using FRINGEFREQ reference file:{self.residual_fringe_filename}")

        # set up regions reference file
        self.regions_filename = self.get_reference_file(input_data, "regions")
        self.log.info(f"Using MRS regions reference file: {self.regions_filename}")

        # Check for a valid reference files. If they are not found skip step
        if self.residual_fringe_filename == "N/A" or self.regions_filename == "N/A":
            if self.residual_fringe_filename == "N/A":
                self.log.warning("No FRINGEFREQ reference file found")
                self.log.warning("Residual Fringe step will be skipped")

            if self.regions_filename == "N/A":
                self.log.warning("No MRS regions reference file found")
                self.log.warning("Residual Fringe step will be skipped")

            input_data.meta.cal_step.residual_fringe = "SKIPPED"
            return input_data

        # Do the correction
        rfc = residual_fringe.ResidualFringeCorrection(
            input_data,
            self.residual_fringe_filename,
            self.regions_filename,
            self.ignore_regions,
            **pars,
        )
        result = rfc.do_correction()
        result.meta.cal_step.residual_fringe = "COMPLETE"
        return result
