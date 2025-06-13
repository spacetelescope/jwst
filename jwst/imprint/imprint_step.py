from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from jwst.background import subtract_images

__all__ = ["ImprintStep"]


class ImprintStep(Step):
    """
    Remove NIRSpec MSA imprint structure from an exposure.

    The imprint structure is removed by subtracting a separate
    imprint (a.k.a. leakcal) exposure.
    """

    class_alias = "imprint"

    spec = """
    """  # noqa: E501

    def process(self, input_data, imprint):
        """
        Subtract an imprint image from the input data.

        If a single imprint image is provided, it is directly subtracted
        without further checks.

        If multiple imprint images are provided, the background target
        flag (`meta.dither.observation.bkgdtarg`) is checked for a match to
        the input data.  If there is a single imprint image matching the
        input data's background flag, then it is directly subtracted without
        further checks.

        If there are multiple imprint images that match the input data's
        background flag, then the imprint is checked for an
        observation ID (`meta.observation.observation_number`) and
        dither position index (`meta.dither.position_number`).  If there is
        a match, then the matching imprint image is subtracted from the
        input data.  If there is no match, then the step is skipped for the
        input data.

        Parameters
        ----------
        input_data : DataModel or str
            Input exposure to be corrected.
        imprint : list of str or DataModel
            Imprint exposures associated with the input.

        Returns
        -------
        DataModel
            The imprint subtracted exposure.
        """
        # Open the input science image and get its dither pattern position number
        input_model = datamodels.open(input_data)
        pos_no = input_model.meta.dither.position_number
        obs_no = input_model.meta.observation.observation_number
        is_bkgd = input_model.meta.observation.bkgdtarg

        # Check all the imprints to see if there is a direct match for the input
        imprint_models = []
        match_model = None
        for imprint_exp in imprint:
            imprint_model = datamodels.open(imprint_exp)
            imprint_pos_no = imprint_model.meta.dither.position_number
            imprint_obs_no = imprint_model.meta.observation.observation_number
            imprint_bkg = imprint_model.meta.observation.bkgdtarg

            if len(imprint) == 1:
                # Exactly one imprint provided: use it.
                match_model = datamodels.open(imprint[0])
                imprint_models.append(match_model)
            elif is_bkgd != imprint_bkg:
                # More than one imprint and imprint does not match input's
                # background status. Don't consider it further.
                imprint_model.close()
            else:
                # Keep the open model for further checks
                imprint_models.append(imprint_model)

                # Check if there is a direct match to the current
                # observation and dither position.
                if pos_no == imprint_pos_no and obs_no == imprint_obs_no:
                    match_model = imprint_model
                    break

        # Check for a single imprint intended for use with all dither positions.
        if match_model is None and len(imprint_models) == 1:
            # No direct match found - use the only one available.
            match_model = imprint_models[0]

        # Copy the input image to the output
        result = input_model.copy()

        if match_model is not None:
            # Subtract the matching imprint image
            self.log.info(
                f"Subtracting imprint image {match_model.meta.filename} "
                f"from {input_model.meta.filename}"
            )
            result = subtract_images.subtract(input_model, match_model)

            # Update the step status and close the imprint model
            result.meta.cal_step.imprint = "COMPLETE"
        else:
            self.log.warning(f"No matching imprint image found for {input_model.meta.filename}")
            self.log.warning("Step will be skipped")
            result.meta.cal_step.imprint = "SKIPPED"

        # Close any open imprint models
        for model in imprint_models:
            model.close()

        # Close the input model
        input_model.close()
        return result
