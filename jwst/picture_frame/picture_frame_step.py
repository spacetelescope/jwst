import logging

from stdatamodels.jwst import datamodels

from jwst.picture_frame.picture_frame import correct_picture_frame
from jwst.stpipe import Step

__all__ = ["PictureFrameStep"]

log = logging.getLogger(__name__)


class PictureFrameStep(Step):
    """Perform picture frame correction."""

    class_alias = "picture_frame"

    spec = """
        n_sigma = float(default=2.0)  # Clipping level for non-background signal.
        save_mask = boolean(default=False)  # Save the created background mask
        save_correction = boolean(default=False)  # Save the thermal correction data
        skip = boolean(default=True)  # By default, skip the step.
    """  # noqa: E501

    reference_file_types = ["pictureframe"]

    def process(self, input_data):
        """
        Scale and subtract the thermal "picture frame" effect from a ramp data set.

        First, a draft rate image is created if needed from the input ramp data.
        A background mask is created from the rate image by masking any known science
        regions, then iteratively sigma-clipping.  Then, the median levels are computed
        from background data in a central region and from the edges of each group image.
        These levels are used to scale and offset the correction image in the
        ``pictureframe`` reference file, then the scaled image is subtracted from the
        input data in each group.

        Input data is expected to be a NIRSpec FULL frame ramp file (RampModel),
        in between jump and ramp fitting steps (before 1/f noise cleaning),
        or a rate file (ImageModel or CubeModel).

        Alternately, a rate image or rateints cube may be provided.  In that case,
        no draft rate image is necessary, and the correction is directly applied
        to the rate image.

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.RampModel` \
                     or `~stdatamodels.jwst.datamodels.ImageModel` \
                     or `~stdatamodels.jwst.datamodels.CubeModel`
            Filename or input datamodel to be corrected. Must be NIRSpec full-frame.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.RampModel` \
                       or `~stdatamodels.jwst.datamodels.ImageModel` \
                       or `~stdatamodels.jwst.datamodels.CubeModel`
            The corrected datamodel, matching the type of the input.
        """
        # Open the input data model
        # Note: this step can be run either in stage 1 on ramp-type
        # data or in stage 2 on rate-type data, so it's important
        # not to provide a specific datamodel to prepare_output.
        output_model = self.prepare_output(input_data)

        pictureframe_file = self.get_reference_file(output_model, "pictureframe")
        exp_type = output_model.meta.exposure.type
        if pictureframe_file == "N/A":
            log.warning(f"Picture frame correction is not available for exposure type {exp_type}.")
            status = "SKIPPED"
            mask_model = None
            correction_model = None
        else:
            log.info(f"Using PICTUREFRAME reference file: {pictureframe_file}")
            pictureframe_model = datamodels.PictureFrameModel(pictureframe_file)

            # Correct for the picture frame effect.
            # The output model is updated in place.
            output_model, mask_model, correction_model, status = correct_picture_frame(
                output_model,
                pictureframe_model,
                n_sigma=self.n_sigma,
                input_dir=self.input_dir,
                save_mask=self.save_mask,
                save_correction=self.save_correction,
            )

        # Save the mask, if requested
        if self.save_mask and mask_model is not None:
            mask_path = self.make_output_path(
                basepath=output_model.meta.filename, suffix="pctfrm_mask"
            )
            log.info(f"Saving mask file {mask_path}")
            mask_model.save(mask_path)
            mask_model.close()

        # Save the correction, if requested
        if self.save_correction and correction_model is not None:
            correction_path = self.make_output_path(
                basepath=output_model.meta.filename, suffix="pctfrm_correction"
            )
            log.info(f"Saving correction file {correction_path}")
            correction_model.save(correction_path)
            correction_model.close()

        # Set the step completion status
        output_model.meta.cal_step.picture_frame = status

        return output_model
