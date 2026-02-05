from jwst.coron import stack_refs
from jwst.stpipe import Step

__all__ = ["StackRefsStep"]


class StackRefsStep(Step):
    """Stack multiple PSF reference exposures into a single CubeModel."""

    class_alias = "stack_refs"

    spec = """
    """  # noqa: E501

    def process(self, input_files):
        """
        Execute the StackRefs calibration step.

        Parameters
        ----------
        input_files : str, `~jwst.datamodels.container.ModelContainer`, or \
                      list of `~stdatamodels.jwst.datamodels.CubeModel`
            Association file or ModelContainer containing input science exposures.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.CubeModel`
            PSF reference exposures stacked into a CubeModel
        """
        # Open the inputs
        input_models = self.prepare_output(input_files)

        # Call the stacking routine
        output_model = stack_refs.make_cube(input_models)
        output_model.meta.cal_step.stack_psfs = "COMPLETE"

        # Output is a new model: close the input if it was opened here
        if input_models is not input_files:
            input_models.close()

        return output_model
