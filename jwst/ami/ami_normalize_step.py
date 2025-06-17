from stdatamodels.jwst import datamodels

from jwst.stpipe import Step

from . import ami_normalize

__all__ = ["AmiNormalizeStep"]


class AmiNormalizeStep(Step):
    """Normalize target LG results using reference LG results."""

    class_alias = "ami_normalize"

    spec = """
    suffix = string(default='aminorm-oi')
    """  # noqa: E501

    def process(self, target, reference):
        """
        Normalize the LG results for science target, using the LG results for reference target.

        Parameters
        ----------
        target : str or model
            Target input
        reference : str or model
            Reference input

        Returns
        -------
        result : AmiOIModel object
            AMI data model that's been normalized
        """
        # Open the target and reference input models
        target_model = datamodels.AmiOIModel(target)
        reference_model = datamodels.AmiOIModel(reference)

        # Call the normalization routine
        result = ami_normalize.normalize_lg(target_model, reference_model)

        result.meta.cal_step.ami_normalize = "COMPLETE"

        # Close the input models
        target_model.close()
        reference_model.close()

        # We're done
        return result
