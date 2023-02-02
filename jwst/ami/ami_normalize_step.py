from stdatamodels.jwst import datamodels

from ..stpipe import Step

from . import ami_normalize

__all__ = ["AmiNormalizeStep"]


class AmiNormalizeStep(Step):
    """
    AmiNormalizeStep: Normalize target LG results using reference LG results
    """

    class_alias = "ami_normalize"

    spec = """
    """

    def process(self, target, reference):
        """
        Normalizes the LG results for a science target, using the LG results
        for a reference target.

        Parameters
        ----------
        target: string or model
            target input

        reference: string or model
            reference input

        Returns
        -------
        result: AmiLgModel object
            AMI data model that's been normalized
        """

        # Open the target and reference input models
        target_model = datamodels.AmiLgModel(target)
        reference_model = datamodels.AmiLgModel(reference)

        # Call the normalization routine
        result = ami_normalize.normalize_LG(target_model, reference_model)

        result.meta.cal_step.ami_normalize = 'COMPLETE'

        # Close the input models
        target_model.close()
        reference_model.close()

        # We're done
        return result
