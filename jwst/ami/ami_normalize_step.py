from stdatamodels.jwst import datamodels

from jwst.ami import ami_normalize
from jwst.stpipe import Step

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
        target : str or `~stdatamodels.jwst.datamodels.AmiOIModel`
            Target input file name or datamodel.
        reference : str or `~stdatamodels.jwst.datamodels.AmiOIModel`
            Reference input file name or datamodel.

        Returns
        -------
        result : `~stdatamodels.jwst.datamodels.AmiOIModel`
            AMI data model that's been normalized
        """
        # Open the target and reference input models
        target_model = self.prepare_output(target, open_as_type=datamodels.AmiOIModel)
        reference_model = self.prepare_output(reference, open_as_type=datamodels.AmiOIModel)

        # Call the normalization routine
        result = ami_normalize.normalize_lg(target_model, reference_model)

        result.meta.cal_step.ami_normalize = "COMPLETE"

        # Close the input models
        if target_model is not target:
            target_model.close()
        if reference_model is not reference:
            reference_model.close()

        # We're done
        return result
