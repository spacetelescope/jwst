import logging

from jwst.datamodels import ModelContainer
from jwst.stpipe import Step

__all__ = ["PFPCStep"]

log = logging.getLogger(__name__)


class PFPCStep(Step):
    """Apply a point fixed pattern correction (PFPC) to dithered point source spectra."""

    class_alias = "pfpc"

    spec = """
    """  # noqa: E501

    def process(self, input_data):
        """
        Extract spectra and apply a point fixed pattern correction (PFPC).

        Parameters
        ----------
        input_data : str or `~stdatamodels.jwst.datamodels.JwstDataModel` \
                     or `~jwst.datamodels.container.ModelContainer`
            The input filename, datamodel, or container of 2D spectral data.
            If the input is a container, each contained model is separately
            processed.

        Returns
        -------
        output_model : `~stdatamodels.jwst.datamodels.JwstDataModel` \
                       or `~jwst.datamodels.container.ModelContainer`
            Extracted spectra with fixed pattern correction applied. If
            no correction can be applied, the input data is returned.
        """
        output_model = self.prepare_output(input_data)
        if isinstance(output_model, ModelContainer):
            models = output_model
        else:
            models = [output_model]

        # Set up output path name to include the ASN ID if available
        self.add_asn_id_to_output_name(models)

        return output_model
