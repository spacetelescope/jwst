from .model_base import JwstDataModel


__all__ = ['MSATargAcqModel']


class MSATargAcqModel(JwstDataModel):
    """
    A data model for MSA Target Acquisition

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/msatargacq.schema"
