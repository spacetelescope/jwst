from .reference import ReferenceFileModel

__all__ = ['FringeFreqModel']


class FringeFreqModel(ReferenceFileModel):
    """
    A data model for 2D fringe correction images.

    Parameters
    __________
    rfc_freq_short : numpy table
    rfc_freq_medium : numpy table
    rfc_freq_long : numpy table
    max_amp : numpy table
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/fringe_freq.schema"
