from .reference import ReferenceFileModel

__all__ = ['MirMrsXArtCorrModel']


class MirMrsXArtCorrModel(ReferenceFileModel):
    """
    A data model for MIRI MRS cross-artifact corrections file.

    Parameters
    __________
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    ch1a_table : numpy table
         Cross artifact correction parameters for Channel 1A

    ch1b_table : numpy table
         Cross artifact correction parameters for Channel 1B

    ch1c_table : numpy table
         Cross artifact correction parameters for Channel 1C

    ch2a_table : numpy table
         Cross artifact correction parameters for Channel 2A

    ch2b_table : numpy table
         Cross artifact correction parameters for Channel 2B

    ch2c_table : numpy table
         Cross artifact correction parameters for Channel 2C

    ch3a_table : numpy table
         Cross artifact correction parameters for Channel 3A

    ch3b_table : numpy table
         Cross artifact correction parameters for Channel 3B

    ch3c_table : numpy table
         Cross artifact correction parameters for Channel 3C

    ch4a_table : numpy table
         Cross artifact correction parameters for Channel 4A

    ch4b_table : numpy table
         Cross artifact correction parameters for Channel 4B

    ch4c_table : numpy table
         Cross artifact correction parameters for Channel 4C

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/miri_mrsxartcorr.schema"
