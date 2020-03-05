from .reference import ReferenceFileModel


__all__ = ['RSCDModel']


class RSCDModel(ReferenceFileModel):
    """
    A data model for the RSCD reference file.

    Parameters
    __________
    rscd_group_skip_table : numpy table
        Reference table for RSCD correction baseline correction
        A table with 3 columms that set the number of group to skip for each
        subarray and readpatt

    rscd_gen_table : numpy table
        Reference table for RSCD correction enhanced correction
        A table with 5 columms that sets up general parameters for the enhanced 

    rscd_int1_table : numpy table
        Reference table for RSCD correction enhanced correction
        A table with 7 columms that sets up correction parameters  for the enhanced 
        correction for integration 1 based on subarray, readpatt, even or odd row

    rscd_int2_table : numpy table
        Reference table for RSCD correction enhanced correction
        A table with 7 columms that sets up correction parameters  for the enhanced 
        correction for integration 2 based on subarray, readpatt, even or odd row

    rscd_int3_table : numpy table
        Reference table for RSCD correction enhanced correction
        A table with 7 columms that sets up correction parameters  for the enhanced 
        correction for integration 3 based on subarray, readpatt, even or odd row

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/rscd.schema"
