from .reference import ReferenceFileModel


__all__ = ['TrapParsModel']


class TrapParsModel(ReferenceFileModel):
    """
    A data model for trap capture and decay parameters.

    Parameters
    __________
    trappars_table : numpy table
         Trap capture and decay parameters
         A table with three columns for trap-capture parameters and one
         column for the trap-decay parameter.  Each row of the table is
         for a different trap family.
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/trappars.schema"
