from .reference import ReferenceFileModel

__all__ = ['IFUCubeParsModel']


class IFUCubeParsModel(ReferenceFileModel):
    """
    A data model for IFU Cube  parameters reference tables.

    Parameters
    __________
    ifucubepars_table : numpy table
         IFU Cube default parameters table
    """
    schema_url = "ifucubepars.schema"


class NirspecIFUCubeParsModel(ReferenceFileModel):
    """
    A data model for Nirspec ifucubepars reference files.
    """
    schema_url = "nirspec_ifucubepars.schema"


class MiriIFUCubeParsModel(ReferenceFileModel):
    """
    A data model for MIRI mrs ifucubepars reference files.

    Parameters
    __________
    ifucubepars_table : numpy table
         default IFU cube  parameters table

    ifucubepars_msn_table : numpy table
         default IFU cube msn parameters table

    ifucubepars_prism_wavetable : numpy table
         default IFU cube prism wavetable

    ifucubepars_med_wavetable : numpy table
         default IFU cube med resolution wavetable

    ifucubepars_high_wavetable : numpy table
         default IFU cube high resolution wavetable
    """
    schema_url = "miri_ifucubepars.schema"
