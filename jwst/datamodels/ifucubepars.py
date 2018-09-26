from .reference import ReferenceFileModel

__all__ = ['IFUCubeParsModel']


class IFUCubeParsModel(ReferenceFileModel):
    """
    A data model for IFU Cube  parameters reference tables.
    """
    schema_url = "ifucubepars.schema.yaml"


class NirspecIFUCubeParsModel(ReferenceFileModel):
    """
    A data model for Nirspec ifucubepars reference files.
    """
    schema_url = "nirspec_ifucubepars.schema.yaml"


class MiriIFUCubeParsModel(ReferenceFileModel):
    """
    A data model for MIRI mrs ifucubepars reference files.
    """
    schema_url = "miri_ifucubepars.schema.yaml"

