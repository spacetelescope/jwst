from .reference import ReferenceFileModel


__all__ = ['DrizParsModel']


class DrizParsModel(ReferenceFileModel):
    """
    A data model for drizzle parameters reference tables.
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/drizpars.schema"
