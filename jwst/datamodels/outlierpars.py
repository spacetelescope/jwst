from .reference import ReferenceFileModel

__all__ = ['OutlierParsModel']


class OutlierParsModel(ReferenceFileModel):
    """
    A data model for outlier detection parameters reference tables.
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/outlierpars.schema"
