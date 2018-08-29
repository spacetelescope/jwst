from .reference import ReferenceFileModel

__all__ = ['DrizParsModel']


class DrizParsModel(ReferenceFileModel):
    """
    A data model for drizzle parameters reference tables.
    """
    schema_url = "drizpars.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        super(DrizParsModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
