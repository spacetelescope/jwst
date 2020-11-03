from stdatamodels import DataModel as _DataModel

# from ..lib.basic_utils import deprecate_class


class JwstDataModel(_DataModel):

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/core.schema"

    def __del__(self):
        pass


# @deprecate_class(JwstDataModel)
class DataModel(_DataModel):

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/core.schema"

    def __del__(self):
        pass
