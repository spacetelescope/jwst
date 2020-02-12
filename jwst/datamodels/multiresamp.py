from . import model_base

__all__ = ['MultiResampModel']


class MultiResampModel(model_base.DataModel):
    """
    A data model for spectroscopic resampled slits.

    This model has a special member `slits` that can be used to
    deal with each DrizProduct at a time.  It behaves like a list::


    Parameters
    __________
    products.items.data : numpy float32 array
         resampled science data

    products.items.wht : numpy float32 array
         drizzle algorithm weight array

    products.items.con : numpy int32 array
         drizzle algorithm context array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/multiresamp.schema"


    def __init__(self, init=None, **kwargs):
        super(MultiResampModel, self).__init__(init=init, **kwargs)
