from . import model_base

__all__ = ['MultiResampModel']


class MultiResampModel(model_base.DataModel):
    """
    A data model for multi-resampled images or resampled multi-slit images.

    This model has a special member `products` that can be used to
    deal with each resampled image/slit data  at a time.  It behaves like a list::

       >>> from . import MultiResampModel
       >>> multiresamp_model = MultiResampModel()
       >>> multiresamp_model.products.append(multiResampModel())
       >>> multiresap_model.products[0] # doctest: +SKIP

    If `init` is a file name or an `MultiResamptModel` instance, an empty
    `MultiResampModel` will be created and assigned to attribute `products[0]`,
    and the `data`, `wht`, and `con` attributes from the
    input file or `DrizProductModel` will be copied to the first element of
    `products`.

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
        super(MultiResampModel, self).__init__(init=None, **kwargs)
        self.update(init)
        #self.products.append(self.products.item())
        #self.products[0].data = init.data
        #self.products[0].wht = init.wht
        #self.products[0].con = init.con
        #return

        super(MultiResampModel, self).__init__(init=init, **kwargs)
