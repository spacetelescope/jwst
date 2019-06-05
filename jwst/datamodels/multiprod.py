from . import model_base
from .drizproduct import DrizProductModel

__all__ = ['MultiProductModel']


class MultiProductModel(model_base.DataModel):
    """
    A data model for multi-DrizProduct images.

    This model has a special member `products` that can be used to
    deal with each DrizProduct at a time.  It behaves like a list::

       >>> from . import DrizProductModel
       >>> multiprod_model = MultiProductModel()
       >>> multiprod_model.products.append(DrizProductModel())
       >>> multiprod_model.products[0] # doctest: +SKIP
       <DrizProductModel>

    If `init` is a file name or an `DrizProductModel` instance, an empty
    `DrizProductModel` will be created and assigned to attribute `products[0]`,
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
    schema_url = "multiproduct.schema"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, DrizProductModel):
            super(MultiProductModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.products.append(self.products.item())
            self.products[0].data = init.data
            self.products[0].wht = init.wht
            self.products[0].con = init.con
            return

        super(MultiProductModel, self).__init__(init=init, **kwargs)
