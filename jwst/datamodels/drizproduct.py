from .model_base import DataModel


__all__ = ['DrizProductModel']


class DrizProductModel(DataModel):
    """
    A data model for drizzle-generated products.

    Attributes
    __________
    data : numpy float32 array
         The science data

    con : numpy int32 array
         Drizzle Context array

    wht : numpy float32 array
         Drizzle Weight array

    relsens : numpy table
         relative sensitivity table
    """
    schema_url = "drizproduct.schema.yaml"

    @property
    def hdrtab(self):
        return self._extra_fits.HDRTAB.data

    @hdrtab.setter
    def hdrtab(self, v):
        self._extra_fits.HDRTAB.data = v
