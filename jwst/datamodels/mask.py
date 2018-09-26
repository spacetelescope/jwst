from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['MaskModel']


class MaskModel(ReferenceFileModel):
    """
    A data model for 2D masks.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    dq : numpy array
        The data quality array.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "mask.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(MaskModel, self).__init__(init=init, **kwargs)

        if self.dq is not None or self.dq_def is not None:
            self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq

    def get_primary_array_name(self):
        """
        Returns the name "primary" array for this model, which
        controls the size of other arrays that are implicitly created.
        This is intended to be overridden in the subclasses if the
        primary array's name is not "data".
        """
        return 'dq'
