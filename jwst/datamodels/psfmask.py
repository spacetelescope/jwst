from . import model_base


__all__ = ['PsfMaskModel']


class PsfMaskModel(model_base.DataModel):
    """
    A data model for coronagraphic 2D PSF mask reference files

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The 2-D mask array
    """
    schema_url = "psfmask.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        super(PsfMaskModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
