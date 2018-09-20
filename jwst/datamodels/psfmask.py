from .reference import ReferenceFileModel


__all__ = ['PsfMaskModel']


class PsfMaskModel(ReferenceFileModel):
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
