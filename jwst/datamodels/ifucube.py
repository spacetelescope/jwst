from . import model_base


__all__ = ['IFUCubeModel']


class IFUCubeModel(model_base.DataModel):
    """
    A data model for 3D IFU  cubes.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.  3-D.

    dq : numpy array
        The data quality array.  3-D.

    err : numpy array
        The error array.  3-D

    weightmap : numpy array
        The weight map array.  3-D
    """
    schema_url = "ifucube.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, 
                 weightmap=None, hdrtab=None,  **kwargs):
        super(IFUCubeModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if weightmap is not None:
            self.weightmap = weightmap

        if hdrtab is not None:
            self.hdrtab = hdrtab
        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
