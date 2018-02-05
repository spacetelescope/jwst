from . import model_base
from .image import ImageModel


__all__ = ['SlitModel', 'SlitDataModel']


class SlitDataModel(model_base.DataModel):
    """
    A data model for 2D images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.

    dq : numpy array
        The data quality array.

    err : numpy array
        The error array.

    relsens : numpy array
        The relative sensitivity table.

    """

    schema_url = "slitdata.schema.yaml"


    def __init__(self, init=None, data=None, dq=None, err=None,
                 wavelength=None, var_poisson=None, var_rnoise=None,
                 relsens=None, area=None, **kwargs):
        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitDataModel, self).__init__(init=None, **kwargs)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.relsens = init.relsens
            self.area = init.area
            if hasattr(init, 'wavelength'):
                self.wavelength = init.wavelength
            if hasattr(init, 'var_poisson'):
                self.var_poisson = init.var_poisson
            if hasattr(init, 'var_rnoise'):
                self.var_rnoise = init.var_rnoise
            for key in kwargs:
                setattr(key, kwargs[key])

            if hasattr(init.meta, 'wcs'):
                self.meta.wcs = init.meta.wcs
            else:
                self.meta.wcs = None
            return
        super(SlitDataModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if wavelength is not None:
            self.wavelength = wavelength

        if var_poisson is not None:
            self.var_poisson = var_poisson

        if var_rnoise is not None:
            self.var_rnoise = var_rnoise

        if relsens is not None:
            self.relsens = relsens

        if area is not None:
            self.area = area

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err


class SlitModel(model_base.DataModel):
    """
    A data model for 2D images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.

    dq : numpy array
        The data quality array.

    err : numpy array
        The error array.

    relsens : numpy array
        The relative sensitivity table.

    """
    schema_url = "slit.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None,
                 wavelength=None, var_poisson=None, var_rnoise=None,
                 bunit_data=None, bunit_err=None, name=None, xstart=None,
                 xsize=None, ystart=None, ysize=None, slitlet_id=None,
                 source_id=None, source_name=None, source_alias=None,
                 stellarity=None, source_type=None, source_xpos=None, source_ypos=None,
                 shutter_state=None, area=None, relsens=None, barshadow=None, **kwargs):

        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.relsens = init.relsens
            self.area = init.area
            if hasattr(init, 'var_poisson'):
                self.var_poisson = init.var_poisson
            if hasattr(init, 'var_rnoise'):
                self.var_rnoise = init.var_rnoise
            if hasattr(init.meta, 'wcs'):
                self.meta.wcs = init.meta.wcs
            else:
                self.meta.wcs = None
            return

        super(SlitModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
        if dq is not None:
            self.dq = dq
        if err is not None:
            self.err = err
        if wavelength is not None:
            self.wavelength = wavelength
        if kwargs:
            for key in kwargs:
                setattr(self, key, kwargs[key])

        if var_poisson is not None:
            self.var_poisson = var_poisson

        if var_rnoise is not None:
            self.var_rnoise = var_rnoise
        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

        if name is not None:
            self.name = name

