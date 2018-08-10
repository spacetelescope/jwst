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
                 relsens=None, area=None,
                 wavelength_pointsource=None, pathloss_pointsource=None,
                 wavelength_uniformsource=None, pathloss_uniformsource=None,
                 wavelength_pointsource2d=None, pathloss_pointsource2d=None,
                 wavelength_uniformsource2d=None, pathloss_uniformsource2d=None,
                 **kwargs):
        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitDataModel, self).__init__(init=None, **kwargs)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.relsens = init.relsens
            self.area = init.area
            if init.hasattr('wavelength'):
                self.wavelength = init.wavelength
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            for key in kwargs:
                setattr(key, kwargs[key])

            if init.meta.hasattr('wcs'):
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

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if wavelength_pointsource is not None:
            self.wavelength_poointsource = wavelength_pointsource

        if pathloss_pointsource is not None:
            self.pathloss_pointsource = pathloss_pointsource

        if wavelength_uniformsource is not None:
            self.wavelength_uniformsource = wavelength_uniformsource

        if pathloss_uniformsource is not None:
            self.pathloss_uniform_source = pathloss_uniformsource

        if pathloss_pointsource2d is not None:
            self.pathloss_pointsource2d = pathloss_pointsource2d

        if pathloss_pointsource2d is not None:
            self.pathloss_pointsource2d = pathloss_pointsource2d

        if wavelength_uniformsource2d is not None:
            self.wavelength_uniformsource2d = wavelength_uniformsource2d

        if pathloss_uniformsource2d is not None:
            self.pathloss_uniform_source2d = pathloss_uniformsource2d


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

    int_times : table
        The int_times table

    """
    schema_url = "slit.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None,
                 wavelength=None, var_poisson=None, var_rnoise=None,
                 bunit_data=None, bunit_err=None, name=None, xstart=None,
                 xsize=None, ystart=None, ysize=None, slitlet_id=None,
                 source_id=None, source_name=None, source_alias=None,
                 stellarity=None, source_type=None, source_xpos=None, source_ypos=None,
                 shutter_state=None, area=None, relsens=None,
                 int_times=None, barshadow=None,
                 wavelength_pointsource=None, pathloss_pointsource=None,
                 wavelength_uniformsource=None, pathloss_uniformsource=None,
                 wavelength_pointsource2d=None, pathloss_pointsource2d=None,
                 wavelength_uniformsource2d=None, pathloss_uniformsource2d=None,
                 **kwargs):

        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.relsens = init.relsens
            self.area = init.area
            if init.hasattr('wavelength'):
                self.wavelength = init.wavelength
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            if init.hasattr('int_times'):
                self.int_times = init.int_times
            if init.meta.hasattr('wcs'):
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

        if int_times is not None:
            self.int_times = int_times

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if bunit_data is not None:
            self.meta.bunit_data = bunit_data

        if bunit_err is not None:
            self.meta.bunit_err = bunit_err

        if name is not None:
            self.name = name

        if wavelength_pointsource is not None:
            self.wavelength_poointsource = wavelength_pointsource

        if pathloss_pointsource is not None:
            self.pathloss_pointsource = pathloss_pointsource

        if wavelength_uniformsource is not None:
            self.wavelength_uniformsource = wavelength_uniformsource

        if pathloss_uniformsource is not None:
            self.pathloss_uniform_source = pathloss_uniformsource

        if pathloss_pointsource2d is not None:
            self.pathloss_pointsource2d = pathloss_pointsource2d

        if pathloss_pointsource2d is not None:
            self.pathloss_pointsource2d = pathloss_pointsource2d

        if wavelength_uniformsource2d is not None:
            self.wavelength_uniformsource2d = wavelength_uniformsource2d

        if pathloss_uniformsource2d is not None:
            self.pathloss_uniform_source2d = pathloss_uniformsource2d
