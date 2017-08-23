from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['SlitModel']


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

    relsens2d: numpy array
        The relative sensitivty 2D array.
    """
    schema_url = "slit.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, relsens=None,
                 relsens2d=None, area=None, pathloss_pointsource=None,
                 wavelength_pathloss=None, pathloss_uniformsource=None, 
                 wavelength_uniformsource=None, wavelength=None, **kwargs):
        super(SlitModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if relsens is not None:
            self.relsens = relsens

        if relsens2d is not None:
            self.relsens2d = relsens2d

        if area is not None:
            self.area = area
        if pathloss_pointsource is not None:
            self.pathlosss_pointsource = pathloss_pointsource
        if wavelength_pathloss is not None:
            self.wavelength_pathlosss = wavelength_pathloss
        # what about zeropoint? Check if ms.slits[0] has zeropoint
        if pathloss_uniformsource is not None:
            self.pathloss_uniformsource = pathloss_uniformsource
        if wavelength_uniformsource is not None:
            self.wavelength_uniformsource = wavelength_uniformsource
        if wavelength is not None:
            self.wavelength = wavelength
        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
