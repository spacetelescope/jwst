from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['RampFitOutputModel']


class RampFitOutputModel(model_base.DataModel):
    """
    A data model for the optional output of the ramp fitting step.

    In the parameter definitions below, `n_int` is the number of
    integrations, `max_seg` is the maximum number of segments that
    were fit, `nreads` is the number of reads in an integration, and
    `ny` and `nx` are the height and width of the image.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    slope : numpy array (n_int, max_seg, ny, nx)

    sigslope : numpy array (n_int, max_seg, ny, nx)

    var_p : numpy array (n_int, max_seg, ny, nx)

    var_r : numpy array (n_int, max_seg, ny, nx)

    yint : numpy array (n_int, max_seg, ny, nx)

    sigyint : numpy array (n_int, max_seg, ny, nx)

    pedestal : numpy array (n_int, max_seg, ny, nx)

    weights : numpy array (n_int, max_seg, ny, nx)

    crmag : numpy array (n_int, max_seg, ny, nx)
    """
    schema_url = "rampfitoutput.schema.yaml"

    def __init__(self, init=None,
                 slope=None,
                 sigslope=None,
                 var_p=None,
                 var_r=None,
                 yint=None,
                 sigyint=None,
                 pedestal=None,
                 weights=None,
                 crmag=None,
                 **kwargs):
        super(RampFitOutputModel, self).__init__(init=init, **kwargs)

        if slope is not None:
            self.slope = slope

        if sigslope is not None:
            self.sigslope = sigslope

        if var_p is not None:
            self.var_p = var_p

        if var_r is not None:
            self.var_r = var_r

        if yint is not None:
            self.yint = yint

        if sigyint is not None:
            self.sigyint = sigyint

        if pedestal is not None:
            self.pedestal = pedestal

        if weights is not None:
            self.weights = weights

        if crmag is not None:
            self.crmag = crmag
