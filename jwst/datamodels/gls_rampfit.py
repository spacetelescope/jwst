from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['GLS_RampFitModel']


class GLS_RampFitModel(model_base.DataModel):
    """
    A data model for the optional output of the ramp fitting step
    for the GLS algorithm.
    """
    schema_url = "gls_rampfit.schema.yaml"

    def __init__(self, init=None,
                 yint=None,
                 sigyint=None,
                 pedestal=None,
                 crmag=None,
                 sigcrmag=None,
                 **kwargs):
        super(GLS_RampFitModel, self).__init__(init=init, **kwargs)

        if yint is not None:
            self.yint = yint

        if sigyint is not None:
            self.sigyint = sigyint

        if pedestal is not None:
            self.pedestal = pedestal

        if crmag is not None:
            self.crmag = crmag

        if sigcrmag is not None:
            self.sigcrmag = sigcrmag
