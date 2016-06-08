from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['AmiLgModel']


class AmiLgModel(model_base.DataModel):
    """
    A data model for AMI LG analysis results.
    """
    schema_url = "amilg.schema.yaml"

    def __init__(self, init=None,
                 fit_image=None,
                 resid_image=None,
                 closure_amp_table=None,
                 closure_phase_table=None,
                 fringe_amp_table=None,
                 fringe_phase_table=None,
                 pupil_phase_table=None,
                 solns_table=None, **kwargs):
        super(AmiLgModel, self).__init__(init=init, **kwargs)

        if fit_image is not None:
            self.fit_image = fit_image

        if resid_image is not None:
            self.resid_image = resid_image

        if closure_amp_table is not None:
            self.closure_amp_table = closure_amp_table

        if closure_phase_table is not None:
            self.closure_phase_table = closure_phase_table

        if fringe_amp_table is not None:
            self.fringe_amp_table = fringe_amp_table

        if fringe_phase_table is not None:
            self.fringe_phase_table = fringe_phase_table

        if pupil_phase_table is not None:
            self.pupil_phase_table = pupil_phase_table

        if solns_table is not None:
            self.solns_table = solns_table
