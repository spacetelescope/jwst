from __future__ import absolute_import, unicode_literals, division, print_function

from . import ramp


__all__ = ['MIRIRampModel']


class MIRIRampModel(ramp.RampModel):
    """
    A data model for MIRI ramps. Includes the ``refout`` array.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    data : numpy array
        The science data.

    pixeldq : numpy array
        2-D data quality array.

    groupdq : numpy array
        3-D or 4-D data quality array.

    err : numpy array
        The error array.

    refout : numpy array
        The array of reference output data.

    group : table
        The group parameters table.

    """
    schema_url = "miri_ramp.schema.yaml"

    def __init__(self, init=None, data=None, pixeldq=None, groupdq=None,
                 err=None, refout=None, zeroframe=None, group=None, **kwargs):
        super(MIRIRampModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if pixeldq is not None:
            self.pixeldq = pixeldq

        if groupdq is not None:
            self.groupdq = groupdq

        if err is not None:
            self.err = err

        if refout is not None:
            self.refout = refout

        if zeroframe is not None:
            self.zeroframe = zeroframe

        if group is not None:
            self.group = group

        # Implicitly create arrays
        self.pixeldq = self.pixeldq
        self.groupdq = self.groupdq
        self.err = self.err
