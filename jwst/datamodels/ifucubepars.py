from .reference import ReferenceFileModel

__all__ = ['IFUCubeParsModel']


class IFUCubeParsModel(ReferenceFileModel):
    """
    A data model for IFU Cube  parameters reference tables.
    """
    schema_url = "ifucubepars.schema.yaml"

    def __init__(self, init=None, ifucubepars_table=None, ifucubepars_msn_table=None,**kwargs):
        super(IFUCubeParsModel, self).__init__(init=init, **kwargs)

        if ifucubepars_table is not None:
            self.ifucubepars_table = ifucubepars_table

        if ifucubepars_msn_table is not None:
            self.ifucubepars_msn_table = ifucubepars_msn_table


class NirspecIFUCubeParsModel(IFUCubeParsModel):
    """
    A data model for Nirspec ifucubepars reference files.
    """
    schema_url = "nirspec_ifucubepars.schema.yaml"

    def __init__(self, init=None, ifucubepars_table=None,
                 ifucubepars_msn_table=None,
                 ifucubepars_prism_wavetable=None,
                 ifucubepars_med_wavetable=None,
                 ifucubepars_high_wavetable=None,
                 **kwargs):
        super(NirspecIFUCubeParsModel, self).__init__(init=init, **kwargs)

        if ifucubepars_table is not None:
            self.ifucubepars_table = ifucubepars_table

        if ifucubepars_msn_table is not None:
            self.ifucubepars_msn_table = ifucubepars_msn_table

        if ifucubepars_prism_wavetable is not None:
            self.ifucubepars_prism_wavetable = ifucubepars_prism_wavetable
        if ifucubepars_med_wavetable is not None:
            self.ifucubepars_med_wavetable = ifucubepars_med_wavetable
        if ifucubepars_high_wavetable is not None:
            self.ifucubepars_high_wavetable = ifucubepars_high_wavetable


class MiriIFUCubeParsModel(IFUCubeParsModel):
    """
    A data model for MIRI mrs ifucubepars reference files.
    """
    schema_url = "miri_ifucubepars.schema.yaml"

    def __init__(self, init=None, ifucubepars_table=None,
                 ifucubepars_msn_table=None,
                 ifucubepars_multichannel_wavetable=None,**kwargs):
        super(MiriIFUCubeParsModel, self).__init__(init=init, **kwargs)

        if ifucubepars_table is not None:
            self.ifucubepars_table = ifucubepars_table

        if ifucubepars_msn_table is not None:
            self.ifucubepars_msn_table = ifucubepars_msn_table

        if ifucubepars_multichannel_wavetable is not None:
            self.ifucubepars_multichannel_wavetable = ifucubepars_multichannel_wavetable
