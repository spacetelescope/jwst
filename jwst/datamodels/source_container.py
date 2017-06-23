from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from . import (
    ModelContainer,
    MultiExposureModel
)

__all__ = ['SourceContainerModel']


class SourceContainerModel(ModelContainer):
    """
    A container to make MultiExposureModel look like ModelContainer
    """
    def __init__(self, init=None, **kwargs):

        if not isinstance(init, MultiExposureModel):
            raise TypeError('Input {0!r} must be a MultiExposureModel')

        super(SourceContainerModel, self).__init__(init=None, **kwargs)

        self._multiexposure = init
        self._models = self._multiexposure.exposures
