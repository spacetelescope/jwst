from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from . import (
    ImageModel,
    ModelContainer,
    MultiExposureModel
)

__all__ = ['SourceModelContainer']


# Models that can initiliaze into a SourceModelContainer
VALID_INITS = (
    MultiExposureModel,
)


class SourceModelContainer(ModelContainer):
    """
    A container to make MultiExposureModel look like ModelContainer
    """
    def __init__(self, init=None, **kwargs):

        if not isinstance(init, (self.__class__, ) + VALID_INITS):
            raise TypeError(
                'Input {0!r} cannot initialize a SourceModelContainer'.format(init)
            )

        if isinstance(init, SourceModelContainer):
            super(SourceModelContainer, self).__init__(init, **kwargs)
            self._multiexposure = init._multiexposure
        elif isinstance(init, MultiExposureModel):

            # Convert exposures to ImageModels to allow
            # iteration over various properties.
            models = [
                ImageModel(exposure.instance)
                for exposure in init.exposures
            ]

            super(SourceModelContainer, self).__init__(init=models, **kwargs)
            self._multiexposure = init

    def save(self, *args, **kwargs):
        """Save out the container as a MultiExposureModel"""
        self._multiexposure.exposures = [
            model.instance
            for model in self._models
        ]
        self._multiexposure.save(*args, **kwargs)
