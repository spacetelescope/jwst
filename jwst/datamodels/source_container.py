from . import ModelContainer, MultiExposureModel

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
            super(SourceModelContainer, self).__init__(init=None, **kwargs)
            self._models = init.exposures
            self._multiexposure = init

    def save(self,
             path=None,
             dir_path=None,
             save_model_func=None,
             *args, **kwargs):
        """Save out the container as a MultiExposureModel"""

        if save_model_func is None:
            self._multiexposure.save(
                path=path,
                dir_path=dir_path,
                *args, **kwargs
            )
        else:
            save_model_func(self._multiexposure, output_file=path)
