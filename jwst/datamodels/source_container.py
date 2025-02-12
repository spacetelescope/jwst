from .container import ModelContainer

from stdatamodels.jwst.datamodels import MultiExposureModel, SlitModel


__all__ = ["SourceModelContainer"]


# Models that can initialize into a SourceModelContainer
VALID_INITS = (MultiExposureModel,)


class SourceModelContainer(ModelContainer):
    """
    A container to make MultiExposureModel look like ModelContainer.

    The `MultiExposureModel.exposures` list contains the data for each exposure
    from a common slit id. Though the information is the same, the structures
    are not true `SlitModel` instances. This container creates a `SlitModel`
    wrapper around each exposure, such that pipeline code can treat each
    as a `DataModel`.
    """

    def __init__(self, init=None, **kwargs):
        """
        Initialize the container from a MultiExposureModel or another SourceModelContainer.

        Parameters
        ----------
        init : MultiExposureModel, SourceModelContainer, or None, optional
            The models to wrap, by default None
        **kwargs : dict
            Additional arguments to pass to the initializer of the parent class,
            e.g. to `MultiExpsoureModel.__init__()`.
        """
        if not isinstance(init, (self.__class__,) + VALID_INITS):
            raise TypeError(f"Input {init} cannot initialize a SourceModelContainer")

        if isinstance(init, SourceModelContainer):
            super(SourceModelContainer, self).__init__(init, **kwargs)
            self._multiexposure = init._multiexposure  # noqa: SLF001
        elif isinstance(init, MultiExposureModel):
            super(SourceModelContainer, self).__init__(init=None, **kwargs)

            # Convert each exposure to an actual SlitModel.
            # Note that the model is not instantiated
            # since a copy is not desired.
            models = []
            for exposure in init.exposures:
                model = SlitModel()
                model._instance.update(exposure._instance)  # noqa: SLF001
                models.append(model)
            self._models = models
            self._multiexposure = init

    @property
    def multiexposure(self):
        """
        Return an updated version of the MultiExposureModel that is being wrapped.

        Returns
        -------
        MultiExposureModel
            The MultiExposureModel being wrapped, be updated with any new data in the container.
        """
        # Reapply models back to the exposures
        for exposure, model in zip(self._multiexposure.exposures, self._models, strict=True):
            exposure._instance.update(model._instance)  # noqa: SLF001

        return self._multiexposure

    def save(self, path=None, dir_path=None, save_model_func=None, *args, **kwargs):
        """
        Save out the container as a MultiExposureModel.

        Parameters
        ----------
        path : str
            The path to the output file.
        dir_path : str
            The path to the output directory.
        save_model_func : callable
            A function to save the model.
        *args, **kwargs : tuple, dict
            Additional arguments to pass to the save function
        """
        if save_model_func is None:
            self.multiexposure.save(*args, path=path, dir_path=dir_path, **kwargs)
        else:
            save_model_func(self.multiexposure, output_file=path)
