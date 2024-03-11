"""exp_to_source: Reformat Level2b multi-source data to be source-based.
"""
import logging

from collections import OrderedDict
from collections.abc import Callable

from stdatamodels.properties import merge_tree
from stdatamodels.jwst.datamodels import MultiExposureModel

from jwst.datamodels import SourceModelContainer

__all__ = ['exp_to_source', 'multislit_to_container']

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def exp_to_source(inputs):
    """Reformat exposure-based MSA data to source-based.

    Parameters
    ----------
    inputs : [~jwst.datamodels.MultiSlitModel, ...]
        List of MultiSlitModel instances to reformat.

    Returns
    -------
    multiexposures : dict
        Returns a dict of MultiExposureModel instances wherein each
        instance contains slits belonging to the same source.
        The key is the ID of each source, i.e. ``source_id``.
    """
    result = DefaultOrderedDict(MultiExposureModel)

    for exposure in inputs:
        log.info(f'Reorganizing data from exposure {exposure.meta.filename}')

        for slit in exposure.slits:
            log.debug(f'Copying source {slit.source_id}')
            result_slit = result[str(slit.source_id)]
            result_slit.exposures.append(slit)
            # store values for later use (after merge_tree)
            # these values are incorrectly getting overwritten by
            # the top model.
            slit_bunit = slit.meta.bunit_data
            slit_bunit_err = slit.meta.bunit_err
            slit_model = slit.meta.model_type
            slit_wcsinfo = slit.meta.wcsinfo.instance
            # exposure.meta.bunit_data and bunit_err does not exist
            # before calling merge_tree save these values
            # Before merge_tree the slits have a model_type of SlitModel.
            # After merge_tree it is overwritten with MultiSlitModel.
            # store the model type to undo overwriting of modeltype.

            merge_tree(result_slit.exposures[-1].meta.instance, exposure.meta.instance)

            result_slit.exposures[-1].meta.bunit_data = slit_bunit
            result_slit.exposures[-1].meta.bunit_err = slit_bunit_err
            result_slit.exposures[-1].meta.model_type = slit_model
            result_slit.exposures[-1].meta.wcsinfo = slit_wcsinfo

            if result_slit.meta.instrument.name is None:
                result_slit.update(exposure)

            result_slit.meta.filename = None  # Resulting merged data doesn't come from one file

        exposure.close()

    # Turn off the default factory
    result.default_factory = None

    return result


def multislit_to_container(inputs):
    """Reformat exposure-based MSA data to source-based containers.

    Parameters
    ----------
    inputs : [~jwst.datamodels.MultiSlitModel, ...]
        List of MultiSlitModel instances to reformat, or just a
        ModelContainer full of MultiSlitModels.

    Returns
    -------
    containers : dict
        Returns a dict of ModelContainer instances wherein each
        instance contains ImageModels of slits belonging to the same source.
        The key is the ID of each slit, i.e. ``source_id``.
    """
    containers = exp_to_source(inputs)
    for id in containers:
        containers[id] = SourceModelContainer(containers[id])

    return containers


class DefaultOrderedDict(OrderedDict):
    # Source http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
                not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))
