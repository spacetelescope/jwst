"""exp_to_source: Reformat Level2b MSA data to be source-based.
"""
from collections import defaultdict

from jwst.datamodels import MultiExposureModel, ModelContainer, ImageModel
from jwst.datamodels.properties import merge_tree

__all__ = ['exp_to_source', 'multislit_to_container']


def exp_to_source(inputs):
    """Reformat exposure-based MSA data to source-based.

    Parameters
    ----------
    inputs: [MultiSlitModel, ...]
        List of MultiSlitModel instances to reformat.

    Returns
    -------
    {str: MultiExposureModel, }
        Returns a dict of MultiExposureModel instances wherein each
        instance contains slits belonging to the same source.
        The key is the name of each source.
    """
    result = defaultdict(MultiExposureModel)
    for exposure in inputs:
        for slit in exposure.slits:
            result[slit.name].exposures.append(slit)
            result[slit.name].exposures[-1].meta = exposure.meta
    return result


def multislit_to_container(inputs):
    """Reformat exposure-based MSA data to source-based containers.

    Parameters
    ----------
    inputs: [MultiSlitModel, ...]
        List of MultiSlitModel instances to reformat, or just a 
        ModelContainer full of MultiSlitModels.

    Returns
    -------
    {str: ModelContainer, }
        Returns a dict of ModelContainer instances wherein each
        instance contains ImageModels of slits belonging to the same source.
        The key is the name of each slit.
    """
    result = defaultdict(ModelContainer)
    for exposure in inputs:
        for slit in exposure.slits:
            result[slit.name].append(ImageModel(slit.instance))
            merge_tree(result[slit.name][-1].meta.instance, 
                exposure.meta.instance)
    return result
