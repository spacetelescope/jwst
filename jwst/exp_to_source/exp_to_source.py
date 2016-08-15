"""exp_to_source: Reformat Level2b MSA data to be source-based.
"""
from collections import defaultdict

from jwst.datamodels import MultiSlitExposureModel


def exp_to_source(inputs):
    """Reformat exposure-based MSA data to source-based.

    Parameters
    ----------
    inputs: [MultiSlitModel, ...]
        List of MultiSlitModel instances to reformat.

    Returns
    -------
    {str: MultiSlitExposureModel, }
        Returns a dict of MultiSlitExposureModel instances wherein each
        instance contains slits belonging to the same source.
        The key is the name of each source.
    """
    result = defaultdict(MultiSlitExposureModel)
    for exposure in inputs:
        for slit in exposure.slits:
            result[slit.name].slits.append(slit)
            try:
                result[slit.name].meta.exposures.append(
                    exposure.meta.to_tree()
                )
            except AttributeError:
                result[slit.name].meta.exposures = [exposure.meta.to_tree()]
    return result
