"""exp_to_source: Reformat Level2b MSA data to be source-based.
"""
from collections import defaultdict

from jwst.datamodels import MultiSlitModel


def exp_to_source(inputs):
    """Reformat exposure-based MSA data to source-based.

    Parameters
    ----------
    inputs: [MultiSlitModel, ...]
        List of MultiSlitModel instances to reformat.

    Returns
    -------
    {str: MultiSlitModel, }
        Returns a dict of MultiSlitModel instances wherein each
        instance contains slits belonging to the same source.
        The key is the name of each source.
    """
    result = defaultdict(MultiSlitModel)
    for exposure in inputs:
        for slit in exposure.slits:
            result[slit.name].slits.append(slit)
    return result
