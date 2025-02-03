"""Submodule for performing outlier detection on coronagraphy data."""

import logging

import numpy as np

from stdatamodels.jwst import datamodels

from jwst.resample.resample_utils import build_mask

from .utils import create_cube_median, flag_model_crs
from ._fileio import save_median

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["detect_outliers"]


def detect_outliers(
    input_model,
    save_intermediate_results,
    good_bits,
    maskpt,
    snr,
    make_output_path,
):
    """
    Flag outliers in coronography data.

    Parameters
    ----------
    input_model : ~jwst.datamodels.CubeModel
        The input cube model.
    save_intermediate_results : bool
        If True, save the median model.
    good_bits : int
        DQ flag bit values indicating good pixels.
    maskpt : float
        The percentage of the mean weight to use as a threshold for masking.
    snr : float
        The signal-to-noise ratio threshold for flagging outliers.
    make_output_path : callable
        A function that generates a path for saving intermediate results.

    Returns
    -------
    ~jwst.datamodels.CubeModel
        The input model with outliers flagged.
    """
    if not isinstance(input_model, datamodels.JwstDataModel):
        input_model = datamodels.open(input_model)

    if not isinstance(input_model, datamodels.CubeModel):
        raise TypeError(f"Input must be a CubeModel: {input_model}")

    # FIXME weight_type could now be used here. Similar to tso data coron
    # data was previously losing var_rnoise due to the conversion from a cube
    # to a ModelContainer (which makes the default ivm weight ignore var_rnoise).
    # Now that it's handled as a cube we could use the var_rnoise.
    input_model.wht = build_mask(input_model.dq, good_bits).astype(np.float32)

    # Perform median combination on set of drizzled mosaics
    median_data = create_cube_median(input_model, maskpt)

    if save_intermediate_results:
        # make a median model
        median_model = datamodels.ImageModel(median_data)
        median_model.update(input_model)
        median_model.meta.wcs = input_model.meta.wcs

        save_median(median_model, make_output_path)
        del median_model

    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    flag_model_crs(
        input_model,
        median_data,
        snr,
    )
    return input_model
