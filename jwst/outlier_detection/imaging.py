"""
Submodule for performing outlier detection on imaging data.
"""

import copy
import logging
import os

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelLibrary
from jwst.resample import resample
from jwst.resample.resample_utils import build_driz_weight
from jwst.stpipe.utilities import record_step_status

from stcal.outlier_detection.utils import compute_weight_threshold

from .utils import nanmedian3D, flag_model_crs, flag_resampled_model_crs, OnDiskMedian
from ._fileio import save_median

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["detect_outliers"]


def detect_outliers(
    input_models,
    save_intermediate_results,
    good_bits,
    maskpt,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    resample_data,
    weight_type,
    pixfrac,
    kernel,
    fillval,
    allowed_memory,
    in_memory,
    asn_id,
    make_output_path,
):
    """
    Flag outliers in imaging data.

    input_models is expected to be a ModelLibrary

    See `OutlierDetectionStep.spec` for documentation of these arguments.
    """
    if not isinstance(input_models, ModelLibrary):
        input_models = ModelLibrary(input_models, on_disk=not in_memory)

    if len(input_models) < 2:
        log.warning(f"Input only contains {len(input_models)} exposures")
        log.warning("Outlier detection will be skipped")
        record_step_status(input_models, "outlier_detection", False)
        return input_models
    
    with input_models:
        example_model = input_models.borrow(0)
        output_path = make_output_path(basepath=example_model.meta.filename,
                    suffix='')
        input_models.shelve(example_model, modify=False)
        del example_model
    output_path = os.path.dirname(output_path)

    if resample_data:
        # Start by creating resampled/mosaic images for
        # each group of exposures
        resamp = resample.ResampleData(
            input_models,
            output=output_path,
            single=True,
            blendheaders=False,
            wht_type=weight_type,
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            good_bits=good_bits,
            in_memory=in_memory,
            asn_id=asn_id,
            allowed_memory=allowed_memory,
        )
        median_wcs = resamp.output_wcs

        j = 0
        for group_id, indices in input_models.group_indices.items():

            # single is hardcoded to True above so it's ok to just call resample_group directly
            drizzled_model = resamp.resample_group(input_models, indices)
            if j == 0:
                ngroups = len(input_models.group_names)
                full_shape = (ngroups,) + drizzled_model.data.shape
                if in_memory:
                    data_frames = np.empty(full_shape, dtype=np.float32)
                else:
                    median_computer = OnDiskMedian(full_shape,
                                                   dtype=drizzled_model.data.dtype,
                                                   buffer_size=None)

            if save_intermediate_results:
                output_name = drizzled_model.meta.filename

                if resamp.output_dir is not None:
                    output_name = os.path.join(resamp.output_dir, output_name)
                drizzled_model.save(output_name)
                log.info(f"Saved model in {output_name}")

            # handle the weights right away, so only data array needs to be saved
            weight_threshold = compute_weight_threshold(drizzled_model.wht, maskpt)
            drizzled_model.data[drizzled_model.wht < weight_threshold] = np.nan

            if in_memory:
                data_frames[j] = drizzled_model.data
            else:
                # write spatial sections to disk
                median_computer.add_image(drizzled_model.data)
            j += 1

    else:
        # for non-dithered data, the resampled image is just the original image
        with input_models:
            for i, drizzled_model in enumerate(input_models):
                drizzled_model.wht = build_driz_weight(
                    drizzled_model,
                    weight_type=weight_type,
                    good_bits=good_bits)
                input_models.shelve(drizzled_model, modify=True)

                if i == 0:
                    # copy for when saving median and input is a filename?
                    median_wcs = copy.deepcopy(drizzled_model.meta.wcs)
                    full_shape = (len(input_models),) + drizzled_model.data.shape
                    if in_memory:
                        data_frames = np.empty(full_shape, dtype=np.float32)
                    else:
                        median_computer = OnDiskMedian(full_shape,
                                                    dtype=drizzled_model.data.dtype,
                                                    buffer_size=None)

                if save_intermediate_results:
                    output_name = drizzled_model.meta.filename
                    if output_path is not None:
                        output_name = os.path.join(output_path, output_name)
                    drizzled_model.save(output_name)
                    log.info(f"Saved model in {output_name}")

                # handle the weights right away, so only data array needs to be saved
                weight_threshold = compute_weight_threshold(drizzled_model.wht, maskpt)
                drizzled_model.data[drizzled_model.wht < weight_threshold] = np.nan

                if in_memory:
                    data_frames[i] = drizzled_model.data
                else:
                    # write spatial sections to disk
                    median_computer.add_image(drizzled_model.data)


    # Perform median combination on set of drizzled mosaics
    if in_memory:
        median_data = nanmedian3D(data_frames)
        del data_frames
    else:
        median_data = median_computer.compute_median()
        median_computer.cleanup()

    if save_intermediate_results:
        # make a median model
        median_model = datamodels.ImageModel(median_data)
        median_model.update(drizzled_model)
        median_model.meta.wcs = median_wcs
        save_median(median_model, make_output_path, asn_id)
        del median_model
    del drizzled_model

    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    with input_models:
        for image in input_models:
            if resample_data:
                flag_resampled_model_crs(image,
                                         median_data,
                                         median_wcs,
                                         snr1,
                                         snr2,
                                         scale1,
                                         scale2,
                                         backg,
                                         save_blot=save_intermediate_results,
                                         make_output_path=make_output_path)
            else:
                flag_model_crs(image, median_data, snr1)
            input_models.shelve(image, modify=True)

    return input_models
