"""
The ever-present utils sub-module. A home for all...
"""
import warnings
import psutil

import numpy as np

from jwst.resample.resample import compute_image_pixel_area
from stcal.outlier_detection.utils import compute_weight_threshold, gwcs_blot, flag_crs, flag_resampled_crs
from stdatamodels.jwst import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
OUTLIER = datamodels.dqflags.pixel['OUTLIER']
_ONE_MB = 1 << 20


def create_cube_median(cube_model, maskpt):
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)

    median = np.nanmedian(
        np.ma.masked_array(cube_model.data, np.less(cube_model.wht, weight_threshold), fill_value=np.nan),
        axis=0,
    )
    return median


def create_median(resampled_models, maskpt, on_disk=True, allowed_memory=None):
    """Create a median image from the singly resampled images.

    Parameters
    ----------
    resampled_models : ModelLibrary
        The singly resampled images.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    on_disk : bool
        If True, the input models are on disk and will be read in chunks.

    allowed_memory : float
        The fraction of available memory to be used by create_median.
        If None, use 50% of the machine's available memory.
        This parameter has no effect if on_disk is False.

    Returns
    -------
    median_image : ndarray
        The median image.
    """
    # Compute the weight threshold for each input model
    weight_thresholds = []
    with resampled_models:
        for resampled in resampled_models:
            weight_threshold = compute_weight_threshold(resampled.wht, maskpt)
            weight_thresholds.append(weight_threshold)
            resampled_models.shelve(resampled, modify=False)

    # compute median over all models
    if not on_disk:
        # easier case: all models in library can be loaded into memory at once
        model_list = []
        with resampled_models:
            for resampled in resampled_models:
                model_list.append(resampled.data)
                resampled_models.shelve(resampled, modify=False)
                del resampled
        return np.nanmedian(np.array(model_list), axis=0)
    else:
        # set up buffered access to all input models
        log.info("Computing median in chunks to save memory")
        if allowed_memory is None:
            allowed_memory = 0.5
        machine_available_memory = psutil.virtual_memory().available + psutil.swap_memory().total
        available_memory = machine_available_memory * allowed_memory

        # buffer_size sets size of chunk in MB that will be read into memory per model
        buffer_size = available_memory / len(resampled_models.group_names) / _ONE_MB
        log.info(f"Chunk size {buffer_size} MB per model for {len(resampled_models)} models (totaling {allowed_memory * 100}% of available memory)")

        with resampled_models:
            example_model = resampled_models.borrow(0)
            shp = example_model.data.shape
            dtype = example_model.data.dtype
            nsections, section_nrows = _compute_buffer_indices(example_model, buffer_size)
            resampled_models.shelve(example_model, modify=False)
            del example_model

        if nsections == 1:
            log.info("Chunk size is large enough to read entire model at once")
        else:
            log.info(f"Data will be read in {nsections} sections of {section_nrows} rows each")

        # get spatial sections of library and compute timewise median, one by one
        resampled_sections = _get_sections(resampled_models, nsections, section_nrows, shp[0])
        median_image_empty = np.empty(shp, dtype) * np.nan
        return _create_median(resampled_sections, resampled_models, weight_thresholds, median_image_empty)


def _get_sections(library, nsections, section_nrows, imrows):
    """Iterator to return sections from a ModelLibrary.

    Parameters
    ----------
    library : ModelLibrary
        The input data models.

    nsections : int
        The number of spatial sections in each model

    section_nrows : int
        The number of rows in each section

    imrows : int
        The total number of rows in the image

    Yields
    ------
    data_subset : ndarray
        array of shape (len(library), section_nrows, ncols) representing a spatial
        subset of all the data arrays in the library

    weight_subset : ndarray
        weights corresponding to data_list

    row_range : tuple
        The range of rows in the image covered by the data arrays
    """
    with library:
        example_model = library.borrow(0)
        dtype = example_model.data.dtype
        dtype_wht = example_model.wht.dtype
        shp = example_model.data.shape
        library.shelve(example_model, 0, modify=False)
        del example_model
    for i in range(nsections):
        row1 = i * section_nrows
        row2 = min(row1 + section_nrows, imrows)
        
        data_subset = np.empty((len(library), row2 - row1, shp[1]), dtype)
        weight_subset = np.empty((len(library), row2 - row1, shp[1]), dtype_wht)
        with library:
            for j, model in enumerate(library):
                data_subset[j] = model.data[row1:row2]
                weight_subset[j] = model.wht[row1:row2]
                library.shelve(model, j, modify=False)
                del model
        yield (data_subset, weight_subset, (row1, row2))


def _compute_buffer_indices(model, buffer_size=None):

    imrows, imcols = model.data.shape
    data_item_size = model.data.itemsize
    min_buffer_size = imcols * data_item_size
    buffer_size = min_buffer_size if buffer_size is None else (buffer_size * _ONE_MB)
    section_nrows = min(imrows, int(buffer_size // min_buffer_size))
    if section_nrows == 0:
        buffer_size = min_buffer_size
        log.warning("WARNING: Buffer size is too small to hold a single row."
                        f"Increasing buffer size to {buffer_size / _ONE_MB}MB")
        section_nrows = 1
    nsections = int(np.ceil(imrows / section_nrows))
    return nsections, section_nrows


def _create_median(resampled_sections, resampled_models, weight_thresholds, median_image_empty):
    median_image = median_image_empty
    for (resampled_sci, resampled_weight, (row1, row2)) in resampled_sections:
        # Create a mask for each input image, masking out areas where there is
        # no data or the data has very low weight
        badmasks = []
        for weight, weight_threshold in zip(resampled_weight, weight_thresholds):
            badmask = np.less(weight, weight_threshold)
            log.debug("Percentage of pixels with low weight: {}".format(
                np.sum(badmask) / len(weight.flat) * 100))
            badmasks.append(badmask)

        # Fill resampled_sci array with nan's where mask values are True
        for f1, f2 in zip(resampled_sci, badmasks):
            for elem1, elem2 in zip(f1, f2):
                elem1[elem2] = np.nan

        del badmasks

        # For a stack of images with "bad" data replaced with Nan
        # use np.nanmedian to compute the median.
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore",
                                    message="All-NaN slice encountered",
                                    category=RuntimeWarning)
            median_image[row1:row2] = np.nanmedian(resampled_sci, axis=0)
        del resampled_sci, resampled_weight

    return median_image


def flag_crs_in_models(
    input_models,
    median_data,
    snr1,
):
    for image in input_models:
        # dq flags will be updated in-place
        flag_model_crs(image, median_data, snr1)
    

def flag_resampled_model_crs(
    input_model,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
):
    if 'SPECTRAL' not in input_model.meta.wcs.output_frame.axes_type:
        input_pixflux_area = input_model.meta.photometry.pixelarea_steradians
        # Set array shape, needed to compute image pixel area
        input_model.meta.wcs.array_shape = input_model.shape
        input_pixel_area = compute_image_pixel_area(input_model.meta.wcs)
        pix_ratio = np.sqrt(input_pixflux_area / input_pixel_area)
    else:
        pix_ratio = 1.0

    blot = gwcs_blot(median_data, median_wcs, input_model.data.shape, input_model.meta.wcs, pix_ratio)
    # dq flags will be updated in-place
    _flag_resampled_model_crs(input_model, blot, snr1, snr2, scale1, scale2, backg)


def _flag_resampled_model_crs(
    input_model,
    blot,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
):
    # If the datamodel has a measured background that has not been subtracted
    # use it instead of the user provided backg.
    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (input_model.meta.background.subtracted is False and
            input_model.meta.background.level is not None):
        backg = input_model.meta.background.level
        log.debug(f"Adding background level {backg} to blotted image")

    cr_mask = flag_resampled_crs(input_model.data, input_model.err, blot, snr1, snr2, scale1, scale2, backg)

    # update the dq flags in-place
    input_model.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)
    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")


def flag_crs_in_models_with_resampling(
    input_models,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
):
    for image in input_models:
        flag_resampled_model_crs(image, median_data, median_wcs, snr1, snr2, scale1, scale2, backg)


def flag_model_crs(image, blot, snr):
    cr_mask = flag_crs(image.data, image.err, blot, snr)
    # update dq array in-place
    image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)
    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")
