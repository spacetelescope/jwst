"""
The ever-present utils sub-module. A home for all...
"""

import warnings
import numpy as np
import tempfile
from pathlib import Path

from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.resample.resample import compute_image_pixel_area
from stcal.outlier_detection.utils import compute_weight_threshold, gwcs_blot, flag_crs, flag_resampled_crs
from stdatamodels.jwst import datamodels

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
OUTLIER = datamodels.dqflags.pixel['OUTLIER']
_ONE_MB = 1 << 20


class DiskAppendableArray:
    
    def __init__(self, slice_shape, dtype='float32', tempdir="", filestem="data"):
        self._temp_dir = tempfile.TemporaryDirectory(dir=tempdir)
        self._temp_path = Path(self._temp_dir.name)
        self._slice_shape = slice_shape
        self._dtype = dtype
        self._filename = self._temp_path / Path(filestem + ".bits")
        self._append_count = 0
        with open(self._filename, 'wb') as f: # Noqa: F841
            pass

    @property
    def shape(self):
        return (self._append_count,) + self._slice_shape
        
    def append(self, data):
        if data.shape != self._slice_shape:
            raise ValueError(f"Data shape {data.shape} does not match slice shape {self._slice_shape}")
        if data.dtype != self._dtype:
            raise ValueError(f"Data dtype {data.dtype} does not match array dtype {self._dtype}")
        with open(self._filename, 'ab') as f:
            f.write(data.tobytes())
        self._append_count += 1

    def read(self):
        shp = (self._append_count,) + self._slice_shape
        return np.fromfile(self._filename, dtype=self._dtype).reshape(shp)


def create_cube_median(cube_model, maskpt):
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)

    median = np.nanmedian(
        np.ma.masked_array(cube_model.data, np.less(cube_model.wht, weight_threshold), fill_value=np.nan),
        axis=0,
    )
    return median


def create_median(resampled_models, maskpt, buffer_size=None):
    """Create a median image from the singly resampled images.

    Parameters
    ----------
    resampled_models : ModelLibrary
        The singly resampled images.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    buffer_size : int
        The buffer size for the median computation, units of bytes.

    Returns
    -------
    median_image : ndarray
        The median image.
    """
    on_disk = resampled_models._on_disk
    if on_disk and buffer_size is None:
        raise ValueError("Buffer size must be provided when resampled models are on disk")
    
    # Compute the weight threshold for each input model
    weight_thresholds = []
    model_list = []
    with resampled_models:
        for resampled in resampled_models:
            weight_threshold = compute_weight_threshold(resampled.wht, maskpt)
            if not on_disk:
                # handle weights right away for in-memory case
                data = resampled.data.copy()
                data[resampled.wht < weight_threshold] = np.nan
                model_list.append(data)
                del data
            else:
                weight_thresholds.append(weight_threshold)
            resampled_models.shelve(resampled, modify=False)
        del resampled
    
    # easier case: all models in library can be loaded into memory at once
    if not on_disk:
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore",
                                    message="All-NaN slice encountered",
                                    category=RuntimeWarning)
            return np.nanmedian(np.array(model_list), axis=0)
    else:
        # set up buffered access to all input models
        # get spatial sections of library and compute timewise median, one by one
        resampled_sections, row_indices = _write_sections(resampled_models, weight_thresholds, buffer_size)
        return _create_median(resampled_sections, row_indices)


def _write_sections(library, weight_thresholds, buffer_size):
    """Write spatial sections from a ModelLibrary into temporary files
    grouped along the time axis.

    Parameters
    ----------
    library : ModelLibrary
        The input data models.

    weight_thresholds : list
        The weight thresholds for masking out low weight pixels.

    buffer_size : int
        The buffer size for the median computation, units of bytes.

    Returns
    -------
    temp_arrays : list
        A list of DiskAppendableArray objects, each holding a spatial section
        of every input model, stacked along the time axis.

    row_indices : list
        A list of tuples, each containing the start and end row indices of the
        spatial section in the original data model.
    """

    with library:

        # get an example model to determine dtype, shape, buffer indices
        example_model = library.borrow(0)
        dtype = example_model.data.dtype
        shp = example_model.data.shape
        itemsize = example_model.data.itemsize
        imrows = shp[0]

        # compute buffer indices
        nsections, section_nrows = _compute_buffer_indices((len(library),)+shp, itemsize, buffer_size)
        library.shelve(example_model, 0, modify=False)
        del example_model

        # set up temp file handlers for each section
        # handle zeroth section separately just to get the tempdir for the rest
        slice_shape = (section_nrows, shp[1])
        arr0 = DiskAppendableArray(slice_shape, filestem="section0", dtype=dtype)
        tempdir = arr0._temp_dir.name
        temp_arrays = [arr0,]
        for i in range(1, nsections-1):
            arr = DiskAppendableArray(slice_shape, tempdir=tempdir, filestem=f"section{i}", dtype=dtype)
            temp_arrays.append(arr)

        # handle the last section separately because it has a different shape
        slice_shape_last = (imrows - (nsections-1) * section_nrows, shp[1])
        arrn = DiskAppendableArray(slice_shape_last, tempdir=tempdir, filestem=f"section{nsections-1}", dtype=dtype)
        temp_arrays.append(arrn)
            
        # now append data from each model to all the sections
        row_indices = []
        for j, model in enumerate(library):
            for i in range(nsections):
                row1 = i * section_nrows
                row2 = min(row1 + section_nrows, imrows)
                if j == 0:
                    row_indices.append((row1, row2))
                arr = temp_arrays[i]
                sci = model.data[row1:row2]

                # handle weight thresholding right here, while array is open
                thresh = weight_thresholds[j]
                sci[model.wht[row1:row2] < thresh] = np.nan
                arr.append(sci)

            library.shelve(model, j, modify=False)
        del model

    return temp_arrays, row_indices


def _compute_buffer_indices(shape, itemsize, buffer_size):
    """
    Parameters
    ----------
    shape : tuple
        The shape of the full input, ie, (n_images, imrows, imcols).

    itemsize : int
        The size of a single array element in bytes.

    buffer_size : int
        The buffer size for the median computation, units of bytes.

    Returns
    -------
    nsections : int
        The number of sections to divide the input data into.

    section_nrows : int
        The number of rows in each section (except the last one).
    """

    nimages, imrows, imcols = shape
    per_model_buffer_size = buffer_size / nimages
    min_buffer_size = imcols * itemsize
    section_nrows = min(imrows, int(per_model_buffer_size // min_buffer_size))

    if section_nrows == 0:
        buffer_size = min_buffer_size * nimages
        log.warning("WARNING: Buffer size is too small to hold a single row."
                        f"Increasing buffer size to {buffer_size / _ONE_MB}MB")
        section_nrows = 1

    nsections = int(np.ceil(imrows / section_nrows))
    log.info(f"Computing median over {nimages} images in {nsections} sections with total memory buffer {buffer_size / _ONE_MB} MB")
    return nsections, section_nrows


def _create_median(resampled_sections, row_indices):
    """
    Parameters
    ----------
    resampled_sections : list
        List of DiskAppendableArray objects, each holding a spatial section
        of every input model, stacked along the time axis.

    row_indices : list
        List of tuples, each containing the start and end row indices of the
        spatial section in the original data model.

    Returns
    -------
    median_image : ndarray
        The median image.
    """

    dtype = resampled_sections[0]._dtype
    output_rows = row_indices[-1][1]
    output_cols = resampled_sections[0].shape[2]
    median_image = np.empty((output_rows, output_cols), dtype) * np.nan

    for i, disk_arr in enumerate(resampled_sections):
        row1, row2 = row_indices[i]
        arr = disk_arr.read()
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore",
                                    message="All-NaN slice encountered",
                                    category=RuntimeWarning)
            median_image[row1:row2] = np.nanmedian(arr, axis=0)
        del arr, disk_arr

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
    save_blot=False,
    make_output_path=None,
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
    if save_blot:
        if make_output_path is None:
            raise ValueError("make_output_path must be provided if save_blot is True")
        model_path = make_output_path(input_model.meta.filename, suffix='blot')
        blot_model = _make_blot_model(input_model, blot)
        blot_model.meta.filename = model_path
        blot_model.save(model_path)
        log.info(f"Saved model in {model_path}")
        del blot_model
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

    # Make sure all data, error, and variance arrays have
    # matching NaNs and DQ flags
    match_nans_and_flags(input_model)


def flag_crs_in_models_with_resampling(
    input_models,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    save_blot=False,
    make_output_path=None,
):
    for image in input_models:
        flag_resampled_model_crs(image,
                                 median_data,
                                 median_wcs,
                                 snr1,
                                 snr2,
                                 scale1,
                                 scale2,
                                 backg,
                                 save_blot=save_blot,
                                 make_output_path=make_output_path)


def flag_model_crs(image, blot, snr):
    cr_mask = flag_crs(image.data, image.err, blot, snr)

    # Update dq array in-place
    image.dq |= cr_mask * np.uint32(DO_NOT_USE | OUTLIER)

    # Make sure all data, error, and variance arrays have
    # matching NaNs and DQ flags
    match_nans_and_flags(image)

    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")


def _make_blot_model(input_model, blot):
    blot_model = type(input_model)()
    blot_model.data = blot
    blot_model.update(input_model)
    return blot_model