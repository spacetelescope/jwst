"""
The ever-present utils sub-module. A home for all...
"""

import copy
from functools import partial
import warnings
import numpy as np
import tempfile
from pathlib import Path

from jwst.lib.pipe_utils import match_nans_and_flags
from jwst.resample.resample import compute_image_pixel_area
from jwst.resample.resample_utils import build_driz_weight
from stcal.outlier_detection.utils import compute_weight_threshold, gwcs_blot, flag_crs, flag_resampled_crs
from stdatamodels.jwst import datamodels
from . import _fileio

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
OUTLIER = datamodels.dqflags.pixel['OUTLIER']
_ONE_MB = 1 << 20


def nanmedian3D(cube, overwrite_input=True):
    """Compute the median of a cube ignoring warnings and with memory efficiency.
    np.nanmedian always uses at least 64-bit precision internally, and this is too
    memory-intensive. Instead, loop over the median calculation to avoid the
    memory usage of the internal upcasting and temporary array allocation.
    The additional runtime of this loop is indistinguishable from zero,
    but this loop cuts overall step memory usage roughly in half for at least one
    test association.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore",
                                message="All-NaN slice encountered",
                                category=RuntimeWarning)
        output_arr = np.empty(cube.shape[1:], dtype=np.float32)
        for i in range(output_arr.shape[0]):
            # this for loop looks silly, but see docstring above
            np.nanmedian(cube[:,i,:], axis=0, overwrite_input=overwrite_input, out=output_arr[i,:])
        return output_arr


def create_cube_median(cube_model, maskpt):
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)

    # not safe to use overwrite_input=True here because we are operating on model.data directly
    return nanmedian3D(
        np.ma.masked_array(cube_model.data, np.less(cube_model.wht, weight_threshold), fill_value=np.nan),
        overwrite_input=False)


def median_without_resampling(input_models,
                              maskpt,
                              weight_type,
                              good_bits,
                              save_intermediate_results=False,
                              make_output_path=None,
                              buffer_size=None):
    """
    Shared code between imaging and spec modes for resampling and median computation

    Parameters
    ----------
    input_models : ModelLibrary
        The input datamodels.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    weight_type : str
        The type of weighting to use when combining images. Options are:
        'ivm' (inverse variance) or 'exptime' (exposure time).

    good_bits : int
        The bit values that are considered good when determining the
        data quality of the input.

    save_intermediate_results : bool
        if True, save the drizzled models and median model to fits.

    make_output_path : function
        The functools.partial instance to pass to save_median. Must be 
        specified if save_intermediate_results is True. Default None.

    buffer_size : int
        The size of chunk in bytes that will be read into memory when computing the median.
        This parameter has no effect if the input library has its on_disk attribute
        set to False.
    """
    in_memory = not input_models._on_disk
    ngroups = len(input_models)

    with input_models:
        for i in range(len(input_models)):

            drizzled_model = input_models.borrow(i)
            drizzled_model.wht = build_driz_weight(drizzled_model,
                                                    weight_type=weight_type,
                                                    good_bits=good_bits)
            median_wcs = copy.deepcopy(drizzled_model.meta.wcs)
            input_models.shelve(drizzled_model, i, modify=True)

            if save_intermediate_results:
                # write the drizzled model to file
                _fileio.save_drizzled(drizzled_model, make_output_path)

            if i == 0:
                input_shape = (ngroups,)+drizzled_model.data.shape
                dtype = drizzled_model.data.dtype
                median_computer = _make_median_computer(input_shape, in_memory, buffer_size, dtype)

            weight_threshold = compute_weight_threshold(drizzled_model.wht, maskpt)
            drizzled_model.data[drizzled_model.wht < weight_threshold] = np.nan
            _append_to_median_computer(median_computer, i, drizzled_model.data, in_memory)

    # Perform median combination on set of drizzled mosaics
    median_data = _evaluate_median_computer(median_computer, in_memory)

    if save_intermediate_results:
        # Save median model to fits
        median_model = datamodels.ImageModel(median_data)
        median_model.update(drizzled_model)
        median_model.meta.wcs = median_wcs
        _fileio.save_median(median_model, make_output_path)
    del drizzled_model

    return median_data, median_wcs


def median_with_resampling(input_models,
                           resamp,
                           maskpt,
                           save_intermediate_results=False,
                           make_output_path=None,
                           buffer_size=None):
    """
    Shared code between imaging and spec modes for resampling and median computation

    Parameters
    ----------
    input_models : ModelLibrary
        The input datamodels.

    resamp : resample.resample.ResampleData object
        The controlling object for the resampling process.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    save_intermediate_results : bool
        if True, save the drizzled models and median model to fits.

    make_output_path : function
        The functools.partial instance to pass to save_median. Must be 
        specified if save_intermediate_results is True. Default None.

    buffer_size : int
        The size of chunk in bytes that will be read into memory when computing the median.
        This parameter has no effect if the input library has its on_disk attribute
        set to False.
    """
    if not resamp.single:
        raise ValueError("median_with_resampling should only be used for resample_many_to_many")
    
    in_memory = not input_models._on_disk
    indices_by_group = list(input_models.group_indices.values())
    ngroups = len(indices_by_group)

    with input_models:
        for i, indices in enumerate(indices_by_group):

            median_wcs = resamp.output_wcs
            drizzled_model = resamp.resample_group(input_models, indices)

            if save_intermediate_results:
                # write the drizzled model to file
                _fileio.save_drizzled(drizzled_model, make_output_path)

            if i == 0:
                input_shape = (ngroups,)+drizzled_model.data.shape
                dtype = drizzled_model.data.dtype
                median_computer = _make_median_computer(input_shape, in_memory, buffer_size, dtype)

            weight_threshold = compute_weight_threshold(drizzled_model.wht, maskpt)
            drizzled_model.data[drizzled_model.wht < weight_threshold] = np.nan
            _append_to_median_computer(median_computer, i, drizzled_model.data, in_memory)


    # Perform median combination on set of drizzled mosaics
    median_data = _evaluate_median_computer(median_computer, in_memory)

    if save_intermediate_results:
        # Save median model to fits
        median_model = datamodels.ImageModel(median_data)
        median_model.update(drizzled_model)
        median_model.meta.wcs = median_wcs
        # drizzled model already contains asn_id
        make_output_path = partial(make_output_path, asn_id=None)
        _fileio.save_median(median_model, make_output_path)
    del drizzled_model

    return median_data, median_wcs


def _make_median_computer(full_shape, in_memory, buffer_size, dtype):

    if in_memory:
        # allocate memory for data arrays that go into median
        median_computer = np.empty(full_shape, dtype=np.float32)
    else:
        # set up temporary storage for data arrays that go into median
        median_computer = OnDiskMedian(full_shape,
                                        dtype=dtype,
                                        buffer_size=buffer_size)
    return median_computer


def _append_to_median_computer(median_computer, idx, data, in_memory):
    if in_memory:
        # populate pre-allocated memory with the drizzled data
        median_computer[idx] = data
    else:
        # distribute the drizzled data into the temporary storage
        median_computer.add_image(data)


def _evaluate_median_computer(median_computer, in_memory):
    if in_memory:
        median_data = nanmedian3D(median_computer)
        del median_computer
    else:
        median_data = median_computer.compute_median()
        median_computer.cleanup()
    return median_data


class DiskAppendableArray:
    """
    Creates a temporary file to which to append data, in order to perform
    timewise operations on a stack of input images without holding all of them
    in memory.
    
    This class is purpose-built for the median computation during outlier detection 
    and is not very flexible. It is assumed that each data array passed to `append`
    represents the same spatial segment of the full dataset. It is also assumed that 
    each data array passed to `append` represents only a single instant in time; 
    the append operation will stack them along a new axis.
    
    The `read` operation is only capable of loading the full array back into memory.
    When working with large datasets that do not fit in memory, the
    required workflow is to create many DiskAppendableArray objects, each holding
    a small spatial segment of the full dataset.
    """
    
    def __init__(self, slice_shape, dtype, filename):
        """
        Parameters
        ----------
        slice_shape : tuple
            The shape of the spatial section of input data to be appended to the array.

        dtype : str
            The data type of the array. Must be a valid numpy array datatype.

        filename : str
            The full file path in which to store the array
        """
        if len(slice_shape) != 2:
            raise ValueError(f"Invalid slice_shape {slice_shape}. Only 2-D arrays are supported.")
        self._filename = Path(filename)
        with open(filename, "wb") as f:   # noqa: F841
            pass
        self._slice_shape = slice_shape
        self._dtype = np.dtype(dtype)
        self._append_count = 0


    @property
    def shape(self):
        return (self._append_count,) + self._slice_shape


    def append(self, data):
        """Add a new slice to the temporary file"""
        if data.shape != self._slice_shape:
            raise ValueError(f"Data shape {data.shape} does not match slice shape {self._slice_shape}")
        if data.dtype != self._dtype:
            raise ValueError(f"Data dtype {data.dtype} does not match array dtype {self._dtype}")
        with open(self._filename, "ab") as f:
            data.tofile(f, sep="")
        self._append_count += 1


    def read(self):
        """Read the 3-D array into memory"""
        shp = (self._append_count,) + self._slice_shape
        with open(self._filename, "rb") as f:
            output = np.fromfile(f, dtype=self._dtype).reshape(shp)
        return output


class OnDiskMedian:

    def __init__(self, shape, dtype='float32', tempdir="", buffer_size=None):
        """
        Set up temporary files to perform operations on a stack of 2-D input arrays
        along the stacking axis (e.g., a time axis) without
        holding all of them in memory. Currently the only supported operation
        is the median.

        Parameters
        ----------
        shape: tuple
            The shape of the entire input, (n_images, imrows, imcols).

        dtype : str
            The data type of the input data.

        tempdir : str
            The parent directory in which to create the temporary directory,
            which itself holds all the DiskAppendableArray tempfiles.
            Default is the current working directory.

        buffer_size : int, optional
            The buffer size, units of bytes.
            Default is the size of one input image.
        """
        if len(shape) != 3:
            raise ValueError(f"Invalid input shape {shape}; only three-dimensional data are supported.")
        self._expected_nframes = shape[0]
        self.frame_shape = shape[1:]
        self.dtype = np.dtype(dtype)
        self.itemsize = self.dtype.itemsize
        self._temp_dir = tempfile.TemporaryDirectory(dir=tempdir)
        self._temp_path = Path(self._temp_dir.name)

        # figure out number of sections and rows per section that are needed
        self.nsections, self.section_nrows = self._get_buffer_indices(buffer_size=buffer_size)
        self.slice_shape = (self.section_nrows, shape[2])
        self._n_adds = 0

        # instantiate a temporary DiskAppendableArray for each section
        self._temp_arrays = self._temparray_setup(dtype)


    def _get_buffer_indices(self, buffer_size=None):
        """
        Parameters
        ----------
        buffer_size : int, optional
            The buffer size for the median computation, units of bytes.

        Returns
        -------
        nsections : int
            The number of sections to divide the input data into.

        section_nrows : int
            The number of rows in each section (except the last one).
        """
        imrows, imcols = self.frame_shape
        if buffer_size is None:
            buffer_size = imrows * imcols * self.itemsize
        per_model_buffer_size = buffer_size / self._expected_nframes
        min_buffer_size = imcols * self.itemsize
        section_nrows = min(imrows, int(per_model_buffer_size // min_buffer_size))

        if section_nrows <= 0:
            buffer_size = min_buffer_size * self._expected_nframes
            log.warning("Buffer size is too small to hold a single row."
                            f"Increasing buffer size to {buffer_size / _ONE_MB}MB")
            section_nrows = 1
        self.buffer_size = buffer_size

        nsections = int(np.ceil(imrows / section_nrows))
        log.info(f"Computing median over {self._expected_nframes} groups in {nsections} "
                    f"sections with total memory buffer {buffer_size / _ONE_MB} MB")
        return nsections, section_nrows


    def _temparray_setup(self, dtype):
        """Set up temp file handlers for each spatial section"""
        temp_arrays = []
        for i in range(0, self.nsections):
            shp = self.slice_shape
            if i == self.nsections - 1:
                # last section has whatever shape is left over
                shp = (self.frame_shape[0] - (self.nsections-1) * self.section_nrows, self.frame_shape[1])
            arr = DiskAppendableArray(shp, dtype, self._temp_path / f"{i}.bin")
            temp_arrays.append(arr)
        return temp_arrays


    def add_image(self, data):
        """Split resampled model data into spatial sections and write to disk."""
        if self._n_adds >= self.nsections:
            raise IndexError(f"Too many calls to add_image. Expected at most {self.nsections} input models.")
        self._validate_data(data)
        self._n_adds += 1
        for i in range(self.nsections):
            row1 = i * self.section_nrows
            row2 = min(row1 + self.section_nrows, self.frame_shape[0])
            arr = self._temp_arrays[i]
            arr.append(data[row1:row2])


    def _validate_data(self, data):
        if data.shape != self.frame_shape:
            raise ValueError(f"Data shape {data.shape} does not match expected shape {self.frame_shape}")
        if data.dtype != self.dtype:
            raise ValueError(f"Data dtype {data.dtype} does not match expected dtype {self.dtype}")
        

    def cleanup(self):
        """Remove the temporary files and directory when finished"""
        self._temp_dir.cleanup()
        return


    def compute_median(self):
        """Read spatial sections from disk and compute the median across groups
        (median over number of exposures on a per-pixel basis)"""
        row_indices = [(i * self.section_nrows, min((i+1) * self.section_nrows, self.frame_shape[0]))
                       for i in range(self.nsections)]

        output_rows = row_indices[-1][1]
        output_cols = self._temp_arrays[0].shape[2]
        median_image = np.full((output_rows, output_cols), np.nan, dtype=self.dtype)

        for i, disk_arr in enumerate(self._temp_arrays):
            row1, row2 = row_indices[i]
            arr = disk_arr.read()
            median_image[row1:row2] = nanmedian3D(arr)
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
        _fileio.save_blot(input_model, blot, make_output_path)
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
