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


def _nanmedian3D(cube, overwrite_input=True):
    """Convenience function to compute the median of a cube ignoring warnings
    and setting memory-efficient flag"""
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore",
                                message="All-NaN slice encountered",
                                category=RuntimeWarning)
        return np.nanmedian(cube, axis=0, overwrite_input=overwrite_input)


def create_cube_median(cube_model, maskpt):
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)

    # not safe to use overwrite_input=True here because we are operating on model.data directly
    return _nanmedian3D(
        np.ma.masked_array(cube_model.data, np.less(cube_model.wht, weight_threshold), fill_value=np.nan),
        overwrite_input=False)


def create_median(resampled_models, maskpt, buffer_size=None):
    """Create a median image from the singly resampled images.

    Parameters
    ----------
    resampled_models : ModelLibrary
        The singly resampled images.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    Returns
    -------
    median_image : ndarray
        The median image.
    """
    on_disk = resampled_models._on_disk
    
    # Compute weight threshold for each input model and set NaN in data where weight is below threshold
    with resampled_models:
        for i, resampled in enumerate(resampled_models):
            weight_threshold = compute_weight_threshold(resampled.wht, maskpt)
            resampled.data[resampled.wht < weight_threshold] = np.nan
            if not on_disk:
                if i == 0:
                    model_cube = np.empty((len(resampled_models),) + resampled.data.shape,
                                         dtype=resampled.data.dtype)
                model_cube[i] = resampled.data
            else:
                if i == 0:
                    # set up temporary storage for spatial sections of all input models
                    shp = (len(resampled_models),) + resampled.data.shape
                    median_computer = OnDiskMedian(shp,
                                                   dtype=resampled.data.dtype,
                                                   buffer_size=buffer_size)
                # write spatial sections to disk
                median_computer.add_image(resampled.data)
            resampled_models.shelve(resampled, i, modify=False)
        del resampled
    
    if not on_disk:
        # easier case: all models in library can be loaded into memory at once
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore",
                                    message="All-NaN slice encountered",
                                    category=RuntimeWarning)
            return _nanmedian3D(model_cube)
    else:
        # retrieve spatial sections from disk one by one and compute timewise median
        med = median_computer.compute_median()
        median_computer.cleanup()
        return med


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
    
    def __init__(self, slice_shape, dtype, tempdir):
        """
        Parameters
        ----------
        slice_shape : tuple
            The shape of the spatial section of input data to be appended to the array.

        dtype : str
            The data type of the array. Must be a valid numpy array datatype.

        tempdir : str
            The directory in which to create the temporary files.
            Default is the current working directory.
        """
        if len(slice_shape) != 2:
            raise ValueError(f"Invalid slice_shape {slice_shape}. Only 2-D arrays are supported.")
        self._temp_file, self._filename = tempfile.mkstemp(dir=tempdir)
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
            f.write(data.tobytes())
        self._append_count += 1


    def read(self):
        """Read the 3-D array into memory, then delete the tempfile"""
        shp = (self._append_count,) + self._slice_shape
        with open(self._filename, "rb") as f:
            output = np.fromfile(f, dtype=self._dtype).reshape(shp)
        return output
    

    def cleanup(self):
        Path.unlink(self._filename)
        return


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
        self.buffer_size = buffer_size
        per_model_buffer_size = buffer_size / self._expected_nframes
        min_buffer_size = imcols * self.itemsize
        section_nrows = min(imrows, int(per_model_buffer_size // min_buffer_size))

        if section_nrows == 0:
            buffer_size = min_buffer_size * self._expected_nframes
            log.warning("Buffer size is too small to hold a single row."
                            f"Increasing buffer size to {buffer_size / _ONE_MB}MB")
            section_nrows = 1

        nsections = int(np.ceil(imrows / section_nrows))
        log.info(f"Computing median over {self._expected_nframes} images in {nsections}"
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
            arr = DiskAppendableArray(shp, dtype, self._temp_path)
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
        [arr.cleanup() for arr in self._temp_arrays]
        Path.rmdir(self._temp_path)
        return


    def compute_median(self):
        """Read spatial sections from disk and compute the median across groups
        (median over number of exposures on a per-pixel basis)"""
        row_indices = [(i * self.section_nrows, min((i+1) * self.section_nrows, self.frame_shape[0]))
                       for i in range(self.nsections)]

        output_rows = row_indices[-1][1]
        output_cols = self._temp_arrays[0].shape[2]
        median_image = np.empty((output_rows, output_cols), self.dtype) * np.nan

        for i, disk_arr in enumerate(self._temp_arrays):
            row1, row2 = row_indices[i]
            arr = disk_arr.read()
            median_image[row1:row2] = _nanmedian3D(arr)
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