"""
The ever-present utils sub-module. A home for all...
"""

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


# not inheriting from MutableSequence here as insert is complicated
class TempArrayHandler:
    def __init__(self, tempdir=""):
        self._temp_dir = tempfile.TemporaryDirectory(dir=tempdir)
        self._temp_path = Path(self._temp_dir.name)
        self._filenames = []
        self._data_shape = None
        self._data_dtype = None

    @property
    def closed(self):
        return not hasattr(self, "_temp_dir")

    def close(self):
        if self.closed:
            return
        self._temp_dir.cleanup()
        del self._temp_dir

    def __del__(self):
        self.close()

    def __len__(self):
        if self.closed:
            raise Exception("use after close")
        return len(self._filenames)

    def __getitem__(self, index):
        if self.closed:
            raise Exception("use after close")
        fn = self._filenames[index]
        return np.load(fn)

    def _validate_input(self, arr):
        if arr.ndim != 2:
            raise Exception(f"Only 2D arrays are supported: {arr.ndim}")
        if self._data_shape is None:
            self._data_shape = arr.shape
        else:
            if arr.shape != self._data_shape:
                raise Exception(
                    f"Input shape mismatch: {arr.shape} != {self._data_shape}"
                )
        if self._data_dtype is None:
            self._data_dtype = arr.dtype
        else:
            if arr.dtype != self._data_dtype:
                raise Exception(
                    f"Input dtype mismatch: {arr.dtype} != {self._data_dtype}"
                )

    def __setitem__(self, index, value):
        self._validate_input(value)
        if self.closed:
            raise Exception("use after close")
        fn = self._filenames[index]
        if fn is None:
            fn = self._temp_path / f"{index}.npy"
        np.save(fn, value, False)
        self._filenames[index] = fn

    def append(self, value):
        if self.closed:
            raise Exception("use after close")
        index = len(self)
        self._filenames.append(None)
        self.__setitem__(index, value)

    def median(self, buffer_size=100 << 20):
        if self.closed:
            raise Exception("use after close")
        if not len(self):
            raise Exception("can't take median of empty list")

        # figure out how big the buffer can be
        n_arrays = len(self)
        allowed_memory_per_array = buffer_size // n_arrays

        n_dim_1 = allowed_memory_per_array // (
            self._data_dtype.itemsize * self._data_shape[0]
        )
        if n_dim_1 < 1:
            # TODO more useful error message
            raise Exception("Not enough memory")
        if n_dim_1 >= self._data_shape[1]:
            return np.nanmedian(self, axis=0)

        buffer = np.empty(
            (n_arrays, self._data_shape[0], n_dim_1), dtype=self._data_dtype
        )
        median = np.empty(self._data_shape, dtype=self._data_dtype)

        e = n_dim_1
        slices = [slice(0, e)]
        while e <= self._data_shape[1]:
            s = e
            e += n_dim_1
            slices.append(slice(s, min(e, self._data_shape[1])))

        for s in slices:
            for i, arr in enumerate(self):
                buffer[i, :, : (s.stop - s.start)] = arr[:, s]
            median[:, s] = np.nanmedian(buffer[:, :, : (s.stop - s.start)], axis=0)
        return median


def create_cube_median(cube_model, maskpt):
    log.info("Computing median")

    weight_threshold = compute_weight_threshold(cube_model.wht, maskpt)

    median = np.nanmedian(
        np.ma.masked_array(cube_model.data, np.less(cube_model.wht, weight_threshold), fill_value=np.nan),
        axis=0,
    )
    return median


def create_median(resampled_models, maskpt, on_disk=True, buffer_size=10.0):
    """Create a median image from the singly resampled images.

    Parameters
    ----------
    resampled_models : ModelLibrary
        The singly resampled images.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    on_disk : bool
        If True, the input models are on disk and will be read in chunks.

    buffer_size : float
        The size of chunk in MB, per input model, that will be read into memory.
        This parameter has no effect if on_disk is False.

    Returns
    -------
    median_image : ndarray
        The median image.
    """
    # initialize storage for median computation
    if on_disk:
        # harder case: need special on-disk numpy array handling
        model_list = TempArrayHandler()
    else:
        # easier case: all models in library can be loaded into memory at once
        model_list = []

    # Compute the weight threshold for each input model
    with resampled_models:
        for resampled in resampled_models:
            weight_threshold = compute_weight_threshold(resampled.wht, maskpt)
            mask = np.less(resampled.wht, weight_threshold)
            resampled.data[mask] = np.nan
            model_list.append(resampled.data)
            # this is still modified if on_disk is False, but doesn't matter because never used again
            resampled_models.shelve(resampled, modify=False) 
            del resampled
    
    if not on_disk:
        return np.nanmedian(np.array(model_list), axis=0)
    else:
        # compute median using on-disk arrays
        median_data = model_list.median(buffer_size=buffer_size*_ONE_MB)
        model_list.close()
        return median_data


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