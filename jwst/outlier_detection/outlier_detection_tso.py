import numpy as np
from .outlier_detection import OutlierDetection, _remove_file, flag_cr
from jwst.resample.resample_utils import build_mask

from jwst import datamodels as dm
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetectionTSO"]


class OutlierDetectionTSO(OutlierDetection):
    """Class to flag outlier pixels in DQ of TSO data. Works similarly to
    imaging outlier detection, but does not resample and uses a rolling median."""

    def __init__(self, input_model: dm.CubeModel | dm.ModelContainer, **pars):
        """Initialize class for TSO data processing.

        Parameters
        ----------
        input_models : ~jwst.datamodels.CubeModel or ~jwst.datamodels.ModelContainer
            The input TSO data cube. if a ModelContainer is passed in, it is assumed
            to be a list of ImageModels, and gets converted to a CubeModel.

        """
        super().__init__(input_model, **pars)
        if isinstance(self.inputs, dm.ModelContainer):
            self._undo_convert_inputs()
        else:
            self.input_models = self.inputs
            self.converted = False


    def _undo_convert_inputs(self):
        """
        Convert input from ModelContainer of ImageModels with identical
        metadata (as is passed into outlier_detection for imaging) to a CubeModel.
        For typical use of calwebb_tso3, the input is already a CubeModel, so this
        is not needed. It is included mainly for backwards compatibility with older
        workflows that may have used a ModelContainer of ImageModels as input.

        This method converts `self.inputs` into a version of
        `self.input_models` suitable for processing by the class.
        """
        num_inputs = len(self.inputs)
        log.debug(f"Converting ModelContainer with {num_inputs} images back to CubeModel")
        shp = (num_inputs, *self.inputs[0].data.shape)
        cube = dm.CubeModel(data=np.empty(shp, dtype=np.float32),
                                    err=np.empty(shp, dtype=np.float32),
                                    dq=np.empty(shp, dtype=np.uint32),)
                                    #var_noise=np.empty(shp, dtype=np.float32),)
        
        for i, im in enumerate(self.inputs):
            cube.data[i] = im.data
            cube.err[i] = im.err
            cube.dq[i] = im.dq
            #cube.var_noise[i] = im.var_noise
            if i == 0:
                cube.meta = im.meta
            del im
        self.converted = True
        self.input_models = cube


    def do_detection(self):
        """Flag outlier pixels in DQ of input images."""
        weighted_cube = self.weight_no_resample()

        maskpt = self.outlierpars.get('maskpt', 0.7)
        weight_threshold = self.compute_weight_threshold(weighted_cube.wht, maskpt)

        rolling_width = self.outlierpars.get('n_ints', 25)
        if (rolling_width > 1) and (rolling_width < weighted_cube.shape[0]):
            medians = self.compute_rolling_median(weighted_cube, weight_threshold, w = rolling_width)

        else:
            medians = np.nanmedian(weighted_cube.data, axis=0)
            # this is a 2-D array, need to repeat it into the time axis
            # for consistent shape with rolling median case
            medians = np.broadcast_to(medians, weighted_cube.shape)

        # turn this into a cube model for saving and passing to flag_cr
        median_model = dm.CubeModel(data=medians)

        # Save median model if pars['save_intermediate_results'] is True
        # this will be a CubeModel with rolling median values.
        if self.outlierpars['save_intermediate_results']:
            with dm.open(weighted_cube[0]) as dm0:
                median_model.update(dm0)
            self.save_median(median_model) 

        # no need for blotting, resample is turned off for TSO
        # go straight to outlier detection
        log.info("Flagging outliers")
        flag_cr(self.input_models, median_model, **self.outlierpars)

        if self.converted:
            # Make sure actual input gets updated with new results
            for i in range(len(self.inputs)):
                self.inputs[i].dq = self.input_models.dq[i]
    
    
    def weight_no_resample(self):
        """
        give weights to model without resampling

        Notes
        -----
        Prior to PR #8473, the `build_driz_weight` function was used to
        create the weights for the input models for TSO data. However, that
        function was simply returning a copy of the DQ array because the 
        var_noise was not being passed in by calwebb_tso3. As of PR #8473,
        a cube model that includes the var_noise is passed into TSO 
        outlier detection, so `build_driz_weight` would weight the cube model
        by the variance. Therefore `build_driz_weight` was removed in order to
        preserve the original behavior. If it is determined later that exposure
        time or inverse variance weighting should be used here, build_driz_weight
        should be re-implemented.
        """
        weighted_cube = self.input_models.copy()
        dqmask = build_mask(self.input_models.dq, self.outlierpars['good_bits'])
        weighted_cube.wht = dqmask.astype(np.float32)
        return weighted_cube
    

    def compute_rolling_median(self, model: dm.CubeModel, weight_threshold: np.ndarray, w:int = 25) -> np.ndarray:
        '''
        Set bad and low-weight data to NaN, then compute the rolling median over the time axis.

        Parameters
        ----------
        model : ~jwst.datamodels.CubeModel
            The input cube model

        weight_threshold : np.ndarray
            The weight thresholds for each integration.

        w : int
            The window size for the rolling median. 

        Returns
        -------
        np.ndarray
            The rolling median of the input data. Same dimensions as input.
        '''

        sci = model.data
        weight = model.wht
        badmask = np.less(weight, weight_threshold)
        log.debug("Percentage of pixels with low weight: {}".format(
                np.sum(badmask) / weight.size * 100))

        # Fill resampled_sci array with nan's where mask values are True
        for f1, f2 in zip(sci, badmask):
            for elem1, elem2 in zip(f1, f2):
                elem1[elem2] = np.nan

        del badmask

        if w > sci.shape[0]:
            raise ValueError("Window size must be less than the number of integrations.")
        meds = moving_median_over_zeroth_axis(sci, w)

        del sci
        return meds
    

def moving_median_over_zeroth_axis(x: np.ndarray, w: int) -> np.ndarray:
    """
    Calculate the median of a moving window over the zeroth axis of an N-d array.
    Algorithm works by expanding the array into an additional dimension 
    where the new axis has the same length as the window size. Each entry in that
    axis is a copy of the original array shifted by 1 with respect to the previous
    entry, such that the rolling median is simply the median over the new axis.
    modified from https://stackoverflow.com/a/71154394, see link for more details.

    Parameters
    ----------
    x : np.ndarray
        The input array.

    w : int
        The window size.

    Returns
    -------
    np.ndarray
        The rolling median of the input array. Same dimensions as input.

    Notes
    -----
    Returns NaN for the first and last w//2 elements.
    """
    if w <= 1:
        raise ValueError("Rolling median window size must be greater than 1.")
    shifted = np.zeros((x.shape[0]+w-1, w, *x.shape[1:]))*np.nan
    for idx in range(w-1):
        shifted[idx:-w+idx+1, idx] = x
    shifted[idx+1:, idx+1] = x
    medians = np.median(shifted, axis=1)
    for idx in range(w-1):
        medians[idx] = np.median(shifted[idx, :idx+1])
        medians[-idx-1] = np.median(shifted[-idx-1, -idx-1:])
    medians = medians[(w-1)//2:-(w-1)//2]

    # Fill in the edges with the nearest valid value
    medians[:w//2] = medians[w//2]
    medians[-w//2:] = medians[-w//2]
    return medians
