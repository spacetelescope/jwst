import numpy as np
from .outlier_detection import OutlierDetection, _remove_file
from jwst.resample.resample_utils import build_driz_weight

from jwst import datamodels as dm
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetectionTSO"]


class OutlierDetectionTSO(OutlierDetection):
    """Class to flag outlier pixels in DQ of TSO data. Works similarly to
    imaging outlier detection, but does not resample and uses a rolling median."""

    def __init__(self, input_models: dm.ModelContainer, **pars):
        """Initialize class for TSO data processing.

        Parameters
        ----------
        input_models : ~jwst.datamodels.ModelContainer, str
            list of data models as ModelContainer or ASN file,
            one data model for each input 2-D ImageModel
        """
        super().__init__(input_models, **pars)
        self._convert_inputs()


    def do_detection(self):

        """Flag outlier pixels in DQ of input images."""
        drizzled_models = self.driz_weight_no_resample()

        maskpt = self.outlierpars.get('maskpt', 0.7)
        weight_thresholds = self.compute_weight_threshold(drizzled_models, maskpt)

        rolling_width = self.outlierpars.get('n_ints', 25)
        if (rolling_width > 1) and (rolling_width < len(drizzled_models)):
            medians = self.compute_rolling_median(drizzled_models, weight_thresholds, w = rolling_width)

        else:
            medians = super().create_median(drizzled_models)
            # this is a 2-D array, need to repeat it into the time axis
            # for consistent shape with rolling median case
            shp = (len(drizzled_models), medians.shape[0], medians.shape[1])
            medians = np.broadcast_to(medians, shp)

        # Save median model if pars['save_intermediate_results'] is True
        # this will be a CubeModel with rolling median values.
        if self.outlierpars['save_intermediate_results']:
            median_model = dm.CubeModel(data=medians)
            with dm.open(drizzled_models[0]) as dm0:
                median_model.update(dm0)
            self.save_median(median_model) 
        else:
            # since we're not saving intermediate results if the drizzled models
            # were written to disk, remove them
            if not self.outlierpars['in_memory']:
                for fn in drizzled_models._models:
                    _remove_file(fn)

        # no need for blotting, resample is turned off for TSO
        blot_models = dm.ModelContainer(open_models=False)
        for i in range(len(self.input_models)):
            median_model = dm.ImageModel(data=medians[i].data)
            blot_models.append(median_model)

        self.detect_outliers(blot_models)

        return
    
    def driz_weight_no_resample(self):
        """
        give weights to drizzled model without resampling
        """
        drizzled_models = self.input_models
        for i in range(len(self.input_models)):
            drizzled_models[i].wht = build_driz_weight(
                self.input_models[i],
                weight_type=self.outlierpars['weight_type'],
                good_bits=self.outlierpars['good_bits'])
        return drizzled_models

    def compute_rolling_median(self, models: dm.ModelContainer, weight_thresholds: list, w:int = 25) -> np.ndarray:
        '''
        Set bad and low-weight data to NaN, then compute the rolling median over the time axis.

        Parameters
        ----------
        models : ~jwst.datamodels.ModelContainer
            The input data models.

        weight_thresholds : list
            The weight thresholds for each integration.

        w : int
            The window size for the rolling median. 

        Returns
        -------
        np.ndarray
            The rolling median of the input data. Same dimensions as input.
        '''

        # Load model into memory as 3-D array
        sci = np.array([model.data for model in models])
        wht = np.array([model.wht for model in models])
        badmasks = []
        for weight, weight_threshold in zip(wht, weight_thresholds):
            badmask = np.less(weight, weight_threshold)
            log.debug("Percentage of pixels with low weight: {}".format(
                np.sum(badmask) / len(weight.flat) * 100))
            badmasks.append(badmask)

        # Fill resampled_sci array with nan's where mask values are True
        for f1, f2 in zip(sci, badmasks):
            for elem1, elem2 in zip(f1, f2):
                elem1[elem2] = np.nan

        del badmasks

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
