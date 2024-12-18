import math
import numpy as np
import warnings

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels.dqflags import pixel

from jwst.assign_wcs.util import create_grism_bbox

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def subtract_wfss_bkg(
    input_model,
    bkg_filename,
    wl_range_name,
    mmag_extract=None,
    rescaler_kwargs={},
):
    """Scale and subtract a background reference image from WFSS/GRISM data.

    Parameters
    ----------
    input_model : JWST data model
        input target exposure data model

    bkg_filename : str
        name of master background file for WFSS/GRISM

    wl_range_name : str
        name of wavelengthrange reference file

    mmag_extract : float, optional, default None
        minimum abmag of grism objects to extract

    rescaler_kwargs : dict, optional, default {}
        Keyword arguments to pass to ScalingFactorComputer

    Returns
    -------
    result : JWST data model
        background-subtracted target data model
    """
    bkg_ref = datamodels.open(bkg_filename)
    
    # get the dispersion axis
    try:
        dispaxis = input_model.meta.wcsinfo.dispersion_direction
    except AttributeError:
        log.warning("Dispersion axis not found in input science image metadata. "
                    "Variance stopping criterion will have no effect for iterative "
                    "outlier rejection (will run until maxiter).")
        dispaxis = None
    rescaler_kwargs["dispersion_axis"] = dispaxis

    # get the source catalog for masking
    if hasattr(input_model.meta, "source_catalog"):
        got_catalog = True
    else:
        log.warning("No source_catalog found in input.meta.")
        got_catalog = False

    # Create a mask from the source catalog, True where there are no sources,
    # i.e. in regions we can use as background.
    if got_catalog:
        bkg_mask = _mask_from_source_cat(input_model, wl_range_name, mmag_extract)
        if not _sufficient_background_pixels(input_model.dq, bkg_mask):
            log.warning("Not enough background pixels to work with.")
            log.warning("Step will be SKIPPED.")
            return None
    else:
        bkg_mask = np.ones(input_model.data.shape, dtype=bool)

    # compute scaling factor for the reference background image
    log.info("Starting iterative outlier rejection for background subtraction.")
    rescaler = _ScalingFactorComputer(**rescaler_kwargs)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                category=RuntimeWarning,
                                message="All-NaN slice encountered")
        # copy to avoid propagating NaNs from iterative clipping into final product
        sci = input_model.data.copy()
        var = input_model.err.copy()**2
        bkg = bkg_ref.data.copy()
        factor, _ = rescaler(sci, bkg, var, mask=~bkg_mask)

    # extract the derived factor and apply it to the unmasked, non-outlier-rejected data
    subtract_this = factor * bkg_ref.data
    result = input_model.copy()
    result.data = input_model.data - subtract_this
    result.dq = np.bitwise_or(input_model.dq, bkg_ref.dq)

    log.info(f"Average of scaled background image = {np.nanmean(subtract_this):.3e}")
    log.info(f"Scaling factor = {factor:.5e}")

    bkg_ref.close()

    return result


class _ScalingFactorComputer:

    def __init__(self, p=1.0, maxiter=5, delta_rms_thresh=0, dispersion_axis=None):
        """
        Parameters
        ----------
        p: float, optional
            Percentile for sigma clipping on both low and high ends per iteration, default 1.0.
            For example, with p=2.0, the middle 96% of the data is kept.
        maxiter: int, optional
            Maximum number of iterations for outlier rejection. Default 5.
        delta_rms_thresh: float, optional
            Stopping criterion for outlier rejection; stops when the rms residuals
            change by less than this fractional threshold in a single iteration.
            For example, assuming delta_rms_thresh=0.1 and a residual RMS of 100
            in iteration 1, the iteration will stop if the RMS residual in iteration
            2 is greater than 90.
            Default 0.0, i.e., ignore this and only stop at maxiter.
        dispersion_axis: int, optional
            The index to select the along-dispersion axis. Used to compute the RMS
            residual, so must be set if rms_thresh > 0. Default None.
        """
        if (delta_rms_thresh > 0) and (dispersion_axis not in [1,2]):
            msg = (f"Unrecognized dispersion axis {dispersion_axis}. "
                    "Dispersion axis must be specified if delta_rms_thresh "
                    "is used as a stopping criterion.")
            raise ValueError(msg)

        self.p = p
        self.maxiter = maxiter
        self.delta_rms_thresh = delta_rms_thresh
        self.dispersion_axis = dispersion_axis


    def __call__(self, sci, bkg, var, mask=None):
        """
        Parameters
        ----------
        sci: ndarray
            The science data.
        bkg: ndarray
            The reference background model.
        var: ndarray
            Total variance (error squared) of the science data.
        mask: ndarray[bool], optional
            Initial mask to be applied to the data, True where bad.
            Typically this would mask out the real sources in the data.

        Returns
        -------
        float
            Scaling factor that minimizes sci - factor*bkg,
            taking into account residuals and outliers.
        ndarray[bool]
            Outlier mask generated by the iterative clipping procedure.
        """
        if mask is None:
            mask = np.zeros(sci.shape, dtype="bool")
        self._update_nans(sci, bkg, var, mask)

        # iteratively reject more and more outliers
        i = 0
        last_rms_resid = np.inf
        while (i < self.maxiter):

            # compute the factor that minimizes the residuals
            factor = self.err_weighted_mean(sci, bkg, var)
            sci_sub = sci-factor*bkg

            # Check fractional improvement stopping criterion before incrementing.
            # Note this never passes in iteration 0 because last_rms_resid is inf.
            if self.delta_rms_thresh > 0:
                rms_resid = self._compute_rms_residual(sci_sub)
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",
                                            category=RuntimeWarning,
                                            message="invalid value encountered in scalar divide")
                    fractional_diff = (last_rms_resid - rms_resid)/last_rms_resid
                if fractional_diff < self.delta_rms_thresh:
                    msg = (f"Stopping at iteration {i}; too little improvement "
                           "since last iteration (hit delta_rms_thresh).")
                    log.info(msg)
                    break
                last_rms_resid = rms_resid

            i += 1

            # Reject outliers based on residual between sci and bkg.
            # Updating the sci, var, and bkg nan values means that
            # they are ignored by nanpercentile in the next iteration
            limits = np.nanpercentile(sci_sub, (self.p, 100-self.p))
            mask += np.logical_or(sci_sub < limits[0], sci_sub > limits[1])
            self._update_nans(sci, bkg, var, mask)

        if i >= self.maxiter:
            log.info(f"Stopped at maxiter ({i}).")

        self._iters_run_last_call = i
        return self.err_weighted_mean(sci, bkg, var), mask


    def err_weighted_mean(self, sci, bkg, var):
        """Remove any var=0 values, which can happen for real data"""
        mask = (var == 0)
        self._update_nans(sci, bkg, var, mask)
        return np.nansum(sci*bkg/var) / np.nansum(bkg*bkg/var)


    def _update_nans(self, sci, bkg, var, mask):
        sci[mask] = np.nan
        bkg[mask] = np.nan
        var[mask] = np.nan


    def _compute_rms_residual(self, sci_sub):
        """
        Calculate the background-subtracted RMS along the dispersion axis, which is found
        by taking the median profile of the image collapsed along the cross-dispersion axis.
        Note meta.wcsinfo.dispersion_axis is 1-indexed coming out of assign_wcs, i.e., in [1,2].
        So we need to """
        collapsing_axis = int(self.dispersion_axis - 1)
        sci_sub_profile = np.nanmedian(sci_sub, axis=collapsing_axis)
        return np.sqrt(np.nanmean(sci_sub_profile**2))


def _sufficient_background_pixels(dq_array, bkg_mask, min_pixels=100):
    """Count number of good pixels for background use.

    Check DQ flags of pixels selected for bkg use - XOR the DQ values with
    the DO_NOT_USE flag to flip the DO_NOT_USE bit. Then count the number
    of pixels that AND with the DO_NOT_USE flag, i.e. initially did not have
    the DO_NOT_USE bit set.
    """
    return np.count_nonzero((dq_array[bkg_mask]
                            ^ pixel['DO_NOT_USE'])
                            & pixel['DO_NOT_USE']
                            ) > min_pixels


def _mask_from_source_cat(input_model, wl_range_name, mmag_extract=None):
    """Create a mask that is False within bounding boxes of sources.

    Parameters
    ----------
    input_model : JWST data model
        input target exposure data model

    wl_range_name : str
        Name of the wavelengthrange reference file

    mmag_extract : float
        Minimum abmag of grism objects to extract

    Returns
    -------
    bkg_mask : ndarray
        Boolean mask:  True for background, False for pixels that are
        inside at least one of the source regions defined in the source
        catalog.
    """

    shape = input_model.data.shape
    bkg_mask = np.ones(shape, dtype=bool)

    reference_files = {"wavelengthrange": wl_range_name}
    grism_obj_list = create_grism_bbox(input_model, reference_files, mmag_extract)

    for obj in grism_obj_list:
        order_bounding = obj.order_bounding
        for order in order_bounding.keys():
            ((ymin, ymax), (xmin, xmax)) = order_bounding[order]
            xmin = int(math.floor(xmin))
            xmax = int(math.ceil(xmax)) + 1  # convert to slice limit
            ymin = int(math.floor(ymin))
            ymax = int(math.ceil(ymax)) + 1
            xmin = max(xmin, 0)
            xmax = min(xmax, shape[-1])
            ymin = max(ymin, 0)
            ymax = min(ymax, shape[-2])
            bkg_mask[..., ymin:ymax, xmin:xmax] = False

    return bkg_mask
