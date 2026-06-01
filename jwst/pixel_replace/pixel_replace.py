import logging
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import nirspec

log = logging.getLogger(__name__)

__all__ = ["PixelReplaceArrays", "PixelReplacement"]


@dataclass
class PixelReplaceArrays:
    """
    Container for data arrays and dispersion direction.

    Algorithms operate on this dataclass rather than on a
    `~stdatamodels.jwst.datamodels.JwstDataModel`.
    This avoids the overhead of constructing intermediate DataModel objects,
    which was slowing runtime for TSO data with thousands of integrations,
    and provides a consistent interface for :meth:`PixelReplacement.mingrad`
    and :meth:`PixelReplacement.fit_profile`.
    """

    data: np.ndarray
    """Science array."""

    dq: np.ndarray
    """Data quality array."""

    err: np.ndarray
    """Total error array."""

    var_poisson: np.ndarray | None
    """Poisson variance array."""

    var_rnoise: np.ndarray | None
    """Read-noise variance array."""

    var_flat: np.ndarray | None
    """Flat-field variance array."""

    dispersion_direction: int
    """Dispersion direction."""


class PixelReplacement:
    """
    Main class for performing pixel replacement.

    This class controls loading the input data model, selecting the
    method for pixel replacement, and executing each step. This class
    should provide modularization to allow for multiple options and possible
    future reference files.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Datamodel with bad pixels to replace. Updated in-place.

    **pars
        Optional parameters to modify how pixel replacement
        will execute.
    """

    # Shortcuts for DQ Flags
    DO_NOT_USE = datamodels.dqflags.pixel["DO_NOT_USE"]
    FLUX_ESTIMATED = datamodels.dqflags.pixel["FLUX_ESTIMATED"]
    NON_SCIENCE = datamodels.dqflags.pixel["NON_SCIENCE"]

    # Shortcuts for dispersion direction for ease of reading
    HORIZONTAL = 1
    VERTICAL = 2
    LOG_SLICE = ["column", "row"]

    def __init__(self, input_model, **pars):
        self.input = input_model
        self.pars = {}
        self.pars.update(pars)
        # Store algorithm options here.
        self.algorithm_dict = {
            "fit_profile": self.fit_profile,
            "mingrad": self.mingrad,
        }

        # Choose algorithm from dict using input par.
        try:
            self.algorithm = self.algorithm_dict[self.pars["algorithm"]]

        except KeyError as err:
            log.critical(
                f"Algorithm name {self.pars['algorithm']} provided does "
                "not match an implemented algorithm!"
            )
            raise KeyError from err

    @staticmethod
    def _arrays_from_model(model):
        """Extract PixelReplaceArrays from DataModel, copying arrays."""  # numpydoc ignore: RT01
        return PixelReplaceArrays(
            data=model.data.copy(),
            dq=model.dq.copy(),
            err=model.err.copy(),
            var_poisson=model.var_poisson.copy(),
            var_rnoise=model.var_rnoise.copy(),
            var_flat=model.var_flat.copy(),
            dispersion_direction=model.meta.wcsinfo.dispersion_direction,
        )

    @staticmethod
    def _model_from_arrays(arrays, model):
        """Write PixelReplaceArrays back into a DataModel in place."""  # numpydoc ignore: RT01
        model.data = arrays.data
        model.dq = arrays.dq
        model.err = arrays.err
        model.var_poisson = arrays.var_poisson
        model.var_rnoise = arrays.var_rnoise
        model.var_flat = arrays.var_flat

    def replace(self):
        """
        Unpack model and apply pixel replacement algorithm.

        Process the input `~stdatamodels.jwst.datamodels.JwstDataModel`,
        unpack any model that holds
        more than one 2D spectrum, then apply selected algorithm
        to each 2D spectrum in input.
        """
        # ImageModel inputs (MIR_LRS-FIXEDSLIT)
        # or 2D SlitModel inputs (e.g. NRS_FIXEDSLIT in spec3)
        if isinstance(self.input, datamodels.ImageModel) or (
            isinstance(self.input, datamodels.SlitModel) and self.input.data.ndim == 2
        ):
            arrays = self._arrays_from_model(self.input)
            arrays = self.algorithm(arrays)
            self._model_from_arrays(arrays, self.input)
            n_replaced = np.count_nonzero(self.input.dq & self.FLUX_ESTIMATED)
            log.info(f"Input model had {n_replaced} pixels replaced.")
        elif isinstance(self.input, datamodels.IFUImageModel):
            # Attempt to run pixel replacement on each throw of the IFU slicer
            # individually.
            xx, yy = np.indices(self.input.data.shape)
            if self.input.meta.exposure.type == "MIR_MRS":
                if self.pars["algorithm"] == "mingrad":
                    # mingrad method
                    arrays = self._arrays_from_model(self.input)
                    arrays = self.algorithm(arrays)
                    self._model_from_arrays(arrays, self.input)
                else:
                    # fit_profile method
                    det2ab = self.input.meta.wcs.get_transform(
                        self.input.meta.wcs.available_frames[0], "alpha_beta"
                    )
                    _, beta_array, _ = det2ab(yy, xx)
                    unique_beta = np.unique(beta_array)
                    unique_beta = unique_beta[~np.isnan(unique_beta)]
                    for i, beta in enumerate(unique_beta):
                        # Define a mask that is True where this trace is located
                        trace_mask = beta_array == beta
                        arrays = self._arrays_from_model(self.input)
                        arrays.dq = np.where(
                            # When not in this trace, set NON_SCIENCE and DO_NOT_USE
                            ~trace_mask,
                            arrays.dq | self.DO_NOT_USE | self.NON_SCIENCE,
                            arrays.dq,
                        )

                        arrays = self.algorithm(arrays)
                        self.input.data = np.where(trace_mask, arrays.data, self.input.data)
                        self.input.dq = np.where(trace_mask, arrays.dq, self.input.dq)
                        self.input.err = np.where(trace_mask, arrays.err, self.input.err)
                        self.input.var_poisson = np.where(
                            trace_mask, arrays.var_poisson, self.input.var_poisson
                        )
                        self.input.var_rnoise = np.where(
                            trace_mask, arrays.var_rnoise, self.input.var_rnoise
                        )
                        self.input.var_flat = np.where(
                            trace_mask, arrays.var_flat, self.input.var_flat
                        )

                        n_replaced = np.count_nonzero(arrays.dq & self.FLUX_ESTIMATED)
                        log.info(
                            f"Input MRS frame had {n_replaced} pixels replaced "
                            f"in IFU slice {i + 1}."
                        )

                n_replaced = np.count_nonzero(self.input.dq & self.FLUX_ESTIMATED)
                log.info(f"Input MRS frame had {n_replaced} total pixels replaced.")
            else:
                if self.pars["algorithm"] == "mingrad":
                    # mingrad method
                    arrays = self._arrays_from_model(self.input)
                    arrays = self.algorithm(arrays)
                    self._model_from_arrays(arrays, self.input)
                else:
                    # fit_profile method - iterate over IFU slices
                    for i in range(30):
                        slice_wcs = nirspec.nrs_wcs_set_input(self.input, i)
                        det2slicer = slice_wcs.get_transform(
                            self.input.meta.wcs.available_frames[0], "slicer"
                        )
                        _, _, wave = det2slicer(yy, xx)

                        # Define a mask that is True where this trace is located
                        trace_mask = wave > 0
                        arrays = self._arrays_from_model(self.input)
                        arrays.dq = np.where(
                            # When not in this trace, set NON_SCIENCE and DO_NOT_USE
                            ~trace_mask,
                            arrays.dq | self.DO_NOT_USE | self.NON_SCIENCE,
                            arrays.dq,
                        )

                        arrays = self.algorithm(arrays)
                        self.input.data = np.where(trace_mask, arrays.data, self.input.data)
                        self.input.dq = np.where(trace_mask, arrays.dq, self.input.dq)
                        self.input.err = np.where(trace_mask, arrays.err, self.input.err)
                        self.input.var_poisson = np.where(
                            trace_mask, arrays.var_poisson, self.input.var_poisson
                        )
                        self.input.var_rnoise = np.where(
                            trace_mask, arrays.var_rnoise, self.input.var_rnoise
                        )
                        self.input.var_flat = np.where(
                            trace_mask, arrays.var_flat, self.input.var_flat
                        )

                        n_replaced = np.count_nonzero(arrays.dq & self.FLUX_ESTIMATED)
                        log.info(
                            f"Input NRS_IFU frame had {n_replaced} pixels "
                            f"replaced in IFU slice {i + 1}."
                        )

                n_replaced = np.count_nonzero(self.input.dq & self.FLUX_ESTIMATED)
                log.info(f"Input NRS_IFU frame had {n_replaced} total pixels replaced.")

        # MultiSlitModel inputs (WFSS, NRS_FIXEDSLIT, ?)
        elif isinstance(self.input, datamodels.MultiSlitModel):
            for i, _slit in enumerate(self.input.slits):
                slit_model = datamodels.SlitModel(self.input.slits[i].instance)
                arrays = self._arrays_from_model(slit_model)
                slit_model.close()

                arrays = self.algorithm(arrays)

                n_replaced = np.count_nonzero(arrays.dq & self.FLUX_ESTIMATED)
                log.info(f"Slit {i} had {n_replaced} pixels replaced.")

                self._model_from_arrays(arrays, self.input.slits[i])

        # CubeModel inputs are TSO (so far?); SlitModel may be NRS_BRIGHTOBJ,
        # also requiring a re-packaging of the data into 2D inputs for the algorithm
        elif isinstance(self.input, datamodels.CubeModel | datamodels.SlitModel):
            dispaxis = self.input.meta.wcsinfo.dispersion_direction

            for i in range(len(self.input.data)):
                # Ensure variance arrays exist
                var_dict = {
                    "var_poisson": None,
                    "var_rnoise": None,
                    "var_flat": None,
                }
                for key in var_dict.keys():
                    if self.input[key] is not None:
                        var_dict[key] = self.input[key][i].copy()

                arrays = PixelReplaceArrays(
                    data=self.input.data[i].copy(),
                    dq=self.input.dq[i].copy(),
                    err=self.input.err[i].copy(),
                    var_poisson=var_dict["var_poisson"],
                    var_rnoise=var_dict["var_rnoise"],
                    var_flat=var_dict["var_flat"],
                    dispersion_direction=dispaxis,
                )
                arrays = self.algorithm(arrays)
                n_replaced = np.count_nonzero(arrays.dq & self.FLUX_ESTIMATED)
                log.info(f"Input TSO integration {i} had {n_replaced} pixels replaced.")

                self.input.data[i] = arrays.data
                self.input.dq[i] = arrays.dq
                self.input.err[i] = arrays.err
                for key in var_dict.keys():
                    if self.input[key] is not None:
                        self.input[key][i] = getattr(arrays, key)

        else:
            # This should never happen, as these should be caught in the step code.
            log.critical(
                "Pixel replacement code did not filter this input correctly - skipping step."
            )
            return

    def fit_profile(self, arrays):
        """
        Replace pixels with the profile fit method.

        Fit a profile to adjacent columns, scale profile to
        column with missing pixel(s), and find flux estimate
        from scaled profile.

        Error and variance values for the replaced pixels
        are similarly estimated, using the scales from the
        profile fit to the data.

        Parameters
        ----------
        arrays : `PixelReplaceArrays`
            Pixel arrays and dispersion direction for the 2D spectrum to process.
            Arrays are modified in place.

        Returns
        -------
        arrays : `PixelReplaceArrays`
            The input with bad pixels now flagged with FLUX_ESTIMATED
            and holding a flux value estimated from the spatial profile.
        """
        # np.nanmedian() entry full of NaN values would produce a numpy
        # warning (despite well-defined behavior - return a NaN)
        # so we suppress that here.
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")

        dispaxis = arrays.dispersion_direction

        # Make a copy of the input DQ, before replacement
        input_dq = arrays.dq.copy()

        # Truncate array to region where good pixels exist
        good_pixels = np.where(~input_dq & self.DO_NOT_USE)
        if np.any(0 in np.shape(good_pixels)):
            log.warning(
                "No good pixels in at least one dimension of "
                "data array - skipping pixel replacement."
            )
            return arrays
        x_range = [np.min(good_pixels[0]), np.max(good_pixels[0]) + 1]
        y_range = [np.min(good_pixels[1]), np.max(good_pixels[1]) + 1]

        valid_shape = [x_range, y_range]
        profile_cut = valid_shape[dispaxis - 1]

        # COMMENTS NOTE:
        # In comments and parameter naming, I will try to be consistent in using
        # "profile" to describe vectors in the spatial, i.e. cross-dispersion direction,
        # and "slice" to describe vectors in the spectral, i.e. dispersion direction.

        # Create set of slice indices which we can later use for profile creation
        valid_profiles = set(range(*valid_shape[2 - dispaxis]))
        profiles_to_replace = set()

        # Loop over axis of data array corresponding to cross-
        # dispersion direction by indexing data shape with
        # strange dispaxis argument. Keep indices in full-frame numbering scheme,
        # but only iterate through slices with valid data.
        for ind in range(*valid_shape[2 - dispaxis]):
            # Exclude regions with no data for dq slice.
            dq_slice = input_dq[self.custom_slice(dispaxis, ind)][profile_cut[0] : profile_cut[1]]
            # Exclude regions with NON_SCIENCE flag
            dq_slice = np.where(dq_slice & self.NON_SCIENCE, self.NON_SCIENCE, dq_slice)
            # Find bad pixels in region containing valid data.
            n_bad = np.count_nonzero(dq_slice & self.DO_NOT_USE)
            n_nonscience = np.count_nonzero(dq_slice & self.NON_SCIENCE)
            if n_bad + n_nonscience == len(dq_slice):
                log.debug(f"Slice {ind} contains no good pixels. Skipping replacement.")
                valid_profiles.discard(ind)
            elif n_bad == 0:
                log.debug(f"Slice {ind} contains no bad pixels.")
            else:
                log.debug(f"Slice {ind} contains {n_bad} bad pixels.")
                profiles_to_replace.add(ind)

        log.debug(f"Number of profiles with at least one bad pixel: {len(profiles_to_replace)}")
        for ind in profiles_to_replace:
            # Use sets for convenient finding of neighboring slices to use in profile creation
            adjacent_inds = set(
                range(ind - self.pars["n_adjacent_cols"], ind + self.pars["n_adjacent_cols"] + 1)
            )
            adjacent_inds.discard(ind)
            valid_adjacent_inds = list(adjacent_inds.intersection(valid_profiles))

            # Cut out valid neighboring profiles
            adjacent_condition = self.custom_slice(dispaxis, valid_adjacent_inds)
            profile_data = arrays.data[adjacent_condition]
            profile_err = arrays.err[adjacent_condition]
            if profile_data.size == 0:
                log.info(
                    f"Profile in {self.LOG_SLICE[dispaxis - 1]} {ind} "
                    f"has no valid adjacent values - skipping."
                )
                continue

            # Mask out bad pixels
            invalid_condition = (input_dq[adjacent_condition] & self.DO_NOT_USE).astype(bool)
            profile_data[invalid_condition] = np.nan
            profile_err[invalid_condition] = np.nan

            # Add additional cut to pull only from region with valid data
            # for convenience (may not be necessary)
            region_condition = self.custom_slice(3 - dispaxis, range(*profile_cut))
            profile_data = profile_data[region_condition]
            profile_snr = np.abs(profile_data / profile_err[region_condition])

            # Normalize profile data
            # TODO: check on signs here - absolute max sometimes picks up
            #  large negative outliers
            profile_norm_scale = np.nanmax(np.abs(profile_data), axis=(dispaxis - 1), keepdims=True)
            # If profile data has SNR < 5 everywhere just use unity scaling
            # (so we don't normalize to noise)
            if np.nanmax(profile_snr) < 5:
                profile_norm_scale[:] = 1.0
            normalized = profile_data / profile_norm_scale

            # Get corresponding error and variance data and scale and mask to match
            # Handle the variance arrays as errors, so the scales match.
            err_names = ["err", "var_poisson", "var_rnoise", "var_flat"]
            norm_errors = {}
            for err_name in err_names:
                if err_name.startswith("var"):
                    if (err_arr := getattr(arrays, err_name)) is None:
                        continue
                    err = np.sqrt(err_arr)
                else:
                    err = getattr(arrays, err_name)
                norm_err = err[adjacent_condition]
                norm_err[invalid_condition] = np.nan
                norm_errors[err_name] = norm_err[region_condition] / profile_norm_scale

            # Pull median for each pixel across profile.
            # Profile entry full of NaN values would produce a numpy
            # warning (despite well-defined behavior - return a NaN)
            # so we suppress that above.
            median_profile = np.nanmedian(normalized, axis=(2 - dispaxis))

            # Do the same for the errors
            for err_name in norm_errors:
                norm_errors[err_name] = np.nanmedian(norm_errors[err_name], axis=(2 - dispaxis))

            # Clean current profile of values flagged as bad
            current_condition = self.custom_slice(dispaxis, ind)
            current_profile = arrays.data[current_condition]
            cleaned_current = np.where(
                input_dq[current_condition] & self.DO_NOT_USE, np.nan, current_profile
            )[range(*profile_cut)]

            replace_mask = np.where(~np.isnan(cleaned_current))[0]
            if len(replace_mask) == 0:
                log.info(
                    f"Profile in {self.LOG_SLICE[dispaxis - 1]} {ind} "
                    f"has no valid values - skipping."
                )
                continue
            min_median = median_profile[replace_mask]
            min_current = cleaned_current[replace_mask]
            norm_current = min_current / np.max(min_current)

            # Scale median profile to current profile with bad pixel - minimize mse?
            # Only do this scaling if we didn't default to all-unity scaling above,
            # and require input values below 1e20 so that we don't overflow the
            # minimization routine with extremely bad noise.
            if (
                (np.nanmedian(profile_norm_scale) != 1.0)
                & (np.nanmax(np.abs(min_median)) < 1e20)
                & (np.nanmax(np.abs(norm_current)) < 1e20)
            ):
                # TODO: check on signs here - absolute max sometimes picks up
                #  large negative outliers
                norm_scale = minimize(
                    self.profile_mse,
                    x0=np.abs(np.nanmax(norm_current)),
                    args=(np.abs(min_median), np.abs(norm_current)),
                    method="Nelder-Mead",
                ).x
                scale = np.max(min_current)
            else:
                norm_scale = 1.0
                scale = 1.0

            # Replace pixels that are do-not-use but not non-science
            current_dq = input_dq[current_condition][range(*profile_cut)]
            replace_condition = (current_dq & self.DO_NOT_USE ^ current_dq & self.NON_SCIENCE) == 1
            replaced_current = np.where(
                replace_condition, median_profile * norm_scale * scale, cleaned_current
            )

            # Change the dq bits where old flag was DO_NOT_USE and new value is not nan
            replaced_dq = np.where(
                replace_condition & ~(np.isnan(replaced_current)),
                current_dq ^ self.DO_NOT_USE ^ self.FLUX_ESTIMATED,
                current_dq,
            )

            # Update data and DQ in the output model
            arrays.data[current_condition][range(*profile_cut)] = replaced_current
            arrays.dq[current_condition][range(*profile_cut)] = replaced_dq

            # Also update the errors and variances
            current_err = arrays.err[current_condition][range(*profile_cut)]
            replaced_err = np.where(
                replace_condition, norm_errors["err"] * norm_scale * scale, current_err
            )
            arrays.err[current_condition][range(*profile_cut)] = replaced_err

            # Some values in NIRSpec variances may overflow in the squares - ignore the warning.
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "overflow encountered", RuntimeWarning)
                for var in ["var_poisson", "var_rnoise", "var_flat"]:
                    if (var_arr := getattr(arrays, var)) is not None:
                        current_var = var_arr[current_condition][range(*profile_cut)]
                        replaced_var = np.where(
                            replace_condition,
                            (norm_errors[var] * norm_scale * scale) ** 2,
                            current_var,
                        )
                        var_arr[current_condition][range(*profile_cut)] = replaced_var
                        setattr(arrays, var, var_arr)

        return arrays

    @staticmethod
    def _interp_neighbors(arr, yindx, xindx):
        """
        Interpolate using neighboring pixels in both horizontal and vertical directions.

        Parameters
        ----------
        arr : ndarray
            2-D input array.
        yindx, xindx : ndarray
            1-D arrays, each length N, of row/column indices of the bad pixels.

        Returns
        -------
        ndarray
            Interpolations with shape of ``(2, N)`` in the horizontal (0th index)
            and vertical (1st index) directions.
        """
        horiz = (arr[yindx, xindx - 1] + arr[yindx, xindx + 1]) / 2.0
        vert = (arr[yindx - 1, xindx] + arr[yindx + 1, xindx]) / 2.0
        return np.array([horiz, vert])

    def mingrad(self, arrays):
        """
        Replace pixels with the minimum gradient replacement method.

        Test the gradient along the spatial and spectral axes using
        immediately adjacent pixels.  Pick whichever dimension has the minimum
        absolute gradient and replace the missing pixel with the average
        of the two adjacent pixels along that dimension.

        This aims to make the process extremely local; near point sources it should do
        the replacement along the spectral axis avoiding sampling issues, while near bright
        extended emission line the replacement should be along the spatial axis.  May still
        be suboptimal near bright emission lines from unresolved point sources.

        Does not attempt any replacement if a NaN value is bordered by another NaN value
        along a given axis.

        Parameters
        ----------
        arrays : `PixelReplaceArrays`
            Pixel arrays and dispersion direction for the 2D spectrum to process.
            Arrays are modified in-place.

        Returns
        -------
        arrays : `PixelReplaceArrays`
            The input with flagged bad pixels now flagged with FLUX_ESTIMATED
            and holding a flux value estimated from adjacent pixels.
        """
        # np.nanmedian() entry full of NaN values would produce a numpy
        # warning (despite well-defined behavior - return a NaN)
        # so we suppress that here.
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")

        log.info("Using minimum gradient method.")

        in_var_dict = {
            "var_poisson": None,
            "var_rnoise": None,
            "var_flat": None,
        }
        interp_rootvar_dict = {
            "var_poisson": None,
            "var_rnoise": None,
            "var_flat": None,
        }
        for key in in_var_dict.keys():
            if getattr(arrays, key) is not None:
                in_var_dict[key] = getattr(arrays, key)

        # Make an array of x/y values on the detector
        (ysize, xsize) = arrays.data.shape
        basex, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))
        pad = 1  # Padding around edge of array to ensure we don't look for neighbors outside array

        # Find NaN-valued pixels
        indx = np.where(
            (~np.isfinite(arrays.data))
            & (basey > pad)
            & (basey < ysize - pad)
            & (basex > pad)
            & (basex < xsize - pad)
        )
        # X and Y indices
        yindx, xindx = indx[0], indx[1]

        # Absolute gradient along each axis from indata, shape (2, N), used to choose direction
        diffs = np.array(
            [
                np.abs(arrays.data[yindx, xindx - 1] - arrays.data[yindx, xindx + 1]),
                np.abs(arrays.data[yindx - 1, xindx] - arrays.data[yindx + 1, xindx]),
            ]
        )

        # Interpolated values for each quantity in both directions, shape (2, N)
        interp_data = self._interp_neighbors(arrays.data, yindx, xindx)
        interp_err = self._interp_neighbors(arrays.err, yindx, xindx)
        # Propagate variance components as errors to get the scales right
        for key in in_var_dict.keys():
            if in_var_dict[key] is not None:
                interp_rootvar_dict[key] = self._interp_neighbors(
                    np.sqrt(in_var_dict[key]), yindx, xindx
                )

        # Replace NaN diffs with inf so argmin naturally prefers the valid direction.
        # Mask is True where at least one valid direction, False elsewhere,
        # such that pixels where both diffs are inf have no usable neighbor pair and are skipped.
        diffs_with_infs = np.where(np.isnan(diffs), np.inf, diffs)  # (2, N)
        mask = ~np.all(np.isinf(diffs_with_infs), axis=0)  # (N,)

        # Per-pixel direction index: 0 = horizontal, 1 = vertical
        indmin = np.argmin(diffs_with_infs, axis=0)  # (N,)
        col_idx = np.arange(len(yindx))

        # Select the minimium-gradient interpolated values and update model with them
        indmin = indmin[mask]
        col_idx = col_idx[mask]
        arrays.data[yindx[mask], xindx[mask]] = interp_data[indmin, col_idx]
        arrays.err[yindx[mask], xindx[mask]] = interp_err[indmin, col_idx]
        # Square the interpolated errors back to variances for insertion
        for key in interp_rootvar_dict.keys():
            if interp_rootvar_dict[key] is not None:
                in_var_dict[key][yindx[mask], xindx[mask]] = (
                    interp_rootvar_dict[key][indmin, col_idx] ** 2
                )
                setattr(arrays, key, in_var_dict[key])

        # Update DQ flags for pixels that were replaced.
        orig_dq = arrays.dq[yindx, xindx]  # (N,)
        remove_dnu = (
            mask
            & (orig_dq & self.DO_NOT_USE).astype(bool)
            & ~(orig_dq & self.NON_SCIENCE).astype(bool)
        )
        arrays.dq[yindx[remove_dnu], xindx[remove_dnu]] -= self.DO_NOT_USE
        arrays.dq[yindx[mask], xindx[mask]] |= self.FLUX_ESTIMATED

        return arrays

    def custom_slice(self, dispaxis, index):
        """
        Construct slice for ease of use with varying dispersion axis.

        Parameters
        ----------
        dispaxis : int
            Using module-defined:

            * 1 = HORIZONTAL
            * 2 = VERTICAL

        index : int or list
            Index or indices of cross-dispersion
            vectors to slice

        Returns
        -------
        tuple
            Slice constructed using numpy
        """
        if dispaxis == self.HORIZONTAL:
            return np.s_[:, index]
        elif dispaxis == self.VERTICAL:
            return np.s_[index, :]
        else:
            raise IndexError("Custom slice requires valid dispersion axis specification!")

    def profile_mse(self, scale, median, current):
        """
        Calculate mean-squared error of fitted profile.

        Parameters
        ----------
        scale : float
            Initial estimate of scale factor to bring
            normalized median profile up to current profile
        median : ndarray
            Median profile constructed from neighboring profile slices
        current : ndarray
            Current profile with bad pixels to be replaced

        Returns
        -------
        float
            Mean-squared error for minimization purposes
        """
        return np.nansum((current - (median * scale)) ** 2.0) / (
            len(median) - np.count_nonzero(np.isnan(current))
        )
