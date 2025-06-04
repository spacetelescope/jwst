import logging
import warnings

import numpy as np
from scipy.optimize import minimize
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import nirspec

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class PixelReplacement:
    """
    Main class for performing pixel replacement.

    This class controls loading the input data model, selecting the
    method for pixel replacement, and executing each step. This class
    should provide modularization to allow for multiple options and possible
    future reference files.
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
        """
        Initialize the class with input data model.

        Parameters
        ----------
        input_model : datamodel, str
            Datamodel or list of data models as ModelContainer
            or ASN file, one datamodel for each input image

        **pars : dict, optional
            Optional parameters to modify how pixel replacement
            will execute.
        """
        self.input = input_model
        self.pars = {}
        self.pars.update(pars)
        self.output = self.input.copy()
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

    def replace(self):
        """
        Unpack model and apply pixel replacement algorithm.

        Process the input DataModel, unpack any model that holds
        more than one 2D spectrum, then apply selected algorithm
        to each 2D spectrum in input.
        """
        # ImageModel inputs (MIR_LRS-FIXEDSLIT)
        # or 2D SlitModel inputs (e.g. NRS_FIXEDSLIT in spec3)
        if isinstance(self.input, datamodels.ImageModel) or (
            isinstance(self.input, datamodels.SlitModel) and self.input.data.ndim == 2
        ):
            self.output = self.algorithm(self.input)
            n_replaced = np.count_nonzero(self.output.dq & self.FLUX_ESTIMATED)
            log.info(f"Input model had {n_replaced} pixels replaced.")
        elif isinstance(self.input, datamodels.IFUImageModel):
            # Attempt to run pixel replacement on each throw of the IFU slicer
            # individually.
            xx, yy = np.indices(self.input.data.shape)
            if self.input.meta.exposure.type == "MIR_MRS":
                if self.pars["algorithm"] == "mingrad":
                    # mingrad method
                    new_model = self.algorithm(self.input)
                    self.output = new_model
                else:
                    # fit_profile method
                    _, beta_array, _ = self.input.meta.wcs.transform(
                        "detector", "alpha_beta", yy, xx
                    )
                    unique_beta = np.unique(beta_array)
                    unique_beta = unique_beta[~np.isnan(unique_beta)]
                    for i, beta in enumerate(unique_beta):
                        # Define a mask that is True where this trace is located
                        trace_mask = beta_array == beta
                        trace_model = self.input.copy()
                        trace_model.dq = np.where(
                            # When not in this trace, set NON_SCIENCE and DO_NOT_USE
                            ~trace_mask,
                            trace_model.dq | self.DO_NOT_USE | self.NON_SCIENCE,
                            trace_model.dq,
                        )

                        trace_model = self.algorithm(trace_model)
                        self.output.data = np.where(
                            # Where trace is located, set replaced values
                            trace_mask,
                            trace_model.data,
                            self.output.data,
                        )

                        # do the same for dq, err, and var
                        self.output.dq = np.where(trace_mask, trace_model.dq, self.output.dq)
                        self.output.err = np.where(trace_mask, trace_model.err, self.output.err)
                        self.output.var_poisson = np.where(
                            trace_mask, trace_model.var_poisson, self.output.var_poisson
                        )
                        self.output.var_rnoise = np.where(
                            trace_mask, trace_model.var_rnoise, self.output.var_rnoise
                        )
                        self.output.var_flat = np.where(
                            trace_mask, trace_model.var_flat, self.output.var_flat
                        )

                        n_replaced = np.count_nonzero(trace_model.dq & self.FLUX_ESTIMATED)
                        log.info(
                            f"Input MRS frame had {n_replaced} pixels replaced "
                            f"in IFU slice {i + 1}."
                        )

                        trace_model.close()

                n_replaced = np.count_nonzero(self.output.dq & self.FLUX_ESTIMATED)
                log.info(f"Input MRS frame had {n_replaced} total pixels replaced.")
            else:
                if self.pars["algorithm"] == "mingrad":
                    # mingrad method
                    new_model = self.algorithm(self.input)
                    self.output = new_model
                else:
                    # fit_profile method - iterate over IFU slices
                    for i in range(30):
                        slice_wcs = nirspec.nrs_wcs_set_input(self.input, i)
                        _, _, wave = slice_wcs.transform("detector", "slicer", yy, xx)
                        # Define a mask that is True where this trace is located
                        trace_mask = wave > 0
                        trace_model = self.input.copy()
                        trace_model.dq = np.where(
                            # When not in this trace, set NON_SCIENCE and DO_NOT_USE
                            ~trace_mask,
                            trace_model.dq | self.DO_NOT_USE | self.NON_SCIENCE,
                            trace_model.dq,
                        )

                        trace_model = self.algorithm(trace_model)
                        self.output.data = np.where(
                            # Where trace is located, set replaced values
                            trace_mask,
                            trace_model.data,
                            self.output.data,
                        )

                        # do the same for dq, err, and var
                        self.output.dq = np.where(trace_mask, trace_model.dq, self.output.dq)
                        self.output.err = np.where(trace_mask, trace_model.err, self.output.err)
                        self.output.var_poisson = np.where(
                            trace_mask, trace_model.var_poisson, self.output.var_poisson
                        )
                        self.output.var_rnoise = np.where(
                            trace_mask, trace_model.var_rnoise, self.output.var_rnoise
                        )
                        self.output.var_flat = np.where(
                            trace_mask, trace_model.var_flat, self.output.var_flat
                        )

                        n_replaced = np.count_nonzero(trace_model.dq & self.FLUX_ESTIMATED)
                        log.info(
                            f"Input NRS_IFU frame had {n_replaced} pixels "
                            f"replaced in IFU slice {i + 1}."
                        )

                        trace_model.close()

                n_replaced = np.count_nonzero(self.output.dq & self.FLUX_ESTIMATED)
                log.info(f"Input NRS_IFU frame had {n_replaced} total pixels replaced.")

        # MultiSlitModel inputs (WFSS, NRS_FIXEDSLIT, ?)
        elif isinstance(self.input, datamodels.MultiSlitModel):
            for i, _slit in enumerate(self.input.slits):
                slit_model = datamodels.SlitModel(self.input.slits[i].instance)

                slit_replaced = self.algorithm(slit_model)

                n_replaced = np.count_nonzero(slit_replaced.dq & self.FLUX_ESTIMATED)
                log.info(f"Slit {i} had {n_replaced} pixels replaced.")

                self.output.slits[i] = slit_replaced

        # CubeModel inputs are TSO (so far?); SlitModel may be NRS_BRIGHTOBJ,
        # also requiring a re-packaging of the data into 2D inputs for the algorithm
        elif isinstance(self.input, datamodels.CubeModel | datamodels.SlitModel):
            for i in range(len(self.input.data)):
                img_model = datamodels.ImageModel(
                    data=self.input.data[i],
                    dq=self.input.dq[i],
                    err=self.input.err[i],
                    var_poisson=self.input.var_poisson[i],
                    var_rnoise=self.input.var_rnoise[i],
                    var_flat=self.input.var_flat[i],
                )
                img_model.update(self.input)
                img_replaced = self.algorithm(img_model)
                n_replaced = np.count_nonzero(img_replaced.dq & self.FLUX_ESTIMATED)
                log.info(f"Input TSO integration {i} had {n_replaced} pixels replaced.")

                self.output.data[i] = img_replaced.data
                self.output.dq[i] = img_replaced.dq
                self.output.err[i] = img_replaced.err
                self.output.var_poisson[i] = img_replaced.var_poisson
                self.output.var_rnoise[i] = img_replaced.var_rnoise
                self.output.var_flat[i] = img_replaced.var_flat
                img_replaced.close()
                img_model.close()

        else:
            # This should never happen, as these should be caught in the step code.
            log.critical(
                "Pixel replacement code did not filter this input correctly - skipping step."
            )
            return

    def fit_profile(self, model):
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
        model : DataModel
            Either the input to the pixel_replace step in the
            case of DataModels containing only one 2D spectrum,
            or a single 2D spectrum from the input DataModel
            containing multiple spectra (i.e. MultiSlitModel).
            Requires data and dq attributes.

        Returns
        -------
        model_replaced : DataModel
            DataModel with flagged bad pixels now flagged with
            FLUX_ESTIMATED and holding a flux value estimated
            from spatial profile, derived from adjacent columns.
        """
        # np.nanmedian() entry full of NaN values would produce a numpy
        # warning (despite well-defined behavior - return a NaN)
        # so we suppress that here.
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")

        dispaxis = model.meta.wcsinfo.dispersion_direction

        model_replaced = model.copy()

        # Truncate array to region where good pixels exist
        good_pixels = np.where(~model.dq & self.DO_NOT_USE)
        if np.any(0 in np.shape(good_pixels)):
            log.warning(
                "No good pixels in at least one dimension of "
                "data array - skipping pixel replacement."
            )
            return model
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
            dq_slice = model.dq[self.custom_slice(dispaxis, ind)][profile_cut[0] : profile_cut[1]]
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
            profile_data = model.data[adjacent_condition]
            profile_err = model.err[adjacent_condition]

            # Mask out bad pixels
            invalid_condition = (model.dq[adjacent_condition] & self.DO_NOT_USE).astype(bool)
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
                    err = np.sqrt(getattr(model, err_name))
                else:
                    err = getattr(model, err_name)
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
            current_profile = model.data[current_condition]
            cleaned_current = np.where(
                model.dq[current_condition] & self.DO_NOT_USE, np.nan, current_profile
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
            current_dq = model.dq[current_condition][range(*profile_cut)]
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
            model_replaced.data[current_condition][range(*profile_cut)] = replaced_current
            model_replaced.dq[current_condition][range(*profile_cut)] = replaced_dq

            # Also update the errors and variances
            current_err = model.err[current_condition][range(*profile_cut)]
            replaced_err = np.where(
                replace_condition, norm_errors["err"] * norm_scale * scale, current_err
            )
            model_replaced.err[current_condition][range(*profile_cut)] = replaced_err

            # Some values in NIRSpec variances may overflow in the squares - ignore the warning.
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "overflow encountered", RuntimeWarning)
                current_var = model.var_poisson[current_condition][range(*profile_cut)]
                replaced_var = np.where(
                    replace_condition,
                    (norm_errors["var_poisson"] * norm_scale * scale) ** 2,
                    current_var,
                )
                model_replaced.var_poisson[current_condition][range(*profile_cut)] = replaced_var

                current_var = model.var_rnoise[current_condition][range(*profile_cut)]
                replaced_var = np.where(
                    replace_condition,
                    (norm_errors["var_rnoise"] * norm_scale * scale) ** 2,
                    current_var,
                )
                model_replaced.var_rnoise[current_condition][range(*profile_cut)] = replaced_var

                current_var = model.var_flat[current_condition][range(*profile_cut)]
                replaced_var = np.where(
                    replace_condition,
                    (norm_errors["var_flat"] * norm_scale * scale) ** 2,
                    current_var,
                )
                model_replaced.var_flat[current_condition][range(*profile_cut)] = replaced_var

        return model_replaced

    def mingrad(self, model):
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
        model : DataModel
            Either the input to the pixel_replace step in the
            case of DataModels containing only one 2D spectrum,
            or a single 2D spectrum from the input DataModel
            containing multiple spectra (i.e. MultiSlitModel).
            Requires data and dq attributes.

        Returns
        -------
        model_replaced : DataModel
            DataModel with flagged bad pixels now flagged with
            FLUX_ESTIMATED and holding a flux value estimated
            from spatial profile, derived from adjacent columns.
        """
        # np.nanmedian() entry full of NaN values would produce a numpy
        # warning (despite well-defined behavior - return a NaN)
        # so we suppress that here.
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")

        log.info("Using minimum gradient method.")

        # Input data, err, and dq values
        indata = model.data
        indq = model.dq
        inerr = model.err

        # Propagate variance components as errors to get the scales right
        in_var_p = np.sqrt(model.var_poisson)
        in_var_r = np.sqrt(model.var_rnoise)
        in_var_f = np.sqrt(model.var_flat)

        # Output data, err, var, and dq values
        model_replaced = model.copy()
        newdata = model_replaced.data
        newdq = model_replaced.dq
        newerr = model_replaced.err
        new_var_p = model_replaced.var_poisson
        new_var_r = model_replaced.var_rnoise
        new_var_f = model_replaced.var_flat

        # Make an array of x/y values on the detector
        (ysize, xsize) = indata.shape
        basex, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))
        pad = 1  # Padding around edge of array to ensure we don't look for neighbors outside array

        # Find NaN-valued pixels
        indx = np.where(
            (~np.isfinite(indata))
            & (basey > pad)
            & (basey < ysize - pad)
            & (basex > pad)
            & (basex < xsize - pad)
        )
        # X and Y indices
        yindx, xindx = indx[0], indx[1]

        # Loop over these NaN-valued pixels
        nreplaced = 0
        for ii in range(0, len(xindx)):
            left_data, right_data = (
                indata[yindx[ii], xindx[ii] - 1],
                indata[yindx[ii], xindx[ii] + 1],
            )
            top_data, bottom_data = (
                indata[yindx[ii] - 1, xindx[ii]],
                indata[yindx[ii] + 1, xindx[ii]],
            )

            left_err, right_err = inerr[yindx[ii], xindx[ii] - 1], inerr[yindx[ii], xindx[ii] + 1]
            top_err, bottom_err = inerr[yindx[ii] - 1, xindx[ii]], inerr[yindx[ii] + 1, xindx[ii]]

            left_var_p, right_var_p = (
                in_var_p[yindx[ii], xindx[ii] - 1],
                in_var_p[yindx[ii], xindx[ii] + 1],
            )
            top_var_p, bottom_var_p = (
                in_var_p[yindx[ii] - 1, xindx[ii]],
                in_var_p[yindx[ii] + 1, xindx[ii]],
            )

            left_var_r, right_var_r = (
                in_var_r[yindx[ii], xindx[ii] - 1],
                in_var_r[yindx[ii], xindx[ii] + 1],
            )
            top_var_r, bottom_var_r = (
                in_var_r[yindx[ii] - 1, xindx[ii]],
                in_var_r[yindx[ii] + 1, xindx[ii]],
            )

            left_var_f, right_var_f = (
                in_var_f[yindx[ii], xindx[ii] - 1],
                in_var_f[yindx[ii], xindx[ii] + 1],
            )
            top_var_f, bottom_var_f = (
                in_var_f[yindx[ii] - 1, xindx[ii]],
                in_var_f[yindx[ii] + 1, xindx[ii]],
            )

            # Compute absolute difference (slope) and average value in each direction (may be NaN)
            diffs = np.array([np.abs(left_data - right_data), np.abs(top_data - bottom_data)])
            interp_data = np.array([(left_data + right_data) / 2.0, (top_data + bottom_data) / 2.0])
            interp_err = np.array([(left_err + right_err) / 2.0, (top_err + bottom_err) / 2.0])
            interp_var_p = np.array(
                [(left_var_p + right_var_p) / 2.0, (top_var_p + bottom_var_p) / 2.0]
            )
            interp_var_r = np.array(
                [(left_var_r + right_var_r) / 2.0, (top_var_r + bottom_var_r) / 2.0]
            )
            interp_var_f = np.array(
                [(left_var_f + right_var_f) / 2.0, (top_var_f + bottom_var_f) / 2.0]
            )

            # Replace with the value from the lowest absolute slope estimator that was not NaN
            try:
                indmin = np.nanargmin(diffs)
                newdata[yindx[ii], xindx[ii]] = interp_data[indmin]
                newerr[yindx[ii], xindx[ii]] = interp_err[indmin]

                # Square the interpolated errors back into variance
                new_var_p[yindx[ii], xindx[ii]] = interp_var_p[indmin] ** 2
                new_var_r[yindx[ii], xindx[ii]] = interp_var_r[indmin] ** 2
                new_var_f[yindx[ii], xindx[ii]] = interp_var_f[indmin] ** 2

                # If original pixel was in the science array, remove
                # the DO_NOT_USE flag
                if (indq[yindx[ii], xindx[ii]] & self.DO_NOT_USE) and not (
                    indq[yindx[ii], xindx[ii]] & self.NON_SCIENCE
                ):
                    newdq[yindx[ii], xindx[ii]] -= self.DO_NOT_USE

                # Either way, add the FLUX_ESTIMATED flag
                newdq[yindx[ii], xindx[ii]] |= self.FLUX_ESTIMATED

                nreplaced += 1

            except (IndexError, ValueError):
                pass

        model_replaced.data = newdata
        model_replaced.err = newerr
        model_replaced.dq = newdq

        return model_replaced

    def custom_slice(self, dispaxis, index):
        """
        Construct slice for ease of use with varying dispersion axis.

        Parameters
        ----------
        dispaxis : int
            Using module-defined HORIZONTAL=1,
            VERTICAL=2

        index : int or list
            Index or indices of cross-dispersion
            vectors to slice

        Returns
        -------
        Tuple
            Slice constructed using np.s_
        """
        if dispaxis == self.HORIZONTAL:
            return np.s_[:, index]
        elif dispaxis == self.VERTICAL:
            return np.s_[index, :]
        else:
            raise IndexError("Custom slice requires valid dispersion axis specification!")

    def profile_mse(self, scale, median, current):
        """
        Calculate mean squared error of fitted profile.

        Parameters
        ----------
        scale : float
            Initial estimate of scale factor to bring
            normalized median profile up to current profile
        median : array
            Median profile constructed from neighboring
            profile slices
        current : array
            Current profile with bad pixels to be
            replaced

        Returns
        -------
        float
            Mean squared error for minimization purposes
        """
        return np.nansum((current - (median * scale)) ** 2.0) / (
            len(median) - np.count_nonzero(np.isnan(current))
        )
