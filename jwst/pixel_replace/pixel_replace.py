import logging
import numpy as np
from stdatamodels.jwst import datamodels
from scipy.optimize import minimize
import warnings
from ..assign_wcs.nirspec import nrs_wcs_set_input

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class PixelReplacement:
    """Main class for performing pixel replacement.

    This class controls loading the input data model, selecting the
    method for pixel replacement, and executing each step. This class
    should provide modularization to allow for multiple options and possible
    future reference files.

    """

    # Shortcuts for DQ Flags
    DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
    FLUX_ESTIMATED = datamodels.dqflags.pixel['FLUX_ESTIMATED']
    NON_SCIENCE = datamodels.dqflags.pixel['NON_SCIENCE']

    # Shortcuts for dispersion direction for ease of reading
    HORIZONTAL = 1
    VERTICAL = 2

    default_suffix = 'pixrep'

    def __init__(self, input_model, **pars):
        """
        Initialize the class with input data model.

        Parameters
        ----------
        input_model : DataModel, str
            list of data models as ModelContainer or ASN file,
            one data model for each input image

        pars : dict, optional
            Optional parameters to modify how pixel replacement
            will execute.
        """
        self.input = input_model
        self.pars = dict()
        self.pars.update(pars)
        self.output = self.input.copy()
        # Store algorithm options here.
        self.algorithm_dict = {
            'fit_profile': self.fit_profile,
            'mingrad': self.mingrad,
        }

        # Choose algorithm from dict using input par.
        try:
            self.algorithm = self.algorithm_dict[self.pars['algorithm']]

        except KeyError:
            log.critical(f"Algorithm name {self.pars['algorithm']} provided does "
                         "not match an implemented algorithm!")
            raise Exception

    def replace(self):
        """
        Process the input DataModel, unpack any model that holds
        more than one 2D spectrum, then apply selected algorithm
        to each 2D spectrum in input.
        """
        # ImageModel inputs (MIR_LRS-FIXEDSLIT)
        if isinstance(self.input, datamodels.ImageModel):
            self.output = self.algorithm(self.input)
            n_replaced = np.count_nonzero(self.output.dq & self.FLUX_ESTIMATED)
            log.info(f"Input model had {n_replaced} pixels replaced.")
        elif isinstance(self.input, datamodels.IFUImageModel):
                # Attempt to run pixel replacement on each throw of the IFU slicer
                # individually.
                xx, yy = np.indices(self.input.data.shape)
                if self.input.meta.exposure.type == 'MIR_MRS':
                    # Mingrad method
                    if (self.pars['algorithm'] == 'mingrad'):
                        new_model = self.algorithm(self.input)
                        self.output = new_model
                    # fit_profile method
                    else:
                        _, beta_array, _ = self.input.meta.wcs.transform('detector', 'alpha_beta', yy, xx)
                        unique_beta = np.unique(beta_array)
                        unique_beta = unique_beta[~np.isnan(unique_beta)]
                        for i, beta in enumerate(unique_beta):
                            # Define a mask that is True where this trace is located
                            trace_mask = (beta_array == beta)
                            trace_model = self.input.copy()
                            trace_model.dq = np.where(
                                # When not in this trace, set NON_SCIENCE and DO_NOT_USE
                                ~trace_mask,
                                trace_model.dq | self.DO_NOT_USE | self.NON_SCIENCE,
                                trace_model.dq
                            )

                            trace_model = self.algorithm(trace_model)
                            self.output.data = np.where(
                                # Where trace is located, set replaced values
                                trace_mask,
                                trace_model.data,
                                self.output.data
                            )

                            self.output.dq = np.where(
                                # Where trace is located, set replaced values
                                trace_mask,
                                trace_model.dq,
                                self.output.dq
                            )

                            n_replaced = np.count_nonzero(trace_model.dq & self.FLUX_ESTIMATED)
                            log.info(f"Input MRS frame had {n_replaced} pixels replaced in IFU slice {i+1}.")

                            trace_model.close()

                    n_replaced = np.count_nonzero(self.output.dq & self.FLUX_ESTIMATED)
                    log.info(f"Input MRS frame had {n_replaced} total pixels replaced.")
                else:
                    # NRS_IFU method - Fixed number of IFU slices to iterate over
                    for i in range(30):
                        slice_wcs = nrs_wcs_set_input(self.input, i)
                        _, _, wave = slice_wcs.transform('detector', 'slicer', yy, xx)

                        # Define a mask that is True where this trace is located
                        trace_mask = (wave > 0)
                        trace_model = self.input.copy()
                        trace_model.dq = np.where(
                            # When not in this trace, set NON_SCIENCE and DO_NOT_USE
                            ~trace_mask,
                            trace_model.dq | self.DO_NOT_USE | self.NON_SCIENCE,
                            trace_model.dq
                        )

                        trace_model = self.algorithm(trace_model)
                        self.output.data = np.where(
                            # Where trace is located, set replaced values
                            trace_mask,
                            trace_model.data,
                            self.output.data
                        )

                        self.output.dq = np.where(
                            # Where trace is located, set replaced values
                            trace_mask,
                            trace_model.dq,
                            self.output.dq
                        )

                        n_replaced = np.count_nonzero(trace_model.dq & self.FLUX_ESTIMATED)
                        log.info(f"Input NRS_IFU frame had {n_replaced} pixels replaced in IFU slice {i + 1}.")

                        trace_model.close()

                    n_replaced = np.count_nonzero(self.output.dq & self.FLUX_ESTIMATED)
                    log.info(f"Input NRS_IFU frame had {n_replaced} total pixels replaced.")

        # MultiSlitModel inputs (WFSS, NRS_FIXEDSLIT, ?)
        elif isinstance(self.input, datamodels.MultiSlitModel):

            for i, slit in enumerate(self.input.slits):
                slit_model = datamodels.SlitModel(self.input.slits[i].instance)

                slit_replaced = self.algorithm(slit_model)

                n_replaced = np.count_nonzero(slit_replaced.dq & self.FLUX_ESTIMATED)
                log.info(f"Slit {i} had {n_replaced} pixels replaced.")

                self.output.slits[i] = slit_replaced

        # CubeModel inputs are TSO (so far?); only SlitModel seen so far is NRS_BRIGHTOBJ, also requiring
        # a re-packaging of the data into 2D inputs for the algorithm.
        elif isinstance(self.input, (datamodels.CubeModel, datamodels.SlitModel)):
            # Initial attempt looped over model.meta.exposure.nints, but test data had mismatch. Could change this.
            for i in range(len(self.input.data)):
                dummy_model = datamodels.ImageModel(data=self.input.data[i], dq=self.input.dq[i])
                dummy_model.update(self.input)
                dummy_replaced = self.algorithm(dummy_model)
                n_replaced = np.count_nonzero(dummy_replaced.dq & self.FLUX_ESTIMATED)
                log.info(f"Input TSO integration {i} had {n_replaced} pixels replaced.")

                self.output.data[i] = dummy_replaced.data
                self.output.dq[i] = dummy_replaced.dq
                dummy_replaced.close()
                dummy_model.close()

        else:
            # This should never happen, as these should be caught in the step code.
            log.critical("Pixel replacement code did not filter this input correctly - skipping step.")
            return

    def fit_profile(self, model):
        """
        Fit a profile to adjacent columns, scale profile to
        column with missing pixel(s), and find flux estimate
        from scaled profile.

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
        warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')

        dispaxis = model.meta.wcsinfo.dispersion_direction

        model_replaced = model.copy()

        # Truncate array to region where good pixels exist
        good_pixels = np.where(~model.dq & 1)
        if np.any(0 in np.shape(good_pixels)):
            log.warning("No good pixels in at least one dimension of "
                        "data array - skipping pixel replacement.")
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
            dq_slice = model.dq[self.custom_slice(dispaxis, ind)][profile_cut[0]: profile_cut[1]]
            # Exclude regions with NON_SCIENCE flag
            dq_slice = np.where(
                dq_slice & self.NON_SCIENCE,
                self.NON_SCIENCE,
                dq_slice
            )
            # Find bad pixels in region containing valid data.
            n_bad = np.count_nonzero(dq_slice & self.DO_NOT_USE)
            if n_bad == len(dq_slice):
                log.debug(f"Slice {ind} contains no good pixels. Skipping replacement.")
                valid_profiles.discard(ind)
            elif n_bad == 0:
                log.debug(f"Slice {ind} contains no bad pixels.")
            else:
                log.debug(f"Slice {ind} contains {n_bad} bad pixels.")
                profiles_to_replace.add(ind)

        log.debug(f"Number of profiles with at least one bad pixel: {len(profiles_to_replace)}")

        for i, ind in enumerate(profiles_to_replace):

            # Use sets for convenient finding of neighboring slices to use in profile creation
            adjacent_inds = set(range(ind - self.pars['n_adjacent_cols'], ind + self.pars['n_adjacent_cols'] + 1))
            adjacent_inds.discard(ind)
            valid_adjacent_inds = list(adjacent_inds.intersection(valid_profiles))
            # Cut out valid neighboring profiles
            profile_data = model.data[self.custom_slice(dispaxis, valid_adjacent_inds)]
            # Mask out bad pixels
            profile_data = np.where(
                model.dq[self.custom_slice(dispaxis, valid_adjacent_inds)] & self.DO_NOT_USE,
                np.nan,
                profile_data
            )
            # Add additional cut to pull only from region with valid data for convenience (may not be necessary)
            profile_data = profile_data[self.custom_slice(3 - dispaxis, list(range(profile_cut[0], profile_cut[1])))]

            # Normalize profile data
            normalized = profile_data / np.nanmax(np.abs(profile_data), axis=(dispaxis - 1), keepdims=True)

            # Pull median for each pixel across profile.
            # Profile entry full of NaN values would produce a numpy
            # warning (despite well-defined behavior - return a NaN)
            # so we suppress that above.
            median_profile = np.nanmedian(normalized, axis=(2 - dispaxis))

            # check_output[:, i] = median_profile

            # Clean current profile of values flagged as bad
            current_profile = model.data[self.custom_slice(dispaxis, ind)]
            cleaned_current = np.where(
                model.dq[self.custom_slice(dispaxis, ind)] & self.DO_NOT_USE,
                np.nan,
                current_profile
            )[range(*profile_cut)]

            replace_mask = np.where(~np.isnan(cleaned_current))[0]
            min_median = median_profile[replace_mask]
            min_current = cleaned_current[replace_mask]
            norm_current = min_current / np.max(min_current)

            # Scale median profile to current profile with bad pixel - minimize mse?
            norm_scale = minimize(self.profile_mse, x0=np.abs(np.nanmax(norm_current)),
                             args=(min_median, norm_current)).x

            scale = np.max(min_current)

            replaced_current = np.where(
                (model.dq[self.custom_slice(dispaxis, ind)][range(*profile_cut)] & self.DO_NOT_USE ^
                 model.dq[self.custom_slice(dispaxis, ind)][range(*profile_cut)] & self.NON_SCIENCE) == 1,
                median_profile * norm_scale * scale,
                cleaned_current
            )

            # Change the dq bits where old flag was DO_NOT_USE and new value is not nan
            current_dq = model.dq[self.custom_slice(dispaxis, ind)][range(*profile_cut)]
            replaced_dq = np.where(
                (current_dq & self.DO_NOT_USE ^ current_dq & self.NON_SCIENCE == 1) &
                ~(np.isnan(replaced_current)),
                current_dq ^ self.DO_NOT_USE ^ self.FLUX_ESTIMATED,
                current_dq
            )

            model_replaced.data[self.custom_slice(dispaxis, ind)][range(*profile_cut)] = replaced_current
            model_replaced.dq[self.custom_slice(dispaxis, ind)][range(*profile_cut)] = replaced_dq

        return model_replaced

    def mingrad(self, model):
        """
        Minimum gradient replacement method; test the gradient along the spatial
        and spectral axes using immediately adjacent pixels.  Pick whichever dimension
        has the minimum absolute gradient and replace the missing pixel with the average
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
        warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')

        log.info("Using minimum gradient method.")

        # Input data, err, and dq values
        indata = model.data
        inerr = model.err
        indq = model.dq

        # Output data, err, and dq values
        model_replaced = model.copy()
        newdata = model_replaced.data
        newerr = model_replaced.err
        newdq = model_replaced.dq

        # Make an array of x/y values on the detector
        xsize, ysize = self.input.meta.subarray.xsize, self.input.meta.subarray.ysize
        basex, basey = np.meshgrid(np.arange(xsize), np.arange(ysize))
        pad = 1 # Padding around edge of array to ensure we don't look for neighbors outside array

        # Find NaN-valued pixels
        indx = np.where((~np.isfinite(indata))
                        & (basey > pad) & (basey < ysize-pad) & (basex > pad) & (basex < xsize-pad))
        # X and Y indices
        yindx, xindx = indx[0], indx[1]

        # Loop over these NaN-valued pixels
        nreplaced = 0
        for ii in range(0, len(xindx)):
            left_data, right_data = indata[yindx[ii], xindx[ii] - 1], indata[yindx[ii], xindx[ii] + 1]
            top_data, bottom_data = indata[yindx[ii] - 1, xindx[ii]], indata[yindx[ii] + 1, xindx[ii]]
            left_err, right_err = inerr[yindx[ii], xindx[ii] - 1], inerr[yindx[ii], xindx[ii] + 1]
            top_err, bottom_err = inerr[yindx[ii] - 1, xindx[ii]], inerr[yindx[ii] + 1, xindx[ii]]

            # Compute absolute difference (slope) and average value in each direction (may be NaN)
            diffs = np.array([np.abs(left_data - right_data), np.abs(top_data - bottom_data)])
            interp_data = np.array([(left_data + right_data) / 2., (top_data + bottom_data) / 2.])
            interp_err = np.array([(left_err + right_err) / 2., (top_err + bottom_err) / 2.])

            # Replace with the value from the lowest absolute slope estimator that was not NaN
            try:
                indmin = np.nanargmin(diffs)
                newdata[yindx[ii], xindx[ii]] = interp_data[indmin]
                newerr[yindx[ii], xindx[ii]] = interp_err[indmin]
                # If original pixel was in the science array set DQ accordingly
                if ((indq[yindx[ii], xindx[ii]] & self.NON_SCIENCE) == 0):
                    newdq[yindx[ii], xindx[ii]] = self.FLUX_ESTIMATED
                # If original pixel was in non-science region set DQ accordingly
                else:
                    newdq[yindx[ii], xindx[ii]] = self.FLUX_ESTIMATED + self.NON_SCIENCE + self.DO_NOT_USE

                nreplaced += 1
            except Exception:
                pass

        model_replaced.data = newdata
        model_replaced.err = newerr
        model_replaced.dq = newdq

        return model_replaced

    def custom_slice(self, dispaxis, index):
        """
        Construct slice for ease of use with varying
        dispersion axis.

        Parameters
        ----------
        dispaxis : int
            Using module-defined HORIZONTAL=1,
            VERTICAL=2

        index : int or slice
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
            raise Exception

    def profile_mse(self, scale, median, current):
        """Function to feed optimization routine
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

        return (np.nansum((current - (median * scale)) ** 2.) /
                (len(median) - np.count_nonzero(np.isnan(current))))
