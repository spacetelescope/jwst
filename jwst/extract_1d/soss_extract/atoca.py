"""
Main classes for the ATOCA (Darveau-Bernier 2021, in prep).

ATOCA: Algorithm to Treat Order ContAmination (English)
       Algorithme de Traitement d’Ordres ContAmines (French)

@authors: Antoine Darveau-Bernier, Geert Jan Talens
"""

# General imports.
import numpy as np
from scipy.sparse import issparse, csr_matrix, diags

# Local imports.
from . import atoca_utils

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class MaskOverlapError(Exception):
    """Exception to raise if there are too few valid pixels in a spectral order."""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ExtractionEngine:
    """
    Run the ATOCA algorithm (Darveau-Bernier 2022, PASP, DOI:10.1088/1538-3873/ac8a77).

    The ExtractionEngine is basically a fitter. On instantiation, it generates a model
    of the detector, including a mapping between the detector pixels and the wavelength
    for each spectral order, the throughput and convolution kernel, and known detector
    bad pixels. This does not require any real data.
    When called, it ingests data and associated errors, then
    generates an output 1-D spectrum that explains the pixel brightnesses in the data
    within the constraints of the model.

    The engine can also run in reverse:
    The `rebuild` method generates a synthetic 2-D detector 'observation'
    from a known or fitted spectrum, and the `compute_likelihood` method
    compares the synthetic data to the real data to generate a likelihood.
    This allows for a likelihood-based optimization of the spectrum.

    This version models the pixels of the detector using an oversampled trapezoidal integration.
    """

    # The desired data-type for computations. 'float64' is recommended.
    dtype = "float64"

    def __init__(
        self,
        wave_map,
        trace_profile,
        throughput,
        kernels,
        wave_grid,
        mask_trace_profile,
        global_mask=None,
        orders=None,
        threshold=1e-3,
    ):
        """
        Initialize the ExtractionEngine with a detector model.

        Parameters
        ----------
        wave_map : list or array of 2-D arrays
            A list or array of the central wavelength position for each
            order on the detector. Has shape (N_ord, N, M).
            It has to have the same (N, M) as `data`.
        trace_profile : list or array of 2-D arrays
            A list or array of the spatial profile for each order
            Has shape (N_ord, N, M).
            on the detector. It has to have the same (N, M) as `data`.
        throughput : list of array or callable
            A list of functions or array of the throughput at each order.
            If callable, the functions depend on the wavelength.
            If array, projected on `wave_grid`. Has shape (N_ord [, N_k]).
        kernels : callable, sparse matrix, or None
            Convolution kernel to be applied on spectrum (f_k) for each orders.
            Can be a callable with the form f(x, x0) where x0 is
            the position of the center of the kernel. In this case, it must
            return a 1D array (len(x)), so a kernel value
            for each pairs of (x, x0). If callable,
            it will be passed to `convolution.get_c_matrix` function
            and the `c_kwargs` can be passed to this function.
            If sparse, the shape has to be (N_k_c, N_k) and it will
            be used directly. N_ker is the length of the effective kernel
            and N_k_c is the length of the spectrum (f_k) convolved.
            If None, the kernel is set to 1, i.e., do not do any convolution.
        wave_grid : array_like, required
            The grid on which f(lambda) will be projected, shape (N_k).
        mask_trace_profile : List or array of 2-D arrays[bool], required
            A list or array of the pixel that need to be used for extraction,
            for each order on the detector.
            It has to have the same shape (N_ord, N, M) as `trace_profile`.
        global_mask : array[bool], optional
            Boolean Mask of the detector pixels to mask for every extraction, e.g. bad pixels.
            Should not be related to a specific order (if so, use `mask_trace_profile` instead).
            Has shape (N, M).
        orders : list, optional
            List of orders considered. Default is orders = [1, 2]
        threshold : float, optional:
            The contribution of any order on a pixel is considered significant if
            its estimated spatial profile is greater than this threshold value.
            If it is not properly modeled (not covered by the wavelength grid),
            it will be masked. Default is 1e-3.
        """
        if orders is None:
            orders = [1, 2]
        # Set the attributes and ensure everything has correct dtype
        self.wave_map = np.array(wave_map).astype(self.dtype)
        self.trace_profile = np.array(trace_profile).astype(self.dtype)
        self.mask_trace_profile = np.array(mask_trace_profile).astype(bool)
        self.threshold = threshold
        self.data_shape = self.wave_map[0].shape

        # Set wave_grid. Ensure it is sorted and strictly increasing.
        is_sorted = (np.diff(wave_grid) > 0).all()
        if not is_sorted:
            log.warning(
                "`wave_grid` is not strictly increasing. It will be sorted and made unique."
            )
            wave_grid = np.unique(wave_grid)
        self.wave_grid = wave_grid.astype(self.dtype).copy()
        self.n_wavepoints = len(wave_grid)

        # Get wavelengths at the boundaries of each pixel for all orders
        wave_p, wave_m = [], []
        for wave in self.wave_map:
            lp, lm = atoca_utils.get_wave_p_or_m(wave)
            wave_p.append(lp)
            wave_m.append(lm)
        self.wave_p = np.array(wave_p, dtype=self.dtype)
        self.wave_m = np.array(wave_m, dtype=self.dtype)

        # Set orders and ensure that the number of orders is consistent with wave_map length
        self.orders = orders
        self.n_orders = len(self.orders)
        if self.n_orders != len(self.wave_map):
            msg = (
                f"The number of orders specified ({self.n_orders}) and the number of "
                f"wavelength maps provided ({len(self.wave_map)}) do not match."
            )
            log.critical(msg.format(self.n_orders, len(self.wave_map)))
            raise ValueError(msg.format(self.n_orders, len(self.wave_map)))

        # Set a first estimate of i_bounds to estimate mask
        self.i_bounds = [[0, len(self.wave_grid)] for _ in range(self.n_orders)]

        # Estimate a global mask and masks for each orders
        self.mask, self.mask_ord = self._get_masks(global_mask)

        # Ensure there are adequate good pixels left in each order
        good_pixels_in_order = np.sum(np.sum(~self.mask_ord, axis=-1), axis=-1)
        min_good_pixels = 25  # hard-code to qualitatively reasonable value
        if np.any(good_pixels_in_order < min_good_pixels):
            msg = (
                f"At least one order has less than {min_good_pixels} valid pixels. "
                "(mask_trace_profile and mask_wave have insufficient overlap)"
            )
            raise MaskOverlapError(msg)

        # Update i_bounds based on masked wavelengths
        self.i_bounds = self._get_i_bnds()

        # if throughput is given as callable, turn it into an array
        # with shape (n_ord, wave_grid.size)
        self.update_throughput(throughput)

        # Re-build global mask and masks for each orders
        self.mask, self.mask_ord = self._get_masks(global_mask)
        # Save mask here as the general mask,  since `mask` attribute can be changed.
        self.general_mask = self.mask.copy()

        # turn kernels into sparse matrix
        self.kernels = self._create_kernels(kernels)

        # Compute integration weights. see method self.get_w() for details.
        self.weights, self.weights_k_idx = self.compute_weights()

        self.pixel_mapping = [None for _ in range(self.n_orders)]
        self.tikho_mat = None
        self.w_t_wave_c = None

    def get_attributes(self, *args, i_order=None):
        """
        Return list of attributes.

        Parameters
        ----------
        *args : str or list[str]
            All attributes to return.
        i_order : None or int, optional
            Index of order to extract. If specified, it will
            be applied to all attributes in args, so it cannot
            be mixed with non-order dependent attributes).

        Returns
        -------
        list
            Result of [getattr(arg) for arg in args], with i_order indexing if provided
        """
        if i_order is None:
            out = [getattr(self, arg) for arg in args]
        else:
            out = [getattr(self, arg)[i_order] for arg in args]

        if len(out) == 1:
            out = out[0]

        return out

    def update_throughput(self, throughput):
        """
        Update internal throughput values.

        Parameters
        ----------
        throughput : array[float] or callable
            Throughput values for each order, given either as an array
            or as a callable function with self.wave_grid as input.
        """
        throughput_new = []
        for throughput_n in throughput:  # Loop over orders.
            if callable(throughput_n):
                throughput_n = throughput_n(self.wave_grid)

            msg = "Throughputs must be given as callable or arrays matching the extraction grid."
            if not isinstance(throughput_n, np.ndarray):
                log.critical(msg)
                raise TypeError(msg)
            if throughput_n.shape != self.wave_grid.shape:
                log.critical(msg)
                raise ValueError(msg)

            throughput_new.append(throughput_n)

        self.throughput = np.array(throughput_new, dtype=self.dtype)

    def _create_kernels(self, kernels):
        """
        Make sparse matrix from input kernels.

        Parameters
        ----------
        kernels : callable, sparse matrix, or None
            Convolution kernel to be applied on the spectrum (f_k) for each order.
            If None, kernel is set to 1, i.e., do not do any convolution.

        Returns
        -------
        kernels_new : list
            List of sparse matrices for each order.
        """
        # Take thresh to be the kernels min_value attribute.
        # It is a way to make sure that the full kernel is used.
        c_kwargs = []
        for ker in kernels:
            try:
                kwargs_ker = {"thresh": ker.min_value}
            except AttributeError:
                # take the get_c_matrix defaults
                kwargs_ker = {}
            c_kwargs.append(kwargs_ker)

        # Define convolution sparse matrix.
        kernels_new = []
        for i_order, kernel_n in enumerate(kernels):
            if kernel_n is None:
                kernel_n = np.array([1.0])
            if not issparse(kernel_n):
                kernel_n = atoca_utils.get_c_matrix(
                    kernel_n, self.wave_grid, i_bounds=self.i_bounds[i_order], **c_kwargs[i_order]
                )

            kernels_new.append(kernel_n)

        return kernels_new

    def _get_masks(self, global_mask):
        """
        Compute a general mask on the detector and for each order.

        Depends on the trace profile and the wavelength grid.

        Parameters
        ----------
        global_mask : array[bool]
            Boolean mask of the detector pixels to mask for every extraction.

        Returns
        -------
        general_mask : array[bool]
            Mask that combines global_mask, wavelength mask, trace_profile mask
        mask_ord : array[bool]
            Mask applied to each order
        """
        # Get needed attributes
        args = ("threshold", "n_orders", "mask_trace_profile", "trace_profile")
        threshold, n_orders, mask_trace_profile, trace_profile = self.get_attributes(*args)

        # Mask pixels not covered by the wavelength grid.
        mask_wave = np.array([self.get_mask_wave(i_order) for i_order in range(n_orders)])

        # combine trace profile mask with wavelength cutoff mask
        # and apply detector bad pixel mask if specified
        if global_mask is None:
            mask_ord = np.any([mask_trace_profile, mask_wave], axis=0)
        else:
            mask = [global_mask for _ in range(n_orders)]  # For each orders
            mask_ord = np.any([mask_trace_profile, mask_wave, mask], axis=0)

        # Find pixels that are masked in each order.
        general_mask = np.all(mask_ord, axis=0)

        # Mask pixels if mask_trace_profile not masked but mask_wave is.
        # This means that an order is contaminated by another
        # order, but the wavelength range does not cover this part
        # of the spectrum. Thus, it cannot be treated correctly.
        is_contaminated = np.array([tr_profile_ord > threshold for tr_profile_ord in trace_profile])
        general_mask |= np.any(mask_wave, axis=0) & np.all(is_contaminated, axis=0)

        # Apply this new general mask to each order.
        mask_ord = mask_wave | general_mask[None, :, :]

        return general_mask, mask_ord

    def _get_i_bnds(self):
        """
        Define wavelength boundaries for each order using the order's mask and wavelength map.

        Returns
        -------
        list[float]
            Wavelength boundaries for each order
        """
        # Figure out boundary wavelengths
        wave_bounds = []
        for i in range(self.n_orders):
            wave = self.wave_map[i][~self.mask_ord[i]]
            wave_bounds.append([wave.min(), wave.max()])

        # Determine the boundary position on the wavelength grid.
        i_bnds_new = []
        for bounds, i_bnds in zip(wave_bounds, self.i_bounds, strict=True):
            a = np.min(np.where(self.wave_grid >= bounds[0])[0])
            b = np.max(np.where(self.wave_grid <= bounds[1])[0]) + 1

            # Take the most restrictive bound
            a = np.maximum(a, i_bnds[0])
            b = np.minimum(b, i_bnds[1])
            i_bnds_new.append([int(a), int(b)])

        return i_bnds_new

    def wave_grid_c(self, i_order):
        """
        Return wave_grid for a given order constrained according to the i_bounds of that order.

        Parameters
        ----------
        i_order : int
            Order to select the wave_grid for.

        Returns
        -------
        array[float]
            Wave_grid for the given order.
        """
        index = slice(*self.i_bounds[i_order])

        return self.wave_grid[index]

    def compute_weights(self):
        """
        Compute integration weights.

        The weights depend on the integration method used to solve
        the integral of the flux over a pixel and are encoded
        in the class method `get_w()`.

        Returns
        -------
        weights, weights_k_idx : list
            Lists of weights and corresponding grid indices
        """
        # Init lists
        weights, weights_k_idx = [], []
        for i_order in range(self.n_orders):
            # Compute weights
            weights_n, k_idx_n = self.get_w(i_order)

            # Convert to sparse matrix
            # First get the dimension of the convolved grid
            n_kc = np.diff(self.i_bounds[i_order]).astype(int)[0]

            # Then convert to sparse
            weights_n = atoca_utils.sparse_k(weights_n, k_idx_n, n_kc)
            weights.append(weights_n), weights_k_idx.append(k_idx_n)

        return weights, weights_k_idx

    def _set_w_t_wave_c(self, i_order, product):
        """
        Save intermediate matrix product for faster repeated computations.

        Saves the matrix product of the weights (w), the throughput (t),
        the wavelength (lam) and the convolution matrix.

        Parameters
        ----------
        i_order : int
            Order index to save the product for.
        product : sparse matrix
            The matrix product to save.
        """
        if self.w_t_wave_c is None:
            self.w_t_wave_c = [[] for _ in range(self.n_orders)]

        self.w_t_wave_c[i_order] = product.copy()

    def grid_from_map(self, i_order=0):
        """
        Return the wavelength grid and the columns for a given order index.

        Parameters
        ----------
        i_order : int, optional
            Order index to get the wavelength grid for. Default is 0.

        Returns
        -------
        wave_grid : array[float]
            Wavelength grid for the given order index.
        icol : array[float]
            Column indices for the wavelength grid.
        """
        attrs = ["wave_map", "trace_profile"]
        wave_map, trace_profile = self.get_attributes(*attrs, i_order=i_order)

        wave_grid, icol = atoca_utils.grid_from_map(wave_map, trace_profile)
        wave_grid = wave_grid.astype(self.dtype)
        return wave_grid, icol

    def get_pixel_mapping(self, i_order, error=None, quick=False):
        """
        Calculate the pixel mapping.

        Compute the matrix `b_n = (P/sig).w.T.lambda.c_n` ,
        where `P` is the spatial profile matrix (diag),
        `w` is the integrations weights matrix,
        `T` is the throughput matrix (diag),
        `lambda` is the convolved wavelength grid matrix (diag),
        `c_n` is the convolution kernel.
        The model of the detector at order n (`model_n`) is given by the system:
        model_n = b_n.c_n.f ,
        where f is the incoming flux projected on the wavelength grid.
        This method updates the `b_n_list` attribute.

        Parameters
        ----------
        i_order : int
            Label of the order (depending on the initiation of the object).
        error : array_like or None, optional
            Estimate of the error on each pixel. Same shape (N, M) as `data`.
            If None, the error is set to 1, which means the method will return
            b_n instead of b_n/sigma. Default is None.
        quick : bool, optional
            If True, only perform one matrix multiplication
            instead of the whole system: (P/sig).(w.T.lambda.c_n)

        Returns
        -------
        array[float]
            Sparse matrix of b_n coefficients
        """
        if (quick) and (self.w_t_wave_c is None):
            msg = "Attribute w_t_wave_c of ExtractionEngine must exist if quick=True"
            raise AttributeError(msg)

        # Special treatment for error map
        # Can be bool or array.
        if error is None:
            # Sigma will have no effect
            error = np.ones(self.data_shape)

        # Get needed attributes ...
        attrs = ["wave_grid", "mask"]
        wave_grid, mask = self.get_attributes(*attrs)

        # ... order dependent attributes
        attrs = ["trace_profile", "throughput", "kernels", "weights", "i_bounds"]
        trace_profile_n, throughput_n, kernel_n, weights_n, i_bnds = self.get_attributes(
            *attrs, i_order=i_order
        )

        # Keep only valid pixels (P and sig are still 2-D)
        # And apply directly 1/sig here (quicker)
        trace_profile_n = trace_profile_n[~mask] / error[~mask]

        # Compute b_n
        # Quick mode if only `p_n` or `sig` has changed
        if quick:
            # Get pre-computed (right) part of the equation
            right = self.w_t_wave_c[i_order]

            # Apply new p_n
            pixel_mapping = diags(trace_profile_n).dot(right)

        else:
            # First (T * lam) for the convolve axis (n_k_c)
            product = (throughput_n * wave_grid)[slice(*i_bnds)]

            # then convolution
            product = diags(product).dot(kernel_n)

            # then weights
            product = weights_n.dot(product)

            # Save this product for quick mode
            self._set_w_t_wave_c(i_order, product)

            # Then spatial profile
            pixel_mapping = diags(trace_profile_n).dot(product)

        # Save new pixel mapping matrix.
        self.pixel_mapping[i_order] = pixel_mapping

        return pixel_mapping

    def build_sys(self, data, error):
        """
        Build linear system arising from the logL maximisation.

        Parameters
        ----------
        data : (N, M) array_like
            A 2-D array of real values representing the detector image.
        error : (N, M) array_like
            Estimate of the error on each pixel.

        Returns
        -------
        scipy.sparse.csr_matrix, array[float]
            A, b from Ax = b being the system to solve.
        """
        # Get the detector model
        b_matrix, data = self.get_detector_model(data, error)

        # (B_T * B) * f = (data/sig)_T * B
        # (matrix ) * f = result
        matrix = b_matrix.T.dot(b_matrix)
        result = data.dot(b_matrix)

        return matrix, result.toarray().squeeze()

    def get_detector_model(self, data, error):
        """
        Get the linear model of the detector pixel, B.dot(flux) = pixels.

        Parameters
        ----------
        data : array_like
            A 2-D array of real values representing the detector image, shape (N, M).
        error : array_like
            Estimate of the error on each pixel, shape (N, M).

        Returns
        -------
        B, pix_array : array[float]
            From the linear equation B.dot(flux) = pix_array
        """
        # Check if `w_t_wave_c` is pre-computed
        quick = self.w_t_wave_c is not None

        # Build matrix B
        # Initiate with empty matrix
        n_i = (~self.mask).sum()  # n good pixels
        b_matrix = csr_matrix((n_i, self.n_wavepoints))

        # Sum over orders
        for i_order in range(self.n_orders):
            # Get sparse pixel mapping matrix.
            b_matrix += self.get_pixel_mapping(i_order, error, quick=quick)

        # Build detector pixels' array
        # Take only valid pixels and apply `error` on data
        data = data[~self.mask] / error[~self.mask]

        return b_matrix, csr_matrix(data)

    @property
    def tikho_mat(self):  # numpydoc ignore=RT01
        """Return the Tikhonov matrix, computing it if needed."""
        if self._tikho_mat is not None:
            return self._tikho_mat

        self._tikho_mat = atoca_utils.finite_first_d(self.wave_grid)
        return self._tikho_mat

    @tikho_mat.setter
    def tikho_mat(self, t_mat):
        self._tikho_mat = t_mat

    def estimate_tikho_factors(self, flux_estimate):
        """
        Estimate an initial guess of the Tikhonov factor.

        The output factor will be used to find the best Tikhonov factor.
        The flux_estimate is used to generate a factor_guess.
        The user should construct a grid with this output in log space,
        e.g. np.logspace(np.log10(flux_estimate)-4, np.log10(flux_estimate)+4, 9).

        Parameters
        ----------
        flux_estimate : callable
            Estimate of the underlying flux (the solution f_k). Must be function
            of wavelengths and it will be projected on self.wave_grid.

        Returns
        -------
        float
            Estimated Tikhonov factor.
        """
        # Get some values from the object
        mask, wave_grid = self.get_attributes("mask", "wave_grid")

        # Find the number of valid pixels
        n_pixels = (~mask).sum()

        # Project the estimate on the wavelength grid
        estimate_on_grid = flux_estimate(wave_grid)

        # Estimate the norm-2 of the regularization term
        reg_estimate = self.tikho_mat.dot(estimate_on_grid)
        reg_estimate = np.nansum(np.array(reg_estimate) ** 2)

        # Estimate of the factor
        factor_guess = (n_pixels / reg_estimate) ** 0.5

        log.info(f"First guess of tikhonov factor: {factor_guess}")

        return factor_guess

    def get_tikho_tests(self, factors, data, error):
        """
        Test different factors for Tikhonov regularization.

        Parameters
        ----------
        factors : 1D list or array-like
            Factors to be tested.
        data : (N, M) array_like
            A 2-D array of real values representing the detector image.
        error : (N, M) array_like
            Estimate of the error on each pixel. Same shape as `data`.

        Returns
        -------
        tests : dict
            Dictionary of the test results
        """
        # Build the system to solve
        b_matrix, pix_array = self.get_detector_model(data, error)

        tikho = atoca_utils.Tikhonov(b_matrix, pix_array, self.tikho_mat)

        # Test all factors
        tests = tikho.test_factors(factors)

        # Save also grid
        tests["grid"] = self.wave_grid

        return tests

    def best_tikho_factor(self, tests, fit_mode):
        """
        Compute the best scale factor for Tikhonov regularization.

        The scale factor is determined by taking the factor giving the lowest reduced chi2 on
        the detector, the highest curvature of the l-curve or when the improvement
        on the chi2 (so the derivative of the chi2, 'd_chi2') reaches a certain threshold.

        Parameters
        ----------
        tests : dictionary
            Results of Tikhonov extraction tests for different factors.
            Must have the keys "factors" and "-logl".
        fit_mode : str
            Which mode is used to find the best Tikhonov factor. Options are
            'all', 'curvature', 'chi2', 'd_chi2'. If 'all' is chosen, the best of the
            three other options will be selected.

        Returns
        -------
        best_fac : float
            The best Tikhonov factor.
        """
        # Modes to be tested
        if fit_mode == "all":
            # Test all modes
            list_mode = ["curvature", "chi2", "d_chi2"]
        else:
            # Single mode
            list_mode = [fit_mode]

        # Evaluate best factor with different methods
        results = {}
        for mode in list_mode:
            best_fac = tests.best_factor(mode=mode)
            results[mode] = best_fac

        if fit_mode == "all":
            # Choose the best factor.
            # In a well-behaved case, the results should be ordered as 'chi2', 'd_chi2', 'curvature'
            # and 'd_chi2' will be the best criterion determine the best factor.
            # 'chi2' usually overfits the solution and 'curvature' may oversmooth the solution
            if results["curvature"] <= results["chi2"] or results["d_chi2"] <= results["chi2"]:
                # In this case, 'chi2' is likely to not overfit the solution, so must be favored
                best_mode = "chi2"
            elif results["curvature"] < results["d_chi2"]:
                # Take the smaller factor between 'd_chi2' and 'curvature'
                best_mode = "curvature"
            elif results["d_chi2"] <= results["curvature"]:
                best_mode = "d_chi2"
            else:
                msg = "Logic error in comparing methods for best Tikhonov factor."
                log.critical(msg)
                raise ValueError(msg)
        else:
            best_mode = fit_mode

        # Get the factor of the chosen mode
        best_fac = results[best_mode]
        log.debug(f"Mode chosen to find regularization factor is {best_mode}")

        return best_fac

    def rebuild(self, spectrum, fill_value=0.0):
        """
        Build current model image of the detector.

        Parameters
        ----------
        spectrum : callable or array-like
            Flux as a function of wavelength if callable
            or array of flux values corresponding to self.wave_grid.
        fill_value : float or np.nan, optional
            Pixel value where the detector is masked. Default is 0.0.

        Returns
        -------
        array[float]
            The modeled detector image.
        """
        # If flux is callable, evaluate on the wavelength grid.
        if callable(spectrum):
            spectrum = spectrum(self.wave_grid)

        # Iterate over all orders
        i_orders = range(self.n_orders)

        # Evaluate the detector model.
        model = np.zeros(self.data_shape)
        for i_order in i_orders:
            # Compute the pixel mapping matrix (b_n) for the current order.
            pixel_mapping = self.get_pixel_mapping(i_order, error=None)

            # Evaluate the model of the current order.
            model[~self.mask] += pixel_mapping.dot(spectrum)

        # Assign masked values
        model[self.mask] = fill_value
        return model

    def compute_likelihood(self, spectrum, data, error):
        """
        Return the log likelihood associated with a particular spectrum.

        Parameters
        ----------
        spectrum : array[float] or callable
            Flux as a function of wavelength if callable
            or array of flux values corresponding to self.wave_grid.
        data : (N, M) array_like
            A 2-D array of real values representing the detector image.
        error : (N, M) array_like
            Estimate of the error on each pixel.
            Same shape as `data`.

        Returns
        -------
        array[float]
            The log-likelihood of the spectrum.
        """
        # Evaluate the model image for the spectrum.
        model = self.rebuild(spectrum)

        # Compute the log-likelihood for the spectrum.
        with np.errstate(divide="ignore"):
            logl = (model - data) / error
        return -np.nansum((logl[~self.mask]) ** 2)

    @staticmethod
    def _solve(matrix, result):
        """
        Solve the linear system using `scipy.spsolve`.

        Parameters
        ----------
        matrix : (N, M) array_like
            A 2-D array of real values representing the detector image.
        result : (N, M) array_like
            The right-hand side of the linear system.

        Returns
        -------
        array[float]
            Solution of the linear system
        """
        # Get valid indices
        idx = np.nonzero(result)[0]

        # Init solution with NaNs.
        sln = np.ones(result.shape[-1]) * np.nan

        # Only solve for valid indices, i.e. wavelengths that are
        # covered by the pixels on the detector.
        # It will be a singular matrix otherwise.
        matrix = matrix[idx, :][:, idx]
        sln[idx] = atoca_utils.try_solve_two_methods(matrix, result[idx])
        return sln

    @staticmethod
    def _solve_tikho(matrix, result, t_mat, **kwargs):
        """
        Solve system using Tikhonov regularization.

        Parameters
        ----------
        matrix : (N, M) array_like
            A 2-D array of real values representing the detector image.
        result : (N, M) array_like
            The right-hand side of the linear system.
        t_mat : (N, M) array_like
            The Tikhonov matrix.
        **kwargs : dict
            Keyword arguments to pass to the solver.

        Returns
        -------
        array[float]
            Solution of the linear system
        """
        # Note that the indexing is applied inside the function
        tikho = atoca_utils.Tikhonov(matrix, result, t_mat)

        return tikho.solve(**kwargs)

    def __call__(self, data, error, tikhonov=False, factor=None):
        """
        Extract underlying flux on the detector.

        Performs an overlapping extraction of the form:
        (B_T * B) * f = (data/sig)_T * B
        where B is a matrix and f is an array.
        The matrix multiplication B * f is the 2d model of the detector.
        We want to solve for the array f.
        The elements of f are labelled by 'k'.
        The pixels are labeled by 'i'.
        Every pixel 'i' is covered by a set of 'k' for each order
        of diffraction.

        TIPS: To be quicker, only specify the psf (`p_list`) in kwargs.
              There will be only one matrix multiplication:
              (P/sig).(w.T.lambda.c_n).

        Parameters
        ----------
        data : (N, M) array_like
            A 2-D array of real values representing the detector image.
        error : (N, M) array_like
            Estimate of the error on each pixel`
            Same shape as `data`.
        tikhonov : bool, optional
            Whether to use Tikhonov extraction
            Default is False.
        factor : float, optional
            The Tikhonov factor to use if tikhonov is True

        Returns
        -------
        spectrum (f_k) : array[float]
            Solution of the linear system
        """
        # Solve with the specified solver.
        if tikhonov:
            if factor is None:
                msg = "Please specify tikhonov `factor`."
                log.critical(msg)
                raise ValueError(msg)

            # Build the system to solve
            b_matrix, pix_array = self.get_detector_model(data, error)

            spectrum = self._solve_tikho(b_matrix, pix_array, self.tikho_mat, factor=factor)

        else:
            # Build the system to solve
            matrix, result = self.build_sys(data, error)

            # Only solve for valid range `i_grid` (on the detector).
            # It will be a singular matrix otherwise.
            spectrum = self._solve(matrix, result)

        return spectrum

    def _get_lo_hi(self, grid, wave_p, wave_m, mask):
        """
        Find the lowest (lo) and highest (hi) index of wave_grid for each pixels and orders.

        Parameters
        ----------
        grid : array[float]
            Wave_grid to check.
        wave_p : array[float]
            Wavelengths on the higher side of each pixel.
        wave_m : array[float]
            Wavelengths on the lower side of each pixel.

        Returns
        -------
        lo, hi : array[float]
            Arrays of indices for lowest and highest values.
        """
        log.debug("Computing lowest and highest indices of wave_grid.")

        # Find lower (lo) index in the pixel
        lo = np.searchsorted(grid, wave_m, side="right") - 1

        # Find higher (hi) index in the pixel
        hi = np.searchsorted(grid, wave_p) - 1

        # Set invalid pixels negative
        lo[mask], hi[mask] = -1, -2

        return lo, hi

    def get_mask_wave(self, i_order):
        """
        Generate mask bounded by limits of wavelength grid.

        Parameters
        ----------
        i_order : int
            Order to select the wave_map on which a mask will be generated

        Returns
        -------
        array[bool]
            A mask with True where wave_map is outside the bounds of wave_grid
        """
        attrs = ["wave_p", "wave_m", "i_bounds"]
        wave_p, wave_m, i_bnds = self.get_attributes(*attrs, i_order=i_order)
        wave_min = self.wave_grid[i_bnds[0]]
        wave_max = self.wave_grid[i_bnds[1] - 1]

        return (wave_m < wave_min) | (wave_p > wave_max)

    def get_w(self, i_order):
        """
        Compute integration weights 'k' for each grid point and pixel 'i'.

        These depend on the type of interpolation used, i.e. the order `n`.

        Parameters
        ----------
        i_order : int
            Order to set the value of n in output arrays.

        Returns
        -------
        w_n : array
            2D array of weights at this specific order `n`. The shape is given by:
            (number of pixels, max number of wavelengths covered by a pixel)
        k_n : array
            2D array of the wavelength grid indices corresponding to the weights.
            Same shape as w_n
        """
        log.debug("Computing weights and k.")

        # get order dependent attributes
        attrs = ["wave_p", "wave_m", "mask_ord", "i_bounds"]
        wave_p, wave_m, mask_ord, i_bnds = self.get_attributes(*attrs, i_order=i_order)

        # Use the convolved grid (depends on the order)
        wave_grid = self.wave_grid[i_bnds[0] : i_bnds[1]]
        # Compute the wavelength coverage of the grid
        d_grid = np.diff(wave_grid)

        # Compute only valid pixels
        wave_p, wave_m = wave_p[~self.mask], wave_m[~self.mask]
        ma = mask_ord[~self.mask]

        # Get lo hi
        lo, hi = self._get_lo_hi(wave_grid, wave_p, wave_m, ma)  # Get indices

        # Number of used pixels
        n_i = len(lo)
        i = np.arange(n_i)

        # Define first and last index of wave_grid for each pixel
        k_first, k_last = -1 * np.ones(n_i), -1 * np.ones(n_i)

        # If lowest value close enough to the exact grid value,
        # NOTE: Could be approximately equal to the exact grid
        # value. It would look like that.
        # >>> lo_dgrid = lo
        # >>> lo_dgrid[lo_dgrid==len(d_grid)] = len(d_grid) - 1
        # >>> cond = (grid[lo]-wave_m)/d_grid[lo_dgrid] <= 1.0e-8
        # But let's stick with the exactly equal
        cond = wave_grid[lo] == wave_m

        # special case (no need for lo_i - 1)
        k_first[cond & ~ma] = lo[cond & ~ma]
        wave_m[cond & ~ma] = wave_grid[lo[cond & ~ma]]

        # else, need lo_i - 1
        k_first[~cond & ~ma] = lo[~cond & ~ma] - 1

        # Same situation for highest value. If we follow the note
        # above (~=), the code could look like
        # >>> cond = (wave_p-grid[hi])/d_grid[hi-1] <= 1.0e-8
        # But let's stick with the exactly equal
        cond = wave_p == wave_grid[hi]

        # special case (no need for hi_i - 1)
        k_last[cond & ~ma] = hi[cond & ~ma]
        wave_p[cond & ~ma] = wave_grid[hi[cond & ~ma]]

        # else, need hi_i
        k_last[~cond & ~ma] = hi[~cond & ~ma]

        # Generate array of all k_i. Set to -1 if not valid
        k_n = atoca_utils.arange_2d(k_first, k_last + 1)
        bad = k_n == -1

        # Number of valid k per pixel
        n_k = np.sum(~bad, axis=-1)

        # Compute array of all w_i. Set to np.nan if not valid
        # Initialize
        w_n = np.zeros(k_n.shape, dtype=float)

        ####################
        # 4 different cases
        ####################

        # Valid for every cases
        w_n[:, 0] = wave_grid[k_n[:, 1]] - wave_m
        w_n[i, n_k - 1] = wave_p - wave_grid[k_n[i, n_k - 2]]

        # Case 1, n_k == 2
        case = (n_k == 2) & ~ma
        if case.any():
            log.debug("n_k = 2 in get_w().")

            # if k_i[0] != lo_i
            cond = case & (k_n[:, 0] != lo)
            w_n[cond, 1] += wave_m[cond] - wave_grid[k_n[cond, 0]]

            # if k_i[-1] != hi_i
            cond = case & (k_n[:, 1] != hi)
            w_n[cond, 0] += wave_grid[k_n[cond, 1]] - wave_p[cond]

            # Finally
            part1 = wave_p[case] - wave_m[case]
            part2 = d_grid[k_n[case, 0]]
            w_n[case, :] *= (part1 / part2)[:, None]

        # Case 2, n_k >= 3
        case = (n_k >= 3) & ~ma
        if case.any():
            log.debug("n_k = 3 in get_w().")
            n_ki = n_k[case]
            w_n[case, 1] = wave_grid[k_n[case, 1]] - wave_m[case]
            w_n[case, n_ki - 2] += wave_p[case] - wave_grid[k_n[case, n_ki - 2]]

            # if k_i[0] != lo_i
            cond = case & (k_n[:, 0] != lo)
            nume1 = wave_grid[k_n[cond, 1]] - wave_m[cond]
            nume2 = wave_m[cond] - wave_grid[k_n[cond, 0]]
            deno = d_grid[k_n[cond, 0]]
            w_n[cond, 0] *= nume1 / deno
            w_n[cond, 1] += nume1 * nume2 / deno

            # if k_i[-1] != hi_i
            cond = case & (k_n[i, n_k - 1] != hi)
            n_ki = n_k[cond]
            nume1 = wave_p[cond] - wave_grid[k_n[cond, n_ki - 2]]
            nume2 = wave_grid[k_n[cond, n_ki - 1]] - wave_p[cond]
            deno = d_grid[k_n[cond, n_ki - 2]]
            w_n[cond, n_ki - 1] *= nume1 / deno
            w_n[cond, n_ki - 2] += nume1 * nume2 / deno

        # Case 3, n_k >= 4
        case = (n_k >= 4) & ~ma
        if case.any():
            log.debug("n_k = 4 in get_w().")
            n_ki = n_k[case]
            w_n[case, 1] += wave_grid[k_n[case, 2]] - wave_grid[k_n[case, 1]]
            w_n[case, n_ki - 2] += wave_grid[k_n[case, n_ki - 2]] - wave_grid[k_n[case, n_ki - 3]]

        # Case 4, n_k > 4
        case = (n_k > 4) & ~ma
        if case.any():
            log.debug("n_k > 4 in get_w().")
            i_k = np.indices(k_n.shape)[-1]
            cond = case[:, None] & (2 <= i_k) & (i_k < n_k[:, None] - 2)
            ind1, ind2 = np.where(cond)
            w_n[ind1, ind2] = d_grid[k_n[ind1, ind2] - 1] + d_grid[k_n[ind1, ind2]]

        # Finally, divide w_n by 2
        w_n /= 2.0

        # Make sure invalid values are masked
        w_n[k_n < 0] = np.nan

        return w_n, k_n
