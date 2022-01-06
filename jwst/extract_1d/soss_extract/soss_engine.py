#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TODO remove use of args and kwargs as much as possible for clearer code.

# General imports.
import numpy as np
from scipy.sparse import issparse, csr_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp1d, Akima1DInterpolator
from scipy.optimize import minimize_scalar

# Local imports.
from SOSS.dms import engine_utils

# Plotting.
import matplotlib.pyplot as plt


class _BaseOverlap:  # TODO Merge with TrpzOverlap?
    """
    Base class for overlaping extraction of the form:
    (B_T * B) * f = (data/sig)_T * B
    where B is a matrix and f is an array.
    The matrix multiplication B * f is the 2d model of the detector.
    We want to solve for the array f.
    The elements of f are labelled by 'k'.
    The pixels are labeled by 'i'.
    Every pixel 'i' is covered by a set of 'k' for each order
    of diffraction.
    The classes inheriting from this class should specify the
    methods get_w which computes the 'k' associated to each pixel 'i'.
    These depends of the type of interpolation used.
    """
    def __init__(self, wave_map, aperture, throughput, kernels,  # TODO rename aperture
                 orders=None, global_mask=None,
                 wave_grid=None, wave_bounds=None, n_os=2,
                 threshold=1e-5, c_kwargs=None,
                 verbose=False):
        """
        Parameters
        ----------
        wave_map : (N_ord, N, M) list or array of 2-D arrays
            A list or array of the central wavelength position for each
            order on the detector.
            It has to have the same (N, M) as `data`.
        aperture : (N_ord, N, M) list or array of 2-D arrays
            A list or array of the spatial profile for each order
            on the detector. It has to have the same (N, M) as `data`.
        throughput : (N_ord [, N_k]) list of array or callable
            A list of functions or array of the throughput at each order.
            If callable, the functions depend on the wavelength.
            If array, projected on `wave_grid`.
        kernels : array, callable or sparse matrix
            Convolution kernel to be applied on the spectrum (f_k) for each orders.
            Can be array of the shape (N_ker, N_k_c).
            Can be a callable with the form f(x, x0) where x0 is
            the position of the center of the kernel. In this case, it must
            return a 1D array (len(x)), so a kernel value
            for each pairs of (x, x0). If array or callable,
            it will be passed to `convolution.get_c_matrix` function
            and the `c_kwargs` can be passed to this function.
            If sparse, the shape has to be (N_k_c, N_k) and it will
            be used directly. N_ker is the length of the effective kernel
            and N_k_c is the length of the spectrum (f_k) convolved.
        global_mask : (N, M) array_like boolean, optional
            Boolean Mask of the detector pixels to mask for every extraction.
        orders: list, optional:
            List of orders considered. Default is orders = [1, 2]
        wave_grid : (N_k) array_like, optional
            The grid on which f(lambda) will be projected.
            Default is a grid from `utils.get_soss_grid`.
            `n_os` will be passed to this function.
        wave_bounds : list or array-like (N_ord, 2), optional
            Boundary wavelengths covered by each orders.
            Default is the wavelength covered by `wave_map`.
        n_os  : int, optional
            if `wave_grid`is None, it will be used to generate
            a grid. Default is 2.
        threshold : float, optional:
            The pixels where the estimated spatial profile is less than
            this value will be masked. Default is 1e-5.
        c_kwargs : list of N_ord dictionnaries or dictionnary, optional
            Inputs keywords arguments to pass to
            `convolution.get_c_matrix` function for each orders.
            If dictionnary, the same c_kwargs will be used for each orders.
        verbose : bool, optional
            Print steps. Default is False.
        """

        # If no orders specified extract on orders 1 and 2.
        if orders is None:
            orders = [1, 2]

        ###########################
        # Save basic parameters
        ###########################

        # Spectral orders and number of orders.
        self.data_shape = wave_map[0].shape
        self.orders = orders
        self.n_orders = len(orders)
        self.threshold = threshold
        self.verbose = verbose

        # Raise error if the number of orders is not consistent.
        if self.n_orders != len(wave_map):
            msg = ("The number of orders specified {} and the number of "
                   "wavelength maps provided {} do not match.")
            raise ValueError(msg.format(self.n_orders, len(wave_map)))

        # Detector image.
        self.data = np.full(self.data_shape, fill_value=np.nan)

        # Error map of each pixels.
        self.error = np.ones(self.data_shape)

        # Set all reference file quantities to None.
        self.wave_map = None
        self.aperture = None
        self.throughput = None
        self.kernels = None

        # Set the wavelength map and aperture for each order.
        self.update_wave_map(wave_map)
        self.update_aperture(aperture)

        # Generate a wavelength grid if none was provided. TODO Requires self.aperture self.wave_map
        if wave_grid is None:

            if self.n_orders == 2:  # TODO should this be mandatory input.
                wave_grid = engine_utils.get_soss_grid(wave_map, aperture, n_os=n_os)  # TODO check difference between get_soss_grid and grid_from_map
            else:
                wave_grid, _ = self.grid_from_map()

        # Set the wavelength grid and its size.
        self.wave_grid = wave_grid.copy()
        self.n_wavepoints = len(wave_grid)

        # Set the throughput for each order.
        self.update_throughput(throughput)  # TODO requires self.wave_grid

        ###################################
        # Build detector mask
        ###################################

        # Assign a first estimate of i_bounds to be able to compute mask.
        self.i_bounds = [[0, len(wave_grid)] for _ in range(self.n_orders)]  # TODO double check how the i_bounds and mask interact.

        # First estimate of a global mask and masks for each orders
        self.mask, self.mask_ord = self._get_masks(global_mask)

        # Correct i_bounds if it was not specified
        self.i_bounds = self._get_i_bnds(wave_bounds)

        # Re-build global mask and masks for each orders
        self.mask, self.mask_ord = self._get_masks(global_mask)

        # Save mask here as the general mask,
        # since `mask` attribute can be changed.
        self.general_mask = self.mask.copy()

        ####################################
        # Build convolution matrix
        ####################################
        self.update_kernels(kernels, c_kwargs)  # TODO requires self.wave_grid self.i_bounds

        #############################
        # Compute integration weights
        #############################
        # The weights depend on the integration method used solve
        # the integral of the flux over a pixel and are encoded
        # in the class method `get_w()`.
        self.weights, self.weights_k_idx = self.compute_weights()  # TODO put shapes in commments, name of indices.

        #########################
        # Save remaining inputs
        #########################

        # Set masked values to zero. TODO may not be necessary.
        self.data[self.mask] = 0

        # Init the pixel mapping (b_n) matrices. Matrices that transforms the 1D spectrum to a the image pixels.
        self.pixel_mapping = [None for _ in range(self.n_orders)]
        self.i_grid = None
        self.tikho = None
        self.tikho_mat = None
        self.w_t_wave_c = None

        return

    def verbose_print(self, *args, **kwargs):
        """Print if verbose is True. Same as `print` function."""

        if self.verbose:
            print(*args, **kwargs)

        return

    def get_attributes(self, *args, i_order=None):
        """Return list of attributes

        Parameters
        ----------
        args: str
            All attributes to return.
        i_order: None or int, optionoal
            Index of order to extract. If specified, it will
            be applied to all attributes in args, so it cannot
            be mixed with non-order dependent attributes).
        """

        if i_order is None:
            out = [getattr(self, arg) for arg in args]
        else:
            out = [getattr(self, arg)[i_order] for arg in args]

        if len(out) == 1:
            out = out[0]

        return out

    def update_wave_map(self, wave_map):

        self.wave_map = [wave_n.copy() for wave_n in wave_map]  # TODO make dict with order number as key.

        return

    def update_aperture(self, aperture):
        """Update the aperture maps."""

        # Update the aperture profile.
        self.aperture = [aperture_n.copy() for aperture_n in aperture]  # TODO make dict with order number as key.

        return

    def update_throughput(self, throughput):
        """Update the throughput values."""

        # Update the throughput values.
        throughput_new = []  # TODO make dict with order number as key.
        for throughput_n in throughput:  # Loop over orders.

            if callable(throughput_n):

                # Througput was given as a callable function.
                throughput_new.append(throughput_n(self.wave_grid))

            elif throughput_n.shape == self.wave_grid.shape:

                # Throughput was given as an array.
                throughput_new.append(throughput_n)

            else:
                msg = 'Throughputs must be given as callable or arrays matching the extraction grid.'
                raise ValueError(msg)

        # Set the attribute to the new values.
        self.throughput = throughput_new

        return

    def update_kernels(self, kernels, c_kwargs):

        # Verify the c_kwargs. TODO Be explict here?
        if c_kwargs is None:
            c_kwargs = [{} for _ in range(self.n_orders)]

        elif isinstance(c_kwargs, dict):
            c_kwargs = [c_kwargs for _ in range(self.n_orders)]

        # Define convolution sparse matrix. TODO make dict with order number as key.
        kernels_new = []
        for i_order, kernel_n in enumerate(kernels):

            if not issparse(kernel_n):
                kernel_n = engine_utils.get_c_matrix(kernel_n, self.wave_grid,
                                                     i_bounds=self.i_bounds[i_order],
                                                     **c_kwargs[i_order])

            kernels_new.append(kernel_n)

        self.kernels = kernels_new

        return

    def get_mask_wave(self, i_order):
        """Mask according to wavelength grid """

        wave = self.wave_map[i_order]
        imin, imax = self.i_bounds[i_order]
        wave_min = self.wave_grid[imin]
        wave_max = self.wave_grid[imax - 1]  # TODO change so -1 not needed?

        mask = (wave <= wave_min) | (wave >= wave_max)

        return mask

    def _get_masks(self, global_mask):
        """
        Compute a general mask on the detector and for each orders.
        Depends on the spatial profile, the wavelength grid
        and the user defined mask (optional). These are all specified
        when initiating the object.
        """

        # Get needed attributes
        threshold, n_orders = self.get_attributes('threshold', 'n_orders')
        throughput, aperture, wave_map = self.get_attributes('throughput', 'aperture', 'wave_map')

        # Mask according to the spatial profile.
        mask_aperture = np.array([aperture_n < threshold for aperture_n in aperture])

        # Mask pixels not covered by the wavelength grid.
        mask_wave = np.array([self.get_mask_wave(i_order) for i_order in range(n_orders)])

        # Apply user defined mask.
        if global_mask is None:
            mask_ord = np.any([mask_aperture, mask_wave], axis=0)
        else:
            mask = [global_mask for _ in range(n_orders)]  # For each orders
            mask_ord = np.any([mask_aperture, mask_wave, mask], axis=0)

        # Find pixels that are masked in each order.
        general_mask = np.all(mask_ord, axis=0)

        # Mask pixels if mask_aperture not masked but mask_wave is.
        # This means that an order is contaminated by another
        # order, but the wavelength range does not cover this part
        # of the spectrum. Thus, it cannot be treated correctly.
        general_mask |= (np.any(mask_wave, axis=0)
                         & np.all(~mask_aperture, axis=0))

        # Apply this new general mask to each orders.
        mask_ord = (mask_wave | general_mask[None, :, :])

        return general_mask, mask_ord

    def update_mask(self, mask):
        """
        Update `mask` attribute by completing the
        `general_mask` attribute with the input `mask`.
        Everytime the mask is changed, the integration weights
        need to be recomputed since the pixels change.
        """

        # Get general mask
        general_mask = self.general_mask

        # Complete with the input mask
        new_mask = (general_mask | mask)

        # Update attribute
        self.mask = new_mask

        # Correct i_bounds if it was not specified
        # self.update_i_bnds()

        # Re-compute weights
        self.weights, self.weights_k_idx = self.compute_weights()

        return

    def _get_i_bnds(self, wave_bounds=None):
        """
        Define wavelength boundaries for each orders using the order's mask.
        """

        wave_grid = self.wave_grid
        i_bounds = self.i_bounds

        # Check if wave_bounds given
        if wave_bounds is None:
            wave_bounds = []
            for i in range(self.n_orders):
                wave = self.wave_map[i][~self.mask_ord[i]]
                wave_bounds.append([wave.min(), wave.max()])

        # What we need is the boundary position
        # on the wavelength grid.
        i_bnds_new = []
        for bounds, i_bnds in zip(wave_bounds, i_bounds):

            a = np.min(np.where(wave_grid >= bounds[0])[0])
            b = np.max(np.where(wave_grid <= bounds[1])[0]) + 1

            # Take the most restrictive bound
            a = np.maximum(a, i_bnds[0])
            b = np.minimum(b, i_bnds[1])

            # Keep value
            i_bnds_new.append([a, b])

        return i_bnds_new

    def update_i_bnds(self):
        """Update the grid limits for the extraction.
        Needs to be done after modification of the mask
        """

        # Get old and new boundaries.
        i_bnds_old = self.i_bounds
        i_bnds_new = self._get_i_bnds()

        for i_order in range(self.n_orders):

            # Take most restrictive lower bound.
            low_bnds = [i_bnds_new[i_order][0], i_bnds_old[i_order][0]]
            i_bnds_new[i_order][0] = np.max(low_bnds)

            # Take most restrictive upper bound.
            up_bnds = [i_bnds_new[i_order][1], i_bnds_old[i_order][1]]
            i_bnds_new[i_order][1] = np.min(up_bnds)

        # Update attribute.
        self.i_bounds = i_bnds_new

        return

    def wave_grid_c(self, i_order):
        """
        Return wave_grid for the convolved flux at a given order.
        """

        index = slice(*self.i_bounds[i_order])

        return self.wave_grid[index]

    def get_w(self, i_order):
        """Dummy method to be able to init this class"""

        return np.array([]), np.array([])

    def compute_weights(self):
        """
        Compute integration weights

        The weights depend on the integration method used solve
        the integral of the flux over a pixel and are encoded
        in the class method `get_w()`.

        Returns the lists of weights and corresponding grid indices
        """

        # Init lists
        weights, weights_k_idx = [], []
        for i_order in range(self.n_orders):  # For each orders

            weights_n, k_idx_n = self.get_w(i_order)  # Compute weigths

            # Convert to sparse matrix
            # First get the dimension of the convolved grid
            n_kc = np.diff(self.i_bounds[i_order]).astype(int)[0]

            # Then convert to sparse
            weights_n = engine_utils.sparse_k(weights_n, k_idx_n, n_kc)
            weights.append(weights_n), weights_k_idx.append(k_idx_n)

        return weights, weights_k_idx

    def _set_w_t_wave_c(self, i_order, product):  # TODO better name? prod_mat tmp_mat, intermediate_mat?
        """
        Save the matrix product of the weighs (w), the throughput (t),
        the wavelength (lam) and the convolution matrix for faster computation.
        """

        if self.w_t_wave_c is None:
            self.w_t_wave_c = [[] for _ in range(self.n_orders)]  # TODO make dict with order number as key.

        # Assign value
        self.w_t_wave_c[i_order] = product.copy()

        return

    def grid_from_map(self, i_order=0):
        """
        Return the wavelength grid and the columns associated
        to a given order index (i_order)
        """

        attrs = ['wave_map', 'aperture']
        wave_map, aperture = self.get_attributes(*attrs, i_order=i_order)

        wave_grid, icol = engine_utils._grid_from_map(wave_map, aperture, out_col=True)

        return wave_grid, icol

    def get_adapt_grid(self, spectrum=None, n_max=3, **kwargs):
        """
        Return an irregular grid needed to reach a
        given precision when integrating over each pixels.

        Parameters (all optional)
        ----------
        spectrum (f_k): 1D array-like
            Input flux in the integral to be optimized.
            f_k is the projection of the flux on self.wave_grid
        n_max: int (n_max > 0)
            Maximum number of nodes in each intervals of self.wave_grid.
            Needs to be greater then zero.

        kwargs (arguments passed to the function get_n_nodes)
        ------
        tol, rtol : float, optional
            The desired absolute and relative tolerances. Defaults are 1.48e-4.
        divmax : int, optional
            Maximum order of extrapolation. Default is 10.

        Returns
        -------
        os_grid  : 1D array
            Oversampled grid which minimizes the integration error based on
            Romberg's method
        See Also
        --------
        utils.get_n_nodes
        scipy.integrate.quadrature.romberg
        References
        ----------
        [1] 'Romberg's method' https://en.wikipedia.org/wiki/Romberg%27s_method

        """
        # Generate the spectrum (f_k) if not given.
        if spectrum is None:
            spectrum = self.extract()

        # Init output oversampled grid
        os_grid = []

        # Iterate starting with the last order
        for i_order in range(self.n_orders - 1, -1, -1):  # TODO easier way of inverse loop?

            # Grid covered by this order
            grid_ord = self.wave_grid_c(i_order)

            # Estimate the flux at this order
            convolved_spectrum = self.kernels[i_order].dot(spectrum)
            # Interpolate with a cubic spline
            fct = interp1d(grid_ord, convolved_spectrum, kind='cubic')

            # Find number of nodes to reach the precision
            n_oversample, _ = engine_utils.get_n_nodes(grid_ord, fct, **kwargs)

            # Make sure n_oversample is not greater than
            # user's define `n_max`
            n_oversample = np.clip(n_oversample, 0, n_max)

            # Generate oversampled grid
            grid_ord = engine_utils.oversample_grid(grid_ord, n_os=n_oversample)

            # Keep only wavelength that are not already
            # covered by os_grid.
            if os_grid:
                # Under or above os_grid
                index = (grid_ord < np.min(os_grid))
                index |= (grid_ord > np.max(os_grid))
            else:
                index = slice(None)

            # Keep these values
            os_grid.append(grid_ord[index])

        # Convert os_grid to 1D array
        os_grid = np.concatenate(os_grid)

        # Return sorted and unique.
        wave_grid = np.unique(os_grid)

        return wave_grid

    def estimate_noise(self, i_order=0, data=None, error=None, mask=None):
        """
        Relative noise estimate over columns.

        Parameters
        ----------
        i_order: int, optional
            index of diffraction order. Default is 0
        data: 2d array, optional
            map of the detector image
            Default is `self.data`.
        error: 2d array, optional
            map of the estimate of the detector noise.
            Default is `self.sig`
        mask: 2d array, optional
            Bool map of the masked pixels for order `i_order`.
            Default is `self.mask_ord[i_order]`

        Returns
        ------
        wave_grid, noise
        """

        # Use object attributes if not given
        if data is None:
            data = self.data

        if error is None:
            error = self.error

        if mask is None:
            mask = self.mask_ord[i_order]

        # Compute noise estimate only on the trace (mask the rest)
        noise = np.ma.array(error, mask=mask)

        # RMS over columns
        noise = np.sqrt((noise**2).sum(axis=0))

        # Relative
        noise /= np.ma.array(data, mask=mask).sum(axis=0)

        # Convert to array with nans
        noise = noise.filled(fill_value=np.nan)

        # Get associated wavelengths
        wave_grid, i_col = self.grid_from_map(i_order)

        # Return sorted according to wavelenghts
        return wave_grid, noise[i_col]

    def get_pixel_mapping(self, i_order, same=False, error=True, quick=False):
        """
        Compute the matrix `b_n = (P/sig).w.T.lambda.c_n` ,
        where `P` is the spatial profile matrix (diag),
        `w` is the integrations weights matrix,
        `T` is the throughput matrix (diag),
        `lambda` is the convolved wavelength grid matrix (diag),
        `c_n` is the convolution kernel.
        The model of the detector at order n (`model_n`)
        is given by the system:
        model_n = b_n.c_n.f ,
        where f is the incoming flux projected on the wavelenght grid.
        This methods updates the `b_n_list` attribute.
        Parameters
        ----------
        i_order: integer
            Label of the order (depending on the initiation of the object).
        same: bool, optional
            Do not recompute, b_n. Take the last b_n computed.
            Useful to speed up code. Default is False.
        error: bool or (N, M) array_like, optional
            If 2-d array, `sig` is the new error estimation map.
            It is the same shape as `sig` initiation input. If bool,
            wheter to apply sigma or not. The method will return
            b_n/sigma if True or array_like and b_n if False. If True,
            the default object attribute `sig` will be use.
        quick: bool, optional
            If True, only perform one matrix multiplication
            instead of the whole system: (P/sig).(w.T.lambda.c_n)

        Returns
        ------
        sparse matrix of b_n coefficients
        """

        # Force to compute if b_n never computed.
        if self.pixel_mapping[i_order] is None:
            same = False

        # Take the last b_n computed if nothing changes
        if same:
            pixel_mapping = self.pixel_mapping[i_order]

        else:
            pixel_mapping = self._get_pixel_mapping(i_order, error=error, quick=quick)

            # Save new pixel mapping matrix.
            self.pixel_mapping[i_order] = pixel_mapping

        return pixel_mapping

    def _get_pixel_mapping(self, i_order, error=True, quick=False):  # TODO merge with get_pixel_mapping?
        """
        Compute the matrix `b_n = (P/sig).w.T.lambda.c_n` ,
        where `P` is the spatial profile matrix (diag),
        `w` is the integrations weights matrix,
        `T` is the throughput matrix (diag),
        `lambda` is the convolved wavelength grid matrix (diag),
        `c_n` is the convolution kernel.
        The model of the detector at order n (`model_n`)
        is given by the system:
        model_n = b_n.c_n.f ,
        where f is the incoming flux projected on the wavelenght grid.
        Parameters
        ----------
        i_order : integer
            Label of the order (depending on the initiation of the object).
        error: bool or (N, M) array_like, optional
            If 2-d array, `sig` is the new error estimation map.
            It is the same shape as `sig` initiation input. If bool,
            wheter to apply sigma or not. The method will return
            b_n/sigma if True or array_like and b_n if False. If True,
            the default object attribute `sig` will be use.
        quick: bool, optional
            If True, only perform one matrix multiplication
            instead of the whole system: (P/sig).(w.T.lambda.c_n)

        Returns
        ------
        sparse matrix of b_n coefficients
        """

        # Special treatment for error map
        # Can be bool or array.
        if error is False:
            # Sigma will have no effect
            error = np.ones(self.data_shape)
        else:
            if error is not True:
                # Sigma must be an array so
                # update object attribute
                self.error = error.copy()

            # Take sigma from object
            error = self.error

        # Get needed attributes ...
        attrs = ['wave_grid', 'mask']
        wave_grid, mask = self.get_attributes(*attrs)

        # ... order dependent attributes
        attrs = ['aperture', 'throughput', 'kernels', 'weights', 'i_bounds']
        aperture_n, throughput_n, kernel_n, weights_n, i_bnds = self.get_attributes(*attrs, i_order=i_order)

        # Keep only valid pixels (P and sig are still 2-D)
        # And apply direcly 1/sig here (quicker)
        aperture_n = aperture_n[~mask] / error[~mask]

        # Compute b_n
        # Quick mode if only `p_n` or `sig` has changed
        if quick:
            # Get pre-computed (right) part of the equation
            right = self.w_t_wave_c[i_order]

            # Apply new p_n
            pixel_mapping = diags(aperture_n).dot(right)

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
            pixel_mapping = diags(aperture_n).dot(product)

        return pixel_mapping

    def get_i_grid(self, d):
        """ Return the index of the grid that are well defined, so d != 0 """

        if self.i_grid is None:  # TODO Shouldn't this update even if the attribute is already set?
            self.i_grid = np.nonzero(d)[0]

        return self.i_grid

    def build_sys(self, data=None, error=True, mask=None, aperture=None, throughput=None):
        """
        Build linear system arising from the logL maximisation.
        TIPS: To be quicker, only specify the psf (`p_list`) in kwargs.
              There will be only one matrix multiplication:
              (P/sig).(w.T.lambda.c_n).
        Parameters
        ----------
        data : (N, M) array_like, optional
            A 2-D array of real values representing the detector image.
            Default is the object attribute `data`.
        error: bool or (N, M) array_like, optional
            Estimate of the error on each pixel.
            If 2-d array, `sig` is the new error estimation map.
            It is the same shape as `sig` initiation input. If bool,
            wheter to apply sigma or not. The method will return
            b_n/sigma if True or array_like and b_n if False. If True,
            the default object attribute `sig` will be use.
        mask : (N, M) array_like boolean, optional
            Additionnal mask for a given exposure. Will be added
            to the object general mask.
        aperture : (N_ord, N, M) list or array of 2-D arrays, optional
            A list or array of the spatial profile for each order
            on the detector. It has to have the same (N, M) as `data`.
            Default is the object attribute `p_list`
        throughput : (N_ord [, N_k]) list or array of functions, optional
            A list or array of the throughput at each order.
            The functions depend on the wavelength
            Default is the object attribute `t_list`

        Returns
        ------
        A and b from Ax = b beeing the system to solve.
        """

        # Check if inputs are suited for quick mode;
        # Quick mode if `t_list` is not specified.
        quick = (throughput is None)

        # and if mask doesn't change
        quick &= (mask is None)
        quick &= (self.w_t_wave_c is not None)  # Pre-computed
        if quick:
            self.verbose_print('Quick mode is on!')

        # Use data from object as default
        if data is None:
            data = self.data
        else:
            # Update data
            self.data = data

        # Update mask if given
        if mask is not None:
            self.update_mask(mask)

        # Take (updated) mask from object
        mask = self.mask

        # Get some dimensions infos
        n_wavepoints, n_orders = self.n_wavepoints, self.n_orders

        # Update aperture maps and throughput values.
        if aperture is not None:
            self.update_aperture(aperture)

        if throughput is not None:
            self.update_throughput(throughput)

        # Calculations

        # Build matrix B
        # Initiate with empty matrix
        n_i = (~mask).sum()  # n good pixels
        b_matrix = csr_matrix((n_i, n_wavepoints))

        # Sum over orders
        for i_order in range(n_orders):

            # Get sparse pixel mapping matrix.
            b_matrix += self.get_pixel_mapping(i_order, error=error, quick=quick)

        # Build system
        # Fisrt get `sig` which have been update`
        # when calling `get_b_n`
        error = self.error

        # Take only valid pixels and apply `error` on data
        data = data[~mask]/error[~mask]

        # (B_T * B) * f = (data/sig)_T * B
        # (matrix ) * f = result
        matrix = b_matrix.T.dot(b_matrix)
        result = csr_matrix(data.T).dot(b_matrix)

        return matrix, result.toarray().squeeze()

    def set_tikho_matrix(self, t_mat=None, t_mat_func=None,
                         fargs=None, fkwargs=None):
        """
        Set the tikhonov matrix attribute.
        The matrix can be directly specified as an input, or
        it can be built using `t_mat_func`

        Parameters
        ----------
        t_mat: matrix-like, optional
            TIkhonov regularisation matrix. scipy.sparse matrix
            are recommended.
        t_mat_func: callable, optional
            Function use to generate `t_mat`is not specified.
            Will take `fargs` and `fkwargs`as imput.
        fargs: tuple, optional
            Arguments passed to `t_mat_func`
        fkwargs: dict, optional
            Keywords arguments passed to `t_mat_func`
        """

        # Generate the matrix with the function
        if t_mat is None:

            # Default function if not specified
            if t_mat_func is None:

                # Use the nyquist sampled gaussian kernel
                t_mat_func = engine_utils.get_nyquist_matrix

            # Default args
            if fargs is None:
                fargs = (self.wave_grid, )
            if fkwargs is None:
                fkwargs = {"integrate": True}

            # Call function
            t_mat = t_mat_func(*fargs, **fkwargs)

        # Set attribute
        self.tikho_mat = t_mat

        return

    def get_tikho_matrix(self, **kwargs):
        """
        Return the tikhonov matrix.
        Generate it with `set_tikho_matrix` method
        if not define yet. If so, all arguments are passed
        to `set_tikho_matrix`. The result is saved as an attribute.
        """

        if self.tikho_mat is None:
            self.set_tikho_matrix(**kwargs)

        return self.tikho_mat

    def get_tikho_tests(self, factors, tikho=None, estimate=None,
                        tikho_kwargs=None, **kwargs):
        """
        Test different factors for Tikhonov regularisation.

        Parameters
        ----------
        factors: 1D list or array-like
            Factors to be tested.
        tikho: Tikhonov object, optional
            Tikhonov regularisation object (see regularisation.Tikhonov).
            If not given, an object will be initiated using the linear system
            from `build_sys` method and kwargs will be passed.
        estimate: 1D array-like, optional
            Estimate of the flux projected on the wavelength grid.
        tikho_kwargs:
            passed to init Tikhonov object. Possible options
            are `t_mat`, `grid` and `verbose`
        data : (N, M) array_like, optional
            A 2-D array of real values representing the detector image.
            Default is the object attribute `data`.
        error : (N, M) array_like, optional
            Estimate of the error on each pixel`
            Same shape as `data`.
            Default is the object attribute `sig`.
        aperture : (N_ord, N, M) list or array of 2-D arrays, optional
            A list or array of the spatial profile for each order
            on the detector. It has to have the same (N, M) as `data`.
            Default is the object attribute `p_list`
        throughput : (N_ord [, N_k]) list or array of functions, optional
            A list or array of the throughput at each order.
            The functions depend on the wavelength
            Default is the object attribute `t_list`

        Returns
        ------
        dictonary of the tests results
        """

        # Build the system to solve
        matrix, result = self.build_sys(**kwargs)

        # Get valid grid index
        i_grid = self.get_i_grid(result)

        if tikho is None:
            t_mat = self.get_tikho_matrix()
            default_kwargs = {'grid': self.wave_grid,
                              'index': i_grid,
                              't_mat': t_mat}
            if tikho_kwargs is None:
                tikho_kwargs = {}
            tikho_kwargs = {**default_kwargs, **tikho_kwargs}
            tikho = engine_utils.Tikhonov(matrix, result, **tikho_kwargs)
            self.tikho = tikho

        # Test all factors
        tests = tikho.test_factors(factors, estimate)

        # Generate logl using solutions for each factors
        logl_list = []

        # Compute b_n only the first iteration. Then
        # use the same value to rebuild the detector.
        same = False
        for sln in tests['solution']:

            # Init the spectrum (f_k) with nan, so it has the adequate shape
            spectrum = np.ones(result.shape[-1]) * np.nan
            spectrum[i_grid] = sln  # Assign valid values
            logl_list.append(self.compute_likelihood(spectrum, same=same))  # log_l
            same = True

        # Save in tikho's tests
        tikho.test['-logl'] = -1 * np.array(logl_list)

        # Save also grid
        tikho.test["grid"] = self.wave_grid[i_grid]
        tikho.test["i_grid"] = i_grid

        return tikho.test

    def best_tikho_factor(self, tests=None, interpolate=True,
                          interp_index=None, i_plot=False):
        """Compute the best scale factor for Tikhonov regularisation.
        It is determine by taking the factor giving the highest logL on
        the detector.

        Parameters
        ----------
        tests: dictionnary, optional
            Results of tikhonov extraction tests
            for different factors.
            Must have the keys "factors" and "-logl".
            If not specified, the tests from self.tikho.tests
            are used.
        interpolate: bool, optional
            If True, use akima spline interpolation
            to find a finer minimum. Default is true.
        interp_index: 2 element list, optional
            Index around the minimum value on the tested factors.
            Will be used for the interpolation.
            For example, if i_min is the position of
            the minimum logL value and [i1, i2] = interp_index,
            then the interpolation will be perform between
            i_min + i1 and i_min + i2 - 1
        i_plot: bool, optional
            Plot the result of the minimization

        Returns
        -------
        Best scale factor (float)
        """

        if interp_index is None:
            interp_index = [-2, 4]

        # Use pre-run tests if not specified
        if tests is None:
            tests = self.tikho.tests

        # Get relevant quantities from tests
        factors = tests["factors"]
        logl = tests["-logl"]

        # Get position of the minimum value
        i_min = np.argmin(logl)

        # Interpolate to get a finer value
        if interpolate:

            # Only around the best value
            i_range = [i_min + d_i for d_i in interp_index]

            # Make sure it's still a valid index
            i_range[0] = np.max([i_range[0], 0])
            i_range[-1] = np.min([i_range[-1], len(logl) - 1])

            # Which index to use
            index = np.arange(*i_range, 1)

            # Akima spline in log space
            x_val, y_val = np.log10(factors[index]), np.log10(logl[index])
            i_sort = np.argsort(x_val)
            x_val, y_val = x_val[i_sort], y_val[i_sort]
            fct = Akima1DInterpolator(x_val, y_val)

            # Find min
            bounds = (x_val.min(), x_val.max())
            opt_args = {"bounds": bounds,
                        "method": "bounded"}
            min_fac = minimize_scalar(fct, **opt_args).x

            # Plot the fit if required
            if i_plot:

                # Original grid
                plt.plot(np.log10(factors), np.log10(logl), ":")

                # Fit sub-grid
                plt.plot(x_val, y_val, ".")

                # Show akima spline
                x_new = np.linspace(*bounds, 100)
                plt.plot(x_new, fct(x_new))

                # Show minimum found
                plt.plot(min_fac, fct(min_fac), "x")

                # Labels
                plt.xlabel(r"$\log_{10}$(factor)")
                plt.ylabel(r"$\log_{10}( - \log L)$")
                plt.tight_layout()

            # Return to linear scale
            min_fac = 10.**min_fac

        # Simply return the minimum value if no interpolation required
        else:
            min_fac = factors[i_min]

        # Return scale factor minimizing the logL
        return min_fac

    def rebuild(self, spectrum=None, i_orders=None, same=False):
        """Build current model image of the detector.

        :param spectrum: flux as a function of wavelength if callable
            or array of flux values corresponding to self.wave_grid.
        :param i_orders: Indices of orders to model. Default is
            all available orders.
        :param same: If True, do not recompute the pixel_mapping matrix (b_n)
            and instead use the most recent pixel_mapping to speed up the computation.
            Default is False.

        :type spectrum: callable or array-like
        :type i_orders: List[int]
        :type same: bool

        :returns: model - the modelled detector image.
        :rtype: array[float]
        """

        # If no spectrum given compute it.
        if spectrum is None:
            spectrum = self.extract()

        # If flux is callable, evaluate on the wavelength grid.
        if callable(spectrum):
            spectrum = spectrum(self.wave_grid)

        # Iterate over all orders by default.
        if i_orders is None:
            i_orders = range(self.n_orders)

        # Get required class attribute.
        mask = self.mask

        # Evaluate the detector model.
        model = np.zeros(self.data_shape)
        for i_order in i_orders:

            # Compute the pixel mapping matrix (b_n) for the current order.
            pixel_mapping = self.get_pixel_mapping(i_order, error=False, same=same)

            # Evaluate the model of the current order.
            model[~mask] += pixel_mapping.dot(spectrum)

        # Set masked values to NaN.
        model[mask] = np.nan

        return model

    def compute_likelihood(self, spectrum=None, same=False):
        """Return the log likelihood asssociated with a particular spectrum.

        :param spectrum: flux as a function of wavelength if callable
            or array of flux values corresponding to self.wave_grid.
            If not given it will be computed by calling self.extract().
        :param same: If True, do not recompute the pixel_mapping matrix (b_n)
            and instead use the most recent pixel_mapping to speed up the computation.
            Default is False.

        :type spectrum: array-like
        :type same: bool

        :return: logl - The log-likelihood of the spectrum.
        :rtype: array[float]

        """

        # If no spectrum given compute it.
        if spectrum is None:
            spectrum = self.extract()

        # Evaluate the model image for the spectrum.
        model = self.rebuild(spectrum, same=same)

        # Get data and error attributes.
        data = self.data
        error = self.error

        # Compute the log-likelihood for the spectrum.
        logl = -np.nansum((model - data)**2/error**2)

        return logl

    @staticmethod
    def _solve(matrix, result, index=slice(None)):
        """
        Simply pass `matrix` and `result`
        to `scipy.spsolve` and apply index.
        """

        return spsolve(matrix[index, :][:, index], result[index])

    @staticmethod
    def _solve_tikho(matrix, result, index=slice(None), **kwargs):
        """Solve system using Tikhonov regularisation"""

        # Note that the indexing is applied inside the function
        return engine_utils.tikho_solve(matrix, result, index=index, **kwargs)

    def extract(self, tikhonov=False, tikho_kwargs=None,  # TODO merge with __call__.
                factor=None, **kwargs):
        """
        Extract underlying flux on the detector.
        All parameters are passed to `build_sys` method.
        TIPS: To be quicker, only specify the psf (`p_list`) in kwargs.
              There will be only one matrix multiplication:
              (P/sig).(w.T.lambda.c_n).

        Parameters
        ----------
        tikhonov : bool, optional
            Wheter to use tikhonov extraction
            (see regularisation.tikho_solve function).
            Default is False.
        tikho_kwargs : dictionnary or None, optional
            Arguments passed to `tikho_solve`.
        factor : the tikhonov factor to use of tikhonov is True
        data : (N, M) array_like, optional
            A 2-D array of real values representing the detector image.
            Default is the object attribute `data`.
        error : (N, M) array_like, optional
            Estimate of the error on each pixel`
            Same shape as `data`.
            Default is the object attribute `sig`.
        mask : (N, M) array_like boolean, optional
            Additionnal mask for a given exposure. Will be added
            to the object general mask.
        aperture : (N_ord, N, M) list or array of 2-D arrays, optional
            A list or array of the spatial profile for each order
            on the detector. It has to have the same (N, M) as `data`.
            Default is the object attribute `p_list`
        throughput : (N_ord [, N_k]) list or array of functions, optional
            A list or array of the throughput at each order.
            The functions depend on the wavelength
            Default is the object attribute `t_list`

        Returns
        -----
        spectrum (f_k): solution of the linear system
        """

        # Build the system to solve
        matrix, result = self.build_sys(**kwargs)

        # Get index of `wave_grid` convered by the pixel.
        # `wave_grid` may cover more then the pixels.
        i_grid = self.get_i_grid(result)

        # Init spectrum with NaNs.
        spectrum = np.ones(result.shape[-1]) * np.nan

        # Solve with the specified solver.
        # Only solve for valid range `i_grid` (on the detector).
        # It will be a singular matrix otherwise.
        if tikhonov:

            if factor is None:
                raise ValueError("Please specify tikhonov `factor`.")

            t_mat = self.get_tikho_matrix()
            default_kwargs = {'grid': self.wave_grid,
                              'index': i_grid,
                              't_mat': t_mat,
                              'factor': factor}

            if tikho_kwargs is None:
                tikho_kwargs = {}

            tikho_kwargs = {**default_kwargs, **tikho_kwargs}
            spectrum[i_grid] = self._solve_tikho(matrix, result, **tikho_kwargs)

        else:
            spectrum[i_grid] = self._solve(matrix, result, index=i_grid)

        return spectrum

    def __call__(self, **kwargs):
        """
        Extract underlying flux on the detector by calling
        the `extract` method.
        All parameters are passed to `build_sys` method.
        TIPS: To be quicker, only specify the psf (`p_list`) in kwargs.
              There will be only one matrix multiplication:
              (P/sig).(w.T.lambda.c_n).
        Parameters
        ----------
        tikhonov : bool, optional
            Wheter to use tikhonov extraction
            (see regularisation.tikho_solve function).
            Default is False.
        tikho_kwargs : dictionnary or None, optional
            Arguments passed to `tikho_solve`.
        data : (N, M) array_like, optional
            A 2-D array of real values representing the detector image.
            Default is the object attribute `data`.
        error : (N, M) array_like, optional
            Estimate of the error on each pixel`
            Same shape as `data`.
            Default is the object attribute `sig`.
        mask : (N, M) array_like boolean, optional
            Additionnal mask for a given exposure. Will be added
            to the object general mask.
        throughput : (N_ord [, N_k]) list or array of functions, optional
            A list or array of the throughput at each order.
            The functions depend on the wavelength
            Default is the object attribute `t_list`
        aperture : (N_ord, N, M) list or array of 2-D arrays, optional
            A list or array of the spatial profile for each order
            on the detector. It has to have the same (N, M) as `data`.
            Default is the object attribute `p_list`

        Returns
        -----
        spectrum (f_k): solution of the linear system
        """

        return self.extract(**kwargs)

    def bin_to_pixel(self, i_order=0, grid_pix=None, grid_f_k=None, convolved_spectrum=None,
                     spectrum=None, bounds_error=False, throughput=None, **kwargs):
        """
        Integrate the convolved_spectrum (f_k_c) over a pixel grid using the trapezoidal rule.
        The concoled spectrum (f_k_c) is interpolated using scipy.interpolate.interp1d and the
        kwargs and bounds_error are passed to interp1d.
        i_order: int, optional
            index of the order to be integrated, default is 0, so
            the first order specified.
        grid_pix: tuple of two 1d-arrays or 1d-array
            If a tuple of 2 arrays is given, assume it is the lower and upper
            integration ranges. If 1d-array, assume it is the center
            of the pixels. If not given, the wavelength map and the psf map
            of `i_order` will be used to compute a pixel grid.
        grid_f_k: 1d array, optional
            grid on which the convolved flux is projected.
            Default is the wavelength grid for `i_order`.
        convolved_spectrum (f_k_c): 1d array, optional
            Convolved flux to be integrated. If not given, `spectrum`
            will be used (and convolved to `i_order` resolution)
        spectrum (f_k): 1d array, optional
            non-convolved flux (result of the `extract` method).
            Not used if `convolved_spectrum` is specified.
        bounds_error and kwargs:
            passed to interp1d function to interpolate the convolved_spectrum.
        throughput: callable, optional
            Spectral throughput for a given order (ì_ord).
            Default is given by the list of throughput saved as
            the attribute `t_list`.
        """
        # Take the value from the order if not given...

        # ... for the flux grid ...
        if grid_f_k is None:
            grid_f_k = self.wave_grid_c(i_order)

        # ... for the convolved flux ...
        if convolved_spectrum is None:
            # Use the spectrum (f_k) if the convolved_spectrum (f_k_c) not given.
            if spectrum is None:
                raise ValueError("`spectrum` or `convolved_spectrum` must be specified.")
            else:
                # Convolve the spectrum (f_k).
                convolved_spectrum = self.kernels[i_order].dot(spectrum)

        # ... and for the pixel bins
        if grid_pix is None:
            pix_center, _ = self.grid_from_map(i_order)

            # Get pixels borders (plus and minus)
            pix_p, pix_m = engine_utils.get_wave_p_or_m(pix_center)

        else:  # Else, unpack grid_pix

            # Could be a scalar or a 2-elements object)
            if len(grid_pix) == 2:

                # 2-elements object, so we have the borders
                pix_m, pix_p = grid_pix

                # Need to compute pixel center
                d_pix = (pix_p - pix_m)
                pix_center = grid_pix[0] + d_pix
            else:

                # 1-element object, so we have the pix centers
                pix_center = grid_pix

                # Need to compute the borders
                pix_p, pix_m = engine_utils.get_wave_p_or_m(pix_center)

        # Set the throughput to object attribute
        # if not given
        if throughput is None:

            # Need to interpolate
            x, y = self.wave_grid, self.throughput[i_order]
            throughput = interp1d(x, y)

        # Apply throughput on flux
        convolved_spectrum = convolved_spectrum * throughput(grid_f_k)

        # Interpolate
        kwargs['bounds_error'] = bounds_error
        fct_f_k = interp1d(grid_f_k, convolved_spectrum, **kwargs)

        # Intergrate over each bins
        bin_val = []
        for x1, x2 in zip(pix_m, pix_p):

            # Grid points that fall inside the pixel range
            i_grid = (x1 < grid_f_k) & (grid_f_k < x2)
            x_grid = grid_f_k[i_grid]

            # Add boundaries values to the integration grid
            x_grid = np.concatenate([[x1], x_grid, [x2]])

            # Integrate
            integrand = fct_f_k(x_grid) * x_grid
            bin_val.append(np.trapz(integrand, x_grid))

        # Convert to array and return with the pixel centers.
        return pix_center, np.array(bin_val)

    @staticmethod
    def _check_plot_inputs(fig, ax):
        """Method to manage inputs for plots methods."""

        # Use ax or fig if given. Else, init the figure
        if (fig is None) and (ax is None):
            fig, ax = plt.subplots(1, 1, sharex=True)
        elif ax is None:
            ax = fig.subplots(1, 1, sharex=True)

        return fig, ax

    def plot_tikho_factors(self):
        """Plot results of tikhonov tests.

        Returns
        ------
        figure and axes (for plot)
        """

        # Use tikhonov extraction from object
        tikho = self.tikho

        # Init figure
        fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

        # logl plot
        tikho.error_plot(ax=ax[0], test_key='-logl')

        # Error plot
        tikho.error_plot(ax=ax[1])

        # Labels
        ax[0].set_ylabel(r'$\log{L}$ on detector')

        # Other details
        fig.tight_layout()

        return fig, ax

    def plot_sln(self, spectrum, fig=None, ax=None, i_order=0,
                 ylabel='Flux', xlabel=r'Wavelength [$\mu$m]', **kwargs):
        """Plot extracted spectrum

        Parameters
        ----------
        spectrum (f_k): array-like
            Flux projected on the wavelength grid
        fig: matplotlib figure, optional
            Figure to use for plot
            If not given and ax is None, new figure is initiated
        ax: matplotlib axis, optional
            axis to use for plot. If not given, a new axis is initiated.
        i_order: int, optional
            index of the order to plot.
            Default is 0 (so the first order given).
        ylabel: str, optional
            Label of y axis
        xlabel: str, optional
            Label of x axis
        kwargs:
            other kwargs to be passed to plt.plot

        Returns
        ------
        fig, ax
        """

        # Manage method's inputs
        fig, ax = self._check_plot_inputs(fig, ax)

        # Set values to plot
        x = self.wave_grid_c(i_order)
        y = self.kernels[i_order].dot(spectrum)

        # Plot
        ax.plot(x, y, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return fig, ax

    def plot_err(self, spectrum, f_th_ord, fig=None, ax=None,
                 i_order=0, error='relative', ylabel='Error',
                 xlabel=r'Wavelength [$\mu$m]', **kwargs):
        """Plot error on extracted spectrum

        Parameters
        ----------
        spectrum (f_k): array-like
            Flux projected on the wavelength grid
        f_th_ord: array-like
            Injected flux projected on the wavelength grid
            and convolved at `i_order` resolution
        fig: matplotlib figure, optional
            Figure to use for plot
            If not given and ax is None, new figure is initiated
        ax: matplotlib axis, optional
            axis to use for plot. If not given, a new axis is initiated.
        i_order: int, optional
            index of the order to plot. Default is 0 (so the first order given)
        error: str, optional
            Which type of error to plot.
            Possibilities: 'relative', 'absolute', 'to_noise'
            Default is 'relative'. To noise is the error relative
            to the expected Poisson noise error
        ylabel: str, optional
            Label of y axis
        xlabel: str, optional
            Label of x axis
        kwargs:
            other kwargs to be passed to plt.plot

        Returns
        ------
        fig, ax
        """

        # Manage method's inputs
        fig, ax = self._check_plot_inputs(fig, ax)

        # Set values to plot
        x = self.wave_grid_c(i_order)
        convolved_spectrum = self.kernels[i_order].dot(spectrum)

        if error == 'relative':
            y = (convolved_spectrum - f_th_ord) / f_th_ord
        elif error == 'absolute':
            y = convolved_spectrum - f_th_ord
        elif error == 'to_noise':
            y = (convolved_spectrum - f_th_ord) / np.sqrt(f_th_ord)
        else:
            raise ValueError('`error` argument is not valid.')

        # Add info to ylabel
        ylabel += ' ({})'.format(error)

        # Plot
        ax.plot(x, y, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return fig, ax


class ExtractionEngine(_BaseOverlap):  # TODO Merge with _BaseOverlap?
    """
    Version of overlaping extraction with oversampled trapezoidal integration
    overlaping extraction solve the equation of the form:
    (B_T * B) * f = (data/sig)_T * B
    where B is a matrix and f is an array.
    The matrix multiplication B * f is the 2d model of the detector.
    We want to solve for the array f.
    The elements of f are labelled by 'k'.
    The pixels are labeled by 'i'.
    Every pixel 'i' is covered by a set of 'k' for each order
    of diffraction.
    """

    def __init__(self, wave_map, aperture, *args, **kwargs):
        """
        Parameters
        ----------
        aperture : (N_ord, N, M) list or array of 2-D arrays
            A list or array of the spatial profile for each order
            on the detector. It has to have the same (N, M) as `data`.
        wave_map : (N_ord, N, M) list or array of 2-D arrays
            A list or array of the central wavelength position for each
            order on the detector.
            It has to have the same (N, M) as `data`.
        throughput : (N_ord [, N_k]) list of array or callable
            A list of functions or array of the throughput at each order.
            If callable, the functions depend on the wavelength.
            If array, projected on `wave_grid`.
        kernels : array, callable or sparse matrix
            Convolution kernel to be applied on spectrum (f_k) for each orders.
            Can be array of the shape (N_ker, N_k_c).
            Can be a callable with the form f(x, x0) where x0 is
            the position of the center of the kernel. In this case, it must
            return a 1D array (len(x)), so a kernel value
            for each pairs of (x, x0). If array or callable,
            it will be passed to `convolution.get_c_matrix` function
            and the `c_kwargs` can be passed to this function.
            If sparse, the shape has to be (N_k_c, N_k) and it will
            be used directly. N_ker is the length of the effective kernel
            and N_k_c is the length of the spectrum (f_k) convolved.
        data : (N, M) array_like, optional
            A 2-D array of real values representing the detector image.
        error : (N, M) array_like, optional
            Estimate of the error on each pixel. Default is one everywhere.
        mask : (N, M) array_like boolean, optional
            Boolean Mask of the bad pixels on the detector.
        orders: list, optional:
            List of orders considered. Default is orders = [1, 2]
        wave_grid : (N_k) array_like, optional
            The grid on which f(lambda) will be projected.
            Default still has to be improved.
        wave_bounds : list or array-like (N_ord, 2), optional
            Boundary wavelengths covered by each orders.
            Default is the wavelength covered by `wave_map`.
        tresh : float, optional:
            The pixels where the estimated spatial profile is less than
            this value will be masked. Default is 1e-5.
        c_kwargs : list of N_ord dictionnaries or dictionnary, optional
            Inputs keywords arguments to pass to
            `convolution.get_c_matrix` function for each orders.
            If dictionnary, the same c_kwargs will be used for each orders.
        verbose : bool, optional
            Print steps. Default is False.
        """

        # Get wavelength at the boundary of each pixel
        # TODO Could also be an input??
        wave_p, wave_m = [], []
        for wave in wave_map:  # For each order
            lp, lm = engine_utils.get_wave_p_or_m(wave)  # Lambda plus or minus
            wave_p.append(lp), wave_m.append(lm)

        self.wave_p, self.wave_m = wave_p, wave_m  # Save values

        # Init upper class
        super().__init__(wave_map, aperture, *args, **kwargs)

    def _get_lo_hi(self, grid, i_order):
        """
        Find the lowest (lo) and highest (hi) index
        of wave_grid for each pixels and orders.

        Returns:
        -------
        1d array of the lowest and 1d array of the highest index.
        the length is the number of non-masked pixels
        """

        self.verbose_print('Compute low high')

        # Get needed attributes
        mask = self.mask

        # ... order dependent attributes
        attrs = ['wave_p', 'wave_m', 'mask_ord']
        wave_p, wave_m, mask_ord = self.get_attributes(*attrs, i_order=i_order)

        # Compute only for valid pixels
        wave_p = wave_p[~mask]
        wave_m = wave_m[~mask]

        # Find lower (lo) index in the pixel
        lo = np.searchsorted(grid, wave_m, side='right')

        # Find higher (hi) index in the pixel
        hi = np.searchsorted(grid, wave_p) - 1

        # Set invalid pixels for this order to lo=-1 and hi=-2
        ma = mask_ord[~mask]
        lo[ma], hi[ma] = -1, -2

        self.verbose_print('Done')

        return lo, hi

    def get_mask_wave(self, i_order):
        """ Mask according to wavelength grid """

        attrs = ['wave_p', 'wave_m', 'i_bounds']
        wave_p, wave_m, i_bnds = self.get_attributes(*attrs, i_order=i_order)
        wave_min = self.wave_grid[i_bnds[0]]
        wave_max = self.wave_grid[i_bnds[1]-1]

        mask = (wave_m < wave_min) | (wave_p > wave_max)

        return mask

    def get_w(self, i_order):
        """
        Compute integration weights for each grid points and each pixels.
        Depends on the order `n`.

        Returns
        ------
        w_n: 2d array
            weights at this specific order `n`. The shape is given by:
            (number of pixels, max number of wavelenghts covered by a pixel)
        k_n: 2d array
            index of the wavelength grid corresponding to the weights.
            Same shape as w_n
        """

        self.verbose_print('Compute weigths and k')

        # Get needed attributes
        wave_grid, mask = self.get_attributes('wave_grid', 'mask')

        # ... order dependent attributes
        attrs = ['wave_p', 'wave_m', 'mask_ord', 'i_bounds']
        wave_p, wave_m, mask_ord, i_bnds = self.get_attributes(*attrs, i_order=i_order)

        # Use the convolved grid (depends on the order)
        wave_grid = wave_grid[i_bnds[0]:i_bnds[1]]

        # Compute the wavelength coverage of the grid
        d_grid = np.diff(wave_grid)

        # Get lo hi
        lo, hi = self._get_lo_hi(wave_grid, i_order)  # Get indexes

        # Compute only valid pixels
        wave_p, wave_m = wave_p[~mask], wave_m[~mask]
        ma = mask_ord[~mask]

        # Number of used pixels
        n_i = len(lo)
        i = np.arange(n_i)

        self.verbose_print('Compute k')

        # Define fisrt and last index of wave_grid
        # for each pixel
        k_first, k_last = -1*np.ones(n_i), -1*np.ones(n_i)

        # If lowest value close enough to the exact grid value,
        # NOTE: Could be approximately equal to the exact grid
        # value. It would look like that.
        # >>> lo_dgrid = lo
        # >>> lo_dgrid[lo_dgrid==len(d_grid)] = len(d_grid) - 1
        # >>> cond = (grid[lo]-wave_m)/d_grid[lo_dgrid] <= 1.0e-8
        # But let's stick with the exactly equal
        cond = (wave_grid[lo] == wave_m)

        # special case (no need for lo_i - 1)
        k_first[cond & ~ma] = lo[cond & ~ma]
        wave_m[cond & ~ma] = wave_grid[lo[cond & ~ma]]

        # else, need lo_i - 1
        k_first[~cond & ~ma] = lo[~cond & ~ma] - 1

        # Same situation for highest value. If we follow the note
        # above (~=), the code could look like
        # >>> cond = (wave_p-grid[hi])/d_grid[hi-1] <= 1.0e-8
        # But let's stick with the exactly equal
        cond = (wave_p == wave_grid[hi])

        # special case (no need for hi_i - 1)
        k_last[cond & ~ma] = hi[cond & ~ma]
        wave_p[cond & ~ma] = wave_grid[hi[cond & ~ma]]

        # else, need hi_i + 1
        k_last[~cond & ~ma] = hi[~cond & ~ma] + 1

        # Generate array of all k_i. Set to -1 if not valid
        k_n, bad = engine_utils.arange_2d(k_first, k_last + 1, dtype=int)
        k_n[bad] = -1

        # Number of valid k per pixel
        n_k = np.sum(~bad, axis=-1)

        # Compute array of all w_i. Set to np.nan if not valid
        # Initialize
        w_n = np.zeros(k_n.shape, dtype=float)
        ####################
        ####################
        # 4 different cases
        ####################
        ####################

        self.verbose_print('compute w')

        # Valid for every cases
        w_n[:, 0] = wave_grid[k_n[:, 1]] - wave_m
        w_n[i, n_k-1] = wave_p - wave_grid[k_n[i, n_k-2]]

        ##################
        # Case 1, n_k == 2
        ##################
        case = (n_k == 2) & ~ma
        if case.any():

            self.verbose_print('n_k = 2')

            # if k_i[0] != lo_i
            cond = case & (k_n[:, 0] != lo)
            w_n[cond, 1] += wave_m[cond] - wave_grid[k_n[cond, 0]]

            # if k_i[-1] != hi_i
            cond = case & (k_n[:, 1] != hi)
            w_n[cond, 0] += wave_grid[k_n[cond, 1]] - wave_p[cond]

            # Finally
            part1 = (wave_p[case] - wave_m[case])
            part2 = d_grid[k_n[case, 0]]
            w_n[case, :] *= (part1 / part2)[:, None]

        ##################
        # Case 2, n_k >= 3
        ##################
        case = (n_k >= 3) & ~ma
        if case.any():

            self.verbose_print('n_k = 3')
            n_ki = n_k[case]
            w_n[case, 1] = wave_grid[k_n[case, 1]] - wave_m[case]
            w_n[case, n_ki-2] += wave_p[case] - wave_grid[k_n[case, n_ki-2]]

            # if k_i[0] != lo_i
            cond = case & (k_n[:, 0] != lo)
            nume1 = wave_grid[k_n[cond, 1]] - wave_m[cond]
            nume2 = wave_m[cond] - wave_grid[k_n[cond, 0]]
            deno = d_grid[k_n[cond, 0]]
            w_n[cond, 0] *= (nume1 / deno)
            w_n[cond, 1] += (nume1 * nume2 / deno)

            # if k_i[-1] != hi_i
            cond = case & (k_n[i, n_k-1] != hi)
            n_ki = n_k[cond]
            nume1 = wave_p[cond] - wave_grid[k_n[cond, n_ki-2]]
            nume2 = wave_grid[k_n[cond, n_ki-1]] - wave_p[cond]
            deno = d_grid[k_n[cond, n_ki-2]]
            w_n[cond, n_ki-1] *= (nume1 / deno)
            w_n[cond, n_ki-2] += (nume1 * nume2 / deno)

        ##################
        # Case 3, n_k >= 4
        ##################
        case = (n_k >= 4) & ~ma
        if case.any():
            self.verbose_print('n_k = 4')
            n_ki = n_k[case]
            w_n[case, 1] += wave_grid[k_n[case, 2]] - wave_grid[k_n[case, 1]]
            w_n[case, n_ki-2] += (wave_grid[k_n[case, n_ki-2]]
                                  - wave_grid[k_n[case, n_ki-3]])

        ##################
        # Case 4, n_k > 4
        ##################
        case = (n_k > 4) & ~ma
        if case.any():
            self.verbose_print('n_k > 4')
            i_k = np.indices(k_n.shape)[-1]
            cond = case[:, None] & (2 <= i_k) & (i_k < n_k[:, None]-2)
            ind1, ind2 = np.where(cond)
            w_n[ind1, ind2] = (d_grid[k_n[ind1, ind2]-1]
                               + d_grid[k_n[ind1, ind2]])

        # Finally, divide w_n by 2
        w_n /= 2.

        # Make sure invalid values are masked
        w_n[k_n < 0] = np.nan

        self.verbose_print('Done')

        return w_n, k_n
