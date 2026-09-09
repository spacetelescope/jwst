"""Pixel dispersion table, which caches the dispersed-pixel offsets for efficient interpolation."""

import numpy as np
from astropy.modeling.mappings import Mapping
from scipy.interpolate import RegularGridInterpolator

__all__ = ["TracePDT", "build_trace_pdt", "get_grism_detector_transform"]


class TracePDT:
    """
    Interpolated replacement for the "detector" to "grism_detector" WCS transform.

    Holds a pixel dispersion table, which is basically a cache of the
    ``(x0, y0, wavelength) -> (dx, dy)`` mapping on a coarse regular grid
    built using the exact transform.
    On init a linear interpolant is built from the sparse regular grid.
    On call the interpolant provides a drop-in replacement for the exact detector-to-grism
    WCS transform to disperse the direct-image pixels.
    """

    def __init__(self, x0_grid, y0_grid, lam_grid, dx_grid, dy_grid, exact_wavelength_grid=False):
        """
        Initialize the lookup table from precomputed grid offsets.

        Parameters
        ----------
        x0_grid, y0_grid : np.ndarray
            1-D, strictly increasing arrays of detector x, y grid positions.
        lam_grid : np.ndarray
            1-D, strictly increasing array of wavelength grid positions.
        dx_grid, dy_grid : np.ndarray
            Arrays of shape ``(len(x0_grid), len(y0_grid), len(lam_grid))`` giving the
            dispersed-pixel offset from ``(x0, y0)`` at each grid point.
        exact_wavelength_grid : bool, optional
            If True, ``lam_grid`` is known to exactly match the wavelength array that
            will be passed to `evaluate_grid`, allowing it to skip wavelength-axis
            interpolation entirely. Defaults to False.
        """
        self._x0_grid = np.asarray(x0_grid)
        self._y0_grid = np.asarray(y0_grid)
        self._lam_grid = np.asarray(lam_grid)
        self._dx_grid = np.asarray(dx_grid)
        self._dy_grid = np.asarray(dy_grid)
        self._exact_wavelength_grid = exact_wavelength_grid

        # Stack dx, dy into a single vector-valued interpolator so the grid-cell
        # lookup is only done once per query point instead of twice
        values = np.stack([self._dx_grid, self._dy_grid], axis=-1)
        self._interp = RegularGridInterpolator(
            (self._x0_grid, self._y0_grid, self._lam_grid),
            values,
            bounds_error=False,
            fill_value=None,
        )

    def __call__(self, x0, y0, wavelength):
        """
        Interpolate the dispersed pixel position for input pixel position(s) and wavelength(s).

        This generic, elementwise path handles arbitrary (independently varying)
        input shapes. For the common case where many pixels share the same
        wavelength array (as in `~jwst.wfss_contam.disperse._disperse_onto_grism`),
        use `evaluate_grid` instead, which is substantially faster.

        Parameters
        ----------
        x0, y0 : np.ndarray
            Detector x, y position(s) of the direct image pixel(s).
        wavelength : np.ndarray
            Wavelength(s) at which to evaluate the dispersion.
            Must have the same shape as x0 and y0.

        Returns
        -------
        x, y : np.ndarray
            Interpolated x, y position(s) in the dispersed image, same shape as input.
        """
        shape = x0.shape
        pts = np.stack([x0.ravel(), y0.ravel(), np.asarray(wavelength).ravel()], axis=-1)
        offsets = self._interp(pts).reshape(*shape, 2)
        return x0 + offsets[..., 0], y0 + offsets[..., 1]

    def evaluate_grid(self, x0, y0, wavelength):
        """
        Evaluate the PDT on a grid.

        Every pixel is dispersed at the same set of wavelengths, so this exploits the
        separability of trilinear interpolation: the ``(x0, y0)`` grid-cell lookup is
        done once per pixel (independent of wavelength), and the wavelength grid-cell
        lookup is done once for the shared wavelength array.

        If this `TracePDT` was built with ``exact_wavelength_grid=True``,
        wavelength-axis interpolation is skipped entirely, reducing this to a 2-D
        (x0, y0) bilinear lookup. Otherwise, wavelength is interpolated linearly.

        Parameters
        ----------
        x0, y0 : np.ndarray
            1-D arrays giving the detector x, y position of each pixel to disperse.
        wavelength : np.ndarray
            1-D array giving the wavelengths to evaluate,
            shared by every pixel. If this `TracePDT` was built with
            ``exact_wavelength_grid=True`` (see `build_trace_pdt`'s ``lam_grid``
            argument), this must be exactly the same array used to build the grid.

        Returns
        -------
        x, y : np.ndarray
            Arrays of shape ``(n_wave, n_pixels)`` giving the interpolated dispersed
            pixel positions.
        """
        ix, wx = self._cell_weights(self._x0_grid, x0)
        iy, wy = self._cell_weights(self._y0_grid, y0)

        if self._exact_wavelength_grid:
            if wavelength.shape != self._lam_grid.shape or not np.array_equal(
                wavelength, self._lam_grid
            ):
                raise ValueError(
                    "This TracePDT was built with exact_wavelength_grid=True, so "
                    "evaluate_grid() must be called with the same wavelength array "
                    "used to build it."
                )
            dx = self._bilinear_only(self._dx_grid, ix, wx, iy, wy).T
            dy = self._bilinear_only(self._dy_grid, ix, wx, iy, wy).T
            return x0[np.newaxis, :] + dx, y0[np.newaxis, :] + dy

        ilam, wlam = self._cell_weights(self._lam_grid, wavelength)
        dx = self._separable_interp(self._dx_grid, ix, wx, iy, wy, ilam, wlam)
        dy = self._separable_interp(self._dy_grid, ix, wx, iy, wy, ilam, wlam)
        return x0[np.newaxis, :] + dx, y0[np.newaxis, :] + dy

    @staticmethod
    def _cell_weights(grid, values):
        """
        Return the lower grid index and fractional weight for linear interpolation.

        Parameters
        ----------
        grid : np.ndarray
            1-D, strictly increasing array of grid points.
        values : np.ndarray
            Values to locate within the grid.

        Returns
        -------
        idx : np.ndarray
            Indices of the lower grid points for each value.
        w : np.ndarray
            Fractional weights for linear interpolation between the lower and upper grid points.
        """
        idx = np.clip(np.searchsorted(grid, values) - 1, 0, len(grid) - 2)
        w = (values - grid[idx]) / (grid[idx + 1] - grid[idx])
        return idx, w

    @staticmethod
    def _bilinear_only(grid_vals, ix, wx, iy, wy):
        """
        Bilinear-in-(x0, y0) interpolation of ``grid_vals``, preserving the wavelength axis.

        Parameters
        ----------
        grid_vals : np.ndarray
            3-D array of shape ``(len(x0_grid), len(y0_grid), len(lam_grid))``
            containing the values to interpolate.
        ix, wx : np.ndarray
            Lower indices and fractional weights for the x0 grid.
        iy, wy : np.ndarray
            Lower indices and fractional weights for the y0 grid.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_pixels, len(lam_grid))``.
        """
        c00 = grid_vals[ix, iy, :]
        c10 = grid_vals[ix + 1, iy, :]
        c01 = grid_vals[ix, iy + 1, :]
        c11 = grid_vals[ix + 1, iy + 1, :]
        return (
            c00 * ((1 - wx) * (1 - wy))[:, np.newaxis]
            + c10 * (wx * (1 - wy))[:, np.newaxis]
            + c01 * ((1 - wx) * wy)[:, np.newaxis]
            + c11 * (wx * wy)[:, np.newaxis]
        )

    @staticmethod
    def _separable_interp(grid_vals, ix, wx, iy, wy, ilam, wlam):
        """
        Bilinear-in-(x0, y0), then linear-in-wavelength interpolation of ``grid_vals``.

        Parameters
        ----------
        grid_vals : np.ndarray
            3-D array of shape ``(len(x0_grid), len(y0_grid), len(lam_grid))``
            containing the values to interpolate.
        ix, wx : np.ndarray
            Lower indices and fractional weights for the x0 grid.
        iy, wy : np.ndarray
            Lower indices and fractional weights for the y0 grid.
        ilam, wlam : np.ndarray
            Lower indices and fractional weights for the wavelength grid.

        Returns
        -------
        np.ndarray
            Interpolated values at the specified pixel and wavelength coordinates.
        """
        curve = TracePDT._bilinear_only(grid_vals, ix, wx, iy, wy)  # (n_pixels, n_wave_grid)

        # Linear interpolation along wavelength, shared across all pixels.
        lo = curve[:, ilam]
        hi = curve[:, ilam + 1]
        return (lo * (1 - wlam)[np.newaxis, :] + hi * wlam[np.newaxis, :]).T


def native_wavelength_grid(imgxy_to_grismxy, order, wmin, wmax, x_ref, y_ref, oversample_factor=1):
    """
    Make a grid with the approximate spacing that disperses each wavelength into a unique pixel.

    Using the direct-image-to-grism transform, disperse the min, max wavelength.
    Compute the distance in the dispersed plane between those wavelength endpoints.
    Use that distance to figure out how many wavelength samples are needed to put roughly one
    unique wavelength into each pixel.
    Finally oversample that by the oversample factor.

    Parameters
    ----------
    imgxy_to_grismxy : `~astropy.modeling.Model`
        The "detector" to "grism_detector" transform, reduced to 2 (x, y) outputs.
    order : int
        Spectral order number.
    wmin, wmax : float
        Minimum, maximum wavelength for the dispersed spectra of this spectral order.
    x_ref, y_ref : float
        Representative detector x, y position at which to evaluate the native spacing.
    oversample_factor : float, optional
        Factor by which to oversample the native wavelength spacing, default 1.

    Returns
    -------
    np.ndarray
        Wavelength grid spanning ``[wmin, wmax]``.
    """
    x_ref = np.atleast_1d(x_ref)
    y_ref = np.atleast_1d(y_ref)
    xwmin, ywmin = imgxy_to_grismxy(x_ref, y_ref, np.atleast_1d(wmin), order)
    xwmax, ywmax = imgxy_to_grismxy(x_ref, y_ref, np.atleast_1d(wmax), order)
    dxw = xwmax - xwmin
    dyw = ywmax - ywmin
    dlam = np.abs((wmax - wmin) / (dyw - dxw))[0] / oversample_factor
    npts = max(int(np.ceil((wmax - wmin) / dlam)), 3)
    return np.linspace(wmin, wmax, npts)


def get_grism_detector_transform(grism_wcs):
    """
    Get the "detector" to "grism_detector" transform, reduced to 2 (x, y) outputs.

    Parameters
    ----------
    grism_wcs : `~gwcs.wcs.WCS`
        The grism WCS object.

    Returns
    -------
    `~astropy.modeling.Model`
        The "detector" to "grism_detector" transform, with only the (x, y) outputs.
    """
    imgxy_to_grismxy = grism_wcs.get_transform("detector", "grism_detector")
    # We only need the x,y outputs, same as in disperse(). Making the number of
    # outputs dynamic handles legacy WCS objects that did not pass x0, y0, and
    # order through the transform unmodified like the current version does.
    n_outputs = len(imgxy_to_grismxy.outputs)
    return imgxy_to_grismxy | Mapping((0, 1), n_inputs=n_outputs)


def build_trace_pdt(
    grism_wcs, order, wmin, wmax, direct_shape, spacing, wave_oversample_factor=1, lam_grid=None
):
    """
    Build a TracePDT for one spectral order by sampling the exact transform on a coarse grid.

    Parameters
    ----------
    grism_wcs : `~gwcs.wcs.WCS`
        The grism WCS object, from which the "detector" to "grism_detector"
        transform is retrieved and evaluated on the coarse grid.
    order : int
        Spectral order number.
    wmin, wmax : float
        Minimum, maximum wavelength for the dispersed spectra of this spectral order.
        Unused if ``lam_grid`` is provided.
    direct_shape : tuple of int
        ``(nx, ny)`` dimensions of the direct image.
    spacing : int
        Spacing of grid points to sample along each of the x and y axes.
    wave_oversample_factor : float, optional
        Factor by which to oversample the native wavelength spacing
        when building the wavelength axis of the grid. Unused if ``lam_grid`` is provided.
    lam_grid : np.ndarray, optional
        Precomputed wavelength grid to use directly instead of computing one from
        ``wmin``, ``wmax``, and ``wave_oversample_factor``. This must exactly match the
        wavelength array that will later be passed to `TracePDT.evaluate_grid`, as
        setting this means interpolation along the wavelength axis will be skipped entirely
        at query time. If None (the default), the grid is computed
        from ``wmin``, ``wmax``, and ``wave_oversample_factor``.

    Returns
    -------
    TracePDT
        Lookup table object that can be called as ``trace_pdt(x0, y0, wavelength)``
        to approximate the exact "detector" to "grism_detector" transform.
    """
    imgxy_to_grismxy = get_grism_detector_transform(grism_wcs)

    nx, ny = direct_shape
    n_grid_x = int(np.ceil(nx / spacing))
    n_grid_y = int(np.ceil(ny / spacing))
    x0_grid = np.linspace(0, nx - 1, n_grid_x)
    y0_grid = np.linspace(0, ny - 1, n_grid_y)

    exact_wavelength_grid = lam_grid is not None
    if exact_wavelength_grid:
        lam_grid = np.asarray(lam_grid)
    else:
        lam_grid = native_wavelength_grid(
            imgxy_to_grismxy,
            order,
            wmin,
            wmax,
            (nx - 1) / 2.0,
            (ny - 1) / 2.0,
            oversample_factor=wave_oversample_factor,
        )

    xx0, yy0 = np.meshgrid(x0_grid, y0_grid, indexing="ij")
    x0_flat = xx0.ravel()
    y0_flat = yy0.ravel()
    n_pix = x0_flat.size
    n_wave = lam_grid.size

    # Match the (n_wave, n_pixels) calling convention. Some backward grism
    # dispersion transforms rely on this broadcasting pattern internally
    # to efficiently invert the wavelength solution for many pixels at once.
    x0_rep = np.repeat(x0_flat[np.newaxis, :], n_wave, axis=0)
    y0_rep = np.repeat(y0_flat[np.newaxis, :], n_wave, axis=0)
    lam_rep = np.repeat(lam_grid[:, np.newaxis], n_pix, axis=1)

    xd, yd = imgxy_to_grismxy(x0_rep, y0_rep, lam_rep, order)

    dx_grid = (xd - x0_rep).reshape(n_wave, n_grid_x, n_grid_y)
    dy_grid = (yd - y0_rep).reshape(n_wave, n_grid_x, n_grid_y)
    dx_grid = np.moveaxis(dx_grid, 0, -1)
    dy_grid = np.moveaxis(dy_grid, 0, -1)

    return TracePDT(
        x0_grid, y0_grid, lam_grid, dx_grid, dy_grid, exact_wavelength_grid=exact_wavelength_grid
    )
