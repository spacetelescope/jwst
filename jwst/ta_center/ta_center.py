import logging
import warnings

import numpy as np
from astropy.modeling import Fittable2DModel, fitting, models
from astropy.stats import sigma_clip
from astropy.utils.exceptions import AstropyUserWarning
from scipy.ndimage import median_filter

import jwst.datamodels as dm
from jwst.pathloss.pathloss import calculate_pathloss_vector

log = logging.getLogger(__name__)

# Constants for MIRI imager and LRS slit dimensions
PIXSCALE = 0.11  # arcsec/pixel for MIRI imager
JWST_DIAMETER = 6.5  # meters


__all__ = ["center_from_ta_image", "PathlossMask"]


class NoFinitePixelsError(Exception):
    """Custom exception raised when no finite pixels are found in the TA image."""

    pass


class BadFitError(Exception):
    """Custom exception raised when the model fit does not meet quality criteria."""

    pass


def center_from_ta_image(
    ta_image, wavelength, ref_center, subarray_origin=(1, 1), pathloss_file=None
):
    """
    Determine the center of a point source from a TA image.

    Parameters
    ----------
    ta_image : ndarray
        2D target acquisition image data.
    wavelength : float
        Wavelength in microns for calculating the Airy disk size.
    ref_center : tuple of float
        (x_ref, y_ref) reference center position in full-frame detector coordinates.
    subarray_origin : tuple of int, optional
        (xstart, ystart) 1-indexed origin of the subarray in full-frame coordinates.
        Default is (1, 1) for full frame.
    pathloss_file : str, optional
        Path to the pathloss reference file for the slit mask model.
        If provided, fit an Airy disk multiplied by the slit mask model.
        Otherwise, fit just the Airy disk.

    Returns
    -------
    x_center : float
        Fitted x center position in full-frame detector coordinates.
    y_center : float
        Fitted y center position in full-frame detector coordinates.
    """
    # Transform reference center from full-frame to subarray coordinates
    # FITS convention: xstart, ystart are 1-indexed
    # Python/array convention: 0-indexed
    # ref_center is in 0-indexed detector coordinates
    ref_center_subarray = (
        ref_center[0] - (subarray_origin[0] - 1),
        ref_center[1] - (subarray_origin[1] - 1),
    )

    log.info(
        f"Reference center (0-indexed): full-frame=({ref_center[0]:.2f}, {ref_center[1]:.2f}), "
        f"subarray=({ref_center_subarray[0]:.2f}, {ref_center_subarray[1]:.2f})"
    )

    # Create cutout around reference center (in subarray coordinates)
    cutout, cutout_origin = _cutout_center(ta_image, ref_center_subarray)

    if pathloss_file is not None:
        # Load pathloss reference file for slit mask model
        with dm.MirLrsPathlossModel(pathloss_file) as pathloss:
            pathloss_table = pathloss.pathloss_table
            pathloss_wcs = pathloss.meta.wcsinfo

        slit_center_cutout = (
            ref_center_subarray[0] - cutout_origin[0],
            ref_center_subarray[1] - cutout_origin[1],
        )
        pathloss_mask = PathlossMask(
            pathloss_table,
            pathloss_wcs,
            xcenter=slit_center_cutout[0],
            ycenter=slit_center_cutout[1],
            wavelength=wavelength,
        )
    else:
        pathloss_mask = None

    # Fit on the cutout
    x_center_cutout, y_center_cutout = _fit_airy_disk(
        cutout,
        wavelength,
        pathloss_model=pathloss_mask,
    )

    # Transform back to subarray coordinates
    x_center_subarray = x_center_cutout + cutout_origin[0]
    y_center_subarray = y_center_cutout + cutout_origin[1]

    # Transform from subarray to full-frame detector coordinates
    x_center = x_center_subarray + (subarray_origin[0] - 1)
    y_center = y_center_subarray + (subarray_origin[1] - 1)

    log.info(
        f"Fitted center (0-indexed): subarray=({x_center_subarray:.2f}, {y_center_subarray:.2f}), "
        f"full-frame=({x_center:.2f}, {y_center:.2f})"
    )

    return x_center, y_center


def _fit_airy_disk(ta_image, wavelength, pathloss_model=None):
    """
    Fit an Airy disk model to target acquisition image data.

    Parameters
    ----------
    ta_image : ndarray
        2D target acquisition image data.
    wavelength : float
        Wavelength in microns for calculating the Airy disk size.
    pathloss_model : Fittable2DModel, optional
        Pathloss model used to compute fitting weights (not included in the
        fitted model itself).

    Returns
    -------
    x_center : float
        Fitted x center position.
    y_center : float
        Fitted y center position.
    """
    # Create coordinate grids
    y, x = np.mgrid[0 : ta_image.shape[0], 0 : ta_image.shape[1]]

    # Filter non-finite values from data
    mask = np.isfinite(ta_image)
    if np.sum(mask) < 10:
        raise NoFinitePixelsError("Most or all pixels contain non-finite values")
    log.debug(f"Excluding {np.sum(~mask)} non-finite values from fit")
    ta_masked = ta_image[mask]

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Add a uniform background model
    clipped_data = sigma_clip(ta_masked, sigma=1, maxiters=5, masked=False)
    clipped_med = np.median(clipped_data)
    background_init = models.Const2D(amplitude=clipped_med)
    # bound background model to be within three sigma of guess
    background_init.amplitude.bounds = (
        clipped_med - 3 * np.std(ta_masked),
        clipped_med + 3 * np.std(ta_masked),
    )

    # Create initial Airy disk model
    # Initial guess for center and amplitude
    # Do a spatial median filter to remove hot pixels if present, then find location of max
    filtered_image = median_filter(ta_image, size=3)
    y_guess, x_guess = np.unravel_index(np.nanargmax(filtered_image), ta_image.shape)
    # Amplitude guess above background
    amplitude_guess = ta_image[int(y_guess), int(x_guess)] - clipped_med

    airy_init = models.AiryDisk2D(
        amplitude=amplitude_guess, x_0=x_guess, y_0=y_guess, radius=radius_pixels
    )
    # Allow radius to vary a bit (0.99 to 1.5 times diffraction limit)
    airy_init.radius.bounds = (0.99 * radius_pixels, 1.5 * radius_pixels)
    # Force amplitude positive
    airy_init.amplitude.bounds = (0.0, None)
    # Keep x, y center within 5 pixels of centroid guess
    airy_init.x_0.bounds = (x_guess - 5.0, x_guess + 5.0)
    airy_init.y_0.bounds = (y_guess - 5.0, y_guess + 5.0)

    # Build model WITHOUT multiplying by the pathloss mask. Use pathloss
    # only to compute weights for the fitter.
    compound_model = airy_init + background_init

    weights = None
    if pathloss_model is not None:
        # Evaluate the pathloss mask at the valid pixel locations to use
        # as weights. Keep as 1D array corresponding to data[mask].
        pw = pathloss_model(x[mask], y[mask])
        weights = np.atleast_1d(pw).astype(float)

    # Fit the compound model to filtered data using optional weights
    fitted_model = _fit_catch_errors(compound_model, x, y, ta_image, mask, weights=weights)

    # Extract parameters from the Airy disk model (first submodel)
    x_center = fitted_model[0].x_0.value
    y_center = fitted_model[0].y_0.value

    return x_center, y_center


def _fit_catch_errors(model_init, x, y, data, mask, weights=None):
    """
    Fit a model to data, catching common fitting errors and raising custom types.

    Parameters
    ----------
    model_init : Fittable2DModel
        Initial model to fit.
    x : ndarray
        X coordinate grid.
    y : ndarray
        Y coordinate grid.
    data : ndarray
        2D data array to fit.
    mask : ndarray
        Boolean mask array indicating valid data points.

    Returns
    -------
    fitted_model : Fittable2DModel
        Fitted model.
    """
    fitter = fitting.TRFLSQFitter()
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("error", category=AstropyUserWarning)
            # astropy fitters accept a `weights` keyword matching the shape of
            # the dependent data vector. If weights is provided it should be
            # a 1-D array corresponding to data[mask]. Otherwise pass None.
            fitted_model = fitter(
                model_init,
                x[mask],
                y[mask],
                data[mask],
                weights=weights,
            )
    except TypeError as e:
        if str(e).startswith("Improper input: func input vector length"):
            raise NoFinitePixelsError("Not enough finite pixels for fitting") from e
    except AstropyUserWarning as e:  # Raised if fit does not converge
        fit_status = fitter.fit_info["ierr"]
        msg = fitter.fit_info["message"]
        raise BadFitError(f"Model fitting failed with status code {fit_status}: {msg}") from e

    # Build the model image from the fitted model and compare with the data
    fitted_image = fitted_model(x, y)
    residuals = data - fitted_image
    residual_std = np.std(residuals[mask])
    log.info(f"Fit residuals standard deviation: {residual_std:.4f}")

    # # Get the fitted center of the Airy disk
    # xcen, ycen = (fitted_model[0].x_0.value, fitted_model[0].y_0.value)
    # import matplotlib.pyplot as plt
    # fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    # axes[0].imshow(data, origin="lower", cmap="viridis")
    # axes[0].set_title("TA Image Data")
    # axes[1].imshow(fitted_image, origin="lower", cmap="viridis")
    # axes[1].set_title("Fitted Model")
    # axes[2].imshow(residuals, origin="lower", cmap="viridis")
    # axes[2].set_title("Residuals")
    # weights_2d = np.zeros_like(data, dtype=float)
    # weights_2d[mask] = weights if weights is not None else 1.0
    # axes[3].imshow(weights_2d, origin="lower", cmap="viridis")
    # axes[3].set_title("Fit weights (pathloss)")
    # for ax in axes:
    #     ax.plot(xcen, ycen, "rx")
    # plt.tight_layout()
    # plt.show()

    # Raise an exception if the fit is bad
    # For now this is extremely generous - need to check what INS wants here
    data_std = np.std(data[mask])
    if residual_std > data_std:
        raise BadFitError("Fitted model residuals are larger than threshold")

    return fitted_model


def _cutout_center(image, center, size=40):
    """
    Cut out a small square region from an image centered on a reference position.

    Parameters
    ----------
    image : ndarray
        2D image array.
    center : tuple of float
        (x_center, y_center) position for the center of the cutout.
    size : int, optional
        Size of the square cutout in pixels. Default is 20.

    Returns
    -------
    cutout : ndarray
        Square cutout of the image.
    cutout_origin : tuple of int
        (x_origin, y_origin) position of the lower-left corner of the cutout
        in the original image coordinates.
    """
    x_center, y_center = center
    ny, nx = image.shape

    # Convert center to integer pixel for cutout boundaries
    x_center_int = int(np.round(x_center))
    y_center_int = int(np.round(y_center))

    # Calculate half-size
    half_size = size // 2

    # Calculate cutout boundaries
    x_min = max(0, x_center_int - half_size)
    x_max = min(nx, x_center_int + half_size)
    y_min = max(0, y_center_int - half_size)
    y_max = min(ny, y_center_int + half_size)

    # Extract cutout
    cutout = image[y_min:y_max, x_min:x_max]

    # Store the origin for coordinate transformation
    cutout_origin = (x_min, y_min)

    return cutout, cutout_origin


class PathlossMask(Fittable2DModel):
    """
    Model for MIRI LRS pathloss correction.

    This model has no fittable parameters - it's a static mask that depends
    on the fixed wavelength and source position.

    Parameters
    ----------
    pathloss_table : table
        Pathloss reference table.
    pathloss_wcs : dict
        WCS information for the pathloss reference data.
    xcenter : float
        X coordinate of the source center in pixel coordinates.
    ycenter : float
        Y coordinate of the source center in pixel coordinates.
    wavelength : float
        Wavelength in microns.
    """

    n_inputs = 2
    n_outputs = 1

    def __init__(self, pathloss_table, pathloss_wcs, xcenter, ycenter, wavelength, **kwargs):
        self.pathloss_table = pathloss_table
        self.pathloss_data = pathloss_table["pathloss"]
        self.pathloss_wave = pathloss_table["wavelength"]
        self.pathloss_wcs = pathloss_wcs

        self.xcenter = xcenter
        self.ycenter = ycenter
        self.wavelength = wavelength
        self.wave_idx = self._wave_to_idx(wavelength)

        super().__init__(**kwargs)

    def evaluate(self, x, y):
        """
        Compute the pathloss correction factor at given coordinates.

        Parameters
        ----------
        x : float or np.ndarray
            X coordinate in pixel coordinates.
        y : float or np.ndarray
            Y coordinate in pixel coordinates.

        Returns
        -------
        pathloss_value : float or np.ndarray
            Pathloss correction factor.
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        # Calculate offsets in arcseconds to pass into calculate_pathloss_vector
        dx_pixels = (x - self.xcenter).flatten()
        dy_pixels = (y - self.ycenter).flatten()
        loss = np.empty_like(dx_pixels)
        # In the future we should vectorize calculate_pathloss_vector to avoid loop.
        # This isn't preventatively slow, so not a high priority
        for i in range(len(dx_pixels)):
            _, pathloss_vector, _is_inside_slit = calculate_pathloss_vector(
                self.pathloss_data, self.pathloss_wcs, dx_pixels[i], dy_pixels[i], calc_wave=False
            )
            loss[i] = pathloss_vector[self.wave_idx]
        loss = loss.reshape(x.shape)

        # Normalize to max of 1
        loss /= np.max(loss)
        return loss

    def _wave_to_idx(self, wavelength):
        """
        Convert wavelength in microns to index in pathloss table.

        Parameters
        ----------
        wavelength : float
            Wavelength in microns.

        Returns
        -------
        wave_index : int
            Index corresponding to the wavelength in the pathloss table.
        """
        wave_diffs = np.abs(self.pathloss_wave - wavelength)
        wave_index = np.argmin(wave_diffs)
        return wave_index
