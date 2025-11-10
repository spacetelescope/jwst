import logging

import numpy as np
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import custom_model, fitting, models

from jwst.stpipe import Step

__all__ = ["TACenterStep"]

log = logging.getLogger(__name__)


# Constants for MIRI imager and LRS slit dimensions
PIXSCALE = 0.11  # arcsec/pixel for MIRI imager
JWST_DIAMETER = 6.5  # meters
SLIT_WIDTH_PIXELS = 4.6  # LRS slit width in pixels (0.51 arcsec / 0.11 arcsec/pix)
SLIT_LENGTH_PIXELS = 42.7  # LRS slit length in pixels (4.7 arcsec / 0.11 arcsec/pix)
SLIT_CENTER = (50, 50)  # Approximate center of the TA image


class TACenterStep(Step):
    """Perform target acquisition centering analysis."""

    class_alias = "ta_center"

    spec = """
    ta_file = string(default=None). # Target acquisition image file name
    """  # noqa: E501

    # reference_file_types = []

    def process(self, step_input):
        """
        Process target acquisition data.

        Parameters
        ----------
        step_input : DataModel
            The input data model.

        Returns
        -------
        result : DataModel
            The output data model with TA centering applied.
        """
        # Open the input data model
        with dm.open(step_input) as input_model:
            # Work on a copy
            result = input_model.copy()

            if self.ta_file is None:
                log.error("No target acquisition file provided. Step will be SKIPPED.")
                result.meta.cal_step.ta_center = "SKIPPED"
                return result

            if input_model.meta.target.srctype != "POINT":
                log.error(
                    "TA centering is only implemented for point sources. Step will be SKIPPED."
                )
                result.meta.cal_step.ta_center = "SKIPPED"
                return result

            with dm.open(self.ta_file) as ta_model:
                log.info(f"Performing TA centering using file: {self.ta_file}")
                ta_image = ta_model.data
                wavelength = _get_wavelength(input_model.meta.instrument.filter)
                if wavelength is None:
                    log.error(
                        "Unknown filter; cannot determine wavelength for TA centering."
                        " Step will be SKIPPED."
                    )
                    result.meta.cal_step.ta_center = "SKIPPED"
                    return result

                # Determine if this is slit or slitless data
                exp_type = input_model.meta.exposure.type
                is_slit = exp_type == "MIR_LRS-FIXEDSLIT"

            if is_slit:
                log.info("Fitting Airy disk with slit mask model for LRS FIXEDSLIT mode")
                # Fit Airy disk multiplied by slit mask
                x_center, y_center = _fit_airy_disk(ta_image, wavelength, use_slit_model=True)
            else:
                log.info("Fitting Airy disk for slitless mode")
                x_center, y_center = _fit_airy_disk(ta_image, wavelength, use_slit_model=False)

            # Set completion status
            result.x_center = x_center
            result.y_center = y_center
            result.meta.cal_step.ta_center = "COMPLETE"

        return result


def _fit_airy_disk(ta_image, wavelength, weights=None, use_slit_model=False):
    """
    Fit an Airy disk model to target acquisition image data.

    Parameters
    ----------
    ta_image : ndarray
        2D target acquisition image data.
    wavelength : float
        Wavelength in microns for calculating the Airy disk size.
    weights : ndarray, optional
        Weights for fitting. If None, all pixels are weighted equally.
        Only used when use_slit_model=False.
    use_slit_model : bool, optional
        If True, fit an Airy disk multiplied by the slit mask model.
        If False, fit a plain Airy disk (optionally with weights).

    Returns
    -------
    x_center : float
        Fitted x center position.
    y_center : float
        Fitted y center position.
    """
    # Create coordinate grids
    y, x = np.mgrid[0 : ta_image.shape[0], 0 : ta_image.shape[1]]

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Initial guess for center (centroid or peak)
    y_guess = np.sum(y * ta_image) / np.sum(ta_image)
    x_guess = np.sum(x * ta_image) / np.sum(ta_image)
    amplitude_guess = np.max(ta_image)

    if use_slit_model:
        # Create compound model: Airy disk * slit mask
        airy_init = models.AiryDisk2D(
            amplitude=amplitude_guess, x_0=x_guess, y_0=y_guess, radius=radius_pixels
        )
        slit_mask_model = SlitMask()
        compound_model = airy_init * slit_mask_model

        # Fit the compound model to the data
        fitter = fitting.LevMarLSQFitter()
        fitted_model = fitter(compound_model, x, y, ta_image)
        # Extract parameters from the left side of the compound model (Airy disk)
        x_center = fitted_model.left.x_0.value
        y_center = fitted_model.left.y_0.value
    else:
        # Create initial Airy disk model
        airy_init = models.AiryDisk2D(
            amplitude=amplitude_guess, x_0=x_guess, y_0=y_guess, radius=radius_pixels
        )

        # Fit the model to the data
        fitter = fitting.LevMarLSQFitter()
        airy_fit = fitter(airy_init, x, y, ta_image, weights=weights)
        x_center = airy_fit.x_0.value
        y_center = airy_fit.y_0.value

    log.info(f"Fitted center: x={x_center:.2f}, y={y_center:.2f}")
    return x_center, y_center


@custom_model
def SlitMask(x, y):  # noqa: N802
    """
    Evaluate the slit mask at given (x, y) coordinates.

    Computes the fractional overlap between a pixel at position (x, y) and the
    slit aperture at subpixel accuracy. The pixel extends from x to x+1 and y to y+1.
    The custom_model decorator allows this function to be used as a
    2-D astropy model in fitting.

    Parameters
    ----------
    x : float or ndarray
        X coordinate(s) in pixel coordinates (pixel indices).
    y : float or ndarray
        Y coordinate(s) in pixel coordinates (pixel indices).

    Returns
    -------
    mask_value : float or ndarray
        Fractional pixel coverage value(s) between 0.0 and 1.0.
    """
    # Define slit boundaries in pixel coordinates
    y_min = SLIT_CENTER[1] - SLIT_WIDTH_PIXELS / 2
    y_max = SLIT_CENTER[1] + SLIT_WIDTH_PIXELS / 2
    x_min = SLIT_CENTER[0] - SLIT_LENGTH_PIXELS / 2
    x_max = SLIT_CENTER[0] + SLIT_LENGTH_PIXELS / 2

    # Pixel extends from (x, y) to (x+1, y+1) in index coordinates
    # Compute overlap in X direction
    pixel_x_min = x
    pixel_x_max = x + 1
    overlap_x_min = np.maximum(pixel_x_min, x_min)
    overlap_x_max = np.minimum(pixel_x_max, x_max)
    overlap_x = np.maximum(0, overlap_x_max - overlap_x_min)

    # Compute overlap in Y direction
    pixel_y_min = y
    pixel_y_max = y + 1
    overlap_y_min = np.maximum(pixel_y_min, y_min)
    overlap_y_max = np.minimum(pixel_y_max, y_max)
    overlap_y = np.maximum(0, overlap_y_max - overlap_y_min)

    # Total fractional coverage
    mask_value = overlap_x * overlap_y

    return mask_value


def _get_wavelength(filter_name):  # numpydoc ignore=RT01
    """Map filter name to central wavelength in microns."""
    filter_wavelengths = {
        "F560W": 5.6,
        "F770W": 7.7,
        "F1000W": 10.0,
        "F1130W": 11.3,
        "F1280W": 12.8,
        "F1500W": 15.0,
        "F1800W": 18.0,
        "F2100W": 21.0,
        "F2550W": 25.5,
    }
    return filter_wavelengths.get(filter_name, None)
