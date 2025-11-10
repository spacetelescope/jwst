import logging

import numpy as np
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import Fittable2DModel, fitting, models

from jwst.stpipe import Step

__all__ = ["TACenterStep"]

log = logging.getLogger(__name__)


# Constants for MIRI imager and LRS slit dimensions
PIXSCALE = 0.11  # arcsec/pixel for MIRI imager
JWST_DIAMETER = 6.5  # meters
SLIT_WIDTH_PIXELS = 4.6  # LRS slit width in pixels (0.51 arcsec / 0.11 arcsec/pix)
SLIT_LENGTH_PIXELS = 42.7  # LRS slit length in pixels (4.7 arcsec / 0.11 arcsec/pix)


class TACenterStep(Step):
    """Determine position of target source from TA verification image."""

    class_alias = "ta_center"

    spec = """
    ta_file = string(default=None). # Target acquisition image file name
    """  # noqa: E501

    reference_file_types = ["specwcs"]

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

            # Ensure TA file is provided
            if self.ta_file is None:
                log.error("No target acquisition file provided. Step will be SKIPPED.")
                result.meta.cal_step.ta_center = "SKIPPED"
                return result

            # Check that this is a point source
            if input_model.meta.target.srctype != "POINT":
                log.error(
                    "TA centering is only implemented for point sources. Step will be SKIPPED."
                )
                result.meta.cal_step.ta_center = "SKIPPED"
                return result

            # Check exposure type
            exp_type = input_model.meta.exposure.type
            if exp_type not in ["MIR_LRS-FIXEDSLIT", "MIR_LRS-SLITLESS"]:
                log.error(
                    "TA centering is only implemented for MIR_LRS-FIXEDSLIT and"
                    " MIR_LRS-SLITLESS modes. Step will be SKIPPED."
                )
                result.meta.cal_step.ta_center = "SKIPPED"
                return result

            # Read the TA image
            with dm.open(self.ta_file) as ta_model:
                log.info(f"Performing TA centering using file: {self.ta_file}")
                ta_image = ta_model.data
                wavelength = _get_wavelength(ta_model.meta.instrument.filter)
                if wavelength is None:
                    log.error(
                        "Unknown filter; cannot determine wavelength for TA centering."
                        " Step will be SKIPPED."
                    )
                    result.meta.cal_step.ta_center = "SKIPPED"
                    return result

            # read specwcs to get necessary reference points on detector
            reffile = self.get_reference_file(input_model, "specwcs")
            refmodel = dm.MiriLRSSpecwcsModel(reffile)

            is_slit = exp_type == "MIR_LRS-FIXEDSLIT"
            if is_slit:
                log.info("Fitting Airy disk with slit mask model for LRS FIXEDSLIT mode")

                ref_center = (refmodel.meta.x_ref, refmodel.meta.y_ref)

                # Create cutout around reference center
                cutout, cutout_origin = _cutout_center(ta_image, ref_center)

                # Fit on the cutout
                x_center_cutout, y_center_cutout = _fit_airy_disk(
                    cutout,
                    wavelength,
                    slit_center=ref_center,
                    cutout_origin=cutout_origin,
                    use_slit_model=True,
                )

                # Transform back to full-frame coordinates
                x_center = x_center_cutout + cutout_origin[0]
                y_center = y_center_cutout + cutout_origin[1]
            else:
                log.info("Fitting Airy disk for slitless mode")

                ref_center = (refmodel.meta.x_ref_slitless, refmodel.meta.y_ref_slitless)

                # Create cutout around reference center
                cutout, cutout_origin = _cutout_center(ta_image, ref_center)

                # Fit on the cutout
                x_center_cutout, y_center_cutout = _fit_airy_disk(
                    cutout, wavelength, use_slit_model=False
                )

                # Transform back to full-frame coordinates
                x_center = x_center_cutout + cutout_origin[0]
                y_center = y_center_cutout + cutout_origin[1]

            # Set completion status
            result.x_center = x_center
            result.y_center = y_center
            result.meta.cal_step.ta_center = "COMPLETE"

        return result


def _fit_airy_disk(
    ta_image, wavelength, use_slit_model=False, slit_center=None, cutout_origin=None
):
    """
    Fit an Airy disk model to target acquisition image data.

    Parameters
    ----------
    ta_image : ndarray
        2D target acquisition image data.
    wavelength : float
        Wavelength in microns for calculating the Airy disk size.
    use_slit_model : bool, optional
        If True, fit an Airy disk multiplied by the slit mask model.
        If False, fit a plain Airy disk.
    slit_center : tuple of float, optional
        (x_center, y_center) position of the slit center in full-frame coordinates.
        Required when use_slit_model=True.
    cutout_origin : tuple of int, optional
        (x_origin, y_origin) position of the cutout origin in full-frame coordinates.
        Required when use_slit_model=True to transform slit_center to cutout coordinates.

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
        # Transform slit center to cutout coordinates
        slit_center_cutout = (slit_center[0] - cutout_origin[0], slit_center[1] - cutout_origin[1])

        # Create compound model: Airy disk * slit mask
        airy_init = models.AiryDisk2D(
            amplitude=amplitude_guess, x_0=x_guess, y_0=y_guess, radius=radius_pixels
        )
        slit_mask_model = SlitMask(x_center=slit_center_cutout[0], y_center=slit_center_cutout[1])
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
        airy_fit = fitter(airy_init, x, y, ta_image)
        x_center = airy_fit.x_0.value
        y_center = airy_fit.y_0.value

    log.info(f"Fitted center: x={x_center:.2f}, y={y_center:.2f}")
    return x_center, y_center


def _cutout_center(image, center, size=20):
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


class SlitMask(Fittable2DModel):
    """
    Model for MIRI LRS slit mask with subpixel accuracy.

    Parameters
    ----------
    x_center : float
        X coordinate of the slit center in pixel coordinates.
    y_center : float
        Y coordinate of the slit center in pixel coordinates.
    """

    n_inputs = 2
    n_outputs = 1

    def __init__(self, x_center=0.0, y_center=0.0, **kwargs):
        self._x_center = x_center
        self._y_center = y_center

        # Define slit boundaries in pixel coordinates
        self.y_min = self._y_center - SLIT_WIDTH_PIXELS / 2
        self.y_max = self._y_center + SLIT_WIDTH_PIXELS / 2
        self.x_min = self._x_center - SLIT_LENGTH_PIXELS / 2
        self.x_max = self._x_center + SLIT_LENGTH_PIXELS / 2

        super().__init__(**kwargs)

    def evaluate(self, x, y):
        """
        Compute the fractional overlap between pixels and the slit aperture at given coordinates.

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
        # Pixel extends from (x, y) to (x+1, y+1) in index coordinates
        # Compute overlap in X direction
        overlap_x_min = np.maximum(x, self.x_min)
        overlap_x_max = np.minimum(x + 1, self.x_max)
        overlap_x = np.maximum(0, overlap_x_max - overlap_x_min)

        # Compute overlap in Y direction
        overlap_y_min = np.maximum(y, self.y_min)
        overlap_y_max = np.minimum(y + 1, self.y_max)
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
