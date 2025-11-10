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
# SLIT_WIDTH_PIXELS = 4.6  # LRS slit width in pixels (0.51 arcsec / 0.11 arcsec/pix)
# SLIT_LENGTH_PIXELS = 42.7  # LRS slit length in pixels (4.7 arcsec / 0.11 arcsec/pix)
SLITMASK_LL = (302, 295)  # approximate lower-left corner of slit mask
SLITMASK_UR = (352, 309)  # approximate upper-right corner of slit mask


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
            if input_model.meta.target.source_type != "POINT":
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

            # Get subarray information from TA image
            subarray_name = ta_model.meta.subarray.name
            xstart = ta_model.meta.subarray.xstart  # 1-indexed in FITS
            ystart = ta_model.meta.subarray.ystart  # 1-indexed in FITS
            log.info(f"TA subarray: {subarray_name}, origin: ({xstart}, {ystart})")

            is_slit = exp_type == "MIR_LRS-FIXEDSLIT"
            if is_slit:
                log.info("Fitting Airy disk with slit mask model for LRS FIXEDSLIT mode")
                ref_center = (refmodel.meta.x_ref, refmodel.meta.y_ref)
                x_center, y_center = center_from_ta_image(
                    ta_image,
                    wavelength,
                    ref_center,
                    subarray_origin=(xstart, ystart),
                    use_slit_model=True,
                )
            else:
                log.info("Fitting Airy disk for slitless mode")
                ref_center = (refmodel.meta.x_ref_slitless, refmodel.meta.y_ref_slitless)
                x_center, y_center = center_from_ta_image(
                    ta_image,
                    wavelength,
                    ref_center,
                    subarray_origin=(xstart, ystart),
                    use_slit_model=False,
                )

            # Set completion status
            log.info(f"Fitted center: x={x_center:.2f}, y={y_center:.2f}")
            result.x_center = x_center
            result.y_center = y_center
            result.meta.cal_step.ta_center = "COMPLETE"

        return result


def center_from_ta_image(
    ta_image, wavelength, ref_center, subarray_origin=(1, 1), use_slit_model=False
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
    use_slit_model : bool, optional
        If True, fit an Airy disk multiplied by the slit mask model.

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
        f"Reference center: full-frame=({ref_center[0]:.2f}, {ref_center[1]:.2f}), "
        f"subarray=({ref_center_subarray[0]:.2f}, {ref_center_subarray[1]:.2f})"
    )

    # Create cutout around reference center (in subarray coordinates)
    cutout, cutout_origin = _cutout_center(ta_image, ref_center_subarray)

    # Fit on the cutout
    x_center_cutout, y_center_cutout = _fit_airy_disk(
        cutout,
        wavelength,
        slit_center=ref_center_subarray,
        cutout_origin=cutout_origin,
        use_slit_model=use_slit_model,
    )

    # Transform back to subarray coordinates
    x_center_subarray = x_center_cutout + cutout_origin[0]
    y_center_subarray = y_center_cutout + cutout_origin[1]

    # Transform from subarray to full-frame detector coordinates
    x_center = x_center_subarray + (subarray_origin[0] - 1)
    y_center = y_center_subarray + (subarray_origin[1] - 1)

    log.info(
        f"Fitted center: subarray=({x_center_subarray:.2f}, {y_center_subarray:.2f}), "
        f"full-frame=({x_center:.2f}, {y_center:.2f})"
    )

    return x_center, y_center


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

    # Filter non-finite values from data
    mask = np.isfinite(ta_image)
    if not np.any(mask):
        log.error("All pixels are non-finite; cannot fit")
        raise ValueError("All pixels contain non-finite values")
    log.debug(f"Excluding {np.sum(~mask)} non-finite values from fit")

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Initial guess for center (centroid of finite values)
    ta_image_masked = np.where(mask, ta_image, 0.0)
    y_guess = np.sum(y * ta_image_masked) / np.sum(ta_image_masked)
    x_guess = np.sum(x * ta_image_masked) / np.sum(ta_image_masked)
    amplitude_guess = np.max(ta_image[mask])

    # Create initial Airy disk model
    airy_init = models.AiryDisk2D(
        amplitude=amplitude_guess, x_0=x_guess, y_0=y_guess, radius=radius_pixels
    )
    # Allow radius to vary a bit (0.99 to 1.5 times diffraction limit)
    airy_init.radius.bounds = (0.99 * radius_pixels, 1.5 * radius_pixels)

    if use_slit_model:
        # Add a model of the slit mask
        slit_center_cutout = (slit_center[0] - cutout_origin[0], slit_center[1] - cutout_origin[1])
        slit_mask_model = SlitMask(x_center=slit_center_cutout[0], y_center=slit_center_cutout[1])
        compound_model = airy_init * slit_mask_model

        # Fit the compound model to filtered data
        fitter = fitting.LevMarLSQFitter()
        fitted_model = fitter(compound_model, x[mask], y[mask], ta_image[mask])
        # Extract parameters from the left side of the compound model (Airy disk)
        x_center = fitted_model.left.x_0.value
        y_center = fitted_model.left.y_0.value

    else:
        # Fit the model to filtered data
        fitter = fitting.LevMarLSQFitter()
        airy_fit = fitter(airy_init, x[mask], y[mask], ta_image[mask])
        x_center = airy_fit.x_0.value
        y_center = airy_fit.y_0.value

    return x_center, y_center


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

        # Calculate offset from reference center to actual slit mask center
        # SLITMASK_LL and SLITMASK_UR are in full-frame coordinates
        slitmask_x_center = (SLITMASK_LL[0] + SLITMASK_UR[0]) / 2.0
        slitmask_y_center = (SLITMASK_LL[1] + SLITMASK_UR[1]) / 2.0

        # Define slit boundaries in pixel coordinates, centered on x_center, y_center
        # Apply offset from slitmask center
        offset_x = x_center - slitmask_x_center
        offset_y = y_center - slitmask_y_center

        self.x_min = SLITMASK_LL[0] + offset_x
        self.x_max = SLITMASK_UR[0] + offset_x
        self.y_min = SLITMASK_LL[1] + offset_y
        self.y_max = SLITMASK_UR[1] + offset_y

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
        mask_value[mask_value == 0.0] = np.nan

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
