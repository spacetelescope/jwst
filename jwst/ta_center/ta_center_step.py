import logging

import numpy as np
from astropy.modeling import Fittable2DModel, fitting, models

import jwst.datamodels as dm
from jwst.assign_wcs.miri import retrieve_filter_offset
from jwst.pathloss.pathloss import calculate_pathloss_vector
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
    skip = boolean(default=True). # Skip this step by default
    """  # noqa: E501

    reference_file_types = ["specwcs", "pathloss", "filteroffset"]

    def process(self, step_input):
        """
        Process target acquisition data.

        Parameters
        ----------
        step_input : file path, ModelContainer, or DataModel
            The input data model or association.

        Returns
        -------
        result : DataModel
            The output data model with TA centering applied.
        """
        result = self.prepare_output(step_input)
        if isinstance(result, dm.ModelContainer):
            # Extract science and TA image from container
            result, ta_model = self._ta_image_from_container(result)
            self.ta_file = ta_model

        # Ensure TA file is provided
        if str(self.ta_file).lower() == "none" or self.ta_file is None:
            log.error("No target acquisition file provided. Step will be SKIPPED.")
            result.meta.cal_step.ta_center = "SKIPPED"
            return result

        # Check that this is a point source
        if result.meta.target.source_type != "POINT":
            log.error("TA centering is only implemented for point sources. Step will be SKIPPED.")
            result.meta.cal_step.ta_center = "SKIPPED"
            return result

        # Check exposure type
        exp_type = result.meta.exposure.type
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
        reffile = self.get_reference_file(result, "specwcs")
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
            pathloss_file = self.get_reference_file(result, "pathloss")
        else:
            log.info("Fitting Airy disk for slitless mode")
            ref_center = (refmodel.meta.x_ref_slitless, refmodel.meta.y_ref_slitless)
            pathloss_file = None

        try:
            x_center, y_center = center_from_ta_image(
                ta_image,
                wavelength,
                ref_center,
                subarray_origin=(xstart, ystart),
                pathloss_file=pathloss_file,
            )
        except Exception as e:
            log.error(f"Error during TA centering: {e}. Step will be SKIPPED.")
            result.meta.cal_step.ta_center = "SKIPPED"
            return result

        # Apply filter offsets
        filteroffset_file = self.get_reference_file(ta_model, "filteroffset")
        with dm.FilteroffsetModel(filteroffset_file) as filteroffset:
            col_offset, row_offset = retrieve_filter_offset(
                filteroffset, ta_model.meta.instrument.filter
            )
            log.info(f"Applying filter offsets: column={col_offset}, row={row_offset}")
            x_center += col_offset
            y_center += row_offset

        # Set completion status
        result.source_xpos = x_center
        result.source_ypos = y_center
        result.meta.cal_step.ta_center = "COMPLETE"

        return result

    def _ta_image_from_container(self, container):
        """
        Extract the TA image from a container (association or ModelContainer).

        Parameters
        ----------
        container : ModelContainer
            The input container.

        Returns
        -------
        sci_model : DataModel
            The science data model.
        ta_model : DataModel
            The TA image data model.
        """
        sci_idx = container.ind_asn_type("science")
        sci_model = container[sci_idx[0]]
        ta_idx = container.ind_asn_type("target_acquisition")
        if not len(ta_idx):
            ta_model = None
        else:
            ta_model = container[ta_idx[0]]
        return sci_model, ta_model


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
        pathloss = dm.MirLrsPathlossModel(pathloss_file)
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
        Pathloss model to multiply the Airy disk by.

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

    if pathloss_model is not None:
        compound_model = airy_init * pathloss_model

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
        dx_arcsec = dx_pixels * PIXSCALE
        dy_arcsec = dy_pixels * PIXSCALE
        loss = np.empty_like(dx_pixels)
        # In the future we should vectorize calculate_pathloss_vector to avoid loop.
        # This isn't preventatively slow, so not a high priority
        for i in range(len(dx_arcsec)):
            _, pathloss_vector, _is_inside_slit = calculate_pathloss_vector(
                self.pathloss_data, self.pathloss_wcs, dx_arcsec[i], dy_arcsec[i], calc_wave=False
            )
            loss[i] = pathloss_vector[self.wave_idx]
        loss = loss.reshape(x.shape)
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
