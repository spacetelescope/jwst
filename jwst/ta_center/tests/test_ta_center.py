import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import models

from jwst.ta_center.ta_center_step import (
    JWST_DIAMETER,
    PIXSCALE,
    SlitMask,
    TACenterStep,
    _get_wavelength,
)

X_REF_SLIT = 326.13
Y_REF_SLIT = 300.7
X_REF_SLITLESS = 38.5
Y_REF_SLITLESS = 829.0
MIRI_DETECTOR_SHAPE = (1024, 1032)  # (ny, nx) for MIRI imager


def _create_slit_mask(image_shape, slit_center):
    """
    Create a slit mask with subpixel accuracy.

    For a rectangular slit aligned with the pixel grid, computes the fractional
    overlap between each pixel and the slit aperture analytically.

    Parameters
    ----------
    image_shape : tuple
        Shape of the image (ny, nx).
    slit_center : tuple
        Center position of the slit (x, y) in pixel coordinates.

    Returns
    -------
    weights : ndarray
        2D array of fractional pixel coverage values (0.0 to 1.0).
    """
    ny, nx = image_shape

    # Create pixel coordinate grids (pixel indices: 0, 1, 2, ...)
    y_indices, x_indices = np.mgrid[0:ny, 0:nx]

    # Use the callable function to evaluate mask at all pixel positions
    slit_mask_model = SlitMask(x_center=slit_center[0], y_center=slit_center[1])
    weights = slit_mask_model(x_indices, y_indices)

    return weights


def make_slitless_data(wavelength=15.0, offset=(0, 0)):
    """
    Make a fake MIRI LRS slitless TAQ verification image.

    Parameters
    ----------
    wavelength : float, optional
        Wavelength in microns.
    offset : tuple, optional
        (x,y) offset of the source from the reference center in pixels.

    Returns
    -------
    data : 2D ndarray
        Simulated slitless TAQ image.
    """
    # Create coordinate grids using full MIRI detector size
    y, x = np.mgrid[0 : MIRI_DETECTOR_SHAPE[0], 0 : MIRI_DETECTOR_SHAPE[1]]

    # Calculate diffraction-limited Airy disk radius
    # First zero of Airy disk: theta = 1.22 * lambda / D (in radians)
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slitless reference position
    x_center = X_REF_SLITLESS + offset[0]
    y_center = Y_REF_SLITLESS + offset[1]

    # Create Airy disk model at diffraction-limited resolution
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data = airy(x, y)

    # # Add some noise
    noise = np.random.normal(0, 1, data.shape)
    data += noise

    return data


def make_taq_image(make_data_func, filt, offset=(0, 0)):
    """
    Generate a TAQ image for testing.

    Parameters
    ----------
    make_data_func : function
        Function to generate the TAQ image data.  Must take in exactly
        two arguments: wavelength and offset.
    filt : str
        Filter name, e.g., 'F560W', 'F770W', etc.
    offset : tuple, optional
        (x,y) offset of the source from the center of the frame in pixels.

    Returns
    -------
    model : ImageModel
        Simulated TAQ image model.
    """
    wavelength = _get_wavelength(filt)
    data = make_data_func(wavelength, offset)

    model = dm.ImageModel()
    model.data = data
    model.meta.instrument.filter = filt

    # Add subarray metadata (use FULL for now)
    model.meta.subarray.name = "FULL"
    model.meta.subarray.xstart = 1  # 1-indexed FITS convention
    model.meta.subarray.ystart = 1
    model.meta.subarray.xsize = MIRI_DETECTOR_SHAPE[1]  # nx
    model.meta.subarray.ysize = MIRI_DETECTOR_SHAPE[0]  # ny

    return model


@pytest.fixture
def slitless_taq_image(tmp_path):
    model = make_taq_image(make_slitless_data, "F1500W", offset=(2, -3))

    # Add NaN values near the slitless source position to test that this still works ok
    model.data[int(Y_REF_SLITLESS) + 5, int(X_REF_SLITLESS) + 3] = np.nan
    model.data[int(Y_REF_SLITLESS) - 4, int(X_REF_SLITLESS) - 2] = np.inf

    filepath = tmp_path / "slitless_taq.fits"
    model.save(filepath)
    return str(filepath)


def make_slit_data(wavelength=15.0, offset=(0, 0)):
    """
    Make a fake MIRI LRS slit TAQ verification image.

    The PSF is truncated by the slit edges, simulating what happens in
    FIXEDSLIT mode observations.

    Parameters
    ----------
    wavelength : float, optional
        Wavelength in microns.
    offset : tuple, optional
        (x,y) offset of the source from the reference center in pixels.

    Returns
    -------
    data : 2D ndarray
        Simulated slit TAQ image with PSF truncated by slit edges.
    """

    # Create coordinate grids using full MIRI detector size
    y, x = np.mgrid[0 : MIRI_DETECTOR_SHAPE[0], 0 : MIRI_DETECTOR_SHAPE[1]]

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slit reference position
    x_center = X_REF_SLIT + offset[0]
    y_center = Y_REF_SLIT + offset[1]

    # Create full Airy disk model
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data_full = airy(x, y)

    # Create slit mask with subpixel accuracy using the same helper as the algorithm
    slit_weights = _create_slit_mask(data_full.shape, slit_center=(X_REF_SLIT, Y_REF_SLIT))

    # Apply slit weights to the data (multiply by fractional coverage)
    data = data_full * slit_weights

    # Add noise only where there's signal (weighted by slit coverage)
    noise = np.random.normal(0, 1, data.shape)
    data += noise * slit_weights

    # import matplotlib.pyplot as plt
    # from matplotlib.colors import LogNorm

    # # Define cutout region around slit center
    # cutout_size = 60  # pixels on each side
    # x_min = int(X_REF_SLIT - cutout_size)
    # x_max = int(X_REF_SLIT + cutout_size)
    # y_min = int(Y_REF_SLIT - cutout_size)
    # y_max = int(Y_REF_SLIT + cutout_size)

    # # Create cutouts
    # data_full_cutout = data_full[y_min:y_max, x_min:x_max]
    # slit_weights_cutout = slit_weights[y_min:y_max, x_min:x_max]
    # data_cutout = data[y_min:y_max, x_min:x_max]

    # fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # # Plot 1: Full Airy disk (before slit truncation)
    # im1 = axes[0].imshow(data_full_cutout, cmap='gray', norm=LogNorm(), origin='lower', extent=[x_min, x_max, y_min, y_max])
    # axes[0].plot(x_center, y_center, 'r+', markersize=15, markeredgewidth=2)
    # axes[0].set_title('Full Airy Disk (Pre-Slit)')
    # axes[0].set_xlabel('X (pixels)')
    # axes[0].set_ylabel('Y (pixels)')
    # plt.colorbar(im1, ax=axes[0], label='Counts')

    # # Plot 2: Slit weights (showing fractional pixel coverage)
    # im2 = axes[1].imshow(slit_weights_cutout, cmap='gray', origin='lower', vmin=0, vmax=1, extent=[x_min, x_max, y_min, y_max])
    # axes[1].plot(x_center, y_center, 'r+', markersize=15, markeredgewidth=2)
    # axes[1].set_title('Slit Weights')
    # axes[1].set_xlabel('X (pixels)')
    # axes[1].set_ylabel('Y (pixels)')
    # plt.colorbar(im2, ax=axes[1], label='Coverage Fraction')

    # # Plot 3: Final slit-truncated data with noise
    # im3 = axes[2].imshow(data_cutout, cmap='gray', norm=LogNorm(vmin=1), origin='lower', extent=[x_min, x_max, y_min, y_max])
    # axes[2].plot(x_center, y_center, 'r+', markersize=15, markeredgewidth=2, label='True Center')
    # axes[2].set_title(f'Slit TAQ Image\nOffset: ({offset[0]}, {offset[1]}) pix')
    # axes[2].set_xlabel('X (pixels)')
    # axes[2].set_ylabel('Y (pixels)')
    # axes[2].legend()
    # plt.colorbar(im3, ax=axes[2], label='Counts')

    # plt.tight_layout()
    # plt.show()

    return data


def _assign_metadata(model):
    """Assign necessary metadata shared between slit and slitless modes for LRS."""
    model.meta.instrument.name = "MIRI"
    model.meta.target.source_type = "POINT"
    model.meta.instrument.filter = "F1500W"
    model.meta.instrument.detector = "MIRIMAGE"

    model.meta.observation.date = "2025-11-10"
    model.meta.observation.time = "12:00:00"


@pytest.fixture
def input_model_slit():
    model = dm.ImageModel()
    _assign_metadata(model)
    model.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    return model


@pytest.fixture
def input_model_slitless():
    model = dm.ImageModel()
    _assign_metadata(model)
    model.meta.exposure.type = "MIR_LRS-SLITLESS"
    return model


def test_ta_center_slitless(input_model_slitless, slitless_taq_image):
    # Run the TA centering algorithm
    result = TACenterStep.call(input_model_slitless, ta_file=slitless_taq_image)
    x_center, y_center = result.source_xpos, result.source_ypos

    # Expected center position (reference position + offset)
    expected_x = X_REF_SLITLESS + 2
    expected_y = Y_REF_SLITLESS - 3

    # Check that the computed center is close to the expected position
    assert np.isclose(x_center, expected_x, atol=0.05), (
        f"X center {x_center} not close to expected {expected_x}"
    )
    assert np.isclose(y_center, expected_y, atol=0.05), (
        f"Y center {y_center} not close to expected {expected_y}"
    )


@pytest.mark.parametrize(
    "offset",
    [
        (0.0, 0.0),  # Perfectly centered
        (4.0, 0.0),  # Offset to the right
        (-3.0, 0.0),  # Offset to the left
        (0.0, 7.0),  # Offset upward
        (0.0, -6.0),  # Offset downward
        (3.0, -6.0),  # Offset right and down
        (-4.0, 3.6),  # Offset left and up
    ],
)
def test_ta_center_slit(input_model_slit, offset, tmp_path):
    """Test TA centering for LRS slit mode with various offsets."""

    # Generate slit data with the specified offset
    taq_image = make_taq_image(make_slit_data, "F1500W", offset=offset)

    # Add NaN values near the slit source position to test handling
    taq_image.data[int(Y_REF_SLIT) + 4, int(X_REF_SLIT) + 2] = np.nan
    taq_image.data[int(Y_REF_SLIT) - 3, int(X_REF_SLIT) - 4] = np.inf

    # Save to file
    filepath = tmp_path / f"slit_taq_{offset[0]}_{offset[1]}.fits"
    taq_image.save(filepath)

    # Run the TA centering algorithm for slit mode
    result = TACenterStep.call(input_model_slit, ta_file=str(filepath))
    x_center, y_center = result.source_xpos, result.source_ypos

    # Expected center position (reference position + offset)
    expected_x = X_REF_SLIT + offset[0]
    expected_y = Y_REF_SLIT + offset[1]

    # Check that the computed center is close to the expected position
    # The slit truncation may cause slightly larger errors, so use slightly larger tolerance
    assert np.isclose(x_center, expected_x, atol=0.05), (
        f"Offset {offset}: X center {x_center:.2f} not close to expected {expected_x:.2f}"
    )
    assert np.isclose(y_center, expected_y, atol=0.05), (
        f"Offset {offset}: Y center {y_center:.2f} not close to expected {expected_y:.2f}"
    )
