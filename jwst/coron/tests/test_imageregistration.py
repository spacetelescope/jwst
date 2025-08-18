import numpy as np
import pytest

from jwst.coron import imageregistration


def test_register_2d_image(target_image, psf_image):
    """Test that image arrays can be registered."""
    reference = target_image.data
    target = psf_image.data
    aligned, shifts = imageregistration.align_array(reference, target)
    assert aligned.shape == target.shape
    assert not np.allclose(aligned, target)
    assert not np.any(shifts == 0)


def test_register_invalid_array():
    reference = np.arange(10)
    target = np.arange(10)
    with pytest.raises(ValueError, match="must be either a 2D or 3D array"):
        imageregistration.align_array(reference, target)


def test_fourier_imshift():
    image = np.full((3, 10, 10), 0.0)
    image[0, 1, 1] = 1.0
    image[1, 2, 2] = 1.0
    image[2, 3, 3] = 1.0
    shift = [[-1, -1], [-2, -2], [-3, -3]]

    shifted = imageregistration.fourier_imshift(image, shift)
    np.testing.assert_allclose(shifted[:, 0, 0], 1.0)
    np.testing.assert_allclose(shifted[:, 1:, :], 0.0, atol=1e-7)
    np.testing.assert_allclose(shifted[:, :, 1:], 0.0, atol=1e-7)


def test_fourier_imshift_wrong_shifts():
    image = np.full((3, 10, 10), 0.0)
    shift = [[-1, -1], [-2, -2]]
    with pytest.raises(ValueError, match="number of provided shifts must be equal"):
        imageregistration.fourier_imshift(image, shift)


def test_fourier_imshift_wrong_array_shape():
    image = np.full(10, 0.0)
    shift = -1
    with pytest.raises(ValueError, match="must be either a 2D or a 3D array"):
        imageregistration.fourier_imshift(image, shift)
