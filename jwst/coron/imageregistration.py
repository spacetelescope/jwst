import numpy as np

from scipy import optimize
from scipy.ndimage import fourier_shift

from stdatamodels.jwst.datamodels import CubeModel

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def align_fourier_lsq(reference, target, mask=None):
    """
    LSQ optimization with Fourier shift alignment.

    Parameters
    ----------
    reference : numpy.ndarray
        A 2D (``NxK``) image to be aligned to

    target : numpy.ndarray
        A 2D (``NxK``) image to align to reference

    mask : numpy.ndarray, None
        A 2D (``NxK``) image indicating pixels to ignore when
        performing the minimization. The masks acts as
        a weighting function in performing the fit.

    Returns
    -------
    results : numpy.ndarray
        A 1D vector containing (`xshift`, `yshift`, `beta`) values from
        LSQ optimization, where `xshift` and `yshift` are the misalignment
        of target from reference and `beta` is the fraction by which the
        target intensity must be reduced to match the intensity
        of the reference.
    """
    init_pars = [0.0, 0.0, 1.0]
    results, _ = optimize.leastsq(
        shift_subtract,
        init_pars,
        args=(reference, target, mask),
        xtol=1e-15,
        ftol=1e-15,
    )
    return results


def shift_subtract(params, reference, target, mask=None):
    """
    Use Fourier Shift theorem for subpixel shifts.

    Parameters
    ----------
    params : tuple
        Tuple of xshift, yshift, beta

    reference : numpy.ndarray
        See :py:func:`align_fourier_lsq`

    target : numpy.ndarray
        See :py:func:`align_fourier_lsq`

    mask : numpy.ndarray, None
        See :py:func:`align_fourier_lsq`

    Returns
    -------
    residual : numpy.ndarray
        1D numpy.ndarray of target-reference residual after
        applying shift and intensity fraction.
    """
    shift = params[:2]
    beta = params[2]

    offset = fourier_imshift(reference, shift)

    if mask is not None:
        return ((target - beta * offset) * mask).ravel()
    else:
        return (target - beta * offset).ravel()


def fourier_imshift(image, shift):
    """
    Shift an image by use of Fourier shift theorem.

    Parameters
    ----------
    image : numpy.ndarray
        A 2D (``NxK``) or 3D (``LxNxK``) image.

    shift : numpy.ndarray
        A 1D or 2D array of shape ``Lx2`` containing pixel values by which
        to shift image slices in the X and Y directions.

    Returns
    -------
    offset : numpy.ndarray
        Shifted image
    """
    ndim = len(image.shape)

    if ndim == 2:
        shift = np.asanyarray(shift)[:2]
        offset_image = fourier_shift(np.fft.fftn(image), shift[::-1])
        offset = np.fft.ifftn(offset_image).real

    elif ndim == 3:
        nslices = image.shape[0]
        shift = np.asanyarray(shift)[:, :2]
        if shift.shape[0] != nslices:
            raise ValueError(
                "The number of provided shifts must be equal "
                "to the number of slices in the input image."
            )

        offset = np.empty_like(image, dtype=float)
        for k in range(nslices):
            offset[k] = fourier_imshift(image[k], shift[k])

    else:
        raise ValueError("Input image must be either a 2D or a 3D array.")

    return offset


def align_array(reference, target, mask=None, return_aligned=True):
    """
    Compute shifts and realign target image to reference.

    Shifts are computed between target image (or image "slices") and the reference
    image and input target images are realigned to the reference image.

    Parameters
    ----------
    reference : numpy.ndarray
        A 2D (``NxK``) reference image to which input images will be aligned.

    target : numpy.ndarray
        A 2D (``NxK``) or 3D (``MxNxK`` first index used to select slices)
        image(s) that need to be aligned to the reference image.

    mask : numpy.ndarray, None
        A 2D (``NxK``) image indicating pixels to ignore when performing the
        minimization. The masks acts as a weighting function in performing
        the fit.
    return_aligned : bool
        Boolean that indicates whether the aligned image is returned

    Returns
    -------
    aligned : numpy.ndarray
        The aligned image (2D or 3D image of the same shape as input target image).  This
        is only returned if the input keyword `returned_aligned` is True.
    shifts : numpy.ndarray
        1D vector of three elements in the case of 2D input `target` image of
        (xshift, yshift, beta) values from LSQ optimization (see :py:func:`align_fourier_lsq`
        for details) for each slice in the `target` array.
    """
    if len(target.shape) == 2:
        shifts = align_fourier_lsq(reference, target, mask=mask)
        if return_aligned:
            aligned = fourier_imshift(target, -shifts)

    elif len(target.shape) == 3:
        nslices = target.shape[0]
        shifts = np.empty((nslices, 3), dtype=np.float64)
        if return_aligned:
            aligned = np.empty_like(target)

        for m in range(nslices):
            sh = align_fourier_lsq(reference, target[m], mask=mask)
            shifts[m, :] = sh
            if return_aligned:
                aligned[m, :, :] = fourier_imshift(target[m], -sh)

    else:
        raise ValueError("Input target image must be either a 2D or 3D array.")

    if not return_aligned:
        return shifts
    return aligned, shifts


def align_models(reference, target, mask):
    """
    Align target image to reference by calculating and applying shifts.

    Parameters
    ----------
    reference : CubeModel
        3D (``LxNxK`` first index used to select 2D slices)
        reference image to which target images will be aligned.
    target : CubeModel
        3D (``MxNxK`` first index used to select 2D slices)
        image(s) that need to be aligned to the reference image.
    mask : ImageModel, None
        A 2D (``NxK``) image indicating pixels to ignore when performing the
        minimization. Mask acts as a weighting function in performing
        the fit.

    Returns
    -------
    output_model : CubeModel
        A CubeModel containing aligned copies of the input ``target``
        cubes aligned to the first slice in the input ``reference`` cube.
    """
    # Create output CubeModel of required dimensions. Since all science integrations
    # are assumed to have the same shift, the output is just a shifted copy of the
    # 3-D PSF cube
    output_model = CubeModel(target.shape)
    output_model.update(target)

    # Compute the shifts of the PSF ("target") images relative to
    # the science ("reference") image in the first integration
    shifts = align_array(reference.data[0], target.data, mask=mask.data, return_aligned=False)

    # Apply the shifts to the PSF images
    output_model.data = fourier_imshift(target.data, -shifts)

    # Apply the same shifts to the PSF error arrays, if they exist
    if target.err is not None:
        output_model.err = fourier_imshift(target.err, -shifts)

    # TODO: in the future we need to add shifts and other info (such as
    # slice ID from the reference image to which target was aligned)
    # to output cube metadata (or property).
    return output_model
