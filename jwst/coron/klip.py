"""Python implementation of the KLIP algorithm."""

import numpy as np

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def klip(target_model, refs_model, truncate):
    """
    Apply KLIP algorithm to science data.

    Parameters
    ----------
    target_model : CubeModel
        The input images of the target (NINTS x NROWS x NCOLS). Multiple integrations
        within a single exposure are stacked along the first (NINTS) axis of the data arrays.
    refs_model : CubeModel
        The input 3D stack of reference images (NINTS_PSF x NROWS x NCOLS). The first
        (NINTS_PSF) axis is the stack of aligned PSF integrations for that
        target image.
    truncate : int
        Indicates how many rows to keep in the Karhunen-Loeve transform.

    Returns
    -------
    output_target : CubeModel
        Science target Cubemodel with PSF subtracted
    output_psf : CubeModel
        CubeModel of PSF fitted to target image
    """
    # Initialize the output models as copies of the input target model
    output_target = target_model.copy()
    output_psf = target_model.copy()

    # Loop over the target integrations
    for i in range(target_model.data.shape[0]):
        # Load the target data array and flatten it from 2-D to 1-D
        target = target_model.data[i].astype(np.float64)
        tshape = target.shape
        target = target.reshape(-1)

        # Load the reference psf arrays and flatten them from 3-D to 2-D
        refs = refs_model.data.astype(np.float64)
        rshape = refs.shape
        nrefs = rshape[0]
        refs = refs.reshape(nrefs, rshape[1] * rshape[2])

        # Make each ref image have zero mean
        for k in range(nrefs):
            refs[k] -= np.mean(refs[k], dtype=np.float64)

        # Compute Karhunen-Loeve transform of ref images and normalize vectors
        klvect, eigval, eigvect = karhunen_loeve_transform(refs, normalize=True)

        # Truncate the Karhunen-Loeve vectors
        klvect = klvect[:truncate]

        # Compute the PSF fit to the target image
        psfimg = np.dot(klvect.T, np.dot(target, klvect.T))

        # Subtract the PSF fit from the target image
        outimg = target - target.mean()
        outimg -= psfimg

        # Unflatten the PSF and subtracted target images from 1-D to 2-D
        # and copy them to the output models
        psfimg = psfimg.reshape(tshape)
        output_psf.data[i] = psfimg
        outimg = outimg.reshape(tshape)
        output_target.data[i] = outimg

        # Compute the ERR for the fitted target image:
        # the ERR is taken as the std-dev of the KLIP results for all of the
        # PSF reference images.
        #
        # First, apply the PSF fit to each PSF reference image
        refs_fit = refs * 0.0
        for k in range(nrefs):
            refs_fit[k] = refs[k] - np.dot(klvect.T, np.dot(refs[k], klvect.T))

        # Now take the standard deviation of the results
        output_target.err[i] = np.std(refs_fit, 0).reshape(tshape)

    return output_target, output_psf


def karhunen_loeve_transform(m, normalize=False):
    """
    Calculate Karhunen-Loeve Transform of the input.

    Parameters
    ----------
    m : numpy.ndarray
        The array of flattened, background subtracted reference arrays
    normalize : bool
        If True, normalize the returned transform

    Returns
    -------
    klvect : numpy.ndarray
        The Karhunen-Loeve Transform of the input arrays
    eigval : numpy.ndarray
        Array of eigenvalues
    eigvect : numpy.ndarray
        Matrix of eigenvectors
    """
    eigval, eigvect = np.linalg.eigh(np.cov(m))

    # Sort eigenvalues (replicate Mathematica's behaviour):
    idx = eigval.argsort()[::-1]
    eigval = eigval[idx]
    eigvect = eigvect[:, idx]

    # Compute Karhunen-Loeve transform:
    klvect = np.dot(eigvect.T, m)

    if normalize:
        for k in range(len(klvect)):
            klvect[k] /= np.linalg.norm(klvect[k])

    return klvect, eigval, eigvect
