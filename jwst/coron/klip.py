"""
    Python implementation of the KLIP algorithm based on the
    Mathematica script from Remi Soummer.

:Authors: Mihai Cara

"""

import numpy as np

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def klip(target_model, refs_model, truncate):
    """
    Parameters
    ----------
    target_model : CubeModel (NINTS x NROWS x NCOLS)
        The input images of the target. Multiple integrations within
        a single exposure are stacked along the first (NINTS) axis of
        the data arrays.

    refs_model : QuadModel (NTARG x NINTS x NROWS x NCOLS)
        The input 4D stack of reference images. The first (NTARG) axis
        corresponds to the index of each target integration. The second
        (NINTS) axis is the stack of aligned PSF integrations for that
        target image. The length of the target_model first (NINTS) axis
        should be equal to the length of the refs_model first (NTARG)
        axis (i.e. one stack of aligned PSF images for each target
        integration).

    truncate : int
        Indicates how many rows to keep in the Karhunen-Loeve transform.
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
        refs = refs_model.data[i].astype(np.float64)
        rshape = refs.shape
        nrefs = rshape[0]
        refs = refs.reshape(nrefs, rshape[1] * rshape[2])

        # Make each ref image have zero mean
        for k in range(nrefs):
            refs[k] -= np.mean(refs[k], dtype=np.float64)

        # Compute Karhunen-Loeve transform of ref images and normalize vectors
        klvect, eigval, eigvect = KarhunenLoeveTransform(refs, normalize=True)

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


def KarhunenLoeveTransform(m, normalize=False):
    """
    Returns Karhunen-Loeve Transform of the input, eigenvalues, and
    a matrix of eigenvectors.

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
