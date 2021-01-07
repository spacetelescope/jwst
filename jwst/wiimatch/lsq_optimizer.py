"""
A module that provides core algorithm for optimal matching of backgrounds of
N-dimensional images using (multi-variate) polynomials.

:Author: Mihai Cara (contact: help@stsci.edu)


"""
import numpy as np

from .utils import create_coordinate_arrays


__all__ = ['build_lsq_eqs', 'pinv_solve', 'rlu_solve']


def build_lsq_eqs(images, masks, sigmas, degree, center=None,
                  image2world=None, center_cs='image'):
    r"""
    Build system of linear equations whose solution would provide image
    intensity matching in the least squares sense.

    Parameters
    ----------
    images : list of numpy.ndarray
        A list of 1D, 2D, etc. `numpy.ndarray` data array whose "intensities"
        must be "matched". All arrays must have identical shapes.

    masks : list of numpy.ndarray
        A list of `numpy.ndarray` arrays of same length as ``images``.
        Non-zero mask elements indicate valid data in the corresponding
        ``images`` array. Mask arrays must have identical shape to that of
        the arrays in input ``images``.

    sigmas : list of numpy.ndarray
        A list of `numpy.ndarray` data array of same length as ``images``
        representing the uncertainties of the data in the corresponding array
        in ``images``. Uncertainty arrays must have identical shape to that of
        the arrays in input ``images``.

    degree : iterable
        A list of polynomial degrees for each dimension of data arrays in
        ``images``. The length of the input list must match the dimensionality
        of the input images.

    center : iterable, None, optional
        An iterable of length equal to the number of dimensions of images in
        ``images`` parameter that indicates the center of the coordinate system
        in **image** coordinates when ``center_cs`` is ``'image'`` otherwise
        center is assumed to be in **world** coordinates (when ``center_cs``
        is ``'world'``). When ``center`` is `None` then ``center`` is
        set to the middle of the "image" as ``center[i]=image.shape[i]//2``.
        If ``image2world`` is not `None` and ``center_cs`` is ``'image'``,
        then supplied center will be converted to world coordinates.

    image2world : function, None, optional
        Image-to-world coordinates transformation function. This function
        must be of the form ``f(x,y,z,...)`` and accept a number of arguments
        `numpy.ndarray` arguments equal to the dimensionality of images.

    center_cs : {'image', 'world'}, optional
        Indicates whether ``center`` is in image coordinates or in world
        coordinates. This parameter is ignored when ``center`` is set to
        `None`: it is assumed to be `False`. ``center_cs`` *cannot be*
        ``'world'`` when ``image2world`` is `None` unless ``center`` is `None`.

    Returns
    -------
    a : numpy.ndarray
        A 2D `numpy.ndarray` that holds the coefficients of the linear system
        of equations.

    b : numpy.ndarray
        A 1D `numpy.ndarray` that holds the free terms of the linear system of
        equations.

    coord_arrays : list
        A list of `numpy.ndarray` coordinate arrays each of ``images[0].shape``
        shape.

    eff_center : tuple
        A tuple of coordinates of the effective center as used in generating
        coordinate arrays.

    coord_system : {'image', 'world'}
        Coordinate system of the coordinate arrays and returned ``center``
        value.

    Notes
    -----
    :py:func:`build_lsq_eqs` builds a system of linear equations

    .. math::
        a \cdot c = b

    whose solution :math:`c` is a set of coefficients of (multivariate)
    polynomials that represent the "background" in each input image (these are
    polynomials that are "corrections" to intensities of input images) such
    that the following sum is minimized:

    .. math::
        L = \sum^N_{n,m=1,n \neq m} \sum_k \frac{\left[I_n(k) - I_m(k) - P_n(k) + P_m(k)\right]^2}{\sigma^2_n(k) + \sigma^2_m(k)}

    In the above equation, index :math:`k=(k_1,k_2,...)` labels a position
    in input image's pixel grid [NOTE: all input images share a common
    pixel grid].

    "Background" polynomials :math:`P_n(k)` are defined through the
    corresponding coefficients as:

    .. math::
        P_n(k_1,k_2,...) = \sum_{d_1=0,d_2=0,...}^{D_1,D_2,...} c_{d_1,d_2,...}^n \cdot k_1^{d_1} \cdot k_2^{d_2}  \cdot \ldots .

    Coefficients :math:`c_{d_1,d_2,...}^n` are arranged in the vector :math:`c`
    in the following order:

    .. math::
        (c_{0,0,\ldots}^1,c_{1,0,\ldots}^1,\ldots,c_{0,0,\ldots}^2,c_{1,0,\ldots}^2,\ldots)

    Examples
    --------
    >>> import numpy as np
    >>> im1 = np.zeros((5, 5, 4), dtype=np.float)
    >>> cbg = 1.32 * np.ones_like(im1)
    >>> ind = np.indices(im1.shape, dtype=np.float)
    >>> im3 = cbg + 0.15 * ind[0] + 0.62 * ind[1] + 0.74 * ind[2]
    >>> mask = np.ones_like(im1, dtype=np.int8)
    >>> sigma = np.ones_like(im1, dtype=np.float)
    >>> a, b, ca, ef, cs = build_lsq_eqs([im1, im3],
    ... [mask, mask], [sigma, sigma], degree=(1,1,1), center=(0,0,0))
    >>> print(a)
    [[   50.   100.   100.   200.    75.   150.   150.   300.   -50.  -100.
        -100.  -200.   -75.  -150.  -150.  -300.]
     [  100.   300.   200.   600.   150.   450.   300.   900.  -100.  -300.
       -200.  -600.  -150.  -450.  -300.  -900.]
     [  100.   200.   300.   600.   150.   300.   450.   900.  -100.  -200.
       -300.  -600.  -150.  -300.  -450.  -900.]
     [  200.   600.   600.  1800.   300.   900.   900.  2700.  -200.  -600.
       -600. -1800.  -300.  -900.  -900. -2700.]
     [   75.   150.   150.   300.   175.   350.   350.   700.   -75.  -150.
       -150.  -300.  -175.  -350.  -350.  -700.]
     [  150.   450.   300.   900.   350.  1050.   700.  2100.  -150.  -450.
       -300.  -900.  -350. -1050.  -700. -2100.]
     [  150.   300.   450.   900.   350.   700.  1050.  2100.  -150.  -300.
       -450.  -900.  -350.  -700. -1050. -2100.]
     [  300.   900.   900.  2700.   700.  2100.  2100.  6300.  -300.  -900.
       -900. -2700.  -700. -2100. -2100. -6300.]
     [  -50.  -100.  -100.  -200.   -75.  -150.  -150.  -300.    50.   100.
        100.   200.    75.   150.   150.   300.]
     [ -100.  -300.  -200.  -600.  -150.  -450.  -300.  -900.   100.   300.
        200.   600.   150.   450.   300.   900.]
     [ -100.  -200.  -300.  -600.  -150.  -300.  -450.  -900.   100.   200.
        300.   600.   150.   300.   450.   900.]
     [ -200.  -600.  -600. -1800.  -300.  -900.  -900. -2700.   200.   600.
        600.  1800.   300.   900.   900.  2700.]
     [  -75.  -150.  -150.  -300.  -175.  -350.  -350.  -700.    75.   150.
        150.   300.   175.   350.   350.   700.]
     [ -150.  -450.  -300.  -900.  -350. -1050.  -700. -2100.   150.   450.
        300.   900.   350.  1050.   700.  2100.]
     [ -150.  -300.  -450.  -900.  -350.  -700. -1050. -2100.   150.   300.
        450.   900.   350.   700.  1050.  2100.]
     [ -300.  -900.  -900. -2700.  -700. -2100. -2100. -6300.   300.   900.
        900.  2700.   700.  2100.  2100.  6300.]]
    >>> print(b)
    [ -198.5  -412.   -459.   -948.   -344.   -710.5  -781.  -1607.    198.5
       412.    459.    948.    344.    710.5   781.   1607. ]
    """
    nimages = len(images)

    if nimages != len(sigmas):
        raise ValueError("Length of sigmas list must match the length of the "
                         "image list.")

    # exclude pixels that have non-positive associated sigmas except the case
    # when all sigmas are non-positive
    for m, s in zip(masks, sigmas):
        ps = (s > 0)
        if not np.all(~ps):
            m &= ps

    # compute squares of sigmas for repeated use later
    sigmas2 = [s**2 for s in sigmas]

    degree1 = tuple([d + 1 for d in degree])

    npolycoeff = 1
    for d in degree1:
        npolycoeff *= d
    sys_eq_array_size = nimages * npolycoeff

    gshape = (nimages,) + degree1

    # pre-compute coordinate arrays:
    coord_arrays, eff_center, coord_system = create_coordinate_arrays(
        images[0].shape,
        center=center,
        image2world=image2world,
        center_cs=center_cs
    )

    # allocate array for the coefficients of the system of equations (a*x=b):
    a = np.zeros((sys_eq_array_size, sys_eq_array_size), dtype=np.float)
    b = np.zeros(sys_eq_array_size, dtype=np.float)

    for i in range(sys_eq_array_size):
        # decompose first (row, or eq) flat index into "original" indices:
        lp = np.unravel_index(i, gshape)
        l = lp[0]
        p = lp[1:]

        # compute known terms:
        for m in range(nimages):
            if m == l:
                continue

            # compute array elements for m!=l:
            b[i] += _image_pixel_sum(
                image_l = images[l],
                image_m = images[m],
                mask_l = masks[l],
                mask_m = masks[m],
                sigma2_l = sigmas2[l],
                sigma2_m = sigmas2[m],
                coord_arrays=coord_arrays,
                p=p
            )

        for j in range(sys_eq_array_size):

            # decompose second (col, or cf) flat index into "original" indices:
            mp = np.unravel_index(j, gshape)
            m = mp[0]
            pp = mp[1:]

            if l == m:  # we will deal with this case in the next iteration
                continue

            a[i, j] = -_sigma_pixel_sum(
                mask_l = masks[l],
                mask_m = masks[m],
                sigma2_l = sigmas2[l],
                sigma2_m = sigmas2[m],
                coord_arrays=coord_arrays,
                p=p,
                pp=pp
            )

    # now compute coefficients of array 'a' for l==m:
    for i in range(sys_eq_array_size):
        # decompose first (row, or eq) flat index into "original" indices:
        lp = np.unravel_index(i, gshape)
        l = lp[0]
        p = lp[1:]

        for ppi in range(npolycoeff):
            pp = np.unravel_index(ppi, degree1)
            j = np.ravel_multi_index((l,) + pp, gshape)

            for m in range(nimages):
                if m == l:
                    continue
                k = np.ravel_multi_index((m,) + pp, gshape)
                a[i, j] -= a[i, k]

    return a, b, coord_arrays, eff_center, coord_system


def pinv_solve(matrix, free_term, nimages, tol=None):
    """
    Solves a system of linear equations

    .. math::
        a \\cdot c = b.

    using Moore-Penrose pseudoinverse.

    Parameters
    ----------
    matrix : numpy.ndarray
        A 2D array containing coefficients of the system.

    free_term : numpy.ndarray
        A 1D array containing free terms of the system of the equations.

    nimages : int
        Number of images for which the system is being solved.

    tol : float, None, optional
        Cutoff for small singular values for Moore-Penrose pseudoinverse.
        When provided, singular values smaller (in modulus) than
        ``tol * |largest_singular_value|`` are set to zero. When
        ``tol`` is `None` (default), cutoff value is determined based on
        the type of the input ``matrix`` argument.

    Returns
    -------
    bkg_poly_coeff : numpy.ndarray
        A 2D `numpy.ndarray` that holds the solution (polynomial coefficients)
        to the system. The solution is grouped by image.

    Examples
    --------
    >>> from jwst.wiimatch.lsq_optimizer import build_lsq_eqs, pinv_solve
    >>> import numpy as np
    >>> im1 = np.zeros((5, 5, 4), dtype=np.float)
    >>> cbg = 1.32 * np.ones_like(im1)
    >>> ind = np.indices(im1.shape, dtype=np.float)
    >>> im3 = cbg + 0.15 * ind[0] + 0.62 * ind[1] + 0.74 * ind[2]
    >>> mask = np.ones_like(im1, dtype=np.int8)
    >>> sigma = np.ones_like(im1, dtype=np.float)
    >>> a, b, _, _, _ = build_lsq_eqs([im1, im3], [mask, mask],
    ... [sigma, sigma], degree=(1,1,1), center=(0,0,0))
    >>> pinv_solve(a, b, 2) # doctest: +FLOAT_CMP
    array([[-6.60000000e-01, -7.50000000e-02, -3.10000000e-01,
             7.10542736e-15, -3.70000000e-01,  8.88178420e-15,
             9.21485110e-15, -2.77555756e-15],
           [ 6.60000000e-01,  7.50000000e-02,  3.10000000e-01,
            -6.43929354e-15,  3.70000000e-01, -7.77156117e-15,
            -9.32587341e-15,  2.99760217e-15]])

    """
    if tol is None:
        tol = np.finfo(matrix.dtype).eps**(2.0/3.0)
    v = np.dot(np.linalg.pinv(matrix, rcond=tol), free_term)
    bkg_poly_coeff = v.reshape((nimages, v.size // nimages))
    return bkg_poly_coeff


def rlu_solve(matrix, free_term, nimages):
    """
    Computes solution of a "reduced" system of linear equations

    .. math::
        a' \\cdot c' = b'.

    using LU-decomposition. If the original system contained a set of
    linearly-dependent equations, then the "reduced" system is formed by
    dropping equations and unknowns related to the first image. The unknowns
    corresponding to the first image initially are assumed to be 0.
    Upon solving the reduced system, these unknowns are recomputed so that
    mean corection coefficients for all images are 0.
    This function uses `~scipy.linalg.lu_solve` and
    `~scipy.linalg.lu_factor` functions.

    Parameters
    ----------
    matrix : numpy.ndarray
        A 2D array containing coefficients of the system.

    free_term : numpy.ndarray
        A 1D array containing free terms of the system of the equations.

    nimages : int
        Number of images for which the system is being solved.

    Returns
    -------
    bkg_poly_coeff : numpy.ndarray
        A 2D `numpy.ndarray` that holds the solution (polynomial coefficients)
        to the system. The solution is grouped by image.

    Examples
    --------
    >>> from jwst.wiimatch.lsq_optimizer import build_lsq_eqs, rlu_solve
    >>> import numpy as np
    >>> im1 = np.zeros((5, 5, 4), dtype=np.float)
    >>> cbg = 1.32 * np.ones_like(im1)
    >>> ind = np.indices(im1.shape, dtype=np.float)
    >>> im3 = cbg + 0.15 * ind[0] + 0.62 * ind[1] + 0.74 * ind[2]
    >>> mask = np.ones_like(im1, dtype=np.int8)
    >>> sigma = np.ones_like(im1, dtype=np.float)
    >>> a, b, _, _, _ = build_lsq_eqs([im1, im3], [mask, mask],
    ... [sigma, sigma], degree=(1, 1, 1), center=(0, 0, 0))
    >>> rlu_solve(a, b, 2)   # doctest: +FLOAT_CMP
    array([[-6.60000000e-01, -7.50000000e-02, -3.10000000e-01,
            -1.19371180e-15, -3.70000000e-01, -1.62003744e-15,
            -1.10844667e-15,  5.11590770e-16],
           [ 6.60000000e-01,  7.50000000e-02,  3.10000000e-01,
             1.19371180e-15,  3.70000000e-01,  1.62003744e-15,
             1.10844667e-15, -5.11590770e-16]])

    """
    drop =  free_term.size // nimages
    if nimages <= 1:
        return np.zeros((1, drop), dtype=np.float)
    from scipy import linalg
    rmat = matrix[drop:, drop:]
    v = linalg.lu_solve(linalg.lu_factor(rmat),
                        free_term[drop:])
    reduced_bkg_poly_coeff = v.reshape((nimages - 1, v.size // (nimages - 1)))
    delta1 = - reduced_bkg_poly_coeff.sum(axis=0) / nimages
    reduced_bkg_poly_coeff += delta1
    bkg_poly_coeff = np.insert(reduced_bkg_poly_coeff, 0, delta1, axis=0)
    return bkg_poly_coeff


def _image_pixel_sum(image_l, image_m, mask_l, mask_m,
                     sigma2_l, sigma2_m, coord_arrays=None, p=None):
    # Compute sum of:
    # coord_arrays^(p) * (image_l - image_m) / (sigma_l**2 + sigma_m**2)
    #
    # If coord_arrays is None, replace it with 1 (this allows code
    # optimization) for the case of constant background (polynomials of zero
    # degree).
    #
    # NOTE: this function does not check that sigma2 arrays have same shapes
    #       as the coord_arrays arrays (for efficiency purpose).

    cmask = np.logical_and(mask_l, mask_m)

    if coord_arrays is None:
        if p is not None:
            raise ValueError("When pixel indices are None then exponent list "
                             "must be None as well.")

        return np.sum((image_l[cmask] - image_m[cmask]) /
                      (sigma2_l[cmask] + sigma2_m[cmask]),
                      dtype=np.float)

    if len(coord_arrays) != len(p):
        raise ValueError("Lengths of the list of pixel index arrays and "
                         "list of the exponents 'p' must be equal.")

    if len(coord_arrays) == 0:
        raise ValueError("There has to be at least one pixel index.")

    i = coord_arrays[0]**p[0]

    for c, ip in zip(coord_arrays[1:], p[1:]):
        i *= c**ip

    return np.sum(i[cmask] * (image_l[cmask] - image_m[cmask]) /
                  (sigma2_l[cmask] + sigma2_m[cmask]),
                  dtype=np.float)


def _sigma_pixel_sum(mask_l, mask_m, sigma2_l, sigma2_m,
                     coord_arrays=None, p=None, pp=None):
    # Compute sum of coord_arrays^(p+pp) ()/ (sigma_l**2 + sigma_m**2)
    #
    # If coord_arrays is None, replace it with 1 (this allows code
    # optimization) for the case of constant background (polynomials of zero
    # degree).
    #
    # NOTE: this function does not check that sigma2 arrays have same shapes
    #       as the coord_arrays arrays (for efficiency purpose).

    cmask = np.logical_and(mask_l, mask_m)

    if coord_arrays is None:
        if p is not None or pp is not None:
            raise ValueError("When pixel indices are None then exponent lists "
                             "must be None as well.")

        return np.sum(1.0 / (sigma2_l[cmask] + sigma2_m[cmask]),
                      dtype=np.float)

    if len(coord_arrays) != len(p) or len(p) != len(pp):
        raise ValueError("Lengths of the list of pixel index arrays and "
                         "lists of the exponents 'p' and 'pp' must be "
                         "equal.")

    if len(coord_arrays) == 0:
        raise ValueError("There has to be at least one pixel index.")

    i = coord_arrays[0]**(p[0] + pp[0])

    for c, ip, ipp in zip(coord_arrays[1:], p[1:], pp[1:]):
        i *= c**(ip + ipp)

    return np.sum(i[cmask] / (sigma2_l[cmask] + sigma2_m[cmask]),
                  dtype=np.float)
