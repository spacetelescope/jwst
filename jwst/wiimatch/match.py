"""
A module that provides main API for optimal (LSQ) "matching" of weighted
N-dimensional image intensity data using (multivariate) polynomials.

:Author: Mihai Cara (contact: help@stsci.edu)


"""

# STDLIB
import numpy as np

from .lsq_optimizer import build_lsq_eqs, pinv_solve, rlu_solve


__all__ = ['match_lsq']

SUPPORTED_SOLVERS = ['RLU', 'PINV']


def match_lsq(images, masks=None, sigmas=None, degree=0,
              center=None, image2world=None, center_cs='image',
              ext_return=False, solver='RLU'):
    """
    Compute coefficients of (multivariate) polynomials that once subtracted
    from input images would provide image intensity matching in the least
    squares sense.

    Parameters
    ----------
    images : list of numpy.ndarray
        A list of 1D, 2D, etc. `numpy.ndarray` data array whose "intensities"
        must be "matched". All arrays must have identical shapes.

    masks : list of numpy.ndarray, None
        A list of `numpy.ndarray` arrays of same length as ``images``.
        Non-zero mask elements indicate valid data in the corresponding
        ``images`` array. Mask arrays must have identical shape to that of
        the arrays in input ``images``. Default value of `None` indicates that
        all pixels in input images are valid.

    sigmas : list of numpy.ndarray, None
        A list of `numpy.ndarray` data array of same length as ``images``
        representing the uncertainties of the data in the corresponding array
        in ``images``. Uncertainty arrays must have identical shape to that of
        the arrays in input ``images``. The default value of `None` indicates
        that all pixels will be assigned equal weights.

    degree : iterable, int
        A list of polynomial degrees for each dimension of data arrays in
        ``images``. The length of the input list must match the dimensionality
        of the input images. When a single integer number is provided, it is
        assumed that the polynomial degree in each dimension is equal to
        that integer.

    center : iterable, None, optional
        An iterable of length equal to the number of dimensions in
        ``image_shape`` that indicates the center of the coordinate system
        in **image** coordinates when ``center_cs`` is ``'image'`` otherwise
        center is assumed to be in **world** coordinates (when ``center_cs``
        is ``'world'``). When ``center`` is `None` then ``center`` is
        set to the middle of the "image" as ``center[i]=image_shape[i]//2``.
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

    ext_return : bool, optional
        Indicates whether this function should return additional values besides
        optimal polynomial coefficients (see ``bkg_poly_coeff`` return value
        below) that match image intensities in the LSQ sense. See **Returns**
        section for more details.

    solver : {'RLU', 'PINV'}, optional
        Specifies method for solving the system of equations.

    Returns
    -------
    bkg_poly_coeff : numpy.ndarray
        When ``nimages`` is `None`, this function returns a 1D `numpy.ndarray`
        that holds the solution (polynomial coefficients) to the system.

        When ``nimages`` is **not** `None`, this function returns a 2D
        `numpy.ndarray` that holds the solution (polynomial coefficients)
        to the system. The solution is grouped by image.

    a : numpy.ndarray
        A 2D `numpy.ndarray` that holds the coefficients of the linear system
        of equations. This value is returned only when ``ext_return``
        is `True`.

    b : numpy.ndarray
        A 1D `numpy.ndarray` that holds the free terms of the linear system of
        equations. This value is returned only when ``ext_return`` is `True`.

    coord_arrays : list
        A list of `numpy.ndarray` coordinate arrays each of ``image_shape``
        shape. This value is returned only when ``ext_return`` is `True`.

    eff_center : tuple
        A tuple of coordinates of the effective center as used in generating
        coordinate arrays. This value is returned only when ``ext_return``
        is `True`.

    coord_system : {'image', 'world'}
        Coordinate system of the coordinate arrays and returned ``center``
        value. This value is returned only when ``ext_return`` is `True`.

    Notes
    -----
    :py:func:`match_lsq` builds a system of linear equations

    .. math::
        a \\cdot c = b

    whose solution :math:`c` is a set of coefficients of (multivariate)
    polynomials that represent the "background" in each input image (these are
    polynomials that are "corrections" to intensities of input images) such
    that the following sum is minimized:

    .. math::
        L = \\sum^N_{n,m=1,n \\neq m} \\sum_k\
    \\frac{\\left[I_n(k) - I_m(k) - P_n(k) + P_m(k)\\right]^2}\
    {\\sigma^2_n(k) + \\sigma^2_m(k)}.

    In the above equation, index :math:`k=(k_1,k_2,...)` labels a position
    in input image's pixel grid [NOTE: all input images share a common
    pixel grid].

    "Background" polynomials :math:`P_n(k)` are defined through the
    corresponding coefficients as:

    .. math::
        P_n(k_1,k_2,...) = \\sum_{d_1=0,d_2=0,...}^{D_1,D_2,...} \
c_{d_1,d_2,...}^n \\cdot k_1^{d_1} \\cdot k_2^{d_2}  \\cdot \\ldots .

    Coefficients :math:`c_{d_1,d_2,...}^n` are arranged in the vector :math:`c`
    in the following order:

    .. math::
        (c_{0,0,\\ldots}^1,c_{1,0,\\ldots}^1,\\ldots,c_{0,0,\\ldots}^2,\
c_{1,0,\\ldots}^2,\\ldots).

    :py:func:`match_lsq` returns coefficients of the polynomials that
    minimize *L*.

    Examples
    --------
    >>> import numpy as np
    >>> im1 = np.zeros((5, 5, 4), dtype=np.float)
    >>> cbg = 1.32 * np.ones_like(im1)
    >>> ind = np.indices(im1.shape, dtype=np.float)
    >>> im3 = cbg + 0.15 * ind[0] + 0.62 * ind[1] + 0.74 * ind[2]
    >>> mask = np.ones_like(im1, dtype=np.int8)
    >>> sigma = np.ones_like(im1, dtype=np.float)
    >>> match_lsq([im1, im3], [mask, mask], [sigma, sigma],   #doctest: +FLOAT_CMP
    ... degree=(1, 1, 1), center=(0, 0, 0))
    array([[-6.60000000e-01, -7.50000000e-02, -3.10000000e-01,
            -1.19371180e-15, -3.70000000e-01, -1.62003744e-15,
            -1.10844667e-15,  5.11590770e-16],
           [ 6.60000000e-01,  7.50000000e-02,  3.10000000e-01,
             1.19371180e-15,  3.70000000e-01,  1.62003744e-15,
             1.10844667e-15, -5.11590770e-16]])
    """
    solver = solver.upper()
    if solver not in SUPPORTED_SOLVERS:
        ns = len(SUPPORTED_SOLVERS)
        raise ValueError("'solver' must be one of the supported solvers: '{}'"
                         .format(SUPPORTED_SOLVERS[0] if ns == 1 else
                                 '\', \''.join(SUPPORTED_SOLVERS[:-1]) +
                                 '\'' + (',' if ns > 2 else '') +
                                 ' or \'{}'.format(SUPPORTED_SOLVERS[-1])))

    # check that all images have the same shape:
    shapes = set([])
    for im in images:
        shapes.add(im.shape)
    if len(shapes) > 1:
        raise ValueError("All images must have identical shapes.")

    nimages = len(images)
    ndim = len(images[0].shape)

    # check that the number of good pixel mask arrays matches the numbers
    # of input images, and if 'masks' is None - set all of them to True:
    if masks is None:
        masks = [np.ones_like(images[0], dtype=np.bool) for i in images]

    else:
        if len(masks) != nimages:
            raise ValueError("Length of masks list must match the length of "
                             "the image list.")

        for m in masks:
            shapes.add(m.shape)
        if len(shapes) > 1:
            raise ValueError("Shape of each mask array must match the shape "
                             "of input images.")

        # make a copy of the masks since we might modify these masks later
        masks = [m.copy() for m in masks]

    # check that the number of sigma arrays matches the numbers
    # of input images, and if 'sigmas' is None - set all of them to 1:
    if sigmas is None:
        sigmas = [np.ones_like(images[0], dtype=np.float) for i in images]

    else:
        if len(sigmas) != nimages:
            raise ValueError("Length of sigmas list must match the length of "
                             "the image list.")

        for s in sigmas:
            shapes.add(s.shape)
        if len(shapes) > 1:
            raise ValueError("Shape of each sigma array must match the shape "
                             "of input images.")

    # check that 'degree' has the same length as the number of dimensions
    # in image arrays:
    if hasattr(degree, '__iter__'):
        if len(degree) != ndim:
            raise ValueError("The length of 'degree' parameter must match "
                             "the number of image dimensions.")
        degree = tuple([int(d) for d in degree])

    else:
        intdeg = int(degree)
        degree = tuple([intdeg for i in range(ndim)])

    # check that 'center' has the same length as the number of dimensions
    # in image arrays:
    if hasattr(center, '__iter__'):
        if len(center) != ndim:
            raise ValueError("The length of 'center' parameter must match "
                             "the number of image dimensions.")

    elif center is not None:
        center = tuple([center for i in range(ndim)])

    # build the system of equations:
    a, b, coord_arrays, eff_center, coord_system = build_lsq_eqs(
        images,
        masks,
        sigmas,
        degree,
        center=center,
        image2world=image2world,
        center_cs=center_cs
    )

    # solve the system:
    if solver == 'RLU':
        bkg_poly_coef = rlu_solve(a, b, nimages)
    else:
        tol = np.finfo(images[0].dtype).eps**(2.0/3.0)
        bkg_poly_coef = pinv_solve(a, b, nimages, tol)

    if ext_return:
        return bkg_poly_coef, a, b, coord_arrays, eff_center, coord_system
    else:
        return bkg_poly_coef
