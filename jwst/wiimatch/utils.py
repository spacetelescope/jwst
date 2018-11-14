"""
This module provides utility functions for use by :py:mod:`wiimatch` module.

:Author: Mihai Cara (contact: help@stsci.edu)

"""

import numpy as np


__all__ = ['create_coordinate_arrays']


def create_coordinate_arrays(image_shape, center=None, image2world=None,
                             center_cs='image'):
    """
    Create a list of coordinate arrays/grids for each dimension in the image
    shape. This function is similar to `numpy.indices` except it returns the
    list of arrays in reversed order. In addition, it can center image
    coordinates to a provided ``center`` and also convert image coordinates to
    world coordinates using provided ``image2world`` function.

    Parameters
    ----------
    image_shape : sequence of int
        The shape of the image/grid.

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

    Returns
    -------
    coord_arrays : list
        A list of `numpy.ndarray` coordinate arrays each of ``image_shape``
        shape.

    eff_center : tuple
        A tuple of coordinates of the effective center as used in generating
        coordinate arrays.

    coord_system : {'image', 'world'}
        Coordinate system of the coordinate arrays and returned ``center``
        value.

    Examples
    --------
    >>> create_coordinate_arrays((3,5,4))
        ((array([[[-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.]],
         [[-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.]],
         [[-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.],
          [-1.,  0.,  1.,  2.]]]),
    array([[[-2., -2., -2., -2.],
          [-1., -1., -1., -1.],
          [ 0.,  0.,  0.,  0.],
          [ 1.,  1.,  1.,  1.],
          [ 2.,  2.,  2.,  2.]],
         [[-2., -2., -2., -2.],
          [-1., -1., -1., -1.],
          [ 0.,  0.,  0.,  0.],
          [ 1.,  1.,  1.,  1.],
          [ 2.,  2.,  2.,  2.]],
         [[-2., -2., -2., -2.],
          [-1., -1., -1., -1.],
          [ 0.,  0.,  0.,  0.],
          [ 1.,  1.,  1.,  1.],
          [ 2.,  2.,  2.,  2.]]]),
    array([[[-2., -2., -2., -2.],
          [-2., -2., -2., -2.],
          [-2., -2., -2., -2.],
          [-2., -2., -2., -2.],
          [-2., -2., -2., -2.]],
         [[-1., -1., -1., -1.],
          [-1., -1., -1., -1.],
          [-1., -1., -1., -1.],
          [-1., -1., -1., -1.],
          [-1., -1., -1., -1.]],
         [[ 0.,  0.,  0.,  0.],
          [ 0.,  0.,  0.,  0.],
          [ 0.,  0.,  0.,  0.],
          [ 0.,  0.,  0.,  0.],
          [ 0.,  0.,  0.,  0.]]])), (1.0, 2.0, 2.0), u'image')

    """
    if center_cs not in ['image', 'world']:
        raise ValueError("Parameter 'center_cs' must take one of the "
                         "following two values: 'image' or 'world'.")

    if center is None:
        # set the center at the center of the image array:
        center = tuple([float(i//2) for i in image_shape])
        center_cs = 'image'

    else:
        if len(center) != len(image_shape):
            raise ValueError("Number of coordinates of the 'center' must "
                             "match the dimentionality of the image.")

        if center_cs == 'world' and image2world is None:
            raise ValueError("'center_cs' cannot be 'world' when 'image2world'"
                             " is not defined.")

    ind = np.indices(image_shape, dtype=np.float)[::-1]

    if image2world is None:
        coord_system = 'image'
        eff_center = tuple([c for c in center])

    else:
        if center_cs == 'world':
            eff_center = tuple(map(float, center))

        else:
            # convert image's center from image to world coordinates:
            eff_center = tuple(map(float, image2world(*center)))

        # convert pixel indices to world coordinates:
        # TODO: get rid of the ravel/reshape dancing once
        #       issue https://github.com/spacetelescope/gwcs/issues/89
        #       is fixed.
        ind = image2world(*[x.ravel() for x in ind])
        ind = [x.reshape(image_shape) for x in ind]
        coord_system = 'world'

    coord_arrays = tuple([i - c for (i, c) in zip(ind, eff_center)])

    return coord_arrays, eff_center, coord_system
