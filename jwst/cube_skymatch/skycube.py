"""
This module provides support for performing sky matching on cubes.
on the sky. Primary use case would use the following
generalized step: ``cube_skymatch``.

:Authors: Mihai Cara (contact: help@stsci.edu)

"""

import numpy as np


__all__ = ['SkyCube']


class SkyCube():
    r"""
    Container that holds information about properties of a *single*
    image such as:

    * image data;
    * data weights;
    * mask associated image data indicating "good" (True) or "bad" (False)
      data;
    * id;
    * sky background value/coefficients;
    * other useful properties.

    This class also provides methods for computing overlaps between two cubes
    using masks created from weight maps and also methods for computing
    background difference between two cubes. Background is fitted to the
    difference between two cubes using a 3D polynomial form:

    .. math::

        B(x,y,z)=\sum_{k,l,m=0}^{Dx,Dy,Dz}a_{k,l,m}(x-x_0)^k(y-y_0)^l(z-z_0)^m

    When a WCS object is provided, the background difference polynomial will
    be computed as a function of world coordinates
    (:math:`\alpha`, :math:`\delta`, :math:`\lambda`) instead of the
    image/cube coordinates (x, y, z).
    """

    def __init__(self, data, wcs=None, wcsinfo=None,
                 weights=None, cube_weight=1.0,
                 bkg_deg=0, bkg_center=None,
                 id=None, meta=None):
        r""" Initializes the SkyCube object.

        Parameters
        ----------
        data : numpy.ndarray
            A 3D array of data (z, y, x).

        wcs : gwcs, None, optional
            When a WCS is provided, the background polynomial will be
            computed as a function of world coordinates. In this case,
            `bkg_center` is expected to be provided in world coordinates.

        wcsinfo : dict-like, optional
            Meta WCS info. When ``wcsnfo`` is provided, then ``wcs``

        weights : numpy.ndarray, None, optional
            A 3D array that indicates the weight of individual voxels in the
            input 3D ``data``. A weight of 0 indicates that the corresponding
            voxel data is invalid. The default value of ``None`` sets the
            weight of all voxels in the input ``data`` to 1.

        cube_weight : float, optional
            Weight of the cube. This weight is used when combining cubes to
            form larger mosaics. Also see `combine_with_other`.

        bkg_deg : int, array-like of shape (3,), optional
            Degree of the polynomial used to approximate the background level.
            Degree of 0 corresponds to a constant background and degree of 1
            corresponds to a tri-linear form. It is also possible to specify
            the maximum degree for each coordinate by providing an iterable
            of length 3 such as (Dx, Dy, Dz).

        bkg_center : array-like of shape (3,), None, optional
            A tuple, list, etc. containing 3 coordinates
            (:math:`z_0`, :math:`y_0`, :math:`x_0`) (or
            (:math:`\alpha_0`, :math:`\delta_0`, :math:`\lambda_0`) when `wcs`
            is not `None`) indicating the center with regard to which
            tri-linear polynomial is defined.

            When the default value is used (`None`) and `wcs` is `None`,
            `bkg_center` is computed as the center of the ``data`` :

            .. math::
                (x_0, y_0, z_0) = ([N_x/2], [N_y/2], [N_z/2])

            where (:math:`N_x`, :math:`N_y`, :math:`N_z`) is the shape of the
            input ``data`` (in reversed order).

            When the default value is used (`None`) and a WCS object is
            provided, `bkg_center` is taken to be the ``CRVAL`` (world
            coordinate of the reference pixel).

        id : anything
            The value of this parameter is simple stored within the `SkyCube`
            object. While it can be of any type, it is preferable that `id` be
            of a type with nice string representation.

        meta : dict, None, optional
            A dictionary of various items to be stored within the `SkyCube`
            object.

        """
        # check 'weights' and 'data' dimensions:
        if len(data.shape) != 3:
            raise ValueError("Input 'data' array must be a 3D cube.")
        if weights is not None and weights.shape != data.shape:
            raise ValueError("'weights' must have same shape as 'data' array.")

        wcsinfo_is_None = (wcsinfo is None or wcsinfo.crval1 is None or
                           wcsinfo.crval2 is None or wcsinfo.crval3 is None)

        if wcsinfo_is_None and wcs is not None:
            raise ValueError("'wcsinfo' cannot be None or have its 'crvaln' "
                             "set to None when 'wcs' is not None.")

        self._data = data
        self.weights = weights

        self.cube_weight = float(cube_weight)
        self.meta = meta
        self._id = id
        self._wcs = wcs
        self._wcsinfo = wcsinfo

        # set background degree:
        if hasattr(bkg_deg, '__iter__'):
            if len(bkg_deg) != 3:
                raise ValueError("When 'bkg_deg' is an iterable it must have "
                                 "three elements.")
            self._bkg_degree = tuple(map(int, bkg_deg))
        else:
            self._bkg_degree = 3 * (int(bkg_deg), )

        # set center of the background polynomials:
        self._set_bkg_center(bkg_center)

        bkg_degree_p1 = tuple((i + 1 for i in self._bkg_degree))
        self._bkg_coeff = np.zeros(bkg_degree_p1, dtype=float)
        self._bkg_status = 1  # not computed
        self._bkg_cube = None
        self._bkg_cube_dirty = True

        # Coordinate mesh grid:
        self._mgx = None
        self._mgy = None
        self._mgz = None

    def _set_bkg_center(self, refpt):
        if refpt is None:
            if self.wcs is None:
                nz, ny, nx = self._data.shape
                self._x0 = int(nx / 2)
                self._y0 = int(ny / 2)
                self._z0 = int(nz / 2)
            else:
                self._x0 = self._wcsinfo.crval1
                self._y0 = self._wcsinfo.crval2
                self._z0 = self._wcsinfo.crval3

        else:
            self._x0, self._y0, self._z0 = refpt

    @property
    def bkg_center(self):
        return self._x0, self._y0, self._z0

    @property
    def weights(self):
        """ Set or get ``weights`` cube. """
        return self._weights

    @weights.setter
    def weights(self, weights):
        if weights is None:
            self._weights = np.ones_like(self.data, dtype=float)
            self._mask = np.ones_like(self.data, dtype=bool)
        else:
            if weights.shape != self.data.shape:
                raise ValueError("Weight map must have the same shape as "
                                 "the image cube")
            self._weights = weights
            self._mask = weights > 0.0

    @property
    def mask(self):
        """ Get mask of valid voxels computed from weights cube. """
        return self._mask

    @property
    def wcs(self):
        """ Get WCS. """
        return self._wcs

    @property
    def data(self):
        """ Get cube data.
        """
        return self._data

    @property
    def cube_weight(self):
        """ Set or get cube's weight.
        """
        return self._cube_weight

    @cube_weight.setter
    def cube_weight(self, cube_weight):
        cube_weight = float(cube_weight)
        if cube_weight < 0.0:
            raise ValueError("Cube's weight must be a non-negative number.")
        self._cube_weight = cube_weight

    @property
    def bkg_coeff(self):
        """ Get background polynomial coefficients [].
        Returns a 3D ``numpy.ndarray``.
        """
        return self._bkg_coeff

    @property
    def bkg_degree(self):
        """ Get background polynomial degree.
        """
        return self._bkg_degree

    @property
    def bkg_status(self):
        """
        Retrieve status of the background computations:
        0 - OK
        1 - not computed
        2 - fit failed

        """
        return self._bkg_status

    @property
    def id(self):
        """
        Set or get `SkyCube`'s `id`.

        While `id` can be of any type, it is preferable that `id` be
        of a type with nice string representation.

        """
        return self._id

    @id.setter
    def id(self, id):
        self._id = id

    @property
    def skysub_data(self):
        """ Retrieve sky background-subtracted data. """
        if self._bkg_status > 0:
            return self._data

        return self._data - self.bkg_cube

    def subtract_sky(self):
        """ Subtract computed sky from cube data. """
        if self._bkg_status > 0:
            return self._data

        self._data -= self.bkg_cube

    @property
    def bkg_cube(self, mx=None, my=None, mz=None):
        self.calc_bkg_cube(mx, my, mz)
        return self._bkg_cube

    def _calc_meshgrid(self):
        # This function should be called ONLY AFTER _data, _wcs, _x0, _y0, _z0,
        # _bkg_degree

        z, y, x = np.indices(self._data.shape, dtype=float)
        if self._wcs is not None:
            # TODO: the need to use ravel/reshape is due to a bug in
            # astropy.modeling that is affecting gwcs. In future, all the code
            # in this "if"-block should be replaced directly with:
            # x, y, z = self.wcs(x, y, z)
            shape = z.shape
            x, y, z = self.wcs(x.ravel(), y.ravel(), z.ravel())
            x = x.reshape(shape)
            y = y.reshape(shape)
            z = z.reshape(shape)
        self._mgx = x - self._x0
        self._mgy = y - self._y0
        self._mgz = z - self._z0

        # compute pseudo-Vandermonde matrix:
        v = np.polynomial.polynomial.polyvander3d(
            self._mgx,
            self._mgy,
            self._mgz,
            self._bkg_degree
        )
        self._vander = v.reshape((-1, v.shape[-1]))

    def free_intermediate_arrays(self):
        """ Releases references to internal intermediate arrays.
        """
        self._mgx = None
        self._mgy = None
        self._mgz = None
        self._vander = None
        self._bkg_cube = None
        self._bkg_cube_dirty = True

    def calc_bkg_cube(self, mx=None, my=None, mz=None):
        if self._bkg_status > 0 or not self._bkg_cube_dirty:
            return

        if mx is None and my is None and mz is None:
            if self._mgx is None:
                self._calc_meshgrid()
            mx = self._mgx
            my = self._mgy
            mz = self._mgz

        elif not all([mx is not None, my is not None, mz is not None]):
            raise ValueError("All mesh grids must simultaneously be either "
                             "None or not None.")

        self._bkg_cube = np.polynomial.polynomial.polyval3d(
            mx,
            my,
            mz,
            self._bkg_coeff
        )

        self._bkg_cube_dirty = False

    def overlap_measure(self, other_cube):
        """
        Get a measure of the overlap between two cubes.

        Currently, this function returns the number of common valid
        (``weights``>0) pixels in both cubes. If ``weights`` is ``None``,
        all pixels of that cube are considered valid.

        """
        return np.count_nonzero(np.logical_and(self.mask, other_cube.mask))

    def fit_background(self, other_cube):
        if self._mgx is None:
            self._calc_meshgrid()

        # find common "good" voxels:
        mask = np.logical_and(self.mask, other_cube.mask)

        # compute difference between this cube and the other cube:
        diff = self._data[mask] - other_cube.data[mask]

        # fit polynomial to the difference
        vectd = diff.ravel()
        c = np.linalg.lstsq(self._vander[mask.ravel()], vectd)[0]
        self._bkg_coeff[:, :, :] = c.reshape(self._bkg_coeff.shape)
        self._bkg_cube_dirty = True
        self._bkg_status = 0

    @property
    def skystat(self):
        """ Stores/retrieves a callable object that takes a either a 2D image
        (2D `numpy.ndarray`) or a list of pixel values (a Nx1 array) and
        returns a tuple of two values: some statistics
        (e.g., mean, median, etc.) and number of pixels/values from the input
        image used in computing that statistics.

        When `skystat` is not set, `SkyCube` will use
        :py:class:`~jwst_pipeline.skymatch.skystatistics.SkyStats` object
        to perform sky statistics on image data.

        """
        return self._skystat

    @skystat.setter
    def skystat(self, skystat):
        self._skystat = skystat

    def combine_with_other(self, other):
        """
        Combine this cube with another cube to create a "mosaic". This cube
        will be updated with data from the other cube.

        """
        m = np.logical_or(self.mask, other.mask)
        w1 = self.cube_weight * self.weights[m]
        w2 = other.cube_weight * other.weights[m]
        tw = w1 + w2
        self._data[m] = (w1 * self.data[m] + w2 * other.skysub_data[m]) / tw
        self.weights[m] = tw
