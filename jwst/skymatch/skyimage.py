"""
The ``skyimage`` module contains algorithms that are used by
``skymatch`` to manage all of the information for footprints (image outlines)
on the sky as well as perform useful operations on these outlines such as
computing intersections and statistics in the overlap regions.

:Authors: Mihai Cara (contact: help@stsci.edu)


"""

# STDLIB
import numpy as np

# THIRD-PARTY
from spherical_geometry.polygon import SphericalPolygon

#LOCAL
from . skystatistics import SkyStats
from . import region


__all__ = ['SkyImage', 'SkyGroup']


class SkyImage:
    """
    Container that holds information about properties of a *single*
    image such as:

    * image data;
    * WCS of the chip image;
    * bounding spherical polygon;
    * id;
    * pixel area;
    * sky background value;
    * sky statistics parameters;
    * mask associated image data indicating "good" (1) data.

    """


    def __init__(self, image, wcs_fwd, wcs_inv, pix_area=1.0, convf=1.0,
                 mask=None, id=None, skystat=None, stepsize=None, meta=None):
        """ Initializes the SkyImage object.

        Parameters
        ----------
        image : numpy.ndarray
            A 2D array of image data.

        wcs_fwd : function
            "forward" pixel-to-world transformation function.

        wcs_inv : function
            "inverse" world-to-pixel transformation function.

        pix_area : float, optional
            Average pixel's sky area.

        convf : float, optional
            Conversion factor that when multiplied to `image` data converts
            the data to "uniform" (across multiple images) surface
            brightness units.

            .. note::

              The functionality to support this conversion is not yet
              implemented and at this moment `convf` is ignored.

        mask : numpy.ndarray
            A 2D array that indicates
            what pixels in the input `image` should be used for sky
            computations (``1``) and which pixels should **not** be used
            for sky computations (``0``).

        id : anything
            The value of this parameter is simple stored within the `SkyImage`
            object. While it can be of any type, it is prefereble that `id` be
            of a type with nice string representation.

        skystat : callable, None, optional
            A callable object that takes a either a 2D image (2D
            `numpy.ndarray`) or a list of pixel values (a Nx1 array) and
            returns a tuple of two values: some statistics (e.g., mean,
            median, etc.) and number of pixels/values from the input image
            used in computing that statistics.

            When `skystat` is not set, `SkyImage` will use
            :py:class:`~jwst_pipeline.skymatch.skystatistics.SkyStats` object
            to perform sky statistics on image data.

        stepsize : int, None, optional
            Spacing between vertices of the image's bounding polygon. Default
            value of `None` creates bounding polygons with four vertices
            corresponding to the corners of the image.

        meta : dict, None, optional
            A dictionary of various items to be stored within the `SkyImage`
            object.

        """
        self.image = image
        self.convf = convf
        self.meta = meta
        self._id = id
        self._pix_area = pix_area

        # WCS
        self.wcs_fwd = wcs_fwd
        self.wcs_inv = wcs_inv

        # initial sky value:
        self._sky = 0.0
        self._sky_is_valid = False

        # check that mask has the same shape as image:
        if mask is None:
            self.mask = None

        else:
            if image is None:
                raise ValueError("'mask' must be None when 'image' is None")

            self.mask = np.asanyarray(mask, dtype=np.bool)

            if self.mask.shape != image.shape:
                raise ValueError("'mask' must have the same shape as 'image'.")

        # create spherical polygon bounding the image
        if image is None or wcs_fwd is None or wcs_inv is None:
            self._radec = [(np.array([]), np.array([]))]
            self._polygon = SphericalPolygon([])
            self._poly_area = 0.0

        else:
            self.calc_bounding_polygon(stepsize)

        # set sky statistics function (NOTE: it must return statistics and
        # the number of pixels used after clipping)
        if skystat is None:
            self.set_builtin_skystat()
        else:
            self.skystat = skystat

    @property
    def id(self):
        """ Set or get `SkyImage`'s `id`.

        While `id` can be of any type, it is prefereble that `id` be
        of a type with nice string representation.

        """
        return self._id

    @id.setter
    def id(self, id):
        self._id = id

    @property
    def pix_area(self):
        """ Set or get mean pixel area.
        """
        return self._pix_area

    @pix_area.setter
    def pix_area(self, pix_area):
        self._pix_area = pix_area

    @property
    def poly_area(self):
        """ Get bounding polygon area in srad units.
        """
        return self._poly_area

    @property
    def sky(self):
        """ Sky background value. See `calc_sky` for more details.
        """
        return self._sky

    @sky.setter
    def sky(self, sky):
        self._sky = sky

    @property
    def is_sky_valid(self):
        """
        Indicates whether sky value was successfully computed.
        Must be set externally.
        """
        return self._sky_is_valid

    @is_sky_valid.setter
    def is_sky_valid(self, valid):
        self._sky_is_valid = valid

    @property
    def radec(self):
        """
        Get RA and DEC of the verteces of the bounding polygon as a
        `~numpy.ndarray` of shape (N, 2) where N is the number of verteces + 1.
        """
        return self._radec

    @property
    def polygon(self):
        """ Get image's bounding polygon.
        """
        return self._polygon

    def intersection(self, skyimage):
        """
        Compute intersection of this `SkyImage` object and another
        `SkyImage`, `SkyGroup`, or
        :py:class:`~spherical_geometry.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        skyimage : SkyImage, SkyGroup, SphericalPolygon
            Another object that should be intersected with this `SkyImage`.

        Returns
        -------
        polygon : SphericalPolygon
            A :py:class:`~spherical_geometry.polygon.SphericalPolygon` that is
            the intersection of this `SkyImage` and `skyimage`.

        """
        if isinstance(skyimage, (SkyImage, SkyGroup)):
            return self._polygon.intersection(skyimage.polygon)
        else:
            return self._polygon.intersection(skyimage)

    def calc_bounding_polygon(self, stepsize=None):
        """ Compute image's bounding polygon.

        Parameters
        ----------
        stepsize : int, None, optional
            Indicates the maximum separation between two adjacent vertices
            of the bounding polygon along each side of the image. Corners
            of the image are included automatically. If `stepsize` is `None`,
            bounding polygon will contain only vertices of the image.

        """
        ny, nx = self.image.shape

        if stepsize is None:
            nintx = 2
            ninty = 2
        else:
            nintx = max(2, int(np.ceil((nx + 1.0) / stepsize)))
            ninty = max(2, int(np.ceil((ny + 1.0) / stepsize)))

        xs = np.linspace(-0.5, nx - 0.5, nintx, dtype=np.float)
        ys = np.linspace(-0.5, ny - 0.5, ninty, dtype=np.float)[1:-1]
        nptx = xs.size
        npty = ys.size

        npts = 2 * (nptx + npty)

        borderx = np.empty((npts + 1,), dtype=np.float)
        bordery = np.empty((npts + 1,), dtype=np.float)

        # "bottom" points:
        borderx[:nptx] = xs
        bordery[:nptx] = -0.5
        # "right"
        sl = np.s_[nptx:nptx + npty]
        borderx[sl] = nx - 0.5
        bordery[sl] = ys
        # "top"
        sl = np.s_[nptx + npty:2 * nptx + npty]
        borderx[sl] = xs[::-1]
        bordery[sl] = ny - 0.5
        # "left"
        sl = np.s_[2 * nptx + npty:-1]
        borderx[sl] = -0.5
        bordery[sl] = ys[::-1]

        # close polygon:
        borderx[-1] = borderx[0]
        bordery[-1] = bordery[0]

        ra, dec = self.wcs_fwd(borderx, bordery, with_bounding_box=False)
        # TODO: for strange reasons, occasionally ra[0] != ra[-1] and/or
        #       dec[0] != dec[-1] (even though we close the polygon in the
        #       previous two lines). Then SphericalPolygon fails because
        #       points are not closed. Threfore we force it to be closed:
        ra[-1] = ra[0]
        dec[-1] = dec[0]

        self._radec = [(ra, dec)]
        self._polygon = SphericalPolygon.from_radec(ra, dec)
        self._poly_area = np.fabs(self._polygon.area())

    @property
    def skystat(self):
        """ Stores/retrieves a callable object that takes a either a 2D image
        (2D `numpy.ndarray`) or a list of pixel values (a Nx1 array) and
        returns a tuple of two values: some statistics
        (e.g., mean, median, etc.) and number of pixels/values from the input
        image used in computing that statistics.

        When `skystat` is not set, `SkyImage` will use
        :py:class:`~jwst_pipeline.skymatch.skystatistics.SkyStats` object
        to perform sky statistics on image data.

        """
        return self._skystat

    @skystat.setter
    def skystat(self, skystat):
        self._skystat = skystat

    def set_builtin_skystat(self, skystat='median', lower=None, upper=None,
                            nclip=5, lsigma=4.0, usigma=4.0, binwidth=0.1):
        """
        Replace already set `skystat` with a "built-in" version of a
        statistics callable object used to measure sky background.

        See :py:class:`~jwst_pipeline.skymatch.skystatistics.SkyStats` for the
        parameter description.

        """
        self._skystat = SkyStats(
            skystat=skystat,
            lower=lower,
            upper=upper,
            nclip=nclip,
            lsig=lsigma,
            usig=usigma,
            binwidth=binwidth
        )

    # TODO: due to a bug in the sphere package, see
    #       https://github.com/spacetelescope/sphere/issues/74
    #       intersections with polygons formed as union does not work.
    #       For this reason I re-implement 'calc_sky' below with
    #       a workaround for the bug.
    #       The original implementation (now called ``_calc_sky_orig``
    #       should replace current 'calc_sky' once the bug is fixed.
    #
    def calc_sky(self, overlap=None, delta=True):
        """
        Compute sky background value.

        Parameters
        ----------
        overlap : SkyImage, SkyGroup, SphericalPolygon, list of tuples, \
None, optional
            Another `SkyImage`, `SkyGroup`,
            :py:class:`spherical_geometry.polygons.SphericalPolygon`, or
            a list of tuples of (RA, DEC) of vertices of a spherical
            polygon. This parameter is used to indicate that sky statistics
            should computed only in the region of intersection of *this*
            image with the polygon indicated by `overlap`. When `overlap` is
            `None`, sky statistics will be computed over the entire image.

        delta : bool, optional
            Should this function return absolute sky value or the difference
            between the computed value and the value of the sky stored in the
            `sky` property.

        Returns
        -------
        skyval : float, None
            Computed sky value (absolute or relative to the `sky` attribute).
            If there are no valid data to perform this computations (e.g.,
            because this image does not overlap with the image indicated by
            `overlap`), `skyval` will be set to `None`.

        npix : int
            Number of pixels used to compute sky statistics.

        polyarea : float
            Area (in srad) of the polygon that bounds data used to compute
            sky statistics.

        """
        if overlap is None:

            if self.mask is None:
                data = self.image
            else:
                data = self.image[self.mask]

            polyarea = self.poly_area

        else:
            fill_mask = np.zeros(self.image.shape, dtype=bool)

            if isinstance(overlap, SkyImage):
                intersection = self.intersection(overlap)
                polyarea = np.fabs(intersection.area())
                radec = list(intersection.to_radec())

            elif isinstance(overlap, SkyGroup):
                radec = []
                polyarea = 0.0
                for im in overlap:
                    intersection = self.intersection(im)
                    polyarea1 = np.fabs(intersection.area())
                    if polyarea1 == 0.0:
                        continue
                    polyarea += polyarea1
                    radec += list(intersection.to_radec())

            elif isinstance(overlap, SphericalPolygon):
                radec = []
                polyarea = 0.0
                for p in overlap._polygons:
                    intersection = self.intersection(SphericalPolygon([p]))
                    polyarea1 = np.fabs(intersection.area())
                    if polyarea1 == 0.0:
                        continue
                    polyarea += polyarea1
                    radec += list(intersection.to_radec())

            else: # assume a list of (ra, dec) tuples:
                radec = []
                polyarea = 0.0
                for r, d in overlap:
                    poly = SphericalPolygon.from_radec(r, d)
                    polyarea1 = np.fabs(poly.area())
                    if polyarea1 == 0.0 or len(r) < 4:
                        continue
                    polyarea += polyarea1
                    radec.append(self.intersection(poly).to_radec())

            if polyarea == 0.0:
                return (None, 0, 0.0)

            for ra, dec in radec:
                if len(ra) < 4:
                    continue

                # set pixels in 'fill_mask' that are inside a polygon to True:
                x, y = self.wcs_inv(ra, dec)
                poly_vert = list(zip(*[x, y]))

                polygon = region.Polygon(True, poly_vert)
                fill_mask = polygon.scan(fill_mask)

            if self.mask is not None:
                fill_mask &= self.mask

            data = self.image[fill_mask]

            if data.size < 1:
                return (None, 0, 0.0)

        # Calculate sky
        try:
            skyval, npix = self._skystat(data)
        except ValueError:
            return (None, 0, 0.0)

        if not np.isfinite(skyval):
            return (None, 0, 0.0)

        if delta:
            skyval -= self._sky

        return skyval, npix, polyarea

    def _calc_sky_orig(self, overlap=None, delta=True):
        """
        Compute sky background value.

        Parameters
        ----------
        overlap : SkyImage, SkyGroup, SphericalPolygon, list of tuples, \
None, optional
            Another `SkyImage`, `SkyGroup`,
            :py:class:`spherical_geometry.polygons.SphericalPolygon`, or
            a list of tuples of (RA, DEC) of vertices of a spherical
            polygon. This parameter is used to indicate that sky statistics
            should computed only in the region of intersection of *this*
            image with the polygon indicated by `overlap`. When `overlap` is
            `None`, sky statistics will be computed over the entire image.

        delta : bool, optional
            Should this function return absolute sky value or the difference
            between the computed value and the value of the sky stored in the
            `sky` property.

        Returns
        -------
        skyval : float, None
            Computed sky value (absolute or relative to the `sky` attribute).
            If there are no valid data to perform this computations (e.g.,
            because this image does not overlap with the image indicated by
            `overlap`), `skyval` will be set to `None`.

        npix : int
            Number of pixels used to compute sky statistics.

        polyarea : float
            Area (in srad) of the polygon that bounds data used to compute
            sky statistics.

        """

        if overlap is None:

            if self.mask is None:
                data = self.image
            else:
                data = self.image[self.mask]

            polyarea = self.poly_area

        else:
            fill_mask = np.zeros(self.image.shape, dtype=bool)

            if isinstance(overlap, (SkyImage, SkyGroup, SphericalPolygon)):
                intersection = self.intersection(overlap)
                polyarea = np.fabs(intersection.area())
                radec = intersection.to_radec()

            else: # assume a list of (ra, dec) tuples:
                radec = []
                polyarea = 0.0
                for r, d in overlap:
                    poly = SphericalPolygon.from_radec(r, d)
                    polyarea1 = np.fabs(poly.area())
                    if polyarea1 == 0.0 or len(r) < 4:
                        continue
                    polyarea += polyarea1
                    radec.append(self.intersection(poly).to_radec())

            if polyarea == 0.0:
                return (None, 0, 0.0)

            for ra, dec in radec:
                if len(ra) < 4:
                    continue

                # set pixels in 'fill_mask' that are inside a polygon to True:
                x, y = self.wcs_inv(ra, dec)
                poly_vert = list(zip(*[x, y]))

                polygon = region.Polygon(True, poly_vert)
                fill_mask = polygon.scan(fill_mask)

            if self.mask is not None:
                fill_mask &= self.mask

            data = self.image[fill_mask]

            if data.size < 1:
                return (None, 0, 0.0)

        # Calculate sky
        try:

            skyval, npix = self._skystat(data)

        except ValueError:

            return (None, 0, 0.0)

        if delta:
            skyval -= self._sky

        return skyval, npix, polyarea

    def copy(self):
        """
        Return a shallow copy of the `SkyImage` object.
        """
        si = SkyImage(
            image=None,
            wcs_fwd=self.wcs_fwd,
            wcs_inv=self.wcs_inv,
            pix_area=self.pix_area,
            convf=self.convf,
            mask=None,
            id=self.id,
            stepsize=None,
            meta=self.meta
        )
        si.image = self.image
        si.mask = self.mask
        si._radec = self._radec
        si._polygon = self._polygon
        si._poly_area = self._poly_area
        si.sky = self.sky
        return si


class SkyGroup:
    """
    Holds multiple :py:class:`SkyImage` objects whose sky background values
    must be adjusted together.

    `SkyGroup` provides methods for obtaining bounding polygon of the group
    of :py:class:`SkyImage` objects and to compute sky value of the group.

    """
    def __init__(self, images, id=None, sky=0.0):

        if isinstance(images, SkyImage):
            self._images = [images]

        elif hasattr(images, '__iter__'):
            self._images = []
            for im in images:
                if not isinstance(im, SkyImage):
                    raise TypeError("Each element of the 'images' parameter "
                                    "must be an 'SkyImage' object.")
                self._images.append(im)

        else:
            raise TypeError("Parameter 'images' must be either a single "
                            "'SkyImage' object or a list of 'SkyImage' objects")

        self._id = id
        self._update_bounding_polygon()
        self._sky = sky
        for im in self._images:
            im.sky += sky

    @property
    def id(self):
        """ Set or get `SkyImage`'s `id`.

            While `id` can be of any type, it is prefereble that `id` be
            of a type with nice string representation.

        """
        return self._id

    @id.setter
    def id(self, id):
        self._id = id

    @property
    def sky(self):
        """ Sky background value. See `calc_sky` for more details.
        """
        return self._sky

    @sky.setter
    def sky(self, sky):
        delta_sky = sky - self._sky
        self._sky = sky
        for im in self._images:
            im.sky += delta_sky

    @property
    def radec(self):
        """
        Get RA and DEC of the verteces of the bounding polygon as a
        `~numpy.ndarray` of shape (N, 2) where N is the number of verteces + 1.

        """
        return self._radec

    @property
    def polygon(self):
        """ Get image's bounding polygon.
        """
        return self._polygon

    def intersection(self, skyimage):
        """
        Compute intersection of this `SkyImage` object and another
        `SkyImage`, `SkyGroup`, or
        :py:class:`~spherical_geometry.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        skyimage : SkyImage, SkyGroup, SphericalPolygon
            Another object that should be intersected with this `SkyImage`.

        Returns
        -------
        intersect_poly : SphericalPolygon
            A :py:class:`~spherical_geometry.polygon.SphericalPolygon` that is
            the intersection of this `SkyImage` and `skyimage`.

        """
        if isinstance(skyimage, (SkyImage, SkyGroup)):
            other = skyimage.polygon
        else:
            other = skyimage

        pts1 = np.sort(list(self._polygon.points)[0], axis=0)
        pts2 = np.sort(list(other.points)[0], axis=0)
        if np.allclose(pts1, pts2, rtol=0, atol=5e-9):
            intersect_poly = self._polygon.copy()
        else:
            intersect_poly = self._polygon.intersection(other)
        return intersect_poly


    def _update_bounding_polygon(self):
        polygons = [im.polygon for im in self._images]
        if len(polygons) == 0:
            self._polygon = SphericalPolygon([])
            self._radec = []
        else:
            self._polygon = SphericalPolygon.multi_union(polygons)
            self._radec = list(self._polygon.to_radec())

    def __len__(self):
        return len(self._images)

    def __getitem__(self, idx):
        return self._images[idx]

    def __setitem__(self, idx, value):
        if not isinstance(value, SkyImage):
            raise TypeError("Item must be of 'SkyImage' type")
        value.sky += self._sky
        self._images[idx] = value
        self._update_bounding_polygon()

    def __delitem__(self, idx):
        del self._images[idx]
        if len(self._images) == 0:
            self._sky = 0.0
            self._id = None
        self._update_bounding_polygon()

    def __iter__(self):
        for image in self._images:
            yield image

    def insert(self, idx, value):
        """Inserts a `SkyImage` into the group.
        """
        if not isinstance(value, SkyImage):
            raise TypeError("Item must be of 'SkyImage' type")
        value.sky += self._sky
        self._images.insert(idx, value)
        self._update_bounding_polygon()

    def append(self, value):
        """Appends a `SkyImage` to the group.
        """
        if not isinstance(value, SkyImage):
            raise TypeError("Item must be of 'SkyImage' type")
        value.sky += self._sky
        self._images.append(value)
        self._update_bounding_polygon()

    def calc_sky(self, overlap=None, delta=True):
        """
        Compute sky background value.

        Parameters
        ----------
        overlap : SkyImage, SkyGroup, SphericalPolygon, list of tuples, \
None, optional
            Another `SkyImage`, `SkyGroup`,
            :py:class:`spherical_geometry.polygons.SphericalPolygon`, or
            a list of tuples of (RA, DEC) of vertices of a spherical
            polygon. This parameter is used to indicate that sky statistics
            should computed only in the region of intersection of *this*
            image with the polygon indicated by `overlap`. When `overlap` is
            `None`, sky statistics will be computed over the entire image.

        delta : bool, optional
            Should this function return absolute sky value or the difference
            between the computed value and the value of the sky stored in the
            `sky` property.

        Returns
        -------
        skyval : float, None
            Computed sky value (absolute or relative to the `sky` attribute).
            If there are no valid data to perform this computations (e.g.,
            because this image does not overlap with the image indicated by
            `overlap`), `skyval` will be set to `None`.

        npix : int
            Number of pixels used to compute sky statistics.

        polyarea : float
            Area (in srad) of the polygon that bounds data used to compute
            sky statistics.

        """

        if len(self._images) == 0:
            return (None, 0, 0.0)

        wght = 0
        area = 0.0

        if overlap is None:
            # compute minimum sky across all images in the group:
            wsky = None

            for image in self._images:
                # make sure all images have the same background:
                image.background = self._sky

                sky, npix, imarea = image.calc_sky(overlap=None, delta=delta)

                if sky is None:
                    continue

                if wsky is None or wsky > sky:
                    wsky = sky
                    wght = npix
                    area = imarea

            return (wsky, wght, area)

        ################################################
        ##  compute weighted sky in various overlaps: ##
        ################################################
        wsky = 0.0

        for image in self._images:
            # make sure all images have the same background:
            image.background = self._sky

            sky, npix, area1 = image.calc_sky(overlap=overlap, delta=delta)

            area += area1

            if sky is not None and npix > 0:
                pix_area = npix * image.pix_area
                wsky += sky * pix_area
                wght += pix_area

        if wght == 0.0 or area == 0.0:
            return (None, wght, area)
        else:
            return (wsky / wght, wght, area)
