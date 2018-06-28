"""
This module provides support for working with image footprints on the sky,
source catalogs, and setting and manipulating tangent-plane corrections
of image WCS.

:Authors: Mihai Cara (contact: help@stsci.edu)


"""

# STDLIB
import logging
import numpy as np
from copy import deepcopy

# THIRD-PARTY
import numpy as np
import gwcs
from astropy import table
from spherical_geometry.polygon import SphericalPolygon
from stsci.stimage import xyxymatch

# LOCAL
from ..transforms.tpcorr import TPCorr, rot_mat3D
from . import linearfit
from . import matchutils
from . import __version__
from . import __vdate__


__all__ = ['convex_hull', 'ImageWCS', 'RefCatalog', 'WCSImageCatalog',
           'WCSGroupCatalog']


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ImageWCS():
    """ A class for holding JWST GWCS information and for managing
    tangent-plane corrections.

    """
    def __init__(self, wcs, v2_ref, v3_ref, roll_ref, ra_ref, dec_ref):
        """
        Parameters
        ----------

        wcs : GWCS
            A `GWCS` object.

        v2ref : float
            V2 position of the reference point in degrees.

        v3ref : float
            V3 position of the reference point in degrees.

        roll : float
            Roll angle in degrees.

        ra_ref : float
            RA of the reference point in degrees.

        dec_ref : float
            DEC of the reference point in degrees.

        """
        if not self._check_wcs_structure(wcs):
            raise ValueError("Unsupported WCS structure.")

        self._ra_ref = ra_ref
        self._dec_ref = dec_ref
        self._v2_ref = v2_ref
        self._v3_ref = v3_ref
        self._roll_ref = roll_ref

        # perform additional check that if tangent plane correction is already
        # present in the WCS pipeline, it is of TPCorr class and that
        # its parameters are consistent with reference angles:
        frms = [f[0] for f in wcs.pipeline]
        if 'v2v3corr' in frms:
            self._v23name = 'v2v3corr'
            self._tpcorr = deepcopy(wcs.pipeline[frms.index('v2v3corr')-1][1])
            self._default_tpcorr = None
            if not isinstance(self._tpcorr, TPCorr):
                raise ValueError("Unsupported tangent-plance correction type.")

            # check that transformation parameters are consistent with
            # reference angles:
            v2ref = self._tpcorr.v2ref.value
            v3ref = self._tpcorr.v3ref.value
            roll = self._tpcorr.roll.value
            eps_v2 = 10.0 * np.finfo(v2_ref).eps
            eps_v3 = 10.0 * np.finfo(v3_ref).eps
            eps_roll = 10.0 * np.finfo(roll_ref).eps
            if not (np.isclose(v2_ref, v2ref, rtol=eps_v2) and
                    np.isclose(v3_ref, v3ref, rtol=eps_v3) and
                    np.isclose(roll_ref, roll, rtol=eps_roll)):
                raise ValueError(
                    "WCS/TPCorr parameters 'v2ref', 'v3ref', and/or 'roll' "
                    "differ from the corresponding reference values."
                )

        else:
            self._v23name = 'v2v3'
            self._tpcorr = None
            self._default_tpcorr = TPCorr(
                v2ref=v2_ref, v3ref=v3_ref,
                roll=roll_ref,
                name='tangent-plane linear correction'
            )


        self._owcs = wcs
        self._wcs = deepcopy(wcs)
        self._update_transformations()

    def _update_transformations(self):
        # define transformations from detector/world coordinates to
        # the tangent plane:
        detname = self._wcs.pipeline[0][0]
        worldname = self._wcs.pipeline[-1][0]

        self._world_to_v23 = self._wcs.get_transform(worldname, self._v23name)
        self._v23_to_world = self._wcs.get_transform(self._v23name, worldname)
        self._det_to_v23 = self._wcs.get_transform(detname, self._v23name)
        self._v23_to_det = self._wcs.get_transform(self._v23name, detname)

        self._det_to_world = self._wcs.__call__
        self._world_to_det = self._wcs.invert

    @property
    def ref_angles(self):
        """ Return a ``wcsinfo``-like dictionary of main WCS parameters. """
        wcsinfo = {
            'ra_ref': self._ra_ref,
            'dec_ref': self._dec_ref,
            'v2_ref': self._v2_ref,
            'v3_ref': self._v3_ref,
            'roll_ref': self._roll_ref
        }
        return wcsinfo

    @property
    def wcs(self):
        """ Get current GWCS object. """
        return self._wcs

    @property
    def original_wcs(self):
        """ Get original GWCS object. """
        return self._owcs

    def copy(self):
        """ Returns a deep copy of this `ImageWCS` object. """
        return deepcopy(self)

    def set_correction(self, matrix=[[1, 0], [0, 1]], shift=[0, 0]):
        """
        Sets a tangent-plane correction of the GWCS object according to
        the provided liniar parameters.

        Parameters
        ----------
        matrix : list, numpy.ndarray
            A ``2x2`` array or list of lists coefficients representing scale,
            rotation, and/or skew transformations.

        shift : list, numpy.ndarray
            A list of two coordinate shifts to be applied to coordinates
            *before* ``matrix`` transformations are applied.

        """
        frms = [f[0] for f in self._wcs.pipeline]

        # if original WCS did not have tangent-plane corrections, create
        # new correction and add it to the WCs pipeline:
        if self._tpcorr is None:
            self._tpcorr = TPCorr(
                v2ref=self._v2_ref, v3ref=self._v3_ref,
                roll=self._roll_ref, matrix=matrix, shift=shift,
                name='tangent-plane linear correction'
            )
            idx_v2v3 = frms.index(self._v23name)
            pipeline = deepcopy(self._wcs.pipeline)
            pf, pt = pipeline[idx_v2v3]
            pipeline[idx_v2v3] = (pf, deepcopy(self._tpcorr))
            pipeline.insert(idx_v2v3 + 1, ('v2v3corr', pt))
            self._wcs = gwcs.WCS(pipeline, name=self._owcs.name)
            self._v23name = 'v2v3corr'

        else:
            # combine old and new corrections into a single one and replace
            # old transformation with the combined correction transformation:
            tpcorr2 = self._tpcorr.__class__(
                v2ref=self._tpcorr.v2ref, v3ref=self._tpcorr.v3ref,
                roll=self._tpcorr.roll, matrix=matrix, shift=shift,
                name='tangent-plane linear correction'
            )

            self._tpcorr = tpcorr2.combine(tpcorr2, self._tpcorr)

            idx_v2v3 = frms.index(self._v23name)
            pipeline = deepcopy(self._wcs.pipeline)
            pipeline[idx_v2v3 - 1] = (pipeline[idx_v2v3 - 1][0],
                                      deepcopy(self._tpcorr))
            self._wcs = gwcs.WCS(pipeline, name=self._owcs.name)

        # reset definitions of the transformations from detector/world
        # coordinates to the tangent plane:
        self._update_transformations()

    def _check_wcs_structure(self, wcs):
        if wcs is None or wcs.pipeline is None:
            return False

        frms = [f[0] for f in wcs.pipeline]
        nframes = len(frms)
        if nframes < 3:
            return False

        if frms.count(frms[0]) > 1 or frms.count(frms[-1]) > 1:
            return False

        if frms.count('v2v3') != 1:
            return False

        idx_v2v3 = frms.index('v2v3')
        if idx_v2v3 == 0 or idx_v2v3 == (nframes - 1):
            return False

        nv2v3corr = frms.count('v2v3corr')
        if nv2v3corr == 0:
            return True
        elif nv2v3corr > 1:
            return False

        idx_v2v3corr = frms.index('v2v3corr')
        if idx_v2v3corr != (idx_v2v3 + 1) or idx_v2v3corr == (nframes - 1):
            return False

        return True

    def det_to_world(self, x, y):
        """
        Convert pixel coordinates to sky coordinates using full
        (i.e., including distortions) transformations.

        """
        ra, dec = self._det_to_world(x, y)
        return ra, dec

    def world_to_det(self, ra, dec):
        """
        Convert sky coordinates to image's pixel coordinates using full
        (i.e., including distortions) transformations.

        """
        x, y = self._world_to_det(ra, dec)
        return x, y

    def det_to_tanp(self, x, y):
        """
        Convert detector (pixel) coordinates to tangent plane coordinates.

        """
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = self._det_to_v23(x, y)
        x, y = tpc.v2v3_to_tanp(v2, v3)
        return x, y

    def tanp_to_det(self, x, y):
        """
        Convert tangent plane coordinates to detector (pixel) coordinates.

        """
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = tpc.tanp_to_v2v3(x, y)
        x, y = self._v23_to_det(v2, v3)
        return x, y

    def world_to_tanp(self, ra, dec):
        """
        Convert tangent plane coordinates to detector (pixel) coordinates.

        """
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = self._world_to_v23(ra, dec)
        x, y = tpc.v2v3_to_tanp(v2, v3)
        return x, y

    def tanp_to_world(self, x, y):
        """
        Convert tangent plane coordinates to world coordinates.

        """
        tpc = self._default_tpcorr if self._tpcorr is None else self._tpcorr
        v2, v3 = tpc.tanp_to_v2v3(x, y)
        ra, dec = self._v23_to_world(v2, v3)
        return ra, dec


class WCSImageCatalog():
    """
    A class that holds information pertinent to an image WCS and a source
    catalog of the sources found in that image.

    """
    def __init__(self, shape, wcs, ref_angles, catalog, name=None, meta={}):
        """
        Parameters
        ----------

        shape : tuple
            A tuple of two integer values indicating the size of the image
            along each axis. Must follow the same convention as the shape of
            a :py:class:`numpy.ndarray` objects. Specifically,
            first size should be indicate the number of rows in the image and
            second size should indicate the number of columns in the image.

        wcs : gwcs.WCS
            WCS associated with the image and the catalog.

        ref_angles : dict
            A Python dictionary providing essential WCS reference angles. This
            parameter must contain at least the following keys:
            ``ra_ref``, ``dec_ref``, ``v2_ref``, ``v3_ref``, and ``roll_ref``.

        catalog : astropy.table.Table
            Source catalog associated with an image. Must contain 'x' and 'y'
            columns which indicate source coordinates (in pixels) in the
            associated image.

        name : str, None, optional
            Image catalog's name.

        meta : dict, optional
            Additional information about image, catalog, and/or WCS to be
            stored (for convenience) within `WCSImageCatalog` object.

        """
        self._shape = shape
        self._name = name
        self._catalog = None
        self._bb_radec = None

        self._imwcs = None
        self.img_bounding_ra = None
        self.img_bounding_dec = None

        self.meta = {}
        self.meta.update(meta)

        self.set_wcs(wcs, ref_angles)
        self.catalog = catalog

    @property
    def imwcs(self):
        """ Get :py:class:`ImageWCS` WCS. """
        return self._imwcs

    @property
    def wcs(self):
        """ Get :py:class:`gwcs.WCS`. """
        if self._imwcs is None:
            return None
        return self._imwcs.wcs

    @property
    def ref_angles(self):
        """ Get ``wcsinfo``. """
        if self._imwcs is None:
            return None
        return self._imwcs.ref_angles

    def set_wcs(self, wcs, ref_angles):
        """ Set :py:class:`gwcs.WCS` and the associated ``wcsinfo```.

        .. note::
            Setting the WCS triggers automatic bounding polygon recalculation.

        Parameters
        ----------

        wcs : gwcs.WCS
            WCS object.

        ref_angles : dict
            A Python dictionary providing essential WCS reference angles. This
            parameter must contain at least the following keys:
            ``ra_ref``, ``dec_ref``, ``v2_ref``, ``v3_ref``, and ``roll_ref``.

        """
        for key in ['dec_ref', 'v2_ref', 'v3_ref', 'roll_ref']:
            if key not in ref_angles:
                raise KeyError("Parameter 'ref_angles' must contain "
                               "'{:s}' key.".format(key))

        self._imwcs = ImageWCS(
            wcs=wcs,
            v2_ref=ref_angles['v2_ref'],
            v3_ref=ref_angles['v3_ref'],
            roll_ref=ref_angles['roll_ref'],
            ra_ref=ref_angles['ra_ref'],
            dec_ref=ref_angles['dec_ref']
        )

        # create spherical polygon bounding the image
        self.calc_bounding_polygon()

    @property
    def name(self):
        """ Get/set :py:class:`WCSImageCatalog` object's name.
        """
        return self._name

    @name.setter
    def name(self, name):
        self._name = name
        if hasattr(self, '_catalog'):
            if self._catalog is not None:
                self._catalog.meta['name'] = name

    @property
    def imshape(self):
        """
        Get/set image's shape. This must be a tuple of two dimensions
        following the same convention as the shape of `numpy.ndarray`.

        """
        return self._shape

    @imshape.setter
    def imshape(self, shape):
        self._shape = shape

    @property
    def catalog(self):
        """ Get/set image's catalog.
        """
        return self._catalog

    @catalog.setter
    def catalog(self, catalog):
        if catalog is None:
            self._catalog = None
            return

        else:
            self._catalog = table.Table(catalog.copy(), masked=True)
            self._catalog.meta['name'] = self._name

        colnames = self._catalog.colnames
        catlen = len(catalog)

        # create spherical polygon bounding the image
        self.calc_bounding_polygon()

    def det_to_world(self, x, y):
        """
        Convert pixel coordinates to sky coordinates using full
        (i.e., including distortions) transformations.

        """
        if self._imwcs is None:
            raise RuntimeError("WCS has not been set")
        return self._imwcs.det_to_world(x, y)

    def world_to_det(self, ra, dec):
        """
        Convert sky coordinates to image's pixel coordinates using full
        (i.e., including distortions) transformations.

        """
        if self._imwcs is None:
            raise RuntimeError("WCS has not been set")
        return self._imwcs.world_to_det(ra, dec)

    def det_to_tanp(self, x, y):
        """
        Convert detector (pixel) coordinates to tangent plane coordinates.

        """
        if self._imwcs is None:
            raise RuntimeError("WCS has not been set")
        return self._imwcs.det_to_tanp(x, y)

    def tanp_to_det(self, x, y):
        """
        Convert tangent plane coordinates to detector (pixel) coordinates.

        """
        if self._imwcs is None:
            raise RuntimeError("WCS has not been set")
        return self._imwcs.tanp_to_det(x, y)

    def tanp_to_world(self, x, y):
        """
        Convert tangent plane coordinates to world coordinates.

        """
        if self._imwcs is None:
            raise RuntimeError("WCS has not been set")
        return self._imwcs.tanp_to_world(x, y)

    def world_to_tanp(self, ra, dec):
        """
        Convert tangent plane coordinates to detector (pixel) coordinates.

        """
        if self._imwcs is None:
            raise RuntimeError("WCS has not been set")
        return self._imwcs.world_to_tanp(ra, dec)

    @property
    def polygon(self):
        """ Get image's (or catalog's) bounding spherical polygon.
        """
        return self._polygon

    def intersection(self, wcsim):
        """
        Compute intersection of this `WCSImageCatalog` object and another
        `WCSImageCatalog`, `WCSGroupCatalog`, or
        :py:class:`~spherical_geometry.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        wcsim : WCSImageCatalog, WCSGroupCatalog, SphericalPolygon
            Another object that should be intersected with this
            `WCSImageCatalog`.

        Returns
        -------
        polygon : SphericalPolygon
            A :py:class:`~spherical_geometry.polygon.SphericalPolygon` that is
            the intersection of this `WCSImageCatalog` and `wcsim`.

        """
        if isinstance(wcsim, (WCSImageCatalog, WCSGroupCatalog)):
            return self._polygon.intersection(wcsim.polygon)
        else:
            return self._polygon.intersection(wcsim)

    # TODO: due to a bug in the sphere package, see
    #       https://github.com/spacetelescope/sphere/issues/74
    #       intersections with polygons formed as union does not work.
    #       For this reason I re-implement 'intersection_area' below with
    #       a workaround for the bug.
    #       The original implementation should be uncommented once the bug
    #       is fixed.
    #
    #def intersection_area(self, wcsim):
        #""" Calculate the area of the intersection polygon.
        #"""
        #return np.fabs(self.intersection(wcsim).area())
    def intersection_area(self, wcsim):
        """ Calculate the area of the intersection polygon.
        """
        if isinstance(wcsim, (WCSImageCatalog, RefCatalog)):
            return np.fabs(self.intersection(wcsim).area())

        else:
            # this is bug workaround for image groups (multi-unions):
            area = 0.0
            for wim in wcsim:
                area += np.fabs(
                    self.polygon.intersection(wim.polygon).area()
                )
            return area

    def _calc_chip_bounding_polygon(self, stepsize=None):
        """
        Compute image's bounding polygon.

        Parameters
        ----------
        stepsize : int, None, optional
            Indicates the maximum separation between two adjacent vertices
            of the bounding polygon along each side of the image. Corners
            of the image are included automatically. If `stepsize` is `None`,
            bounding polygon will contain only vertices of the image.

        """
        if self.wcs is None:
            return

        ny, nx = self.imshape

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

        ra, dec = self.det_to_world(borderx, bordery)
        # TODO: for strange reasons, occasionally ra[0] != ra[-1] and/or
        #       dec[0] != dec[-1] (even though we close the polygon in the
        #       previous two lines). Then SphericalPolygon fails because
        #       points are not closed. Threfore we force it to be closed:
        ra[-1] = ra[0]
        dec[-1] = dec[0]

        self.img_bounding_ra = ra
        self.img_bounding_dec = dec
        self._polygon = SphericalPolygon.from_radec(ra, dec)

    def _calc_cat_convex_hull(self):
        """
        Compute convex hull that bounds the sources in the catalog.

        """
        if self.wcs is None or self.catalog is None:
            return

        x = self.catalog['x']
        y = self.catalog['y']

        if len(x) == 0:
            # no points
            raise RuntimeError("Unexpected error: Contact sofware developer")

        elif len(x) == 1:
            # one point. build a small box around it:
            x, y = convex_hull(x, y, wcs=self.det_to_tanp)

            xv = [x[0] - 0.5, x[0] - 0.5, x[0] + 0.5, x[0] + 0.5, x[0] - 0.5]
            yv = [y[0] - 0.5, y[0] + 0.5, y[0] + 0.5, y[0] - 0.5, y[0] - 0.5]

            ra, dec = self.tanp_to_world(xv, yv)

        elif len(x) == 3:
            # two points. build a small box around them:
            x, y = convex_hull(x, y, wcs=self.det_to_tanp)

            vx = -(y[1] - y[0])
            vy = x[1] - x[0]
            norm = 2.0 * np.sqrt(vx * vx + vy * vy)
            vx /= norm
            vy /= norm

            xv = [x[0] + vx, x[1] + vx, x[1] + vx, x[0] - vx, x[0] + vx]
            yv = [y[0] + vy, y[1] + vy, y[1] + vy, x[0] - vx, x[0] + vx]

            ra, dec = self.tanp_to_world(xv, yv)

        else:
            ra, dec = convex_hull(x, y, wcs=self.det_to_world)

        # TODO: for strange reasons, occasionally ra[0] != ra[-1] and/or
        #       dec[0] != dec[-1] (even though we close the polygon in the
        #       previous two lines). Then SphericalPolygon fails because
        #       points are not closed. Threfore we force it to be closed:
        ra[-1] = ra[0]
        dec[-1] = dec[0]

        self._bb_radec = (ra, dec)
        self._polygon = SphericalPolygon.from_radec(ra, dec)
        self._poly_area = np.fabs(self._polygon.area())

    def calc_bounding_polygon(self):
        """
        Calculate bounding polygon of the image or of the sources in the
        catalog (if catalog was set).

        """
        # we need image's footprint for later:
        self._calc_chip_bounding_polygon()

        # create smallest convex spherical polygon bounding all sources:
        if self._catalog is not None and len(self.catalog) > 0:
            self._calc_cat_convex_hull()

    @property
    def bb_radec(self):
        """
        Get a 2xN `numpy.ndarray` of RA and DEC of the vertices of the
        bounding polygon.

        """
        return self._bb_radec

class WCSGroupCatalog():
    """
    A class that holds together `WCSImageCatalog` image catalog objects
    whose relative positions are fixed and whose source catalogs should be
    fitted together to a reference catalog.

    """
    def __init__(self, images, name=None):
        """
        Parameters
        ----------
        images : list of WCSImageCatalog
            A list of `WCSImageCatalog` image catalogs.

        name : str, None, optional
            Name of the group.

        """

        self._catalog = None

        if isinstance(images, WCSImageCatalog):
            self._images = [images]

        elif hasattr(images, '__iter__'):
            self._images = []
            for im in images:
                if not isinstance(im, WCSImageCatalog):
                    raise TypeError("Each element of the 'images' parameter "
                                    "must be an 'WCSImageCatalog' object.")
                self._images.append(im)

        else:
            raise TypeError("Parameter 'images' must be either a single "
                            "'WCSImageCatalog' object or a list of "
                            "'WCSImageCatalog' objects")

        self._name = name
        self.update_bounding_polygon()
        self._catalog = self.create_group_catalog()

    @property
    def name(self):
        """ Get/set :py:class:`WCSImageCatalog` object's name.
        """
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def polygon(self):
        """ Get image's (or catalog's) bounding spherical polygon.
        """
        return self._polygon

    def intersection(self, wcsim):
        """
        Compute intersection of this `WCSGroupCatalog` object and another
        `WCSImageCatalog`, `WCSGroupCatalog`, or
        :py:class:`~spherical_geometry.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        wcsim : WCSImageCatalog, WCSGroupCatalog, SphericalPolygon
            Another object that should be intersected with this
            `WCSGroupCatalog`.

        Returns
        -------
        polygon : SphericalPolygon
            A :py:class:`~spherical_geometry.polygon.SphericalPolygon` that is
            the intersection of this `WCSGroupCatalog` and `wcsim`.

        """
        if isinstance(wcsim, (WCSGroupCatalog, WCSImageCatalog)):
            return self._polygon.intersection(wcsim.polygon)
        else:
            return self._polygon.intersection(wcsim)

    # TODO: due to a bug in the sphere package, see
    #       https://github.com/spacetelescope/sphere/issues/74
    #       intersections with polygons formed as union does not work.
    #       For this reason I re-implement 'intersection_area' below with
    #       a workaround for the bug.
    #       The original implementation should be uncommented once the bug
    #       is fixed.
    #
    #def intersection_area(self, wcsim):
        #""" Calculate the area of the intersection polygon.
        #"""
        #return np.fabs(self.intersection(wcsim).area())
    def intersection_area(self, wcsim):
        """ Calculate the area of the intersection polygon.
        """
        area = 0.0
        for im in self._images:
            area += im.intersection_area(wcsim)
        return area

    def update_bounding_polygon(self):
        """ Recompute bounding polygons of the member images.
        """
        polygons = [im.polygon for im in self._images]
        if len(polygons) == 0:
            self._polygon = SphericalPolygon([])
        else:
            self._polygon = SphericalPolygon.multi_union(polygons)

    def __len__(self):
        return len(self._images)

    def __getitem__(self, idx):
        return self._images[idx]

    def __iter__(self):
        for image in self._images:
            yield image

    @property
    def catalog(self):
        """ Get/set image's catalog.
        """
        return self._catalog

    def create_group_catalog(self):
        """
        Combine member's image catalogs into a single group's catalog.

        Returns
        -------

        group_catalog : astropy.table.Table
            Combined group catalog.

        """
        catalogs = []
        catno = 0
        for image in self._images:
            catlen = len(image.catalog)

            if image.name is None:
                catname = 'Catalog #{:d}'.format(catno)
            else:
                catname = image.name

            col_catname = table.MaskedColumn(catlen * [catname],
                                             name='cat_name')
            col_imcatidx = table.MaskedColumn(catlen * [catno],
                                              name='_imcat_idx')
            col_id = table.MaskedColumn(image.catalog['id'])
            col_x = table.MaskedColumn(image.catalog['x'], dtype=np.float64)
            col_y = table.MaskedColumn(image.catalog['y'], dtype=np.float64)
            ra, dec = image.det_to_world(
                image.catalog['x'], image.catalog['y']
            )
            col_ra = table.MaskedColumn(ra, dtype=np.float64, name='RA')
            col_dec = table.MaskedColumn(dec, dtype=np.float64, name='DEC')

            cat = table.Table(
                [col_imcatidx, col_catname, col_id, col_x,
                 col_y, col_ra, col_dec],
                masked=True
            )

            catalogs.append(cat)
            catno += 1

        return table.vstack(catalogs, join_type='exact')

    def get_unmatched_cat(self):
        """
        Retrieve only those sources from the catalog that have **not** been
        matched to the sources in the reference catalog.

        """
        mask = self._catalog['matched_ref_id'].mask
        return self._catalog[mask]

    def get_matched_cat(self):
        """
        Retrieve only those sources from the catalog that **have been**
        matched to the sources in the reference catalog.

        """
        mask = np.logical_not(self._catalog['matched_ref_id'].mask)
        return self._catalog[mask]

    def recalc_catalog_radec(self):
        """ Recalculate RA and DEC of the sources in the catalog.
        """
        for k, image in enumerate(self._images):

            idx = (self._catalog['_imcat_idx'] == k)
            if not np.any(idx):
                continue

            ra, dec = image.det_to_world(
                self._catalog['x'][idx], self._catalog['y'][idx]
            )
            self._catalog['RA'][idx] = ra
            self._catalog['DEC'][idx] = dec

    def calc_tanp_xy(self, tanplane_wcs):
        """
        Compute x- and y-positions of the sources from the image catalog
        in the tangent plane. This creates the following
        columns in the catalog's table: ``'xtanp'`` and ``'ytanp'``.

        Parameters
        ----------
        tanplane_wcs : ImageWCS
            A `ImageWCS` object that will provide transformations to
            the tangent plane to which sources of this catalog a should be
            "projected".

        """
        if 'RA' not in self._catalog.colnames or \
           'DEC' not in self._catalog.colnames:
            raise RuntimeError("'recalc_catalog_radec()' should have been run "
                               "prior to calc_tanp_xy()")

        # compute x & y in the reference WCS:
        xtp, ytp = tanplane_wcs.world_to_tanp(self.catalog['RA'],
                                              self.catalog['DEC'])
        self._catalog['xtanp'] = table.MaskedColumn(
            xtp, name='xtanp', dtype=np.float64, mask=False
        )
        self._catalog['ytanp'] = table.MaskedColumn(
            ytp, name='ytanp', dtype=np.float64, mask=False
        )

    def match2ref(self, refcat, minobj=15, searchrad=1.0,
                  separation=0.5, use2dhist=True, xoffset=0.0, yoffset=0.0,
                  tolerance=1.0):
        """ Uses xyxymatch to cross-match sources between this catalog and
            a reference catalog.

        Parameters
        ----------
        refcat : RefCatalog
            A `RefCatalog` object that contains a catalog of reference sources
            as well as a valid reference WCS.

        minobj : int, None, optional
            Minimum number of identified objects from each input image to use
            in matching objects from other images. If the default `None` value
            is used then `align` will automatically deternmine the minimum
            number of sources from the value of the `fitgeom` parameter.

        searchrad : float, optional
            The search radius for a match.

        separation : float, optional
            The  minimum  separation for sources in the input and reference
            catalogs in order to be considered to be disctinct sources.
            Objects closer together than 'separation' pixels
            are removed from the input and reference coordinate lists prior
            to matching. This parameter gets passed directly to
            :py:func:`~stsci.stimage.xyxymatch` for use in matching the object
            lists from each image with the reference image's object list.

        use2dhist : bool, optional
            Use 2D histogram to find initial offset?

        xoffset : float, optional
            Initial estimate for the offset in X between the images and the
            reference frame. This offset will be used for all input images
            provided. This parameter is ignored when `use2dhist` is `True`.

        yoffset : float (Default = 0.0)
            Initial estimate for the offset in Y between the images and the
            reference frame. This offset will be used for all input images
            provided. This parameter is ignored when `use2dhist` is `True`.

        tolerance : float, optional
            The matching tolerance in pixels after applying an initial solution
            derived from the 'triangles' algorithm.  This parameter gets passed
            directly to :py:func:`~stsci.stimage.xyxymatch` for use in
            matching the object lists from each image with the reference
            image's object list.

        """

        colnames = self._catalog.colnames

        if 'xtanp' not in colnames or 'ytanp' not in colnames:
            raise RuntimeError("'calc_tanp_xy()' should have been run prior "
                               "to match2ref()")

        im_xyref = np.asanyarray([self._catalog['xtanp'],
                                  self._catalog['ytanp']]).T
        refxy = np.asanyarray([refcat.catalog['xtanp'],
                               refcat.catalog['ytanp']]).T

        log.info("Matching sources from '{}' with sources from reference "
                 "{:s} '{}'".format(self.name, 'image', refcat.name))

        xyoff = (xoffset, yoffset)

        if use2dhist:
            # Determine xyoff (X,Y offset) and tolerance
            # to be used with xyxymatch:
            zpxoff, zpyoff, flux, zpqual = matchutils.build_xy_zeropoint(
                im_xyref,
                refxy,
                searchrad=searchrad
            )

            if zpqual is not None:
                xyoff = (zpxoff, zpyoff)
                # set tolerance as well
                # This value allows initial guess to be off by 1 in both and
                # still pick up the identified matches
                tolerance = 1.5

        matches = xyxymatch(
            im_xyref,
            refxy,
            origin=xyoff,
            tolerance=tolerance,
            separation=separation
        )

        nmatches = len(matches)
        self._catalog.meta['nmatches'] = nmatches
        minput_idx = matches['input_idx']

        catlen = len(self._catalog)

        # matched_ref_id:
        if 'matched_ref_id' not in colnames:
            c = table.MaskedColumn(name='matched_ref_id', dtype=int,
                                   length=catlen, mask=True)
            self._catalog.add_column(c)
        else:
            self._catalog['matched_ref_id'].mask = True
        self._catalog['matched_ref_id'][minput_idx] = \
            self._catalog['id'][minput_idx]
        self._catalog['matched_ref_id'].mask[minput_idx] = False

        # this is needed to index reference catalog directly without using
        # astropy table indexing which at this moment is experimental:
        if '_raw_matched_ref_idx' not in colnames:
            c = table.MaskedColumn(name='_raw_matched_ref_idx',
                                   dtype=int, length=catlen, mask=True)
            self._catalog.add_column(c)
        else:
            self._catalog['_raw_matched_ref_idx'].mask = True
        self._catalog['_raw_matched_ref_idx'][minput_idx] = \
            matches['ref_idx']
        self._catalog['_raw_matched_ref_idx'].mask[minput_idx] = False

        log.info("Found {:d} matches for '{}'...".format(nmatches, self.name))

        return matches

    def fit2ref(self, refcat, tanplane_wcs, fitgeom='general', nclip=3,
                sigma=3.0):
        """
        Perform linear fit of this group's combined catalog to the reference
        catalog.


        Parameters
        ----------

        refcat : RefCatalog
            A `RefCatalog` object that contains a catalog of reference sources.

        tanplane_wcs : ImageWCS
            A `ImageWCS` object that will provide transformations to
            the tangent plane to which sources of this catalog a should be
            "projected".

        fitgeom : {'shift', 'rscale', 'general'}, optional
            The fitting geometry to be used in fitting the matched object
            lists. This parameter is used in fitting the offsets, rotations
            and/or scale changes from the matched object lists. The 'general'
            fit geometry allows for independent scale and rotation for
            each axis.

        nclip : int, optional
            Number (a non-negative integer) of clipping iterations in fit.

        sigma : float, optional
            Clipping limit in sigma units.

        """
        im_xyref = np.asanyarray([self._catalog['xtanp'],
                                  self._catalog['ytanp']]).T
        refxy = np.asanyarray([refcat.catalog['xtanp'],
                               refcat.catalog['ytanp']]).T

        mask = np.logical_not(self._catalog['matched_ref_id'].mask)
        im_xyref = im_xyref[mask]
        ref_idx = self._catalog['_raw_matched_ref_idx'][mask]
        refxy = refxy[ref_idx]

        fit = linearfit.iter_linear_fit(
            refxy, im_xyref, fitgeom=fitgeom,
            nclip=nclip, sigma=sigma, center=(0, 0)
        )

        xy_fit = fit['img_coords'] + fit['resids']
        fit['fit_xy'] = xy_fit
        fit['fit_RA'], fit['fit_DEC'] = tanplane_wcs.tanp_to_world(*(xy_fit.T))

        log.info("Computed '{:s}' fit for {}:".format(fitgeom, self.name))
        if fitgeom == 'shift':
            log.info("XSH: {:.6g}  YSH: {:.6g}"
                     .format(fit['offset'][0], fit['offset'][1]))
        elif fitgeom == 'rscale' and fit['proper']:
            log.info("XSH: {:.6g}  YSH: {:.6g}    ROT: {:.6g}    SCALE: {:.6g}"
                     .format(fit['offset'][0], fit['offset'][1],
                             fit['rot'], fit['scale'][0]))
        elif fitgeom == 'general' or (fitgeom == 'rscale' and not
                                      fit['proper']):
            log.info("XSH: {:.6g}  YSH: {:.6g}    PROPER ROT: {:.6g}    "
                     .format(fit['offset'][0], fit['offset'][1], fit['rot']))
            log.info("<ROT>: {:.6g}  SKEW: {:.6g}    ROT_X: {:.6g}  "
                     "ROT_Y: {:.6g}".format(fit['rotxy'][2], fit['skew'],
                                            fit['rotxy'][0], fit['rotxy'][1]))
            log.info("<SCALE>: {:.6g}  SCALE_X: {:.6g}  SCALE_Y: {:.6g}"
                     .format(fit['scale'][0], fit['scale'][1],
                             fit['scale'][2]))
        else:
            raise ValueError("Unsupported fit geometry.")

        log.info("")
        log.info("XRMS: {:.6g}    YRMS: {:.6g}".format(fit['rms'][0],
                                                       fit['rms'][1]))
        log.info("Final solution based on {:d} objects."
                 .format(fit['resids'].shape[0]))

        return fit


    def apply_affine_to_wcs(self, tanplane_wcs, matrix, shift):
        """ Applies a general affine transformation to the WCS.
        """
        # compute the matrix for the scale and rotation correction
        matrix = matrix.T
        shift = -np.dot(np.linalg.inv(matrix), shift)

        for imcat in self:
            # compute linear transformation from the tangent plane used for
            # alignment to the tangent plane of another image in the group:
            if imcat.imwcs == tanplane_wcs:
                m = matrix.copy()
                s = shift.copy()
            else:
                r1, t1 = _tp2tp(imcat.imwcs, tanplane_wcs)
                r2, t2 = _tp2tp(tanplane_wcs, imcat.imwcs)
                m = np.linalg.multi_dot([r2, matrix, r1])
                s = t1 + np.dot(np.linalg.inv(r1), shift) + \
                    np.dot(np.linalg.inv(np.dot(matrix, r1)), t2)

            imcat.imwcs.set_correction(m, s)
            imcat.meta['image_model'].meta.wcs = imcat.wcs

    def align_to_ref(self, refcat, minobj=15, searchrad=1.0, separation=0.5,
                     use2dhist=True, xoffset=0.0, yoffset=0.0, tolerance=1.0,
                     fitgeom='rscale', nclip=3, sigma=3.0):
        """
        Matches sources from the image catalog to the sources in the
        reference catalog, finds the affine transformation between matched
        sources, and adjusts images' WCS according to this fit.

        Parameters
        ----------

        refcat : RefCatalog
            A `RefCatalog` object that contains a catalog of reference sources
            as well as a valid reference WCS.

        minobj : int, None, optional
            Minimum number of identified objects from each input image to use
            in matching objects from other images. If the default `None` value
            is used then `align` will automatically deternmine the minimum
            number of sources from the value of the `fitgeom` parameter.

        searchrad : float, optional
            The search radius for a match.

        separation : float, optional
            The  minimum  separation for sources in the input and reference
            catalogs in order to be considered to be disctinct sources.
            Objects closer together than 'separation' pixels
            are removed from the input and reference coordinate lists prior
            to matching. This parameter gets passed directly to
            :py:func:`~stsci.stimage.xyxymatch` for use in matching the object
            lists from each image with the reference image's object list.

        use2dhist : bool, optional
            Use 2D histogram to find initial offset?

        xoffset : float, optional
            Initial estimate for the offset in X between the images and the
            reference frame. This offset will be used for all input images
            provided. This parameter is ignored when `use2dhist` is `True`.

        yoffset : float (Default = 0.0)
            Initial estimate for the offset in Y between the images and the
            reference frame. This offset will be used for all input images
            provided. This parameter is ignored when `use2dhist` is `True`.

        tolerance : float, optional
            The matching tolerance in pixels after applying an initial solution
            derived from the 'triangles' algorithm.  This parameter gets passed
            directly to :py:func:`~stsci.stimage.xyxymatch` for use in
            matching the object lists from each image with the reference
            image's object list.

        fitgeom : {'shift', 'rscale', 'general'}, optional
            The fitting geometry to be used in fitting the matched object
            lists. This parameter is used in fitting the offsets, rotations
            and/or scale changes from the matched object lists. The 'general'
            fit geometry allows for independent scale and rotation for each
            axis.

        nclip : int, optional
            Number (a non-negative integer) of clipping iterations in fit.

        sigma : float, optional
            Clipping limit in sigma units.

        """
        if len(self._images) == 0:
            return

        tanplane_wcs = deepcopy(self._images[0].imwcs)

        self.calc_tanp_xy(tanplane_wcs=tanplane_wcs)
        refcat.calc_tanp_xy(tanplane_wcs=tanplane_wcs)
        self.match2ref(refcat=refcat, minobj=minobj, searchrad=searchrad,
                       separation=separation,
                       use2dhist=use2dhist, xoffset=xoffset, yoffset=yoffset,
                       tolerance=tolerance)
        fit = self.fit2ref(refcat=refcat, tanplane_wcs=tanplane_wcs,
                           fitgeom=fitgeom, nclip=nclip, sigma=sigma)
        self.apply_affine_to_wcs(
            tanplane_wcs=tanplane_wcs,
            matrix=fit['fit_matrix'],
            shift=fit['offset']
        )


def _tp2tp(imwcs1, imwcs2):
    x = np.array([0.0, 1.0, 0.0], dtype=np.float)
    y = np.array([0.0, 0.0, 1.0], dtype=np.float)
    xrp, yrp = imwcs2.world_to_tanp(*imwcs1.tanp_to_world(x, y))

    matrix = np.array([(xrp[1:] - xrp[0]), (yrp[1:] - yrp[0])])
    shift = -np.dot(np.linalg.inv(matrix), [xrp[0], yrp[0]])

    return matrix, shift


class RefCatalog():
    """
    An object that holds a reference catalog and provides
    tools for coordinate convertions using reference WCS as well as
    catalog manipulation and expansion.

    """
    def __init__(self, catalog, name=None):
        """
        Parameters
        ----------
        catalog : astropy.table.Table
            Reference catalog.

            ..note::
                Reference catalogs (:py:class:`~astropy.table.Table`)
                *must* contain *both* ``'RA'`` and ``'DEC'`` columns.

        name : str, None, optional
            Name of the reference catalog.

        """
        self._name = name
        self._catalog = None

        # make sure catalog has RA & DEC
        if catalog is not None:
            self.catalog = catalog

    def _check_catalog(self, catalog):
        if catalog is None:
            raise ValueError("Reference catalogs cannot be None")

        if 'RA' not in catalog.colnames or 'DEC' not in catalog.colnames:
            raise KeyError("Reference catalogs *must* contain *both* 'RA' "
                           "and 'DEC' columns.")

    @property
    def name(self):
        """ Get/set :py:class:`WCSImageCatalog` object's name.
        """
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def catalog(self):
        """ Get/set image's catalog.
        """
        return self._catalog

    @catalog.setter
    def catalog(self, catalog):
        self._check_catalog(catalog)

        if len(catalog) == 0:
            raise ValueError("Catalog must contain at least one source.")

        self._catalog = catalog.copy()

        # create spherical polygon bounding the sources
        self.calc_bounding_polygon()

    @property
    def poly_area(self):
        """ Area of the bounding polygon (in srad).
        """
        return self._poly_area

    @property
    def polygon(self):
        """ Get image's (or catalog's) bounding spherical polygon.
        """
        return self._polygon

    def intersection(self, wcsim):
        """
        Compute intersection of this `WCSImageCatalog` object and another
        `WCSImageCatalog`, `WCSGroupCatalog`, `RefCatalog`, or
        :py:class:`~spherical_geometry.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        wcsim : WCSImageCatalog, WCSGroupCatalog, RefCatalog, SphericalPolygon
            Another object that should be intersected with this
            `WCSImageCatalog`.

        Returns
        -------
        polygon : SphericalPolygon
            A :py:class:`~spherical_geometry.polygon.SphericalPolygon` that is
            the intersection of this `WCSImageCatalog` and `wcsim`.

        """
        if isinstance(wcsim, (WCSImageCatalog, WCSGroupCatalog, RefCatalog)):
            return self._polygon.intersection(wcsim.polygon)
        else:
            return self._polygon.intersection(wcsim)

    # TODO: due to a bug in the sphere package, see
    #       https://github.com/spacetelescope/sphere/issues/74
    #       intersections with polygons formed as union does not work.
    #       For this reason I re-implement 'intersection_area' below with
    #       a workaround for the bug.
    #       The original implementation should be uncommented once the bug
    #       is fixed.
    #
    #def intersection_area(self, wcsim):
        #""" Calculate the area of the intersection polygon.
        #"""
        #return np.fabs(self.intersection(wcsim).area())
    def intersection_area(self, wcsim):
        """ Calculate the area of the intersection polygon.
        """
        if isinstance(wcsim, (WCSImageCatalog, RefCatalog)):
            return np.fabs(self.intersection(wcsim).area())

        else:
            # this is bug workaround:
            area = 0.0
            for wim in wcsim:
                area += np.fabs(
                    self.polygon.intersection(wim.polygon).area()
                )
            return area

    def _calc_cat_convex_hull(self):
        """
        Calculate spherical polygon corresponding to the convex hull bounding
        the sources in the catalog.

        """
        if self.catalog is None:
            return

        # Find an "optimal" tangent plane to the catalog points based on the
        # mean point and then construct a WCS based on the mean point.
        # Compute x, y coordinates in this tangent plane based on the
        # previously computed WCS and return the set of x, y coordinates and
        # "reference WCS".
        x, y, z = TPCorr.spherical2cartesian(
            self.catalog['RA'], self.catalog['DEC']
        )
        ra_ref, dec_ref = TPCorr.cartesian2spherical(
            x.mean(dtype=np.float64),
            y.mean(dtype=np.float64),
            z.mean(dtype=np.float64)
        )
        rotm = [rot_mat3D(np.deg2rad(alpha), 2 - axis)
                for axis, alpha in enumerate([ra_ref, dec_ref])]
        euler_rot = np.linalg.multi_dot(rotm)
        inv_euler_rot = np.linalg.inv(euler_rot)
        xr, yr, zr = np.dot(euler_rot, (x, y, z))
        r0 = TPCorr.r0
        x = r0 * yr / xr
        y = r0 * zr / xr

        xv, yv = convex_hull(x, y)

        if len(xv) == 0:
            # no points
            raise RuntimeError("Unexpected error: Contact sofware developer")

        elif len(xv) == 1:
            # one point. build a small box around it:
            x, y = convex_hull(x, y, wcs=None)

            xv = [x[0] - 0.5, x[0] - 0.5, x[0] + 0.5, x[0] + 0.5, x[0] - 0.5]
            yv = [y[0] - 0.5, y[0] + 0.5, y[0] + 0.5, y[0] - 0.5, y[0] - 0.5]

        elif len(xv) == 3:
            # two points. build a small box around them:
            x, y = convex_hull(x, y, wcs=None)

            vx = -(y[1] - y[0])
            vy = x[1] - x[0]
            norm = 2.0 * np.sqrt(vx * vx + vy * vy)
            vx /= norm
            vy /= norm

            xv = [x[0] + vx, x[1] + vx, x[1] + vx, x[0] - vx, x[0] + vx]
            yv = [y[0] + vy, y[1] + vy, y[1] + vy, x[0] - vx, x[0] + vx]

        # "unrotate" cartezian coordinates back to their original
        # ra_ref and dec_ref "positions":
        xt = np.full_like(xv, r0)
        xcr, ycr, zcr = np.dot(inv_euler_rot, (xt, xv, yv))
        # convert cartesian to spherical coordinates:
        ra, dec = TPCorr.cartesian2spherical(xcr, ycr, zcr)

        # TODO: for strange reasons, occasionally ra[0] != ra[-1] and/or
        #       dec[0] != dec[-1] (even though we close the polygon in the
        #       previous two lines). Then SphericalPolygon fails because
        #       points are not closed. Threfore we force it to be closed:
        ra[-1] = ra[0]
        dec[-1] = dec[0]

        self._radec = [(ra, dec)]
        self._polygon = SphericalPolygon.from_radec(ra, dec)
        self._poly_area = np.fabs(self._polygon.area())

    def calc_bounding_polygon(self):
        """ Calculate bounding polygon of the sources in the catalog.
        """
        # create spherical polygon bounding the sources
        self._calc_cat_convex_hull()

    def expand_catalog(self, catalog):
        """
        Expand current reference catalog with sources from another catalog.

        Parameters
        ----------
        catalog : astropy.table.Table
            A catalog of new sources to be added to the existing reference
            catalog. `catalog` *must* contain *both* ``'RA'`` and ``'DEC'``
            columns.

        """
        self._check_catalog(catalog)
        cat = catalog.copy()
        if self._catalog is None:
            self._catalog = cat
        else:
            self._catalog = table.vstack([self.catalog, cat],
                                         join_type='outer')
        self.calc_bounding_polygon()

    def calc_tanp_xy(self, tanplane_wcs):
        """
        Compute x- and y-positions of the sources from the reference catalog
        in the tangent plane provided by the `tanplane_wcs`.
        This creates the following columns in the catalog's table:
        ``'xtanp'`` and ``'ytanp'``.

        Parameters
        ----------
        tanplane_wcs : ImageWCS
            A `ImageWCS` object that will provide transformations to
            the tangent plane to which sources of this catalog a should be
            "projected".

        """
        # compute x & y in the reference WCS:
        xtp, ytp = tanplane_wcs.world_to_tanp(self.catalog['RA'],
                                              self.catalog['DEC'])
        self._catalog['xtanp'] = table.MaskedColumn(
            xtp, name='xtanp', dtype=np.float64, mask=False
        )
        self._catalog['ytanp'] = table.MaskedColumn(
            ytp, name='ytanp', dtype=np.float64, mask=False
        )


def convex_hull(x, y, wcs=None):
    """Computes the convex hull of a set of 2D points.

    Implements `Andrew's monotone chain algorithm <http://en.wikibooks.org\
/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain>`_.
    The algorithm has O(n log n) complexity.

    Credit: `<http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/\
Convex_hull/Monotone_chain>`_

    Parameters
    ----------

    points : list of tuples
        An iterable sequence of (x, y) pairs representing the points.

    Returns
    -------
    Output : list
        A list of vertices of the convex hull in counter-clockwise order,
        starting from the vertex with the lexicographically smallest
        coordinates.

    """
    ndarray = isinstance(x, np.ndarray) or isinstance(y, np.ndarray)

    # Sort the points lexicographically (tuples are compared
    # lexicographically). Remove duplicates to detect the case we have
    # just one unique point.
    points = sorted(set(zip(x, y)))

    # Boring case: no points or a single point,
    # possibly repeated multiple times.
    if len(points) == 0:
        if not ndarray:
            return (np.array([]), np.array([]))
        else:
            return ([], [])
    elif len(points) == 1:
        if not ndarray:
            return (np.array([points[0][0]]), np.array([points[0][1]]))
        else:
            return ([points[0][0]], [points[0][1]])

    # 2D cross product of OA and OB vectors, i.e. z-component of their
    # 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the
    # beginning of the other list.
    total_hull = np.asanyarray(lower[:-1] + upper)

    ptx = total_hull[:, 0]
    pty = total_hull[:, 1]

    if wcs is None:
        if not ndarray:
            return (ptx.tolist(), pty.tolist())
        else:
            return (ptx, pty)

    # convert x, y vertex coordinates to RA & DEC:
    ra, dec = wcs(ptx, pty)
    ra[-1] = ra[0]
    dec[-1] = dec[0]

    if not ndarray:
        return (ra.tolist(), dec.tolist())
    else:
        return (ra, dec)
