"""
This module provides support for working with image footprints on the sky and
source catalogs.

:Authors: Mihai Cara (contact: help@stsci.edu)

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

# STDLIB
import logging
import numpy as np

# THIRD-PARTY
from astropy import table
from stsci.sphere.polygon import SphericalPolygon
from stsci.stimage import xyxymatch

# LOCAL
from .tpcorr import ImageWCS
from . import linearfit
from . import matchutils
from . import wcsutils
from . import tpcorr
from . import __version__
from . import __vdate__


__all__ = ['WCSImageCatalog', 'WCSGroupCatalog', 'RefCatalog', 'convex_hull']


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class WCSImageCatalog(object):
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
        :py:class:`~stsci.sphere.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        wcsim : WCSImageCatalog, WCSGroupCatalog, SphericalPolygon
            Another object that should be intersected with this
            `WCSImageCatalog`.

        Returns
        -------
        polygon : SphericalPolygon
            A :py:class:`~stsci.sphere.polygon.SphericalPolygon` that is
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


class WCSGroupCatalog(object):
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
        :py:class:`~stsci.sphere.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        wcsim : WCSImageCatalog, WCSGroupCatalog, SphericalPolygon
            Another object that should be intersected with this
            `WCSGroupCatalog`.

        Returns
        -------
        polygon : SphericalPolygon
            A :py:class:`~stsci.sphere.polygon.SphericalPolygon` that is
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
        print("here")
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
        tanplane_wcs : tpcorr.ImageWCS
            A `tpcorr.ImageWCS` object that will provide transformations to
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

        print("\nmatches:\n{}\n".format(matches))

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

    def fit2ref(self, refcat, fitgeom='general', nclip=3, sigma=3.0):
        """
        Perform linear fit of this group's combined catalog to the reference
        catalog.


        Parameters
        ----------

        refcat : RefCatalog
            A `RefCatalog` object that contains a catalog of reference sources
            as well as a valid reference WCS.

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

        print("\nim_xyref:\n")
        for x, y in im_xyref:
            print("({:.16g}, {:.16g})".format(x, y))
        print("\nrefxy:\n")
        for x, y in refxy:
            print("({:.16g}, {:.16g})".format(x, y))
        print("\n")

        print("nclip = {}, fitgeom = '{}', sigma = {}\n".format(nclip, fitgeom, sigma))

        fit = linearfit.iter_linear_fit(
            refxy, im_xyref, fitgeom=fitgeom,
            nclip=nclip, sigma=sigma, center=(0, 0)#refcat.crpix
        )

        print(fit,"\n")
        print("\nm11 = {:.16g} m12 = {:.16g}\nm21 = {:.16g} m22 = {:.16g}"
              .format(*(fit['fit_matrix'].ravel())))
        print("sh1 = {:.16g} sh2 = {:.16g}\n"
              .format(*fit['offset']))

        fit['rms_keys'] = {'RMS_RA': 0, 'RMS_DEC': 0, 'NMATCH': 0} #self._compute_fit_rms(fit, refcat)
        xy_fit = fit['img_coords'] + fit['resids']
        fit['fit_xy'] = xy_fit
        #ra_fit, dec_fit = refcat.wcs_pix2world(xy_fit[0], xy_fit[1])
        fit['fit_RA'] = xy_fit[0] #ra_fit
        fit['fit_DEC'] = xy_fit[1] #dec_fit

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
        log.info("RMS_RA: {:g} (deg)   RMS_DEC: {:g} (deg)"
                 .format(fit['rms_keys']['RMS_RA'],
                         fit['rms_keys']['RMS_DEC']))
        log.info("Final solution based on {:d} objects."
                 .format(fit['rms_keys']['NMATCH']))

        return fit

    def _compute_fit_rms(self, fit, refwcs):
        # start by interpreting the fit to get the RMS values
        crpix = refwcs.crpix + fit['rms']
        rms_ra, rms_dec = refwcs.wcs_pix2world(crpix[0], crpix[1])
        crval1, crval2 = refwcs.crval
        rms_ra = np.abs(rms_ra - crval1)
        rms_dec = np.abs(rms_dec - crval2)
        nmatch = fit['resids'].shape[0]
        return {'RMS_RA': rms_ra, 'RMS_DEC': rms_dec, 'NMATCH': nmatch}

    def apply_affine_to_wcs(self, refcat, matrix, xsh, ysh):
        """ Applies a general affine transformation to the WCS.
        """
        for imcat in self:
            wcsutils.apply_affine_to_wcs(
                imwcs=imcat.imwcs,
                xsh=xsh,
                ysh=ysh,
                matrix=matrix
            )
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

        self.calc_tanp_xy(tanplane_wcs=self._images[0])
        refcat.calc_tanp_xy(tanplane_wcs=self._images[0])
        self.match2ref(refcat=refcat, minobj=minobj, searchrad=searchrad,
                       separation=separation,
                       use2dhist=use2dhist, xoffset=xoffset, yoffset=yoffset,
                       tolerance=tolerance)
        print(self._catalog)
        print(refcat._catalog)
        fit = self.fit2ref(refcat=refcat, fitgeom=fitgeom,
                           nclip=nclip, sigma=sigma)
        print(fit)
        self.apply_affine_to_wcs(refcat=refcat, matrix=fit['fit_matrix'],
                                 xsh=fit['offset'][0], ysh=fit['offset'][1])


class RefCatalog(object):
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
        :py:class:`~stsci.sphere.polygon.SphericalPolygon`
        object.

        Parameters
        ----------
        wcsim : WCSImageCatalog, WCSGroupCatalog, RefCatalog, SphericalPolygon
            Another object that should be intersected with this
            `WCSImageCatalog`.

        Returns
        -------
        polygon : SphericalPolygon
            A :py:class:`~stsci.sphere.polygon.SphericalPolygon` that is
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
        x, y, z = tpcorr.TPCorr.spherical2cartesian(
            self.catalog['RA'], self.catalog['DEC']
        )
        ra_ref, dec_ref = tpcorr.TPCorr.cartesian2spherical(
            x.mean(dtype=np.float64),
            y.mean(dtype=np.float64),
            z.mean(dtype=np.float64)
        )
        rotm = [tpcorr.rot_mat3D(np.deg2rad(alpha), 2 - axis)
                for axis, alpha in enumerate([ra_ref, dec_ref])]
        euler_rot = np.linalg.multi_dot(rotm)
        inv_euler_rot = np.linalg.inv(euler_rot)
        xr, yr, zr = np.dot(euler_rot, (x, y, z))
        r0 = tpcorr.TPCorr.r0
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
        ra, dec = tpcorr.TPCorr.cartesian2spherical(xcr, ycr, zcr)

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
        tanplane_wcs : tpcorr.ImageWCS
            A `tpcorr.ImageWCS` object that will provide transformations to
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
