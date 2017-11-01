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
from .simplewcs import SimpleWCS
from . import linearfit
from . import matchutils
from . import wcsutils


__all__ = ['WCSImageCatalog', 'WCSGroupCatalog', 'RefCatalog', 'convex_hull']
__version__ = '0.7.1'
__vdate__ = '18-April-2016'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class WCSImageCatalog(object):
    """
    A class that holds information pertinent to an image WCS and a source
    catalog of the sources found in that image.

    """
    def __init__(self, shape, wcs, catalog, name=None, meta={}):
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

        self._swcs = None
        self.img_bounding_ra = None
        self.img_bounding_dec = None

        self.meta = {}
        self.meta.update(meta)

        self.wcs = wcs
        self.catalog = catalog

    @property
    def swcs(self):
        """ Get :py:class:`SimpleWCS` WCS.
        """
        return self._swcs

    @property
    def wcs(self):
        """ Get or set :py:class:`gwcs.WCS`.

        .. note::
            Setting the WCS triggers automatic bounding polygon recalculation.

        """
        if self._swcs is None:
            return None
        return self._swcs.wcs

    @wcs.setter
    def wcs(self, wcs):
        self._swcs = SimpleWCS(wcs, copy=False)

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

    @property
    def pscale(self):
        """ Pixel scale in arcsec assuming CD is in deg/pix.
        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs._pscale

    def all_pix2world(self, x, y):
        """
        Convert pixel coordinates to sky coordinates using full
        (i.e., including distortions) transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.all_pix2world(x, y)

    def all_world2pix(self, ra, dec):
        """
        Convert sky coordinates to image's pixel coordinates using full
        (i.e., including distortions) transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.all_world2pix(ra, dec)

    def wcs_pix2world(self, x, y):
        """
        Convert pixel coordinates to sky coordinates using standard WCS
        transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.wcs_pix2world(x, y)

    def wcs_world2pix(self, ra, dec):
        """
        Convert sky coordinates to image's pixel coordinates using
        standard WCS transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.wcs_world2pix(ra, dec)

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

        ra, dec = self.all_pix2world(borderx, bordery)
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

        ra, dec = convex_hull(x, y, wcs=self.wcs_pix2world)

        if len(ra) == 0:
            # no points
            raise RuntimeError("Unexpected error: Contact sofware developer")

        elif len(ra) == 1:
            # one point. build a small box around it:
            x, y = convex_hull(x, y, wcs=None)

            xv = [x[0] - 0.5, x[0] - 0.5, x[0] + 0.5, x[0] + 0.5, x[0] - 0.5]
            yv = [y[0] - 0.5, y[0] + 0.5, y[0] + 0.5, y[0] - 0.5, y[0] - 0.5]

            ra, dec = self.wcs_pix2world(xv, yv)

        elif len(ra) == 3:
            # two points. build a small box around them:
            x, y = convex_hull(x, y, wcs=None)

            vx = -(y[1] - y[0])
            vy = x[1] - x[0]
            norm = 2.0 * np.sqrt(vx * vx + vy * vy)
            vx /= norm
            vy /= norm

            xv = [x[0] + vx, x[1] + vx, x[1] + vx, x[0] - vx, x[0] + vx]
            yv = [y[0] + vy, y[1] + vy, y[1] + vy, x[0] - vx, x[0] + vx]

            ra, dec = self.wcs_pix2world(xv, yv)

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
            ra, dec = image.all_pix2world(
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

            ra, dec = image.all_pix2world(
                self._catalog['x'][idx], self._catalog['y'][idx]
            )
            self._catalog['RA'][idx] = ra
            self._catalog['DEC'][idx] = dec

    def calc_xyref(self, refcat):
        """
        Compute x- and y-positions of the sources from the image catalog
        in the reference image plane. This create the following
        columns in the catalog's table: ``'xref'`` and ``'yref'``.

        """
        if 'RA' not in self._catalog.colnames or \
           'DEC' not in self._catalog.colnames:
            raise RuntimeError("'recalc_catalog_radec()' should have been run "
                               "prior to calc_xyref()")

        # compute x & y in the reference WCS:
        xref, yref = refcat.all_world2pix(self.catalog['RA'],
                                          self.catalog['DEC'])
        self._catalog['xref'] = table.MaskedColumn(
            xref, name='xref', dtype=np.float64, mask=False
        )
        self._catalog['yref'] = table.MaskedColumn(
            yref, name='yref', dtype=np.float64, mask=False
        )

    def match2ref(self, refcat, minobj=15, searchrad=1.0,
                  searchunits='arcseconds', separation=0.5,
                  use2dhist=True, xoffset=0.0, yoffset=0.0, tolerance=1.0):
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

        searchunits : str, optional
            Units for search radius.

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

        if 'xref' not in colnames or 'yref' not in colnames:
            raise RuntimeError("'calc_xyref()' should have been run prior "
                               "to match2ref()")

        im_xyref = np.asanyarray([self._catalog['xref'],
                                  self._catalog['yref']]).T
        refxy = np.asanyarray([refcat.catalog['xref'],
                               refcat.catalog['yref']]).T

        log.info("Matching sources from '{}' with sources from reference "
                 "{:s} '{}'".format(self.name, 'image', refcat.name))

        # convert tolerance from units of arcseconds to pixels, as needed
        if searchunits == 'arcseconds':
            searchrad /= refcat.pscale

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

        im_xyref = np.asanyarray([self._catalog['xref'],
                                  self._catalog['yref']]).T
        refxy = np.asanyarray([refcat.catalog['xref'],
                               refcat.catalog['yref']]).T

        mask = np.logical_not(self._catalog['matched_ref_id'].mask)
        im_xyref = im_xyref[mask]
        ref_idx = self._catalog['_raw_matched_ref_idx'][mask]
        refxy = refxy[ref_idx]

        fit = linearfit.iter_linear_fit(
            im_xyref, refxy, fitgeom=fitgeom,
            nclip=nclip, sigma=sigma, center=refcat.crpix
        )

        fit['rms_keys'] = self._compute_fit_rms(fit, refcat)
        xy_fit = fit['img_coords'] + fit['resids']
        fit['fit_xy'] = xy_fit
        ra_fit, dec_fit = refcat.wcs_pix2world(xy_fit[0], xy_fit[1])
        fit['fit_RA'] = ra_fit
        fit['fit_DEC'] = dec_fit

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
                imwcs=imcat.wcs,
                imshape=imcat.imshape,
                refwcs=refcat.wcs,
                xsh=xsh,
                ysh=ysh,
                matrix=matrix
            )

    def align_to_ref(self, refcat, minobj=15, searchrad=1.0,
                     searchunits='arcseconds', separation=0.5,
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

        searchunits : str, optional
            Units for search radius.

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
        self.calc_xyref(refcat=refcat)
        self.match2ref(refcat=refcat, minobj=minobj, searchrad=searchrad,
                       searchunits=searchunits, separation=separation,
                       use2dhist=use2dhist, xoffset=xoffset, yoffset=yoffset,
                       tolerance=tolerance)
        fit = self.fit2ref(refcat=refcat, fitgeom=fitgeom,
                           nclip=nclip, sigma=sigma)
        self.apply_affine_to_wcs(refcat=refcat, matrix=fit['fit_matrix'],
                                 xsh=fit['offset'][0], ysh=fit['offset'][1])


class RefCatalog(object):
    """
    An object that holds a reference catalog and an associated WCS and provides
    tools for coordinate convertions using reference WCS as well as
    catalog manipulation and expansion.

    """
    def __init__(self, wcs, catalog, name=None, max_cat_name_len=256):
        # TODO: max_cat_name_len is not implemented yet
        """
        Parameters
        ----------

        wcs : gwcs.WCS
            WCS of the reference frame.

        catalog : astropy.table.Table
            Reference catalog.

            ..note::
                Reference catalogs (:py:class:`~astropy.table.Table`)
                *must* contain *both* ``'RA'`` and ``'DEC'`` columns.

        name : str, None, optional
            Name of the reference catalog.

        max_cat_name_len : int, optional
            Not implemented yet. Reserved for future use.

        """
        self._name = name
        self._max_cat_name_len = max_cat_name_len

        self._catalog = None

        # WCS
        if wcs is None:
            raise ValueError("Reference catalog's wcs' cannot be 'None'")

        self._swcs = SimpleWCS(wcs, copy=False)

        # make sure catalog has RA & DEC
        if catalog is not None:
            self.catalog = catalog

    def _rd2xyref(self, catalog):
        return self.wcs_world2pix(catalog['RA'], catalog['DEC'])

    def _recalc_catalog_xy(self, catalog):
        colnames = catalog.colnames

        if 'xref' not in colnames and 'yref' not in colnames:
            x, y = self._rd2xyref(catalog)
            txy = table.Table([x, y], names=('xref', 'yref'),
                              dtype=(np.float64, np.float64))
            catalog = table.hstack([catalog, txy], join_type='exact')

        elif 'xref' in colnames and 'yref' in colnames:
            x, y = self._rd2xyref(catalog)
            catalog['xref'] = x
            catalog['yref'] = y

        else:
            raise ValueError("Reference catalog must either have or not "
                             "have *both* 'x' and 'y' columns.")

        return catalog

    def recalc_catalog_xy(self):
        """
        Recalculate catalog's ``'xref'`` and ``'yref'`` from catalog's
        ``'RA'`` and ``'DEC'`` values.

        """
        self._catalog = self._recalc_catalog_xy(self.catalog)

    def _check_catalog(self, catalog):
        if catalog is None:
            raise ValueError("Reference catalogs cannot be None")

        if 'RA' not in catalog.colnames or 'DEC' not in catalog.colnames:
            raise KeyError("Reference catalogs *must* contain *both* 'RA' "
                           "and 'DEC' columns.")

    @property
    def wcs(self):
        """ Get :py:class:`gwcs.WCS`.
        """
        if self._swcs is None:
            return None
        return self._swcs.wcs

    @property
    def pscale(self):
        """ Pixel scale in arcsec assuming CD is in deg/pix.
        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.pscale

    @property
    def crpix(self):
        """ Retrieve "CRPIX" of the WCS.
        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.crpix

    @property
    def crval(self):
        """ Retrieve "CRVAL" of the WCS.
        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.crval

    def all_pix2world(self, x, y):
        """
        Convert pixel coordinates to sky coordinates using full
        (i.e., including distortions) transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.all_pix2world(x, y)

    def all_world2pix(self, ra, dec):
        """
        Convert sky coordinates to image's pixel coordinates using full
        (i.e., including distortions) transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.all_world2pix(ra, dec)

    def wcs_pix2world(self, x, y):
        """
        Convert pixel coordinates to sky coordinates using standard WCS
        transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.wcs_pix2world(x, y)

    def wcs_world2pix(self, ra, dec):
        """
        Convert sky coordinates to image's pixel coordinates using
        standard WCS transformations.

        """
        if self._swcs is None:
            raise RuntimeError("WCS has not been set")
        return self._swcs.wcs_world2pix(ra, dec)

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

        self.recalc_catalog_xy()

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
        if self.wcs is None or self.catalog is None:
            return

        x = self.catalog['x']
        y = self.catalog['y']

        ra, dec = convex_hull(x, y, wcs=self.wcs_pix2world)

        if len(ra) == 0:
            # no points
            raise RuntimeError("Unexpected error: Contact sofware developer")

        elif len(ra) == 1:
            # one point. build a small box around it:
            x, y = convex_hull(x, y, wcs=None)

            xv = [x[0] - 0.5, x[0] - 0.5, x[0] + 0.5, x[0] + 0.5, x[0] - 0.5]
            yv = [y[0] - 0.5, y[0] + 0.5, y[0] + 0.5, y[0] - 0.5, y[0] - 0.5]

            ra, dec = self.wcs_pix2world(xv, yv)

        elif len(ra) == 3:
            # two points. build a small box around them:
            x, y = convex_hull(x, y, wcs=None)

            vx = -(y[1] - y[0])
            vy = x[1] - x[0]
            norm = 2.0 * np.sqrt(vx * vx + vy * vy)
            vx /= norm
            vy /= norm

            xv = [x[0] + vx, x[1] + vx, x[1] + vx, x[0] - vx, x[0] + vx]
            yv = [y[0] + vy, y[1] + vy, y[1] + vy, x[0] - vx, x[0] + vx]

            ra, dec = self.wcs_pix2world(xv, yv)

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
        cat = self._recalc_catalog_xy(cat)
        if self._catalog is None:
            self._catalog = cat
        else:
            self._catalog = table.vstack([self.catalog, cat],
                                         join_type='outer')
        self._calc_cat_convex_hull()


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
