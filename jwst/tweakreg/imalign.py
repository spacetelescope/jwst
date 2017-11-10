"""
A module that provides functions for matching sky in overlapping images.

:Authors: Mihai Cara (contact: help@stsci.edu)

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

# STDLIB
import logging
from datetime import datetime

# THIRD PARTY
import numpy as np

# LOCAL
from . wcsimage import *
from . wcsutils import create_ref_wcs


__all__ = ['align', 'overlap_matrix', 'max_overlap_pair', 'max_overlap_image']

__version__ = '0.8.0'
__vdate__ = '17-April-2016'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def align(imcat, refcat=None, enforce_user_order=True,
          expand_refcat=False, minobj=None, searchrad=1.0,
          searchunits='arcseconds',
          use2dhist=True, separation=0.5, tolerance=1.0,
          xoffset=0.0, yoffset=0.0,
          fitgeom='general', nclip=3, sigma=3.0):
    """
    Align (groups of) images by adjusting the parameters of their WCS based on
    fits between matched sources in these images and a reference catalog which
    may be automatically created from one of the input images.


    Parameters
    ----------
    imcat : list of WCSImageCatalog or WCSGroupCatalog
        A list of `WCSImageCatalog` or `WCSGroupCatalog` objects whose WCS
        should be adjusted. The difference between `WCSImageCatalog` and
        `WCSGroupCatalog` is that the later is used to find a *single* fit
        to all sources in all component images *simultaneously*. This fit is
        later applied to each component image's WCS.

        .. warning::
            This function modifies the WCS of the input images provided
            through the `imcat` parameter. On return, each input image WCS
            will be updated with an "aligned" version.

    refcat : RefCatalog, optional
        A `RefCatalog` object that contains a catalog of reference sources
        as well as (optionally) a valid reference WCS. When `refcat` is `None`,
        a reference catalog will be created from one of the input (group of)
        images.

    enforce_user_order : bool, optional
        Specifies whether images should be aligned in the order specified in
        the `file` input parameter or `align` should optimize the order
        of alignment by intersection area of the images. Default value (`True`)
        will align images in the user specified order, except when some images
        cannot be aligned in which case `align` will optimize the image
        alignment order. Alignment order optimization is available *only*
        when `expand_refcat` = `True`.

    expand_refcat : bool, optional
        Specifies whether to add new sources from just matched images to
        the reference catalog to allow next image to be matched against an
        expanded reference catalog. By delault, the reference catalog is not
        being expanded.

    minobj : int, None, optional
        Minimum number of identified objects from each input image to use
        in matching objects from other images. If the default `None` value is
        used then `align` will automatically deternmine the minimum number
        of sources from the value of the `fitgeom` parameter.

    searchrad : float, optional
        The search radius for a match.

    searchunits : str, optional
        Units for search radius.

    use2dhist : bool, optional
        Use 2D histogram to find initial offset?

    separation : float, optional
        The  minimum  separation for sources in the input and reference
        catalogs in order to be considered to be disctinct sources.
        Objects closer together than 'separation' pixels
        are removed from the input and reference coordinate lists prior
        to matching. This parameter gets passed directly to
        :py:func:`~stsci.stimage.xyxymatch` for use in matching the object
        lists from each image with the reference image's object list.

    tolerance : float, optional
        The matching tolerance in pixels after applying an initial solution
        derived from the 'triangles' algorithm.  This parameter gets passed
        directly to :py:func:`~stsci.stimage.xyxymatch` for use in
        matching the object lists from each image with the reference image's
        object list.

    xoffset : float, optional
        Initial estimate for the offset in X between the images and the
        reference frame. This offset will be used for all input images
        provided. This parameter is ignored when `use2dhist` is `True`.

    yoffset : float (Default = 0.0)
        Initial estimate for the offset in Y between the images and the
        reference frame. This offset will be used for all input images
        provided. This parameter is ignored when `use2dhist` is `True`.

    fitgeom : {'shift', 'rscale', 'general'}, optional
        The fitting geometry to be used in fitting the matched object lists.
        This parameter is used in fitting the offsets, rotations and/or scale
        changes from the matched object lists. The 'general' fit geometry
        allows for independent scale and rotation for each axis.

    nclip : int, optional
        Number (a non-negative integer) of clipping iterations in fit.

    sigma : float, optional
        Clipping limit in sigma units.

    """

    function_name = align.__name__

    # Time it
    runtime_begin = datetime.now()

    log.info(" ")
    log.info("***** {:s}.{:s}() started on {}"
             .format(__name__, function_name, runtime_begin))
    log.info("      Version {} ({})".format(__version__, __vdate__))
    log.info(" ")

    # check fitgeom:
    fitgeom = fitgeom.lower()
    if fitgeom not in ['shift', 'rscale', 'general']:
        raise ValueError("Unsupported 'fitgeom'. Valid values are: "
                         "'shift', 'rscale', or 'general'")

    # check searchunits:
    searchunits = searchunits.lower()
    if searchunits not in ['arcseconds', 'pixel']:
        raise ValueError("Unsupported 'searchunits'. Valid values are: "
                         "'arcseconds', or 'pixel'")

    if minobj is None:
        if fitgeom == 'general':
            minobj = 3
        elif fitgeom == 'rscale':
            minobj = 2
        else:
            minobj = 1

    # check that we have enough input images:
    if refcat is None or refcat.catalog is None:
        if len(imcat) < 2:
            raise ValueError("Too few input catalogs")
    else:
        if len(imcat) == 0:
            raise ValueError("Too few input catalogs")

    # make sure each input item is a WCSGroupCatalog:
    old_imcat = imcat
    imcat = []
    for img in old_imcat:
        if img.catalog is None:
            raise ValueError("Each image/catalog must have a valid catalog")
        if isinstance(img, WCSImageCatalog):
            imcat.append(WCSGroupCatalog(img, name=img.name))
        elif isinstance(img, WCSGroupCatalog):
            imcat.append(img)
        else:
            raise TypeError("Each input element of 'images' must be a "
                            "'WCSGroupCatalog'")

    # create reference WCS if needed:
    if refcat is None:
        refwcs, refshape = create_ref_wcs(imcat)
        refcat = RefCatalog(wcs=refwcs, catalog=None)

    elif refcat.wcs is None:  # TODO: currently this is not supported
        refwcs, refshape = create_ref_wcs(imcat)
        refcat.wcs = refwcs

    # get the first image to be aligned and
    # create reference catalog if needed:
    if refcat.catalog is None:
        # create reference catalog:
        ref_imcat, current_imcat = max_overlap_pair(
            images=imcat,
            expand_refcat=expand_refcat,
            enforce_user_order=enforce_user_order
        )
        log.info("Selected image '{}' as reference image"
                 .format(ref_imcat.name))
        refcat.expand_catalog(ref_imcat.catalog)

    else:
        # find the first image to be aligned:
        current_imcat = max_overlap_image(
            refimage=refcat,
            images=imcat,
            expand_refcat=expand_refcat,
            enforce_user_order=enforce_user_order
        )

    aligned_imcat = []

    while current_imcat is not None:
        log.info("Aligning image catalog '{}' to the reference catalog."
                 .format(current_imcat.name))

        current_imcat.align_to_ref(
            refcat=refcat,
            minobj=minobj,
            searchrad=searchrad,
            searchunits=searchunits,
            separation=separation,
            use2dhist=use2dhist,
            xoffset=xoffset,
            yoffset=yoffset,
            tolerance=tolerance,
            fitgeom=fitgeom,
            nclip=nclip,
            sigma=sigma
        )

        aligned_imcat.append(current_imcat)

        # add unmatched sources to the reference catalog:
        if expand_refcat:
            refcat.expand_catalog(current_imcat.get_unmatched_cat())
            log.info("Added unmatched sources from '{}' to the reference "
                     "catalog.".format(current_imcat.name))

        # find the next image to be aligned:
        current_imcat = max_overlap_image(
            refimage=refcat,
            images=imcat,
            expand_refcat=expand_refcat,
            enforce_user_order=enforce_user_order
        )

    # log running time:
    runtime_end = datetime.now()
    log.info(" ")
    log.info("***** {:s}.{:s}() ended on {}"
             .format(__name__, function_name, runtime_end))
    log.info("***** {:s}.{:s}() TOTAL RUN TIME: {}"
             .format(__name__, function_name, runtime_end - runtime_begin))
    log.info(" ")

    return aligned_imcat


def overlap_matrix(images):
    nimg = len(images)
    m = np.zeros((nimg, nimg), dtype=np.float)
    for i in range(nimg):
        for j in range(i + 1, nimg):
            area = images[i].intersection_area(images[j])
            m[j, i] = area
            m[i, j] = area
    return m


def max_overlap_pair(images, expand_refcat, enforce_user_order):
    nimg = len(images)

    if nimg == 0:
        return None, None

    elif nimg == 1:
        return images[0], None

    elif nimg == 2 or not expand_refcat or enforce_user_order:
        # for the special case when only two images are provided
        # return (refimage, image) in the same order as provided in 'images'.
        # Also, when ref. catalog is static - revert to old tweakreg behavior
        im1 = images.pop(0)  # reference image
        im2 = images.pop(0)
        return (im1, im2)

    m = overlap_matrix(images)
    imgs = [f.name for f in images]
    n = m.shape[0]
    index = m.argmax()
    i = index / n
    j = index % n
    si = np.sum(m[i])
    sj = np.sum(m[:, j])

    if si < sj:
        c = j
        j = i
        i = c

    if i < j:
        j -= 1

    im1 = images.pop(i)  # reference image
    im2 = images.pop(j)

    # Sort the remaining of the input list of images by overlap area
    # with the reference image (in decreasing order):
    row = m[i]
    row = np.delete(row, i)
    row = np.delete(row, j)
    sorting_indices = np.argsort(row)[::-1]
    images_arr = np.asarray(images)[sorting_indices]
    while len(images) > 0:
        del images[0]
    for k in range(images_arr.shape[0]):
        images.append(images_arr[k])

    return (im1, im2)


def max_overlap_image(refimage, images, expand_refcat, enforce_user_order):
    nimg = len(images)
    if len(images) < 1:
        return None

    if not expand_refcat or enforce_user_order:
        # revert to old tweakreg behavior
        return images.pop(0)

    area = np.zeros(nimg, dtype=np.float)
    for i in range(nimg):
        area[i] = refimage.intersection_area(images[i])

    # Sort the remaining of the input list of images by overlap area
    # with the reference image (in decreasing order):
    sorting_indices = np.argsort(area)[::-1]
    images_arr = np.asarray(images)[sorting_indices]
    while len(images) > 0:
        del images[0]
    for k in range(images_arr.shape[0]):
        images.append(images_arr[k])

    return images.pop(0)
