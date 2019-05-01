"""
A module that provides functions for "aligning" images: specifically, it
provides functions for computing corrections to image WCS so that
images catalogs "align" to the reference catalog *on the sky*.

:Authors: Mihai Cara (contact: help@stsci.edu)


"""

# STDLIB
import logging
from datetime import datetime

# THIRD PARTY
import numpy as np

# LOCAL
from . wcsimage import (WCSGroupCatalog, WCSImageCatalog, RefCatalog)
from . linearfit import SingularMatrixError, NotEnoughPointsError


__all__ = ['align', 'overlap_matrix', 'max_overlap_pair', 'max_overlap_image',
           'NotEnoughCatalogsError']

__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class NotEnoughCatalogsError(Exception):
    """ An error class used to report when there are too few catalogs. """
    pass


def align(imcat, refcat=None, enforce_user_order=True,
          expand_refcat=False, minobj=None, searchrad=1.0,
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

    use2dhist : bool, optional
        Use 2D histogram to find initial offset?

    separation : float, optional
        The  minimum  separation for sources in the input and reference
        catalogs in order to be considered to be disctinct sources.
        Objects closer together than 'separation' *in arcseconds*
        are removed from the input and reference coordinate lists prior
        to matching. This parameter gets passed directly to
        :py:func:`~stsci.stimage.xyxymatch` for use in matching the object
        lists from each image with the reference image's object list.

    tolerance : float, optional
        The matching tolerance *in arcseconds* after applying an initial
        solution derived from the 'triangles' algorithm.  This parameter gets
        passed directly to :py:func:`~stsci.stimage.xyxymatch` for use in
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
    log.info(" ")

    # check fitgeom:
    fitgeom = fitgeom.lower()
    if fitgeom not in ['shift', 'rscale', 'general']:
        raise ValueError("Unsupported 'fitgeom'. Valid values are: "
                         "'shift', 'rscale', or 'general'")

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
            raise NotEnoughCatalogsError(
                "At least two catalogs are needed for alignment."
            )
    else:
        if len(imcat) == 0:
            raise NotEnoughCatalogsError(
                "No catalogs provided for alignment."
            )

    # make sure each input item is a WCSGroupCatalog:
    old_imcat = imcat
    imcat = []
    skipped_imcat = []
    for img in old_imcat:
        if img.catalog is None:
            raise ValueError("Each image/catalog must have a valid catalog")

        if isinstance(img, WCSImageCatalog):
            group = WCSGroupCatalog(img, name=img.name)
        elif isinstance(img, WCSGroupCatalog):
            group = img
        else:
            raise TypeError("Each input element of 'images' must be either a"
                            "'WCSImageCatalog' or a 'WCSGroupCatalog'.")

        if len(group.catalog) < minobj:
            for i in group:
                log.warning("Not enough sources to align image {:s}. "
                            "This image will not be aligned."
                            .format(i.meta['image_model'].meta.filename))
                i.meta['image_model'].meta.cal_step.tweakreg = "SKIPPED"
            skipped_imcat.append(group)

        else:
            imcat.append(group)

    if len(imcat) < 2:
        log.warning(
            "Not enough catalogs containing minimum required number of "
            "sources. Stopping image alignment."
        )
        if len(imcat) == 1:
            for i in imcat[0]:
                i.meta['image_model'].meta.cal_step.tweakreg = "SKIPPED"
        skipped_imcat += imcat
        # log running time:
        runtime_end = datetime.now()
        log.info(" ")
        log.info("***** {:s}.{:s}() ended on {}"
                 .format(__name__, function_name, runtime_end))
        log.info("***** {:s}.{:s}() TOTAL RUN TIME: {}"
                 .format(__name__, function_name, runtime_end - runtime_begin))
        log.info(" ")

        return [], skipped_imcat

    # get the first image to be aligned and
    # create reference catalog if needed:
    if refcat is None or refcat.catalog is None:
        # create reference catalog:
        ref_imcat, current_imcat = max_overlap_pair(
            images=imcat,
            enforce_user_order=enforce_user_order or not expand_refcat
        )
        log.info("Selected image '{}' as reference image"
                 .format(ref_imcat.name))
        refcat = RefCatalog(ref_imcat.catalog, name=ref_imcat.name,
                            footprint_tol=tolerance)

    else:
        # find the first image to be aligned:
        current_imcat = max_overlap_image(
            refimage=refcat,
            images=imcat,
            enforce_user_order=enforce_user_order or not expand_refcat
        )

    aligned_imcat = []

    while current_imcat is not None:
        log.info("Aligning image catalog '{}' to the reference catalog."
                 .format(current_imcat.name))

        try:
            current_imcat.align_to_ref(
                refcat=refcat,
                minobj=minobj,
                searchrad=searchrad,
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

            for i in current_imcat:
                i.meta['image_model'].meta.cal_step.tweakreg = "COMPLETE"

            # add unmatched sources to the reference catalog:
            if expand_refcat:
                refcat.expand_catalog(current_imcat.get_unmatched_cat())
                log.info("Added unmatched sources from '{}' to the reference "
                         "catalog.".format(current_imcat.name))

        except (SingularMatrixError, NotEnoughPointsError) as e:
            for i in current_imcat:
                log.warning(
                    "Image {:s} was not aligned. Reported exception: {:s}"
                    .format(i.meta['image_model'].meta.filename, e.args[0])
                )
                i.meta['image_model'].meta.cal_step.tweakreg = "SKIPPED"
            skipped_imcat.append(current_imcat)

        # find the next image to be aligned:
        current_imcat = max_overlap_image(
            refimage=refcat,
            images=imcat,
            enforce_user_order=enforce_user_order or not expand_refcat
        )

    # log running time:
    runtime_end = datetime.now()
    log.info(" ")
    log.info("***** {:s}.{:s}() ended on {}"
             .format(__name__, function_name, runtime_end))
    log.info("***** {:s}.{:s}() TOTAL RUN TIME: {}"
             .format(__name__, function_name, runtime_end - runtime_begin))
    log.info(" ")

    return aligned_imcat, skipped_imcat


def overlap_matrix(images):
    """
    Compute overlap matrix: non-diagonal elements (i,j) of this matrix are
    absolute value of the area of overlap on the sky between i-th input image
    and j-th input image.

    .. note::
        The diagonal of the returned overlap matrix is set to ``0.0``, i.e.,
        this function does not compute the area of the footprint of a single
        image on the sky.

    Parameters
    ----------

    images : list of WCSImageCatalog, WCSGroupCatalog, or RefCatalog
        A list of catalogs that implement :py:meth:`intersection_area` method.

    Returns
    -------
    m : numpy.ndarray
        A `numpy.ndarray` of shape ``NxN`` where ``N`` is equal to the
        number of input images. Each non-diagonal element (i,j) of this matrix
        is the absolute value of the area of overlap on the sky between i-th
        input image and j-th input image. Diagonal elements are set to ``0.0``.

    """
    nimg = len(images)
    m = np.zeros((nimg, nimg), dtype=np.float)
    for i in range(nimg):
        for j in range(i + 1, nimg):
            area = images[i].intersection_area(images[j])
            m[j, i] = area
            m[i, j] = area
    return m


def max_overlap_pair(images, enforce_user_order):
    """
    Return a pair of images with the largest overlap.

    .. warning::
        Returned pair of images is "poped" from input ``images`` list and
        therefore on return ``images`` will contain a smaller number of
        elements.

    Parameters
    ----------

    images : list of WCSImageCatalog, WCSGroupCatalog, or RefCatalog
        A list of catalogs that implement :py:meth:`intersection_area` method.

    enforce_user_order : bool
        When ``enforce_user_order`` is `True`, a pair of images will be
        returned **in the same order** as they were arranged in the ``images``
        input list. That is, image overlaps will be ignored.

    Returns
    -------
    (im1, im2)
        Returns a tuple of two images - elements of input ``images`` list.
        When ``enforce_user_order`` is `True`, images are returned in the
        order in which they appear in the input ``images`` list. When the
        number of input images is smaller than two, ``im1`` and ``im2`` may
        be `None`.


    """
    nimg = len(images)

    if nimg == 0:
        return None, None

    elif nimg == 1:
        return images[0], None

    elif nimg == 2 or enforce_user_order:
        # for the special case when only two images are provided
        # return (refimage, image) in the same order as provided in 'images'.
        # Also, when ref. catalog is static - revert to old tweakreg behavior
        im1 = images.pop(0)  # reference image
        im2 = images.pop(0)
        return (im1, im2)

    m = overlap_matrix(images)
    n = m.shape[0]
    index = m.argmax()
    i = index // n
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


def max_overlap_image(refimage, images, enforce_user_order):
    """
    Return the image from the input ``images`` list that has the largest
    overlap with the ``refimage`` image.

    .. warning::
        Returned image of images is "poped" from input ``images`` list and
        therefore on return ``images`` will contain a smaller number of
        elements.

    Parameters
    ----------

    images : list of WCSImageCatalog, or WCSGroupCatalog
        A list of catalogs that implement :py:meth:`intersection_area` method.

    enforce_user_order : bool
        When ``enforce_user_order`` is `True`, returned image is the first
        image from the ``images`` input list regardless ofimage overlaps.

    Returns
    -------
    image: WCSImageCatalog, WCSGroupCatalog, or None
        Returns an element of input ``images`` list. When input list is
        empty - `None` is returned.

    """
    if len(images) < 1:
        return None

    if enforce_user_order:
        # revert to old tweakreg behavior
        return images.pop(0)

    area = [refimage.intersection_area(im) for im in images]
    idx = np.argmax(area)
    return images.pop(idx)
