"""
A module that provides functions for matching sky in overlapping images.

:Authors: Mihai Cara

"""
import logging
from datetime import datetime
import numpy as np

# LOCAL
from . skyimage import SkyImage, SkyGroup


__all__ = ['match']


__author__ = 'Mihai Cara'


#DEBUG
__local_debug__ = True

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def match(images, skymethod='global+match', match_down=True, subtract=False):
    """
    A function to compute and/or "equalize" sky background in input images.

    .. note::
       Sky matching ("equalization") is possible only for **overlapping**
       images.

    Parameters
    ----------
    images : list of SkyImage or SkyGroup
        A list of of :py:class:`~jwst.skymatch.skyimage.SkyImage` or
        :py:class:`~jwst.skymatch.skyimage.SkyGroup` objects.

    skymethod : {'local', 'global+match', 'global', 'match'}, optional
        Select the algorithm for sky computation:

        * **'local'** : compute sky background values of each input image or
          group of images (members of the same "exposure"). A single sky value
          is computed for each group of images.

          .. note::
            This setting is recommended when regions of overlap between images
            are dominated by "pure" sky (as opposite to extended, diffuse
            sources).

        * **'global'** : compute a common sky value for all input image and
          groups of images. In this setting `match` will compute
          sky values for each input image/group, find the minimum sky value,
          and then it will set (and/or subtract) sky value of each input image
          to this minimum value. This method *may* be
          useful when input images have been already matched.

        * **'match'** : compute differences in sky values between images
          and/or groups in (pair-wise) common sky regions. In this case
          computed sky values will be relative (delta) to the sky computed
          in one of the input images whose sky value will be set to
          (reported to be) 0. This setting will "equalize" sky values between
          the images in large mosaics. However, this method is not recommended
          when used in conjunction with
          `astrodrizzle <http://stsdas.stsci.edu/stsci_python_sphinxdocs_2.13/drizzlepac/astrodrizzle.html>`_
          because it computes relative sky values while `astrodrizzle` needs
          "measured" sky values for median image generation and CR rejection.

        * **'global+match'** : first use **'match'** method to
          equalize sky values between images and then find a minimum
          "global" sky value in all input images.

          .. note::
            This is the *recommended* setting for images
            containing diffuse sources (e.g., galaxies, nebulae)
            covering significant parts of the image.

    match_down : bool, optional
        Specifies whether the sky *differences* should be subtracted from
        images with higher sky values (`match_down` = `True`) to match the
        image with the lowest sky or sky differences should be added to the
        images with lower sky values to match the sky of the image with the
        highest sky value (`match_down` = `False`).

        .. note::
          This setting applies *only* when `skymethod` parameter is
          either `'match'` or `'global+match'`.

    subtract : bool (Default = False)
        Subtract computed sky value from image data.


    Raises
    ------

    TypeError
        The `images` argument must be a Python list of
        :py:class:`~jwst.skymatch.skyimage.SkyImage` and/or
        :py:class:`~jwst.skymatch.skyimage.SkyGroup` objects.


    Notes
    -----

    :py:func:`match` provides new algorithms for sky value computations
    and enhances previously available algorithms used by, e.g.,
    `astrodrizzle <http://stsdas.stsci.edu/stsci_python_sphinxdocs_2.13/drizzlepac/astrodrizzle.html>`_.

    Two new methods of sky subtraction have been introduced (compared to the
    standard ``'local'``): ``'global'`` and ``'match'``, as well as a
    combination of the two -- ``'global+match'``.

    - The ``'global'`` method computes the minimum sky value across *all*
      input images and/or groups. That sky value is then considered to be
      the background in all input images.

    - The ``'match'`` algorithm is somewhat similar to the traditional sky
      subtraction method (`skymethod` = `'local'`) in the sense that it
      measures the sky indipendently in input images (or groups). The major
      differences are that, unlike the traditional method,

        #. ``'match'`` algorithm computes *relative* (delta) sky values with
           regard to the sky in a reference image chosen from the input list
           of images; *and*

        #. Sky statistics is computed only in the part of the image
           that intersects other images.

      This makes ``'match'`` sky computation algorithm particularly useful
      for "equalizing" sky values in large mosaics in which one may have
      only (at least) pair-wise intersection of images without having
      a common intersection region (on the sky) in all images.

      The `'match'` method works in the following way: for each pair
      of intersecting images, an equation is written that
      requires that average surface brightness in the overlapping part of
      the sky be equal in both images. The final system of equations is then
      solved for unknown background levels.

      .. warning::

        Current algorithm is not capable of detecting cases when some subsets
        of intersecting images (from the input list of images) do not intersect
        at all other subsets of intersecting images (except for the simple
        case when *single* images do not intersect any other images). In these
        cases the algorithm will find equalizing sky values for each
        intersecting subset of images and/or groups of images.
        However since these subsets of images do not intersect each other,
        sky will be matched only within each subset and the "inter-subset"
        sky mismatch could be significant.

        Users are responsible for detecting such cases and adjusting processing
        accordingly.

    - The ``'global+match'`` algorithm combines ``'match'`` and
      ``'global'`` methods in order to overcome the limitation of the
      ``'match'`` method described in the note above: it uses ``'global'``
      algorithm to find a baseline sky value common to all input images
      and the ``'match'`` algorithm to "equalize" sky values in the mosaic.
      Thus, the sky value of the "reference" image will be equal to the
      baseline sky value (instead of 0 in ``'match'`` algorithm alone).

    **Remarks:**
      * :py:func:`match` works directly on *geometrically distorted*
        flat-fielded images thus avoiding the need to perform distortion
        correction of input images.

        Initially, the footprint of a chip in an image is aproximated by a
        2D planar rectangle representing the borders of chip's distorted
        image. After applying distortion model to this rectangle and
        progecting it onto the celestial sphere, it is approximated by
        spherical polygons. Footprints of exposures and mosaics are
        computed as unions of such spherical polygons while overlaps
        of image pairs are found by intersecting these spherical polygons.

    **Limitations and Discussions:**
      Primary reason for introducing "sky match" algorithm was to try to
      equalize the sky in large mosaics in which computation of the
      "absolute" sky is difficult due to the presence of large diffuse
      sources in the image. As discussed above, :py:func:`match`
      accomplishes this by comparing "sky values" in a pair of images in the
      overlap region (that is common to both images). Quite obviously the
      quality of sky "matching" will depend on how well these "sky values"
      can be estimated. We use quotation marks around *sky values* because
      for some image "true" background may not be present at all and the
      measured sky may be the surface brightness of large galaxy, nebula, etc.

      In the discussion below we will refer to parameter names in
      :py:class:`~jwst.skymatch.skystatistics.SkyStats` and these
      parameter names may differ from the parameters of the actual `skystat`
      object passed to initializer of the
      :py:class:`~jwst.skymatch.skyimage.SkyImage`.

      Here is a brief list of possible limitations/factors that can affect
      the outcome of the matching (sky subtraction in general) algorithm:

      * Since sky subtraction is performed on *flat-fielded* but
        *not distortion corrected* images, it is important to keep in mind
        that flat-fielding is performed to obtain uniform surface brightness
        and not flux. This distinction is important for images that have
        not been distortion corrected. As a consequence, it is advisable that
        point-like sources be masked through the user-supplied mask files.
        Values different from zero in user-supplied masks indicate "good" data
        pixels. Alternatively, one can use `upper` parameter to limit the use
        of bright objects in sky computations.

      * Normally, distorted flat-fielded images contain cosmic rays. This
        algorithm does not perform CR cleaning. A possible way of minimizing
        the effect of the cosmic rays on sky computations is to use
        clipping (`nclip` > 0) and/or set `upper` parameter to a value
        larger than most of the sky background (or extended source) but
        lower than the values of most CR pixels.

      * In general, clipping is a good way of eliminating "bad" pixels:
        pixels affected by CR, hot/dead pixels, etc. However, for
        images with complicated backgrounds (extended galaxies, nebulae,
        etc.), affected by CR and noise, clipping process may mask different
        pixels in different images. If variations in the background are
        too strong, clipping may converge to different sky values in
        different images even when factoring in the "true" difference
        in the sky background between the two images.

      * In general images can have different "true" background values
        (we could measure it if images were not affected by large diffuse
        sources). However, arguments such as `lower` and `upper` will
        apply to all images regardless of the intrinsic differences
        in sky levels.

    """
    function_name = match.__name__

    # Time it
    runtime_begin = datetime.now()

    log.info(" ")
    log.info("***** {:s}.{:s}() started on {}"
             .format(__name__, function_name, runtime_begin))
    log.info(" ")

    # check sky method:
    skymethod = skymethod.lower()
    if skymethod not in ['local', 'global', 'match', 'global+match']:
        raise ValueError("Unsupported 'skymethod'. Valid values are: "
                         "'local', 'global', 'match', or 'global+match'")
    do_match = 'match' in skymethod
    do_global = 'global' in skymethod
    show_old = subtract

    log.info("Sky computation method: '{}'".format(skymethod))
    if do_match:
        log.info("Sky matching direction: {:s}"
                 .format("DOWN" if match_down else "UP"))

    log.info("Sky subtraction from image data: {:s}"
             .format("ON" if subtract else "OFF"))

    # check that input file name is a list of either SkyImage or SkyGroup:
    nimages = 0
    for img in images:
        if isinstance(img, SkyImage):
            nimages += 1
        elif isinstance(img, SkyGroup):
            nimages += len(img)
        else:
            raise TypeError("Each element of the 'images' must be either a "
                            "'SkyImage' or a 'SkyGroup'")

    if nimages == 0:
        raise ValueError("Argument 'images' must contain at least one image")

    log.debug(
        "Total number of images to be sky-subtracted and/or matched: {:d}"
        .format(nimages)
    )

    ###################################
    ##   Print conversion factors:   ##
    ###################################
    log.debug(" ")
    log.debug("----  Image data conversion factors:")

    for img in images:
        img_type = 'Image' if isinstance(img, SkyImage) else 'Group'

        if img_type == 'Group':
            log.debug("   *  Group ID={}. Conversion factors:".format(img.id))
            for im in img:
                log.debug("      - Image ID={}. Conversion factor = {:G}"
                          .format(im.id, im.convf))
        else:
            log.debug("   *  Image ID={}. Conversion factor = {:G}"
                      .format(img.id, img.convf))

    ###############################################################
    ## 1. Method: "match" (or "global+match").                   ##
    ##    Find sky "deltas" that will match sky across all       ##
    ##    (intersecting) images.                                 ##
    ###############################################################
    if do_match:

        log.info(" ")
        log.info("----  Computing differences in sky values in "
                 "overlapping regions.")

        # find "optimum" sky changes:
        sky_deltas = _find_optimum_sky_deltas(images, apply_sky=not subtract)
        sky_good = np.isfinite(sky_deltas)

        if np.any(sky_good):
            # match sky "Up" or "Down":
            if match_down:
                refsky = np.amin(sky_deltas[sky_good])
            else:
                refsky = np.amax(sky_deltas[sky_good])
            sky_deltas[sky_good] -= refsky

        # convert to Python list and replace numpy.nan with None
        sky_deltas = [
            float(skd) if np.isfinite(skd) else None for skd in sky_deltas
        ]

        _apply_sky(images, sky_deltas, False, subtract, show_old)
        show_old = True

    ###############################################################
    ## 2. Method: "local". Compute the minimum sky background    ##
    ##    value in each sky group/image.                         ##
    ##    This is an improved (use of masks) replacement         ##
    ##    for the classical 'subtract' used by astrodrizzle.  ##
    ##                                                           ##
    ##    NOTE: incompatible with "match"-containing             ##
    ##          'skymethod' modes.                               ##
    ##                                                           ##
    ## 3. Method: "global". Compute the minimum sky background   ##
    ##    value *across* *all* sky line members.                 ##
    ###############################################################
    if do_global or not do_match:

        log.info(" ")
        if do_global:
            minsky = None
            log.info("----  Computing \"global\" sky - smallest sky value "
                     "across *all* input images.")
        else:
            log.info("----  Sky values computed per image and/or image "
                     "groups.")

        sky_deltas = []
        for img in images:
            sky = img.calc_sky(delta=not subtract)[0]
            sky_deltas.append(sky)
            if do_global and (minsky is None or sky < minsky):
                minsky = sky

        if do_global:
            log.info(" ")
            if minsky is None:
                log.warning("   Unable to compute \"global\" sky value")
            sky_deltas = len(sky_deltas) * [minsky]
            log.info("   \"Global\" sky value correction: {} "
                     "[not converted]".format(minsky))

        if do_match:
            log.info(" ")
            log.info("----  Final (match+global) sky for:")

        _apply_sky(images, sky_deltas, do_global, subtract, show_old)

    # log running time:
    runtime_end = datetime.now()
    log.info(" ")
    log.info("***** {:s}.{:s}() ended on {}"
             .format(__name__, function_name, runtime_end))
    log.info("***** {:s}.{:s}() TOTAL RUN TIME: {}"
             .format(__name__, function_name, runtime_end - runtime_begin))
    log.info(" ")


def _apply_sky(images, sky_deltas, do_global, do_skysub, show_old):
    for img, sky in zip(images, sky_deltas):
        is_group = not isinstance(img, SkyImage)

        if do_global:
            if sky is None:
                valid = img[0].is_sky_valid if is_group else img.is_sky_valid
                sky = 0.0
            else:
                valid = True

        else:
            valid = sky is not None
            if not valid:
                log.warning("   *  {:s} ID={}: Unable to compute sky value"
                            .format('Group' if is_group else 'Image', img.id))
                sky = 0.0

        if is_group:
            # apply sky change:
            old_img_sky = [im.sky for im in img]
            if do_skysub:
                for im in img:
                    im.image -= sky
            img.sky += sky
            new_img_sky = [im.sky for im in img]

            # log sky values:
            log.info("   *  Group ID={}. Sky background of "
                     "component images:".format(img.id))

            for im, old_sky, new_sky in zip(img, old_img_sky, new_img_sky):
                c = 1.0 / im.convf
                if show_old:
                    log.info("      - Image ID={}. Sky background: {:G} "
                             "(old={:G}, delta={:G})"
                             .format(im.id, c * new_sky, c * old_sky, c * sky))
                else:
                    log.info("      - Image ID={}. Sky background: {:G}"
                             .format(im.id, c * new_sky))

                im.is_sky_valid = valid

        else:
            # apply sky change:
            old_sky = img.sky
            if do_skysub:
                img.image -= sky
            img.sky += sky
            new_sky = img.sky

            # log sky values:
            c = 1.0 / img.convf
            if show_old:
                log.info("   *  Image ID={}. Sky background: {:G} "
                         "(old={:G}, delta={:G})"
                         .format(img.id, c * new_sky, c * old_sky, c * sky))
            else:
                log.info("   *  Image ID={}. Sky background: {:G}"
                         .format(img.id, c * new_sky))

            img.is_sky_valid = valid


# TODO: due to a bug in the sphere package, see
#       https://github.com/spacetelescope/sphere/issues/74
#       intersections with polygons formed as union does not work.
#       For this reason I re-implement '_overlap_matrix' below with
#       a workaround for the bug.
#       The original implementation should be uncommented once the bug
#       is fixed.
#

# Original version:
#def _overlap_matrix(images, apply_sky=True):
    ##TODO: to improve performance, the nested loops could be parallelized
    ## since _calc_sky() here can be called independently from previous steps.
    #ns = len(images)
    #A = np.zeros((ns, ns), dtype=float)
    #W = np.zeros((ns, ns), dtype=float)
    #for i in range(ns):
        #for j in range(i+1, ns):
            #overlap = images[i].intersection(images[j])
            #s1, w1, area1 = images[i].calc_sky(overlap=overlap, delta=apply_sky)
            #s2, w2, area2 = images[j].calc_sky(overlap=overlap, delta=apply_sky)
            #if area1 == 0.0 or area2 == 0.0 or s1 is None or s2 is None:
                #continue
            #A[j,i] = s1
            #W[j,i] = w1
            #A[i,j] = s2
            #W[i,j] = w2
    #return A, W

# bug workaround version:
def _overlap_matrix(images, apply_sky=True):
    #TODO: to improve performance, the nested loops could be parallelized
    # since _calc_sky() here can be called independently from previous steps.
    ns = len(images)
    A = np.zeros((ns, ns), dtype=float)
    W = np.zeros((ns, ns), dtype=float)

    for i in range(ns):
        for j in range(i + 1, ns):
            s1, w1, area1 = images[i].calc_sky(
                overlap=images[j], delta=apply_sky
            )

            s2, w2, area2 = images[j].calc_sky(
                overlap=images[i], delta=apply_sky
            )

            if area1 == 0.0 or area2 == 0.0 or s1 is None or s2 is None:
                continue

            A[j, i] = s1
            W[j, i] = w1
            A[i, j] = s2
            W[i, j] = w2

    return A, W


def _find_optimum_sky_deltas(images, apply_sky=True):
    ns = len(images)
    A, W = _overlap_matrix(images, apply_sky=apply_sky)

    def is_valid(i, j):
        return (W[i, j] > 0 and W[j, i] > 0)

    # We need to know how many "non-trivial" (at least for now... - we will
    # compute rank later) equations can be built so that we know the
    # shape of the arrays that need to be created...
    # NOTE: for now use only pairs that *both* have weights > 0 (but a
    # different scenario when only one image has a valid weight can be
    # considered):
    neq = 0
    for i in range(ns):
        for j in range(i + 1, ns):
            if is_valid(i, j):
                neq += 1

    # average weights:
    Wm = 0.5 * (W + W.T)

    # create arrays for coefficients and free terms:
    K = np.zeros((neq, ns), dtype=float)
    F = np.zeros(neq, dtype=float)
    invalid = (ns) * [True]

    # now process intersections between the rest of the images:
    ieq = 0
    for i in range(0, ns):
        for j in range(i + 1, ns):
            if is_valid(i, j):
                K[ieq, i] = Wm[i, j]
                K[ieq, j] = -Wm[i, j]
                F[ieq] = Wm[i, j] * (A[j, i] - A[i, j])
                invalid[i] = False
                invalid[j] = False
                ieq += 1

    try:
        rank = np.linalg.matrix_rank(K, 1.0e-12)
    except np.linalg.LinAlgError:
        log.warning("Unable to compute sky: No valid data in common "
                    "image areas")
        deltas = np.full(ns, np.nan, dtype=np.float)
        return deltas

    if rank < ns - 1:
        log.warning("There are more unknown sky values ({}) to be solved for"
                    .format(ns))
        log.warning("than there are independent equations available "
                    "(matrix rank={}).".format(rank))
        log.warning("Sky matching (delta) values will be computed only for")
        log.warning("a subset (or more independent subsets) of input images.")
    invK = np.linalg.pinv(K, rcond=1.0e-12)

    deltas = np.dot(invK, F)
    deltas[np.asarray(invalid, dtype=bool)] = np.nan
    return deltas
