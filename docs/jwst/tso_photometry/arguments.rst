Step Arguments
==============

The ``tso_photometry`` step has the following optional arguments:

*  ``--save_catalog`` (boolean, default=False)

    If ``save_catalog`` is set to True, the output table of times and photometry
    will be written to an ecsv file with suffix "phot".

    Note that when this step is run as part of the
    :ref:`calwebb_tso3 <calwebb_tso3>` pipeline,
    the ``save_catalog`` argument should *not* be set, because the output
    catalog will always be saved by the pipeline module itself.  The
    ``save_catalog`` argument is useful only when the ``tso_photometry`` step
    is run standalone.

* ``--radius`` (float, default=3.0)

    A floating point value for the photometric aperture radius in pixels.

* ``--radius_inner`` (float, default=4.0)

   A floating point value for the background annulus inner radius in pixels.

* ``--radius_outer`` (float, default=5.0)

   A floating point value for the background annulus outer radius in pixels.

* ``--centroid_source`` (boolean, default=True)

    If True, the source position will be derived from a centroid near the
    planned position before computing the photometry.  If False, the
    planned source position will be used directly.

* ``--search_box_width`` (integer, default=41)

    If ``centroid_source`` is True, then this value is used as the full width
    of the box used for an initial search for the source.  The value must be
    an odd integer.

* ``--fit_box_width`` (integer, default=11)

    If ``centroid_source`` is True, then this value is used as the full width
    of the box used for a final fit to the source.  The value must be
    an odd integer.

* ``--moving_centroid`` (boolean, default=False)

    If ``centroid_source`` is True, a centroid is fit to each integration.
    If ``moving_centroid`` is True, these centroids are used directly to compute
    the photometry for each integration, such that the source position is
    allowed to vary over the exposure.

    If ``moving_centroid`` is False, the median of the computed centroid positions
    is taken to be the "true" source position, and is used for every integration in the exposure.
    That is, the source position is NOT allowed to vary over the exposure.
