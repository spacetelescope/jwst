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

    When ``centroid_source`` is True, a centroid is fit to each integration.
    If ``moving_centroid`` is True, these fit values are used directly to compute
    the photometry for each integration. In this case, the source position is
    allowed to vary over the exposure.

    If ``moving_centroid`` is False, the source position used for all integrations
    is the median of all the computed centroid positions.  In this case, the same
    source position is used for the whole exposure.
