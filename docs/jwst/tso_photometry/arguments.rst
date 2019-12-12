Step Arguments
==============

The ``tso_photometry`` step has one step-specific argument:

*  ``--save_catalog``

If ``save_catalog`` is set to True (the default is False),
the output table of times and photometry will be written to an ecsv file
with suffix "phot".

Note that when this step is run as part of the
:ref:`calwebb_tso3 <calwebb_tso3>` pipeline,
the ``save_catalog`` argument should *not* be set, because the output
catalog will always be saved by the pipeline module itself.  The
``save_catalog`` argument is useful only when the ``tso_photometry`` step
is run standalone.
