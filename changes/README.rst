Writing news fragments for the change log
#########################################

This ``changes/`` directory contains "news fragments": small reStructuredText (``.rst``) files describing a change in a few sentences.
When making a release, run ``towncrier build --version <VERSION>`` to consume existing fragments in ``changes/`` and insert them as a full change log entry at the top of ``CHANGES.rst`` for the released version.

News fragment filenames consist of the pull request number and the change log category (see below).
A single change can have more than one news fragment, if it spans multiple categories.
For example, `pull request #10042 <https://github.com/spacetelescope/jwst/pull/10042>`_ has two change notes: ``tso_photometry`` and ``breaking``:

.. code-block::

  10028.other.rst
  10042.breaking.rst
  10042.tso_photometry.rst
  10043.tweakreg.rst
  10139.docs.rst

Change log categories
*********************

Typically, changes to the JWST calibration pipeline will be to one or more steps and / or shared pipeline modules.
Make a news fragment for every relevant category affected by your change.

If a change breaks **step-level or public API** (`as defined in the docs <https://jwst.readthedocs.io/en/latest/jwst/user_documentation/more_information.html#api-public-vs-private>`_), also add a ``<PR#>.breaking.rst`` fragment describing what changes user may need to make to their code:

- ``<PR#>.breaking.rst``: Also add this fragment if your change breaks **step-level or public API**

Shared Pipeline Modules
=======================

- ``<PR#>.stpipe.rst``
- ``<PR#>.datamodels.rst``
- ``<PR#>.scripts.rst``
- ``<PR#>.set_telescope_pointing.rst``
- ``<PR#>.pipeline.rst``
- ``<PR#>.associations.rst``

Stage 1 Steps
=============

- ``<PR#>.group_scale.rst``
- ``<PR#>.dq_init.rst``
- ``<PR#>.emicorr.rst``
- ``<PR#>.saturation.rst``
- ``<PR#>.ipc.rst``
- ``<PR#>.firstframe.rst``
- ``<PR#>.lastframe.rst``
- ``<PR#>.reset.rst``
- ``<PR#>.superbias.rst``
- ``<PR#>.refpix.rst``
- ``<PR#>.linearity.rst``
- ``<PR#>.rscd.rst``
- ``<PR#>.persistence.rst``
- ``<PR#>.dark_current.rst``
- ``<PR#>.charge_migration.rst``
- ``<PR#>.jump.rst``
- ``<PR#>.picture_frame.rst``
- ``<PR#>.clean_flicker_noise.rst``
- ``<PR#>.ramp_fitting.rst``
- ``<PR#>.gain_scale.rst``

Stage 2 Steps
=============

- ``<PR#>.adaptive_trace_model.rst``
- ``<PR#>.assign_wcs.rst``
- ``<PR#>.badpix_selfcal.rst``
- ``<PR#>.msaflagopen.rst``
- ``<PR#>.imprint.rst``
- ``<PR#>.background.rst``
- ``<PR#>.extract_2d.rst``
- ``<PR#>.master_background.rst``
- ``<PR#>.wavecorr.rst``
- ``<PR#>.srctype.rst``
- ``<PR#>.targ_centroid.rst``
- ``<PR#>.straylight.rst``
- ``<PR#>.wfss_contam.rst``
- ``<PR#>.flatfield.rst``
- ``<PR#>.fringe.rst``
- ``<PR#>.pathloss.rst``
- ``<PR#>.barshadow.rst``
- ``<PR#>.photom.rst``
- ``<PR#>.pixel_replace.rst``
- ``<PR#>.resample_spec.rst``
- ``<PR#>.residual_fringe.rst``
- ``<PR#>.cube_build.rst``
- ``<PR#>.extract_1d.rst``
- ``<PR#>.resample.rst``

Stage 3 Steps
=============

- ``<PR#>.assign_mtwcs.rst``
- ``<PR#>.tweakreg.rst``
- ``<PR#>.skymatch.rst``
- ``<PR#>.exp_to_source.rst``
- ``<PR#>.outlier_detection.rst``
- ``<PR#>.tso_photometry.rst``
- ``<PR#>.stack_refs.rst``
- ``<PR#>.align_refs.rst``
- ``<PR#>.klip.rst``
- ``<PR#>.spectral_leak.rst``
- ``<PR#>.source_catalog.rst``
- ``<PR#>.combine_1d.rst``
- ``<PR#>.ami.rst``

Other Changes
=============

- ``<PR#>.wfs_combine.rst``
- ``<PR#>.white_light.rst``
- ``<PR#>.engdb_tools.rst``
- ``<PR#>.guider_cds.rst``

- ``<PR#>.docs.rst``
- ``<PR#>.other.rst``: Infrastructure or miscellaneous change

Note
----

This README was adapted from the Astropy changelog readme under the terms
of BSD license, which in turn adapted from the Numpy changelog readme under
the terms of the MIT licence.
