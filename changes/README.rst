Changelog
=========

This directory contains "news fragments" which are short files that contain a
small **ReST**-formatted text that will be added to the full changelog.

Make sure to use full sentences with correct case and punctuation.

Consuming news fragments in `changes/` into a new change log entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running `towncrier build` will read all existing fragment files in `changes/`
and create a new entry at the top of `CHANGES.rst` with the specified version number.

```shell
pip install towncrier
towncrier build --version <VERSION>
```

News fragment change types
--------------------------

- ``<PR#>.breaking.rst``: Also add a fragment of this type if your change breaks if your change breaks **step-level or public API** ([as defined in the docs](https://jwst.readthedocs.io/en/latest/jwst/user_documentation/more_information.html#api-public-vs-private))

General Pipeline Changes
""""""""""""""""""""""""

- ``<PR#>.stpipe.rst``
- ``<PR#>.datamodels.rst``
- ``<PR#>.scripts.rst``
- ``<PR#>.set_telescope_pointing.rst``
- ``<PR#>.pipeline.rst``
- ``<PR#>.associations.rst``

Stage 1
"""""""

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

Stage 2
"""""""

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

Stage 3
"""""""

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

Other
"""""

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
