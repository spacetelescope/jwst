0.11.0(Unreleased)
==================

ami
---

assign_wcs
----------

- Fixed a bug in ``get_msa_open_slits`` which prevented the code
  from finding the msa metafile.                                 [#2322]

- Fixed a bug in computing the slit_y locations for Nirspec MSA
  slitlets with more than one shutter.                           [#2325]

- Added a wavelength correction for the effective velocity of JWST
  relative to the barycenter.                                  [#2359, #2406]

- NRC_TSGRISM assigns source location to set pixel [#1235]
associations
------------

- Implemented Rule for Level 2 Nirspec Fixed Slit background. [#2307]
- Handle both numeric and named slits for Level3 products. [#2330]
- Remove MIR_LRS-SLITLESS and NIS_SOSS from the permanent TSO list. [#2330]
- Implement new Level2a rule `Asn_Lv2NRSLAMP`. [#2177]
- Allow "N/A" as a valid, but False, value in association pools. [#2334]
- Sync code version with jwst package version. [#2458]

background
----------

barshadow
---------

- Fixed a bug in ``bar_shadow.py`` interpolate() that caused
  array index to be nan                                        [#2384]

combine_1d
----------

coron
-----

csv_tools
---------

cube_build
----------

cube_skymatch
-------------

dark_current
------------

datamodels
----------

- The ``DataModel`` ``__hasattr__`` method has been replaced by ``hasattr``.
  The former created the attribute when it was accessed. [#2275]

- Improved error messaging when loading fits files into data models. [#2298]

- New warning message when opening a file without DATAMODL keyword. [#2248]

- New info method, similar to the method in astropy fits [#2268]

- Removed BaseExtension class, it was not being used [#2003]

dq_init
-------

emission
--------

engdblog
--------

exp_to_source
-------------

extract_1d
----------

extract_2d
----------
- NRC_TSGRISM implemented with set source location and extraction options [#1710, #1235]

firstframe
----------

fits_generator
--------------

- updated pyparsing to v 2.2.0 [#2382]

flatfield
---------

fringe
------

gain_scale
----------

group_scale
-----------

guider_cds
----------

imprint
-------

ipc
---

jump
----

jwpsf
-----

lastframe
---------

lib
---

linearity
---------

model_blender
-------------

mrs_imatch
----------

msaflagopen
-----------

outlier_detection
-----------------

pathloss
--------

persistence
-----------

photom
------

pipeline
--------

- Fixed a typo in calspec2 which prevented the srctype
  step from running. [#2318]

- Enable resample_spec to run on MIRI fixed slit data in calspec2 [#2424]

- Implement new `Spec2Pipeline` configuration for NIRSpec LAMP exposures [#2174]

- Implement specific exit status for "no science on detector" [#2336]

ramp_fitting
------------

refpix
------

resample
--------

- Fixed spectral resampling so the 2D output for MIRI LRS and NIRSpec MSA
  has the correct orientation and a dispersion that matches the input, i.e.
  non-linear if a prism is in the optical path. [#2348]

- Fixed bug in spectral resampling of MIRI LRS where the interpolation of the
  dispersion was failing. [#2422]

reset
-----

rscd
----

saturation
----------

skymatch
--------

source_catalog
--------------

srctype
-------

stpipe
------

straylight
----------

superbias
---------

timeconversion
--------------

transforms
----------

tso_photometry
--------------

tweakreg
--------

- Fixed the coordinate frames in the output of tweakreg. [#2404]

wfs_combine
-----------

white_light
-----------

wiimatch
--------

