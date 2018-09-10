0.12.0 (unreleased)
===================

ami
---

assign_wcs
----------

associations
------------

- Updated Level2 product naming to use pipeline's remove_suffix. [#2481]


background
----------

barshadow
---------


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


firstframe
----------


fits_generator
--------------


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

ramp_fitting
------------

refpix
------

resample
--------

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

scripts
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

wfs_combine
-----------

white_light
-----------

wiimatch
--------

0.11.0 (2018-09-10)
===================
=======
0.11.0
======
>>>>>>> 516140f8986e2f848b3e224e51ce56d592937491

The 0.11.0 release is highlighted by the inclusion of steps for resampling
spectral images and time series grism observations.   In addition, this
release had 39 issues closed and a number of pull requests to improve PEP8
compliance, improve performance, and enhance the testing.  The release also
included updated documentation for acessing CRDS when running the JWST 
pipeline and updates to the reference file documentation. 

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

- Updated NRC_TSGRISM to assign source location to set pixel [#2286]

- Fixed bug in assign_wcs for ordering of slits for NIRSPEC MSA data [#2366]

- Implemented support for reading and writing WCS information in the 
  WAVE-TAB format [#2350]

- Fixed bug in the ording of cube footprint [#2371]

associations
------------

- Implemented Rule for Level 2 Nirspec Fixed Slit background. [#2307]

- Included Handling of both numeric and named slits for Level3 products. [#2330]

- Removed MIR_LRS-SLITLESS and NIS_SOSS from the permanent TSO list. [#2330]

- Implemented new Level2a rule `Asn_Lv2NRSLAMP`. [#2177]

- Allowed "N/A" as a valid, but False, value in association pools. [#2334]

- Implemented new association types tso_image2 and tso_spec2. [#2431]

- Synced code version with jwst package version. [#2458]

- Implemented source naming for NIRISS WFSS Level3 associations [#2443]

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

- Fixed bug in cube_build.blot_images that was failing for  NIRSPEC IFU images
  with the slide position defined in the WCS [#2345]

- Updated the construction of cube footprint [#2371, #2364, #2327]

cube_skymatch
-------------

dark_current
------------

datamodels
----------

- Added a new info method, similar to the method in astropy fits [#2268]

- The ``DataModel`` ``__hasattr__`` method has been replaced by ``hasattr``.
  The former created the attribute when it was accessed. [#2275]

- Improved error messaging when loading fits files into data models. [#2298]

- New warning message when opening a file without DATAMODL keyword. [#2248]

- Included the ability to handle 'allOf' when reading in  schemas [#2407]

- Removed BaseExtension class, it was not being used [#2430]


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

- Fixed step calling error for unreferenced attribute [#2463]

- Fixed type specification for optional grism mode inputs [#2467]

firstframe
----------

- Unit tests added to the first frame step [#2365]

fits_generator
--------------

- Updated pyparsing to v 2.2.0 [#2382]

- Updated fits_generator to ignore files begining with '.' [#2333]

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

- Unit tests added for lastframe [#2412]

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

- Enabled resample_spec to run on MIRI fixed slit data in calspec2 [#2424]

- Implemented new `Spec2Pipeline` configuration for NIRSpec LAMP exposures [#2174]

- Implemented specific exit status for "no science on detector" [#2336]

- Enabled `extract_2d` for NRC_TSGRISM [#2460]

- Turn off `resample` in `Spec2Pipeline` for multi-integration cube data [#2456]

ramp_fitting
------------

refpix
------

<<<<<<< HEAD
- The memory performance of refpix was improved [#2315]
=======
* The memory performance of refpix was improved [#2315]
>>>>>>> 516140f8986e2f848b3e224e51ce56d592937491

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

scripts
-------

- Added a new script for adding or removing files from an association [#2468]

stpipe
------

- Fixed bug to allow not being able to find a default input file name [#2461]

straylight
----------

superbias
---------

timeconversion
--------------

- Updated the utc_to_tdb module to compute the radial velocity (m / s) of JWST with respect to the solar-system barycenter, and to assign that value to keyword VELOSYS in the SCI header of the specified FITS file. [#2359]

transforms
----------

tso_photometry
--------------

- Updated tso_photometry step for SUB64P/WLP8 mode #2358


tweakreg
--------

- Fixed the coordinate frames in the output of tweakreg. [#2404]

- Updated TPCorr to work with V2, V3 in arcseconds instead of degrees [#2342]

wfs_combine
-----------

white_light
-----------

wiimatch
--------

