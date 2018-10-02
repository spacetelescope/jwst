0.12.0 (unreleased)
===================

ami
---

assign_wcs
----------

- The bounding box for Nirspec WCS objects was modified to include the
  edges of the pixels. [#2491]

- Updated assign_wcs to compute the sky footprint of MIRI MRS and Nirspec
  IFU observations. [#2472]

- Fix minor bug in catalog.utl.get_object_info()[#2550]

- Make GrismObject.partial_order a lookup dict on order and fix partial_order logic

associations
------------

- Updated Level2 product naming to use pipeline's remove_suffix. [#2481]

- Added rule Asn_Lv2NRSIFUNod to handle nod backgrounds for NIRspec IFU [#2532]

- Changed deprecated logger.warn to logger.warning. [#2519]

- Made NIRISS WFSS Level2 associations exclusive. [#2555]

- Added new rule Asn_Lv2WFSC and new association type wfs-image2 [#2599]

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

- Initialize arrays and tables from function args in model_base [#2502]


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

- WFSS modes updated to allow specific order extraction, updated wavelengthrange reference file delivered as part of these changes [#1801]

- add bounding box to wfss output SlitModel


firstframe
----------


fits_generator
--------------

- NIRSpec data now automatically sanitizes the GWA_TILT keyword. [#2494]


flatfield
---------

- Modified the code to find the dispersion direction. [#2492]

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

- Added new Image2Pipeline configuration calwebb_wfs-image2.cfg for WFS&C processing [#2599]

ramp_fitting
------------


refpix
------

resample
--------

- Make finding dispersion axis more robust in resample [#2644]

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

- Fixed bug in logging configuration for `set_telescope_pointing.py`. [#2521]

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

- NIRISS models updated to allow for negative filter wheel rotations [#1801]

- make partial_order attribute of GrismObject as lookup dict on order

tso_photometry
--------------

tweakreg
--------

- Updated tweakreg to use `wcs.available_frames` to get the names of the frames in 
  a WCS pipeline. [#2590]
wfs_combine
-----------

white_light
-----------

wiimatch
--------

0.11.0 (2018-09-10)
===================

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
---------
- Removed spaxel.py and replace class with set of arrays [#2472]
- reworked in mapping of the detector pixel to the sky spaxel so that consistent
  code can be used for both MIRI and NIRSPEC data [#2472]
- Removed some loops in cube_cloud.py for finding which pixels fall in roi
  of spaxels [#2472]
- In a test with MIRI data there was a 13% improvement in the speed of making IFUcubes. In the
  NIRSPEC case there was a 40% improvment in the speed of creating IFUCubes.

- Fixed bug in cube_build.blot_images that was failing for  NIRSPEC IFU images
  with the slide position defined in the WCS [#2345]

- Updated the construction of cube footprint [#2371, #2364, #2327]

cube_skymatch
-------------

dark_current
------------

datamodels
----------

- Initialize arrays and tables from function args in model_base [#2351]

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

- NRC_TSGRISM extract_height honored, bounding box fixed

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

- added support for NIRISS SOSS [#2588]

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

- The memory performance of refpix was improved [#2315]

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

0.10.0 (2018-07-30)
===================

The 0.10.0 release is a snapshot release for DMS testing.   The release
is highlighted by the inclusion of steps for time series observations.
This release had 39 closed issues included a number of improvements
to the wavelength calibration for NIRSPEC observations. 


ami
---

assign_wcs
----------

- Improved the error handling for missing entries in the wavelengthrange reference file [#2213]

- Fix to correctly calculate the wavelength for NIRSPEC Prism observations [#2163]

- process NRS_AUTOFLAT as a MOS observation [#2166]

- fix wavelength units of inverse transform [#2158]

- fix input units to meters when filter=OPAQUE [#2134]

associations
------------

- Implement NIRSpec MSA Background Nod rules #2249


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
---------


cube_skymatch
-------------

dark_current
------------

datamodels
----------

- When reference files are validated, they can either throw a warning or an
  error if strict validation is set. [#2210]

- Update schema enum lists for keywords FILTER, PUPIL, READPATT, and EXP_TYPE [#2226]

- Enable and improved tests for datamodel schemas using the ASDF schema checker [#2240, #2241]

- Update IRS2 data model and add regredssion tests [#2295] 


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
- An example has been added to the model_blener documentation for how to blend meta information [#2206]

mrs_imatch
----------

msaflagopen
-----------

- Added documentation for the msaflagopen step [#2283]

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

-Update the calwebb_tso1 cfg file to skip the firstframe step and save the corrected ramp product. [#2280]


- Implement TSO-specific Level2 configurations [#2297]

ramp_fitting
------------
- Corrected handling of ramps whose first differences are all NaNs (such as ramps with all groups saturated) [#2289]

refpix
------

- Refpix has been updated to handle subarray exposures [#2207]


resample
--------
- Fixed update_fits_wcs() to work on DrizProductModels [#2222]

- A major re-factoring of the resampling code to allow for spectroscopic resampling [#2245]

reset
-----

rscd
----

- The performance of the RSCD step was improved by a factor of 20 [#2247]

- Update to the RSCD documentation [#2211]


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

- A script was written to read the UTC columns (at the start, middle, and end of each integration) from the INT_TIMES table, call the timeconversion module to compute the corresponding times at the solar-system barycenter (TDB), and update the columns in the INT_TIMES table.  [#2285] 

- Fix the problem in timeconversion that was caused by a recent addition of a new field to the ephemeris by retrieving only the fields needed. [#2296]


transforms
----------

tso_photometry
--------------

- MIRI aperture photometry was added to the TSO photometry [#2215]

- Added a new model for setting parameters for TSO photometry [#2239]

- Add a  reference file for use with tso_photometry [#2254, #2264]



tweakreg
--------

wfs_combine
-----------

white_light
-----------

wiimatch
--------

