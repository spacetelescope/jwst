0.13.2 (Unreleased)
===================

background
----------

- Verify the exposures to be used as background have the same NIRSpec GWA
  tilt values as the science exposures. If the background and science
  exposures do not have matching GWA tilt values, then skip the background
  subtraction step in calspec2. [#3252]

calwebb_spec3
-------------

- Add the ``master_background`` subtraction step to the pipeline. [#3296]

combine_1d
----------

- Fix call to wcs.invert, and don't weight flux by sensitivity if the net
  column is all zeros. [#3274]

datamodels
----------

- Fix ``url_mapper`` for fits-schema to allow URLs with of the format
  http://stsci.edu/schemas/fits-schema/ to map to the correct location
  in the ``jwst`` package. [#3239]

- Change ``ModelContainer`` to load and instantiate datamodels from an
  association on init.  This reverts #1027. [#3264]

- Keyword updates to data model schemas, including OBSFOLDR, MIRNGRPS,
  MIRNFRMS, and new PATTTYPE values. [#3266]

extract_1d
----------

- This step can now use a reference image for IFU data.  The reference
  image (for IFU) may be either 2-D or 3-D.  When using a reference image
  for non-IFU data, background smoothing is now done after scaling the
  background count rate. [#3258]

- Unit tests were added for IFU data. [#3285]

master_background
-----------------

- Modified the unit tests for ``expand_to_2d``. [#3242]

- Modified ``MasterBackgroundStep`` to be skipped if ``BackgroundStep``
  was already run on the data.  A new ``force_subtract`` parameter is
  added to override this logic.  [#3263]

outlier_detection
-----------------

- Fixed a bug that was causing the step to crash when calling the
  ``cube_build`` step for MIRI MRS data. [#3296]

reffile_utils
-------------

- Improved error messages when problems are encountered in extracting
  subarrays from reference files. [#3268]

set_telescope_pointing
----------------------

- Fix ``populate_model_from_siaf`` to convert SIAF pixel scale from
  arcsec to degress for CDELTn keywords. [#3248]

tweakreg
--------

- Bug fix: Improved 2D Histogram (pre-match shift) algorithm in Python. [#3281]

- Fixed a bug in handling situations when no useable sources are
  detected in any of the input images. [#3286]

- Enhanced source catalog extraction algorithm to filter out sources outside
  the WCS domain of definition (when available). [#3292]

- Changed the type of exception raised when input has incorrect type. [#3297]

0.13.1 (2019-03-07)
===================

combine_1d
----------

- Added parameter ``background``; for background data, scale the flux,
  error, and net by 1 / NPIXELS, and include NPIXELS in the weight;
  changed the default for ``exptime_key`` to "exposure_time". [#3180]

- There is now a direct interface for calling the step.  This function,
  ``combine_1d_spectra``, may be passed either a ModelContainer or a
  MultiSpecModel object.  Previously this function expected the name of
  an association file. [#3220]

datamodels
----------

- Add back BaseExtension class so url-to-schema mapping works again [#3227]

extract_1d
----------

- If flux conversion is done, the FLUX is now set to zero (instead of
  copying the NET) if the wavelength of a pixel is outside the range of
  the RELSENS array. [#3190]

- Added a parameter ``subtract_background`` to ``extract_1d`` indicating
  whether the local background should be subtracted. If None, the value
  in the extract_1d reference file is used. [#3157, #3186]

- ``extract_1d`` can be run by calling ``extract.do_extract1d`` and
  passing a dictionary of reference file information. [#3202]

- ``ref_dict`` was None in ``run_extract1d``, and a check for that was
  missing. [#3233]

master_background
-----------------

- Added unit tests for expand_to_2d.  Support CombinedSpecModel data
  for the 1-D user-supplied background spectrum. [#3188]

set_bary_helio_times
--------------------

- Raise an exception when unable to compute converted times. [#3197]

set_telescope_pointing
----------------------

- Added population of CDELTn keywords based on SIAF values and fixed bug in calculation
  of S_REGION corners. [#3184]

0.13.0 (2019-02-15)
===================

ami
---

assign_wcs
----------

- Removed ``transform_bbox_from_datamodels`` in favor of
  ``transform_bbox_from_shape`` which now works by using last two dimensions
  in the ``shape``. [#3040]

- Added velocity correction model to the WFSS and TSGRISM wcs pipelines. [#2801]

- Refactored how the pipeline handles subarrays in the WCS. Fixed a bug
  where the bounding box was overwritten in full frame mode. [#2980]

- Rename several functions dealing with calculating bounding boxes for clarity. [#3014]

- The bounding box of the MIRI LRS WCS is now in "image" coordinates, not full frame. [#3063]

- FITS WCS keywords are written out only if the observation is one of the IMAGING_MODES. [#3066]

associations
------------

- Updated docstrings and written documentation. [#2856, #2862]

background
----------

barshadow
---------

combine_1d
----------

coron
-----

- Updated the `stack_refs` routine to update the output data model with metadata
  from the first input model. [#3111]

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
- Updated to recognize NRC_TSGRISM as WFSS data.  SlitDataModel schema now
  specifies that the wavelength attribute should be 2-D, with a default
  value of 0. [#2911]

- Reverse order of RELSENS wavelength and response if the wavelengths are
  not increasing. [#3005]

- Add a test for constant wavelengths (or constant slope). [#3032]

- Fix issue regarding mixing of the syntax for Boolean arrays and for
  integer index arrays. [#3045]

- Changed the names of time-related keywords for extracted spectra. [#3058]

- A new NPIXELS column has been added to the output table. [#3108]

extract_2d
----------
- Moved the update of meta information to the MultiSlitModel instead of the
  SlitModels that compose it. [#2988]

firstframe
----------

fits_generator
--------------

flatfield
---------
- Updated to not extrapolate for wavelengths that are out of bounds,
  either due to the WCS, or the wavelengths for a flat-field image cube,
  or the wavelengths for the fast-variation component. [#2775]

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
- Updated the docstrings [#2822]

jump
----

jwpsf
-----

lastframe
---------

lib
---

- ``set_telescope_pointing`` now populates WCS keywords from the SIAF file. [#3066]

linearity
---------

master_background
-----------------

- Implement the basic step scaffolding for `MasterBackgroundStep`. [#3090]

- Record user-supplied master background in MSTRBKGD keyword [#3101]

- Add step documentation for master background subtraction [#3102]

- Make master background step actually work [#3110]

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
- Ramp-fitting returning zero for all background pixels; Issue #2848, JP-453.

- MIRI ramps with jumps flagged at group 2 result in slopes of 0 in the rate
  image; Issue #2233,

- Processing pixels in ramp fitting in which all groups are saturated; Issue
  #2885.

- Ramp Fit fails when only two groups are in a segment after cosmic ray hits.;
  Issue #2832, JP-450.

- Fixed a bug in which the keywords from the input were not included in the OPT
  output header.

- Simplified and clarified classification of segment types based on DQ flags.

- Added handling of ramps ending in 2 saturated groups.

- Fix units for Read Noise Variance in ramp_fit (PR #2767). This may needed to
  revised based on Mike Regan's comment when he closed this PR.

- Added check to handle integration-specific variances for too short segments.

- More robust handling of ramps flagged as DO_NOT_USE (PR #3016)

refpix
------

- Added a description of processing for IRS2 readout mode data. [#2889]
- Fixed a mistake in the time to read one pixel. [#2923]

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

- Add `Step.record_step_status()` method for use by this step (and any other
  pipeline or pipeline step) [#3110]

straylight
----------

superbias
---------

timeconversion
--------------
- Updated the docstrings [#3020]

transforms
----------

tso_photometry
--------------

tweakreg
--------

- Use a more numerically stable ``numpy.linalg.inv`` instead of own matrix
  inversion. [#3033]

- Bug fix: Use integer division in Python 3. [#3072]


wfs_combine
-----------

white_light
-----------

wiimatch
--------

0.12.3 (2019-01-10)
===================

scripts
-------

- ``set_telescope_pointing.py``: Update method of choosing pointing parameters. [#2900, #3008, #3022]

- ``set_telescope_pointing.py``: Allow undefined SIAF. [#3002, #3006]

0.12.2 (2018-11-15)
===================

associations
------------

- Updated rules based on actual OTB phasing data. [#2831]

wfs_combine
-----------

- Renamed the configuration from `wfs_combine` to `calwebb_wfs-image3`. [#2831]


0.12.1 (2018-10-30)
===================

The 0.12.0 release is highlighted by the completion of updates for level-2b WFSS
processing, support for non-linear wavelength sampling in IFU cubes, and several
Associations updates to support WFS&C observations and background nodding.
This release had 53 issues closed and a number of pull requests to improve PEP8
compliance, improve performance, enhance the testing, and remove all python2
dependencies.  The release also included updated documentation of CRDS reference files.

ami
---

assign_wcs
----------

- The bounding box for NIRSpec WCS objects was modified to include the
  edges of the pixels. [#2491]

- Updated assign_wcs to compute the sky footprint of MIRI MRS and NIRSpec
  IFU observations. [#2474]

- Fixed minor bug in catalog.utl.get_object_info [#2550]

- Fixed bug in bounding_box_from_shape function [#2558]

- Make GrismObject.partial_order a lookup dict on order and fix partial_order logic [#2643]

- Added unit tests for grism modes [#2649]

- Augmented the logic for choosing a Nirspec WCS mode to include a check for the value
  of ``GRATING``. If ``GRATING=MIRROR`` imaging mode is chosen reegardless of ``EXP_TYPE``. [#2761]

- Added new NIRSpec target acq exposure types NRS_WATA and NRS_MSATA to be
  assigned an imaging WCS. Removed NRS_BOTA. [#2781]

associations
------------

- Updated Level2 product naming to use pipeline's remove_suffix. [#2481]

- Added rule Asn_Lv2NRSIFUNod to handle nod backgrounds for NIRSpec IFU [#2532]

- Changed deprecated logger.warn to logger.warning. [#2519]

- Made NIRISS WFSS Level2 associations exclusive. [#2555]

- Added new rule Asn_Lv2WFSC and new association type wfs-image2, including a new
  configuration file "calwebb_wfs-image2.cfg" [#2599]

- Added new rule Asn_Lv2MIRLRSFixedSlitNod to handle LRS Fixed-slit nodding. [#2663]

- Updated MIRI Dark and Flat exposure keywords. [#2698, #2710]

- Updated coronagraphy associations to be integrations-based. [#2773]

- Updated NIRSpec Lamp calibrations to be grating-specific. [#2780]

- Added new NIRSpec target acq exposure types NRS_WATA and NRS_MSATA. [#2780]

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

- Added support for creating IFU Cubes with non-linear wavelength sampling,
  including use of FITS WCS "WAVE-TAB" standard. [#2598]
- Correctly writing TDIM2 to WCS-TABLE extension [#2719]
- Fixed error when making IFUCubes with weighting='miripsf' [#2719]

cube_skymatch
-------------

dark_current
------------

datamodels
----------

- Initialize arrays and tables from function args in model_base [#2502]

- Updated guidestar centroid table column data type [#2526]

- Updated BAND keyword allowed values to include cross-dichroic combinations [#2530]

- Truncate long schema validation error messages to 2000 characters [#2657]

- Various keyword changes, including new EXP_ONLY keyword [#2414]

- Added validate_required_fields to datamodels base, so that "fits_required" is
  checked when writing a model to a file [#2589]

- Added new keywords PWFSEET, NWFSEST, DATE-BEG and made updates to conform to
  FITS convention for units included in keyword comments [#2595]

- Updated allowed SUBARRAY names for FGS and NIRCam [#2667]

- Fixed bug in default value when schema contains combiner [#2668]

- Updates for python 2 to 3 conversion [#2678]

- Updated EXP_TYPE allowed values to include "MIR_DARKALL", "MIR_DARKIMG",
  "MIR_DARKMRS", "MIR_FLATALL", "MIR_FLATIMAGE-EXT", and "MIR_FLATMRS-EXT" [#2709]

- Updated the MiriResolutionModel schema to have column names match the actual
  reference files [#2757]

- Updated EXP_TYPE allowed values to remove NRS_BOTA and replace with NRS_MSATA
  and NRS_WATA [#2772]

documentation
-------------

- Clarifications of input and output file naming. [#2727]


dq_init
-------

- Added ValueError check when loading the input into a data model [#2543]

emission
--------

engdblog
--------

exp_to_source
-------------

extract_1d
----------

- Added or modified docstrings [#2769]

extract_2d
----------

- WFSS modes updated to only extract specific orders, including delivery of updated
  wavelengthrange reference file [#1801]

- Fixed NIRSpec cutout size bug related to FITS 1-indexing [#2541]

- Added bounding box to WFSS output SlitModel [#2643]

- Added unit tests for grism modes [#2649]

- Bounding box sizes in extracted WFSS exposures now correctly cover entire extraction [#2799]

firstframe
----------


fits_generator
--------------

- NIRSpec data now automatically sanitizes the GWA_TILT keyword. [#2494]


flatfield
---------

- Modified the code to find the dispersion direction. [#2492]

- Changed the handling of zero wavelengths for NIRSpec data. [#2659]

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

- Updated step docs, as well as gain and readnoise reference file docs [#2689]

jwpsf
-----

lastframe
---------

lib
---

- Updated reffiles_utils to no longer issue warnings about mismatch in
  data array size params for NIRSpec IRS2 readouts. [#2664]

- Updated reffiles_utils to regard IRS2 science exposures as a match with normal
  sized reference files. [#2755]

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

- Renamed calwebb_tso_image2, calwebb_tso_spec2, and calwebb_nrslamp_spec2 configuration files to
  calwebb_tso-image2.cfg, calwebb_tso-spec2.cfg, and calwebb_nrslamp-spec2.cfg [#2639]

- Updated the order of MIRI steps in calwebb_detector1 and calwebb_dark. [#2669]

- Updated Image2Pipeline and Spec2Pipeline to properly return "cal" results. [#2676]


ramp_fitting
------------

- Improved memory management; Corrected handling of groups in which all pixels have
  insufficient data for a first difference; Corrected handling of ramps whose initial group
  is saturated; Corrected handling of ramps whose single good segment is a single group. [#2464]

- Updated gain and readnoise reference file docs [#2689]

- Fixed bug so that an integration-specific (_rateints) product is only created when
  NINTS>1; Skip MIRI first and/or last groups when flagged as DO_NOT_USE. [#2760]

- Fixed bug in which the number of segments returned exceeds the number
  of groups, which had occurred for a MIRI dataset in which the first or last
  group was flagged as DO_NOT_USE and also flagged as a jump. [#2834]

refpix
------

resample
--------

- Made finding the dispersion axis more robust [#2644]

reset
-----

rscd
----

saturation
----------

- Updated step docs, as well as saturation reference file docs [#2689]

skymatch
--------

- Made skymatch to not fail in 'match' mode when images do not overlap [#2803]

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

- Made partial_order attribute of GrismObject as lookup dict on order [#2643]

tso_photometry
--------------

tweakreg
--------

- Modified default configuration settings: increased "kernel_fwhm" from 2.0
  to 2.5, increased "snr_threshold" from 3 to 10,
  and changed "enforce_user_order" from True to False. [#2510]

- Updated tweakreg to use ``wcs.available_frames`` to get the names of the
  frames in a WCS pipeline. [#2590, #2594, #2629]

- Made the code more robust with images without sources [#2796]

- Made the logic for computations of footprints more reliable for the
  case of 1 or 2 sources in a catalog. [#2797]


- Added two new parameters: ``brightest`` to keep the top ``brightest``
  (based on the flux) objects in the object catalog *after all other
  filtering has been applied* and ``peakmax`` to exclude sources with
  peak pixel values larger or equal to ``peakmax``. ``brightest`` can be used
  to eliminate false detections and ``peakmax`` can be used to filter out
  saturated sources (instrument-specific value).[#2706]

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

- NRC_TSGRISM extract_height honored, bounding box fixed [#2643]

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

- Added support for correcting NIRISS SOSS mode exposures [#2588]

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

- Removed python2-3 dependency in crds_client [#2593]

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

- Update the calwebb_tso1 cfg file to skip the firstframe step and save the corrected ramp product. [#2280]

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

