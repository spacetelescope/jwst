1.5.2 (2022-05-20)
==================

extract_1d
----------

- In SOSS ATOCA, catch negative infinite values in centroid finder;
  catch spline-fit errors in first order flux estimate [#6854]

align_refs
----------

- Change median filter warning message to debug level [#6853]

linearity
---------

- Correct bug when using ZEROFRAME data. [#6851]

ramp_fitting
------------

- Remove the logic that only copied the INT_TIMES table content when processing
  TSO exposures, so that it shows up in all ``rateints`` products [#6852]

- Updated the one good group ramp suppression handler. [spacetelescope/stcal#92]

1.5.1 (2022-05-17)
==================

cube_build
----------

- Fix for residual spectral tearing in MIRI MRS multiband cubes [#6786]

dark_current
------------

- Eliminated extra copying of input model when step gets skipped [#6841]

datamodels
----------

- Update keyword comments/titles for V2_REF, V3_REF, FITXOFFS, FITYOFFS [#6822]

extract_1d
----------

- Fix bug in SOSS algorithm for bad data by replacing source of possible
  infinite values with NaNs, caused by zero division [#6836]

- Exit gracefully if data is with F277W filter; avoid masking entire wavemap
  if subarray is SUBSTRIP96 [#6840]

jump
----
- Enable multiprocessing in jump detection [#6845]

lib
---

- Update ``test_siafdb`` unit test due to recent SIAF DB update [#6842]

linearity
---------

- Adding feature to process ZEROFRAME data with the linearity step. [#6782]

ramp_fitting
----------

- Adding feature to use ZEROFRAME for ramps that are fully saturated, but
  the ZEROFRAME data for that ramp is good. [#6782]

refpix
------

- Adding feature to process ZEROFRAME data with the refpix step. [#6782]

saturation
----------

- Adding feature to process ZEROFRAME data with the saturation step. [#6782]

stpipe
------

- Log the CRDS context for pipeline and standalone step processing [#6835]

superbias
---------

- Adding feature to process ZEROFRAME data with the superbias step. [#6782]

tweakreg
--------

- Changed default value of ``fitgeom`` from ``'general'`` to ``'rshift'``
  at the request of CalWG. [#6838]

1.5.0 (2022-05-05)
==================

associations
------------

- Implement PoolRow to avoid deep copy of the AssociationPool table [#6787]

- Added valid optical paths for NRS_LAMP observations to generate
  or exclude associations using lamp, disperser and detector [#6695]

- Include filename extension for `asn_pool` entry, to maintain consistency
  with `asntable` entry [#6699]

- Add constraint on NIRCam TSGRISM exposures, preventing level 2 and 3
  associations for detector NRCBLONG [#6709]

- Add fgsid option to set_telescope_pointing [#6717]

- Formalize the mkpool utility [#6746]

align_refs
----------

- Fixed behavior generating many unnecessary and slow logging warnings on
  MIRI coronagraphy data, due to large contiguous regions of NON_SCIENCE
  pixels [#6722]

ami
---

- Allow AmiAverageStep to be run on list in command line interface [#6797]

assign_wcs
----------

- Corrected computation of crpix by backward transform of fiducial, allow
  for reference outside of detector frame [#6789]

- Fixed parsing the ``filteroffset`` file which resulted in the offset
  not being used by the WCS. [#6831]

- Fixed assignment of ``wcs.bounding_box`` in MIRI, NIRISS and NIRCAM imaging mode. [#6831]

background
----------

- Added the step parameter ``wfss_mmag_extract`` to allow for setting the
  minimum magnitude of source catalog objects to be used in the WFSS
  background subtraction process [#6788]

- Added a check to make sure that a sufficient number of background
  (source-free) pixels are available in a WFSS image before attempting
  to compute statistics and scale the WFSS background reference image
  [#6788]

cube_build
----------

- Fixed a bug in how the DQ plane of NIRSpec data is set [#6718]

- Use drizzle weight function by default instead of EMSM. [#6820]

- Fix bug for internal_cal cubes produces by move to drizzle default. [#6826]

- Fix bug for Single type cubes called by mrs_imatch using drizzle. [#6827]

cube_skymatch
-------------

- Enabled support for mnemonic DQ codes in the ``cube_skymatch`` step.
  [#6733, #6736]

datamodels
----------

- Added the new keyword "BKGMETH" for use in the ``skymatch`` step.
  [#6736]

- Drop references to transform-1.2.0 from datamodel schemas to prevent
  issues with schema features not supported by stdatamodels. [#6752]

- Remove FILETYPE keyword from core schema, and all assignments to it [#6772]

- Update rscd model to increase the size of group_skip_table to allow FASTR1, SLOWR1, FASTR100 [#6776]

- Correcting the default ZEROFRAME allocation. [#6791]

- Add the new MIRI MRS point source correction reference file data model
  MirMrsPtCorrModel. [#6762]

- Add new datamodel and schema for MIRI MRS cross-artifact reference file
  MirMrsXArtCorrModel [#6800]

- Create MSA_TARG_ACQ table extension schema [#6757]

- Added selector keywords ``readpatt`` and ``preadpatt`` to MIRI flat schema. [#6825]

documentation
-------------

- Added documentation for processing NIRSpec lamp mode data in Spec2Pipeline
  description [#6812]

- Document parameter reference files in the same manor as other references [#6806]

extract_1d
----------

- Clean the logging statements made by `extract_1d` to make the log
  more useful [#6696]

- Check for non-zero array size before computing sigma-clipped
  statistics in IFU mode [#6728]

- Propagate non-differentiated errors for IFU mode observations [#6732]

- Remove temporary `soss_atoca` parameter and make ATOCA the default
  algorithm for SOSS data [#6734]

- Add separate behavior for 2D vs (3D data with only one image)
  by passing appropriate integ value [#6745]

- Allow reference files to specify extraction region for extended
  sources, modify `bkg_fit` default to None while retaining `poly`
  as default mode [#6793]

flatfield
---------

- Change DQ flags for NIRSpec flatfield where one or more component flats
  (fflat, dflat, sflat) is bad (#6794)

general
-------

- Added aliases to all steps, following step_defs naming conventions [#6740]

- Require scikit-image as a dependency (for source catalog deblending).
  [#6816]

lib
---

- Updated default suffix names for RampFit and GuiderCDS steps to
  'ramp_fit' and 'guider_cds' to match alias convention [#6740]

mrs_imatch
----------

- Use drizzle weight function by default instead of EMSM. [#6820]

photom
------

- Allow SOSS input as MultiSpecModel, and do correction on extracted 1d
  spectra [#6734]

pipeline
--------

- Improve memory performance of `calwebb_detector1` pipeline [#6758]

- Update the `calwebb_spec2` pipeline to allow for the creation of an
  optional WFSS product that's in units of e-/sec [#6783]

- Updated `calwebb_spec2`, `calwebb_spec3`, and `calwebb_tso3` to reorder
  step processing for SOSS data - `photom` now comes after `extract_1d` [#6734]

- Added ResetStep back into `calwebb_dark` for MIRI exposures [#6798]

ramp_fitting
------------

- Updated step docs to clarify exactly what calculations are used for
  the various flavors of variance and ERR stored in the output
  products [#6715]

- Adding feature to turn off calculations of ramps with good 0th group,
  but all other groups are saturated. [#6737]

- Fix for handling jumps in the first good group following dropped groups.
  [spacetelescope/stcal#84]

regtest
-------

- Added a residual fringe correction test [#6771]

resample
--------

- Fixed ``resample_spec`` output spectrum centering issue for MIRI LRS
  fixed-slit. [#6777]

- Re-designed algorithm for computation of the output WCS for the
  ``resemple_spec`` step for ``NIRSpec`` data. [#6747, #6780]

- Fixed handling of user-supplied ``weight_type`` parameter value for
  ``resample_spec``. [#6796]

- Fixed an issue with axis number for the spectral axis in ``resample_spec``. [#6802]

reset
-----

- Fix bug in how segemented data is corrected [#6784]

residual_fringe
---------------

- Replaced fitting the background with an astropy fitting package [#6739]

saturation
----------

- Updated to allow the step to flag neighbors of saturated pixels, which is
  controlled by the new step param ``n_pix_grow_sat``, to account for charge
  migration. [spacetelescope/stcal#83] [#6818] [#6830]

skymatch
--------

- Updated the step docs to clarify the details of the various global,
  match, and global+match methods. [#6726]

- Enabled support for mnemonic DQ codes in the ``skymatch`` step. Also
  changed default value for ``dqbits`` from 0 (exclude ALL flagged in DQ
  pixels) to ``'~DO_NOT_USE+NON_SCIENCE'``. [#6733, #6736]

- Updated to populate the "BKGMETH" keyword in output files. [#6736]

- Increased tolerance value for considering two sky polygons identical. [#6805]

source_catalog
--------------

- Fixed the KDTree calculation to use only finite source positions to
  prevent memory issues on Linux systems. [#6765]

- Updated the roundness and sharpness properties to use the source
  centroid position instead of the peak position. [#6766]

- Updated the catalog metadata. [#6813]

srctype
-------

- Add command line option to override source type [#6720]

tweakreg
--------

- Make ``fit_quality_is_good()`` member private and rename it to
  ``_is_wcs_correction_small()``. [#6781]

- Change default settings for ``searchrad``, ``tolerance``, and ``separation``
  parameters for the ``tweakreg`` step. [#6809]

- Change default value of ``brightest`` parameter in the ``tweakreg`` step. [#6810]


1.4.6 (2022-03-25)
==================

set_telescope_pointing
----------------------

- Add option --force-level1bmodel. [#6778]

1.4.5 (2022-03-23)
==================

datamodels
----------

- Updated reset model to include NINTS, NGROUPS keywords and the subarray.schema [#6749]

- Update reset model to include keyword_preadpatt.schema [#6769]

- Update rscd model to increase the size of group_skip_table to allow FASTR1, SLOWR1, FASTR100 [#6776]

reset
-----

- Read NINTS and NGROUPS from model.meta for reset reference file and data instead of using the
  shape of the data to define these values [#6749]

1.4.4 (2022-03-16)
==================

set_telescope_pointing
----------------------
- Set CRVAL* from GS_* for guider exposures. [#6751]

- Add fgsid option to set_telescope_pointing [#6717]

- Further restrict default models that can be updated. [#6767]

- Update COARSE handling of FGS, pysiaf importing, model opening,
  and removal of stale code. [#6735]


1.4.3 (2022-02-03)
==================

set_telescope_pointing
----------------------

- JP-2509 Update COARSE algorithm to use FGS1 exclusively. [#6700]


1.4.2 (2022-01-20)
==================

assign_wcs
----------

- Modified requirements for grism bounding box location to have
  width greater than one pixel [#6579]

associations
------------

- Changed restriction on Level2b creation for ``NRC_TACQ`` exposures
  to ``NRC_IMAGE`` to allow asn creation for tacq but not science [#6681]

extract_1d
----------

- Removed unnecessary verbose parameter being passed amongst
  extract_1d functions, but not user-accessible [#6579]

outlier_detection
-----------------

- Added MIRI MRS cross bands to options for the type of IFU cubes being created [#6666]

skymatch
--------

- Changed default value of ``skymethod`` step parameter to 'match' [#6580]

1.4.1 (2022-01-15)
==================

dark_current
------------

- Added docs mention of external algorithm in ``stcal`` [#6566]

datamodels
----------

- Update names of MIRI coronagraphic flats dither patterns [#6573]

jump
----

- Added docs mention of external algorithm in ``stcal`` [#6566]

- Fix issue in jump detection that occurred when there were only 2 usable
  differences with no other groups flagged. This PR also added tests and
  fixed some of the logging statements in two-point difference. [spacetelescope/stcal#74]

linearity
---------

- Added docs mention of external algorithm in ``stcal`` [#6566]

ramp_fitting
------------

- Added docs mention of external algorithm in ``stcal`` [#6566]

saturation
----------

- Added docs mention of external algorithm in ``stcal`` [#6566]

1.4.0 (2022-01-10)
==================

ami_analyze
-----------

- Call ``img_median_replace`` to replace NaN's and pixels flagged with
  DO_NOT_USE in the input image. [#6334]

assign_wcs
----------

- Open the specwcs reference file for WFSS modes using the ``with`` context
  manager. [#6160]

- Fix bug in NIRspec where ``bounding_box`` can be oversized in height for
  some of the slits. [#6257]

- Updated ``create_grism_bbox`` to be more robust against failures caused by
  bad input data. [#6309]

- Added a function that, when given RA, Dec, lambda, computes which ones project
  into a given NIRSpec IFU slice. [#6316]

- Changed in_ifu_slice in util.py to return the indices of elements in slice.
  Also the x tolerance on finding slice elements was increased. [#6326]

- Fix a bug due to which, under certain circumstances, ``PC``-matrix and
  ``CDELT`` coefficients may be removed from FITS WCS of data products. [#6453]

- Fixed bug in NIRSpec MOS slitlet meta data calculations for background slits
  consisting of multiple shutters. [#6454]

- Enabled spectral order as input to WCS for NIRISS SOSS mode. [#6496]

associations
------------

- Remove MIR_FLATMRS from Asn_Lv3MIRMRS rule. [#6548]

- Fix bug causing ``pytest`` to encounter an error in test collection when
  running with recent commits to ``astropy`` main (``5.0.dev``). [#6176]

- Enhanced level-2b ASN rules for NIRSpec internal lamp exposures to
  handle certain opmode/grating/lamp combinations that result in no data
  on one of the detectors. [#6304]

- Removed Constraint_ExtCal from Asn_Lv2WFSC constraints, as it was
  redundant with Constraint_Image_Science present. [#6384]

- Added constraint to Asn_Lv2ImageNonScience to prevent creation of asns
  for NRC_TACQ exposures with WFSC_LOS_JITTER in the DMS_NOTE. Also added
  new reduce method, Constraint.notall [#6404]

barshadow
---------

- Modify computation of correction array to remove dependencies on the
  centering of the source within the slitlet, because for extended/uniform
  sources the centering is irrelevant. This fixes bugs encountered when
  computing the correction for background signal contained within slits with
  an off-center point source. [#6454, #6459]

cube_build
----------

- Fix bug when creating cubes using output_type=channel. [#6138]

- Move computationally intensive routines to C extensions and
  removed miri psf weight function. [#6093]

- Moved cube blotting to a C extension [#6256]

- Moved variable definitions to top of code in C extension to
  support changes in #6093. [#6255]
- Added weighting option driz (3D drizzling) [#6297]

- Using assign_wsc.utils.in_ifu_slice function to determine which NIRSpec
  sky values mapped to each detector slice. [#6326]

- Fixed error in final exposure times calculated by blend headers. Only the input models
  used in the IFU cube are passed to blend headers. [#6360]

- Update of documentation to explain 3d drizzling and remove miri psf weighting [#6371]

- Fix a bug when creating internal_cal type cubes [#6398]

- Fix incorrect spatial footprint for single band MRS IFU cubes [#6478]

dark_current
----------

- Refactored the code in preparation for moving the code to STCAL. [#6336]

- Moved dark current code from JWST to STCAL. [spacetelecope/stcal#63] [#6444]

- Updated step docs to explain unique form of MIRI dark reference data [#6529]

datamodels
----------

- Remove astropy.io registration of JwstDataModel. [#6179]

- Update VELOSYS keyword comment [#6298]

- Added new keywords FPE_SIDE and ICE_SIDE to core schema [#6314]

- Fix bug preventing extra arguments when calling ``datamodels.open``
  on an ASDF file. [#6327]

- Implement memmap argument when calling ``datamodels.open`` on an ASDF
  file. [#6327]

- Fix bug in schema that disallowed valid p_grating values. [#6333]

- Add ``NDArrayType`` to list of valid types for ``RegionsModel.regions``. [#6333]

- Fix a bug in wcs_ref_models where SpecwcsModel was failing the SimpleModel
  validation as it contains a list of models rather than one simple model.
  Also add some missing allowed BAND values for MIRI MRS distortion
  and regions files.  Fix an incorrect comment on
  FilteroffsetModel. [#6362]

- Changed reference file model name from ``ResidualFringeModel`` to
  ``FringeFreq`` [#6385]

- Updated data products documentation to indicate that variance and error arrays
  are now included in resampled products. [#6420]

- Added SOSS-specific extraction parameters to core schema; add new
  datamodel to store SOSS model traces and aperture weights [#6422]

- Added the ``MirLrsPathlossModel`` for use in the ``pathloss` step. [#6435]

- Added new column 'reference_order' to 'planned_star_table' in
  guider_raw and guider_cal schemas [#6368]

- Moved new column 'reference_order' in guider schemas' planned
  star table to second in order, after 'guide_star_order' [#6465]

- Updated moving_target schema changing mt_detector_x/y to mt_sci_x/y [#6485]

- Fixed names of NIRISS SOSS extract_1d parameter keywords to be legal FITS [#6499]

- Update PATTTYPE enum values to match spellings used in keyword dictionary [#6501]

- Updated documentation to point to stdatamodels.util
  for calls to create_history_entry [#6537]

- Added keyword EXP_TYPE to PsfMaskModel schema [#6540]

- Updated FILTEROFFSET reference file docs to add NIRCam information. [#6541]

dark_current
------------

- Fixed bug during save of optional averaged darks output, bug with
  providing step a file instead of a datamodel, added regression test [#6450]

documentation
-------------

- Update text to point to the JWST CRDS website. [#6549]

- Update to calwebb_detector documentation to include the reset step as one of the steps applied
  to MIRI data [#6785]

extract_1d
----------

- Updated to propagate SRCTYPE keyword during extraction of MIRI LRS
  fixed-slit inputs that are in `SlitModel` form. [#6212]

- Assign 0-indexed integration number to INT_NUM if input
  INT_TIMES table is empty. [#6369]

- Change INT_NUM assignment to 1-indexed. [#6388]

- Added NRS_LAMP as an exp_type that has the extract1d ref file in asdf format [#6460]

- Added the ``center_xy`` step argument to allow user-specified x/y
  center of IFU extraction apertures [#6503]

- Delivery of new algorithm `ATOCA` for SOSS extraction, along with four new reference
  files: speckernel, specprofile, spectrace and wavemap. [#6467]

- Added step parameter `soss_atoca` to turn ATOCA algorithm on, with box extraction
  the default algorithm [#6551]

flatfield
---------

- Updated flatfield step docs to include complete details on how the
  variance and error arrays are updated. [#6245]

- Fixed a bug in flatfield for NIRSpec BrightObj mode where the S-flat cutout
  was calculated incorrectly by not accounting for the slit offset [#6332]

- Added check to NRS_LAMP exposures that routes imaging exposures to the imaging
  half of flatfield, where they will skip the step as expected [#6462]

jump
----

- Updated jump detection step to use common code moved to stcal [#6089]

- In stcal (pr #72), several changes were made to fix existing bugs in the
  twopoint difference routine for jump detection. Some of these issues
  resulted in jumps erroneously being flagged for pixels with only two
  usable groups (i.e one usable difference). This PR on the JWST side
  fixes one of the unit tests to account for this. [#6552]

lib
---

- Implement the MAST AUI interface to the Engineering Database. [#6288]

- Fix ROLL_REF and angle_to_vector calculations [#6452]

- Fix bad implementation of ``angle_to_vector`` in ``set_telescope_pointing``. [#6452]

- Use TRACK algorithms for moving target exposures. [#6452]

- Move setting of the default method to calc_transforms. [#6482]

linearity
--------

- Use the common code in STCAL for linearity correction. [#6386]

- Update of linearity test to support STCAL PR65 [#6509]

outlier_detection
-----------------

- Revert back to using 'linear' interpolation method as default for ``blot``.
  The bug in the implementation of the bilinear interpolator in the ``drizzle``
  package is now fixed. [#6146]

- Log number of flagged outliers in ``outlier_detection`` [#6260]

pathloss
--------

- Updated the ``pathloss`` step and documentation to include processing of
  MIRI LRS fixed-slit exposures. [#6435]

persistence
-----------

- Changed logger from root to `__name__` [#6389]

pipeline
--------

- Added wfss_contam step to `calwebb_spec2` pipeline flow for WFSS modes [#6207]

- Changed logger from root to `__name__` for Ami3, Detector1, Dark, and Guider
  Pipelines [#6389]

- Updated the ``calwebb_spec2`` pipeline to apply the ``pathloss`` step to
  MIRI LRS fixed-slit exposures. [#6435]

ramp_fitting
------------

- Fix special handling for 2 group ramp. [spacetelescope/stcal#70]

- Fix issue with inappropriately including a flagged group at the beginning
  of a ramp segment. [spacetelescope/stcal#68]

- Pixels with negative median rates will have VAR_POISSON set to zero.
  [spacetelescope/stcal#59]

- Update ``RampFitStep`` to pass DQ flags as a parameter to the ``ramp_fit``
  algorithm code in stcal.  Bump version requirement for stcal.  [#6072]

refpix
------

- Refactored the ``subtract_reference`` routine for NRS IRS2 data to reduce
  memory usage. [#6356]

- Updated bad link to JDox in the step documentation. [#6372]

regtest
-------

- Update okifying to handle full folder updates for associations [#6218]

- Remove default cfg usage from all relevant regtests; replaced with
  either pipeline alias or Step instance [#6391]

resample
--------

- Refactor ``resample_spec`` to use a separate function for computing the output
  rectified WCS for lamp data.  [#6296]

- Fix a crash in ``resample_spec`` due to undefined variance arrays. [#6305]

- Fix handling of ``weight_type`` parameter to allow for user override. [#6406]

- Add support for specifying custom output WCS parameters to the resample
  step. [#6364]

- Make ``output_shape`` to be in the "normal" (``nx, ny``) order. [#6417]

- Updated ``drizzle`` version to ``1.13.4`` which contains a fix for the
  bug due to which some 0-weight input pixels may contribute to the output
  image. [#6517]

- Updated step docs to indicate that the default weighting type is
  now "ivm" [#6529]

- Fixed a bug in the ``ResampleSpecData.build_interpolated_output_wcs()``
  due to which, under cerain circumstances, computed output image shape
  could be very large resulting in (very) large memory usage and/or
  incorrect output WCS. [#6533]

residual_fringe
---------------
 - Added documentation on step [#6387]
 - Fixed incorrect data model name [#6487]
 - Added user option to give wavelength range that no correction will be applied [#6545]

skymatch
--------

- Improved reliability when matching sky in images with very close sky
  footprints. [#6421]

- Updated code in ``skymatch.region.py`` with latest improvements and bug fixes
  from ``stsci.skypac.regions.py``. [#6451]

- Updated documentation to clarify details of flat-fielding versus distortion
  corrections [#6546]

source_catalog
--------------

- Fixed issue with non-finite positions for aperture photometry. [#6206]

- Fixed the documentation for ``bkg_boxsize`` to reflect that its data
  type should be integer. [#6300]

- Renamed ``filter_kernel`` to ``kernel`` in the call to ``detect_sources``
  to match the new name of the argument in photutils. [#6527]

wavecorr
--------

- Location of source in NIRSpec fixed slit updated
  (keywords ``SCRCXPOS``, ``SRCYPOS``). [#6243, #6261]

- Fixed the computation of ``model.slits[i].source_xpos``
  for Nirspec fixed slit data. [#6457]

wfs_combine
-----------

- Changed method of loading input association from datamodels.load() to
  Step.load_as_level3_asn() to prevent error when target acq exposure
  not present [#6464]

wfss_contam
-----------

- Updated to process all defined spectral orders for the grism mode [#6175]

- Added step documentation [#6210]

- Fixed handling of filter/pupil names for NIRISS WFSS mode [#6233]


1.3.3 (2021-10-05)
==================

- Avoid using photutils 1.2.0 [#6378]


1.3.2 (2021-09-03)
==================

associations
------------

- Enhanced level-2b ASN rules for NIRSpec internal lamp exposures to
  handle certain opmode/grating/lamp combinations that result in no data
  on one of the detectors. [#6304]

cube_build
----------

- Fix bug when creating cubes using output_type=channel. [#6138]

- Move computationally intensive routines to c extensions and
  removed miri psf weight function. [#6093]

- Moved variable definitions to top of code in c extension to
  support changes in #6093. [#6255]

- Moved cube blotting to a c extension [#6256]

pipeline
--------

- Updated calwebb_tso3 to be more robust in handling null results from
  the ``tso_photometry`` step. [#6265]


1.3.1 (2021-08-09)
==================

lib
---

- Fixed a bug in set_telescope_pointing that was setting wrong meta for the pointing quality [#6264]


1.3.0 (2021-07-31)
==================

associations
------------

- Ensure no Lv3_WFSC associations created on group candidates [#6131]

datamodels
----------

- Add new PATTTYPE values for MIRI Coronagraphic flats:
  4QPM_LFLAT, 4QPM_PFLAT, LYOT_LFLAT, LYOT_PFLAT. [#6232]

- Update ``DarkModel`` to use uint32 for DQ array. [#6228]

- Add NOUTPUTS keyword to the `DarkModel` schema. [#6213]

lib
---

- Add overriding of the matrix calculations to ``set_telescope_pointing.py`` [#5843]

- Add guide star-based pointing algorithm to ``set_telescope_pointing.py`` [#5843]

resample
--------

- Fix the extreme memory consumption seen in resampling of variance arrays. [#6251]

tweakreg
--------

- Add an upper tweak threshold of 10 arcsec to tweakreg [#6252]

wfs_combine
-----------

- Add option to flip the dither locations so that images with different
  filters will have the same pixel locations [#6101]

- Fixed the refine option to correctly use the cross correlation to align
  the images if the WCS is off [#6101]


1.2.3 (2021-06-08)
==================

datamodels
----------

- Add back and use "CALCULATED" for ENGQLPTG. [#6135]

- Convert incoming Path objects to strings in datamodels.open [#6130]


1.2.2 (2021-06-08)
==================

ami_analyze
-----------

- Fix to AMI pupil phases sign error [#6128]

datamodels
----------

- Update moving target schema to match b7.8 keyword schema. [#6129]


1.2.1 (2021-06-07)
==================

associations
------------

- Asn_Lv2WFSS: Add instrument constraint. [#6114]

- Asn_Lv2NRSLAMPSpectral: Allow msaspec only if msametfl is available. [#6085]

combine_1d
----------

- Added SRCTYPE to COMBINE1D output extension headers, propagated from
  EXTRACT1D inputs [#6079]

cube_build
----------

- Fix some typos in the the arguments documentation. [#6077]

datamodels
----------

- Updated enum lists for ENGQLPTG and PATTTYPE keywords [#6081]

- Removed obsolete keyword NDITHPTS and updated attributes for NRIMDTPT [#6083]

- Added units to CombinedSpecModel table output [#6082]

- Added keywords OSS_VER, DETMODE, CMD_TSEL, NOD_TYPE, and GS_V3_PA to
  the core schema [#6086]

- Remove ``ModelContainer`` schema and refactor use of association table
  metadata within. [#6094]

general
-------

- Make CRDS context reporting pytest plugin disabled by default. [#6070]

- Removed all usage of sys.path, in associations and jwst.stpipe [#6098]

lib
---

- Updated set_telescope_pointing to populate ENGQLPTG keyword with new
  allowed values [#6088]

outlier_detection
-----------------

- Avoid using 'linear' interpolation method as default for ``blot`` due to
  a bug in the implementation of the bilinear interpolator in the ``drizzle``
  package. Now the default value will be 'poly5'. [#6116]

ramp_fitting
------------

- Re-enable multiprocessing in ``RampFitStep`` by moving code back from
  stcal package. [#6119]

scripts
-------

- Add migrate_data command with support for migrating spec_table in
  x1d files produced with <= 1.1.0 of this package. [#6055]

tweakreg
--------

- Remove attached tweakreg catalog from datamodel before exiting step [#6102]


1.2.0 (2021-05-24)
==================

ami_analyze
-----------

- Create copy of input datamodel to avoid overwriting input. [#5828]

assign_wcs
----------
- Convert the ra values to array in util.wrap_ra, but if input is a list return
  a list [#6031]

- Moved the routine wrap_ra from cube_build to assign_wcs.util. The s_region is
  now correct for data that cross ra boundary. [#6026]

- Changed evaluation of grism bounding box center from averaged extrema of
  transformed bounding box to transformed centroid of source_cat object [#5809]

- Added pixel shift to MSA slits due to 0-indexing in NIRSpec slit validation
  code, fixing difference between bounding box locations during the separate
  halves of assign_wcs runs [#5927]

- Added logic to prevent the sending of an empty list of slits to the
  validate_open_slits function, so a proper error message is provided to
  the user [#5939]

- Added computed ``spectral_region`` to ``model.meta.wcsinfo``. [#5969]

associations
------------

- Add rule Asn_MIRMRSBackground to treat background as science. [#6046]

- Updated level2b WFSS rules to only consider exposures from the same
  instrument channel when matching direct images with grism images in
  NIRCam WFSS observations. [#5786]

- Removed PATTTYPE='None' constraint from Lv3MIRMRS association rule to
  generate spec3 associations for undithered MRS observations. [#5804]

- Updated level2b WFSS rules to only consider exposures using the same
  PUPIL value (cross filter) when matching direct images with grism images
  in NIRISS WFSS observations. [#5896]

- Updated level2b and level3 TSO rules to exclude exposures with
  EXP_TYPE=NRC_TSGRISM and PUPIL=CLEAR, which can result from NIRCam
  engineering template observations. [#5946]

- Updated level2b NIRSpec FS rules to exclude exposures sharing a primary
  dither location from the list of background exposures [#5994]

background
----------

- Remove unused ``SubtractImagesStep`` [#5919]

- Added new step parameter to optionally save the combined, average
  background image: ``save_combined_background``. [#5954]

calwebb_spec2
-------------

- Updated documentation to indicate that master_background is applied to
  NIRSpec MOS exposures in the calwebb_spec2 pipeline [#5913]

calwebb_spec3
-------------

- Updated documentation to indicate that master_background is applied to
  NIRSpec MOS exposures in the calwebb_spec2 pipeline [#5913]

csv_tools
---------

- The ``csv_tools`` subpackage was removed [#6006]

cube_build
----------

- Fixed typo in ``CubeBuildStep`` spec for grating [#5839]

- Update code to read in spectral and spatial size of exposure on the sky [#5991]

- For calspec2 pipeline skip determining the dq plane in ``cube_build`` [#5991]

- Remove certain WCS keywords that are irrelevant after ``cube_build``. [#6032]

datamodels
----------

- Added ``is_star`` to ``slitmeta`` [#5788]

- Update keyword comments for NIRSpec grating wheel (GWA) keywords [#5844]

- Moved functions in ``dqflags`` and ``dynamic_mask`` to ``stcal`` [#5898]

- API change - ``stcal.dqflags.interpret_bit_flags`` and ``stcal.dynamicdq.dynamic_mask``
  now require the ``mnemonic_map`` as input. [#5898, #5914]

- Implemented new data models ``SpecKernelModel``, ``SpecProfileModel``,
  ``SpecTraceModel``, and ``WaveMapModel`` for use by new NIRISS SOSS
  reference files in optimized 1D extraction [#5925]

- Added ``FULLP`` to SUBARRAY enum list in core, subarray,
  and keyword_psubarray schemas [#5947]

- Moved JWST_[XYZ] and JWST_[DXDYDZ] keywords from primary to SCI extension
  header and updated their comment fields to indicate they'll now be in the
  barycentric frame. Also added the new OBSGEO[XYZ] keywords to the SCI
  extension header, which are in the geocentric frame. [#6050]

- Added a new datamodel, ``SegmentationMapModel`` that has an uint32 data array
  for storing the segmentation map output from ``source_catalog``. [#6051]

documentation
-------------

- Update documentation, deprecating primary use of CFG files [#5901]

- Update pipeline introduction document to include segmentation map (``segm``)
  in list of data products [#5956]

- Update ``assign_mtwcs`` step docs and reference the ``assign_mtwcs`` step in the
  ``calwebb_image3`` and ``calwebb_spec3`` pipeline docs [#6024]

extract_1d
----------

- Implemented error and variance propagation for all modes but those
  utilizing IFU cubes [#6014]

extract_2d
----------

- For WFSS removed setting srctype to UNKNOWN; added setting ``is_star`` in slitmeta [#5788]

- In NRC_TSGRISM mode replaced FITS WCS keywords with JWST specific ones. [#6005]

- Added ``specsys`` to slits. [#6005]

- Added the step parameter ``wfss_nbright`` to allow for only the N brightest
  objects to be extracted from WFSS exposures. Also changed the name of the
  ``mmag_extract`` param to ``wfss_mmag_extract``, for consistency with other
  WFSS-specific params. [#6788]

general
-------

- Update file naming conventions documentation to clarify when optional components
  will be used. [#5796]

- Update DQFLAGS table in RTD docs with new definitions for persistence and
  ad_floor in bits five and six [#5815]

- Update data products, ``calwebb_image3``, and ``source_catalog`` docs to include
  information about the segmentation map product [#5949]

- Replace documentation references to ambiguous class names with full
  paths. [#6017]

jump
-----------------

- Update the step to detect jumps in three and four group integrations [#5915].

- Change the default S/N ratio for not flagging neighbors to be a higher value to
  better reflect the correct IPC.

lib
---

- Update ``update_mt_kwds`` function in ``set_telescope_pointing.py`` to  populate the TARG_RA/TARG_DEC [#5808]

- moved ``basic_utils.multiple_replace`` to stcal. [#5898]

- Implemented window clipping algorithm for WFSS contamination corrections. [#5978]

- Updated ``set_velocity_aberration`` and ``utc_to_tdb`` to access the JWST
  position and velocity keywords from the SCI extension header, rather than the
  primary header. [#6050]

master_background
-----------------

- Updated documentation to more fully describe the various ways in which the
  step is applied [#5913]

outlier_detection
-----------------

- Outlier detection on non-dithered images is implemented with a simple sigma
  clipping, dithered outlier detection cleaned up and HST specific steps removed
  and additional tests added. [#5822]

ramp_fitting
------------

- Refactoring OLS code for ramp fitting to improve readability and maintenance.
  Also, reference to ``nreads`` is being removed and replaced with ``ngroups``
  to remove and confusion on functionality. [#5872]

- Refactoring ramp fit code separating OLS and GLS code into their own file. [#5951]

- Refactoring ramp fit code in preparation for moving code to STCAL. [#6010]

- Moved ramp fit code to STCAL. [#6023]

- Now that ramp fitting has been moved to STCAL, for the JWST unit tests to
  pass need to use STCAL 0.2.1 or greater.  The bug fix for JP-1920 were made
  in STCAL, which affected JWST unit tests for ramp fitting. [#6038]

refpix
------

- Added code to handle NIR subarrays that use 4 readout amplifiers.  Uses and
  applies reference pixel signal from available amplifiers and side reference
  pixel regions, including odd-even column separation if requested [#5926]

- Fixed a bug introduced in #5926 that affected refpix calibration of 1-amp NIR
  subarrays [#5937]

- Added regression test and unit test for NIR 4-amp subarray correction [#5967]

resample
--------

- Fix ``resample_spec`` output size from input images crossing RA=0 [#5929]

- Propagate variance arrays into ``SlitModel`` used as input for ``ResampleSpecStep`` [#5941]

- Remove certain WCS keywords that are irrelevant after resampling. [#5971]

- Propagate error and variance arrays in ``ResampleStep`` for imaging data. [#6036]

- Propagate error and variance arrays in ``ResampleSpecStep`` for 2D spectral data [#6041]

- Record ``pixel_scale_ratio`` and ``pixfrac`` from ``ResampleStep`` in header
  keywords PXSCLRT and PIXFRAC, respectively, or ``meta.resample.pixel_scale_ratio``
  and ``meta.resample.pixfrac``. [#6044]

source_catalog
--------------

- Updated the concentration indices to be calculated as flux ratios
  instead of magnitude differences. The CI column names have also been
  renamed to list the larger EE first, e.g. ``CI_50_30``. [#5810]

- Aperture-corrected total fluxes and magnitudes are now computed for
  all sources. [#5996]

- Photometric errors are now computed using the new resampled total
  error array. [#5997]

- The ``nn_dist`` column was replaced by a ``nn_label`` column
  indicating the label number of the nearest neighbor. [#5998]

- The ``is_star`` column was replaced by a ``is_extended`` column with
  inverted boolean values. [#6018]

- Circular aperture sizes now scale in the case of non-native pixel
  scales in the resampled image. [#6045]

- Segmentation map output dtype is now ``uint32`` [#6051]

srctype
-------

- Added section for WFSS mode data to set srctype based on ``is_star`` value [#5788]

transforms
----------

- Added ``is_star`` to GrismObject [#5788]

tweakreg
--------

- Fix a bug due to ``models_grouped`` now returning ``odict_values`` instead
  of lists. [#6022]

- Updated documentation to include the new "rshift" option for fit geometry [#5899]

wfss_contam
-----------

- Implemented basic step structure to apply WFSS contamination corrections, along with
  the necessary grism library modules [#5508]


1.1.0 (2021-02-26)
==================

assign_wcs
----------

- Added spectral frames to the output WCS frame of TSO and WFSS observations. [#5771]

associations
------------

- Ignore duplicate product names while handling Level 2 associations [#5780]

- Constraint added to Asn_Lv3Coron to remove background exposures [#5781]

extract_1d
----------

- Determine the background using sigma clipping of entire extended region for
  extended source IFU data [#5743]

resample
--------

- Make inverse variance ``weight_type="ivm"`` the default weighting scheme for
  multiple exposures resampled into a single output. [#5738]


1.0.0 (2021-02-22)
==================

assign_mtwcs
------------

- Fixed a bug which caused the step to fail with ``MultiSlitModel`` input. [#5758]

assign_wcs
----------

- Added velocity aberration-corrected frame ``'v2v3vacorr'`` to the WCS
  pipeline which takes into account DVA effects. [#5602]

- Renamed MIRI frame ``'V2_V3_spatial'`` to ``'v2v3_spatial'`` and
  ``'V2_V3_vacorr_spatial'`` to ``'v2v3vacorr_spatial'``. Added axes names
  to the ``'v2v3'`` frame for ``nircam``, ``niriss``, ``miri``, and ``fgs``.
  Renamed axes for ``nirspec`` from ``V2`` and ``V3`` to
  ``v2`` and ``v3``. [#5765]

- Changed units of the ``'v2v3'`` frame for ``nircam`` from ``u.deg`` to
  ``u.arcsec`` [#5765]

associations
------------

- Warn about duplicate product names and do not write duplicate associations [#5721]

- Added new Lvl2 rule, Asn_Lv2NRSLAMPImage, to run Image2 pipeline for NRSLAMP
  exposures with OPMODE=image [#5740]


combine_1d
----------

- Pull source_id from input x1d headers (from source_catalog) to populate
  c1d output headers [#5759]

cube_build
----------

- Added support for cross-dichroic configurations [#5722]

- Added infrastructure to support NIRSpec opaque + grating options to build lamp mode data [#5757]

- When building MIRI internal_cal type cubes removed the requirement that cdelt1=cdelt2 [#5757]


datamodels
----------

- Updated keyword_readpatt, core, preadpatt schemas for new MIRI detector
  readout patterns 'FASTR1', 'FASTR100' and 'SLOWR1' [#5670]

- Added extr_x and extr_y to multispec datamodel. These values are center
  of extraction region for IFU data [#5685]

- Added segmentation map output file name to core schema keywords, under
  keyword 'SEGMFILE' [#5730]

- Added '1LOS' to PATTTYPE enum list in core.schema datamodel [#5728]

- Added 'IMAGE' to OPMODE enum list [#5745]

- Added source_id to combinedspec and multicombinedspec schemas to populate
  combine1d output headers [#5759]

extract_1d
----------

- Adding writing SRCTYPE, EXTR_X, and EXTR_Y to extracted spec for IFU data [#5685]

- Only update the output x1d data using the PRIMARY input data. Prevents SCI data in x1d data [#5694]

- Fixed bug in background region fitting for image columns/rows that have zero weight
  for all pixels [#5696]

group_scale
-----------

- Fix premature model closing in group_scale_step [#5692]


lib
---

- Make EngDB_Value public for JSDP use [#5669]

- Update code in ``set_velocity_aberration.py`` functions based on Colin Cox
  suggestions: simplify DVA scale computation and improve apparent ``RA`` and
  ``DEC`` aberrated position computation. Also, attributes ``ra_offset`` and
  ``dec_offset`` of ``datamodel.meta.velocity_aberration`` have been renamed to
  ``va_ra_ref`` and ``va_dec_ref`` and their corresponding FITS keywords
  have been renamed from ``DVA_RA`` and ``DVA_DEC`` to
  ``VA_RA`` and ``VA_DEC``. [#5666]

- Make get_wcs_values_from_siaf public for JSDP use [#5669]


outlier_detection
-----------------

- Remove hard-coded MRS outlier detection values now that a parameter reference
  file exists. [#5753]

photom
------

- Fixed handling of NIRSpec IFU extended source data, so that the flux
  calibration gets converted to surface brightness [#5761]


pipeline
--------

- Empty remaining cfg files of any content [#5766]

- Remove references to Numpy globals ``np.int``, ``np.float``, ``np.bool`` and
  ``np.str`` in the package. [#5769]


ramp_fitting
------------

- Fixed bug in handling NGROUPS=2 exposures for pixels that saturate in group 2.
  Proper slope, err, and other quantities are now computed from the good data
  in group 1. [#5700]

- Update documentation to define optimal weighting algorithm [#5682]

source_catalog
--------------

- Added the segmentation map as an output data file, with
  suffix "segm". [#5730]

srctype
-------

- Changed default SRCTYPE for non-primary NIRSpec slits in a FIXEDSLIT
  exposure to 'EXTENDED' rather than 'POINT' [#5671]

- Changed logic for handling NIRSpec MOS exposures to blank out the "global"
  value of SRCTYPE, to ensure that only the individual slit-specific values
  of SRCTYPE get used downstream. [#5754]

stpipe
------

- Make jwst.stpipe independent of the rest of the jwst package and move
  core code to spacetelescope/stpipe. [#5695, #5720, #5752]

0.18.3 (2021-01-25)
===================

- Update documentation introduction to include installation and CRDS setup
  instructions. [#5659]

combine1d
---------

- Fixed code error in combine1d, creating extensions per spectral order
  with the same input data [#5644]

ramp_fitting
------------

- Fix a bug in estimating the max number of segments that will be needed
  to fit any pixel [#5653]

set_telescope_pointing
----------------------

- Update the check in set_telescope_pointing that determines whether an
  exposure is TSO mode to always consider hardwired TSO EXP_TYPEs as TSO,
  regardless of TSOVISIT and NINTS settings. [#5657]

white_light
-----------

- Fixed error causing multi-segment data to reject int_times
  for MJDs [#5566]


0.18.2 (2021-01-19)
===================

associations
------------

- JWSTDMS-410 Asn_Lv2NRSLAMPSpectral: Break out the negative cases [#5635]

- Update MIRI LRS-Fixedslit ALONG-SLIT-NOD backgrounds strategies [#5620]

cube_build
----------

- Do not allow variables defined in spec (part of the cube_build_step class) to
  be changed, to allow calspec2 to loop over a list of files and run the
  pipeline. [#5603]

datamodels
----------

- Updated schemas for new keywords CROWDFLD, PRIDTYPE, PRIDTPTS, PATTNPTS, SMGRDPAT,
  changed name of SUBPXPNS to SUBPXPTS, and new allowed values for PATTTYPE. [#5618]

flat_field
----------

- Added DO_NOT_USE to pixels flagged as NON_SCIENCE for non-NIRSpec data [#5601]

outlier_detection
-----------------

- Account for the background subtracted data in the blot image for determining
  the noise image used in flagging outliers [#5601]

set_telescope_pointing
----------------------

- Updated to populate XREF_SCI, YREF_SCI keywords for all exposures with
  TSOVISIT=True, not just NRC_TSGRISM mode. [#5616]

0.18.1 (2021-01-08)
===================

combine1d
---------

- Output FITS now contains separate combine1d extensions for each spectral
  order present in the data [#5204]

cube_build
----------

- Tweaked pixel wavelength preselection range to avoid truncation at the ends
  of the cubes. [#5598]

datamodels
----------

- Fix missing CHANNEL entry in distortion reffile schema. [#5553]

extract_1d
----------

- For IFU data (NIRSpec and MIRI) the extraction radius is now a varying size
  based on wavelength. The apcorr correction is a function of wavelength and
  radius size. Fixes a bug in units conversion for applying the apcorr correction.
  The units are now correctly converted from arcseconds to pixels. Added an
  new method to apply the apcorr correction for IFU data. [#5506]

pipeline
--------

- Removed all unnecessary parameter settings from cfg files for all steps
  and pipelines, and removed references to step config files from most
  pipeline modules (only kept those that are necessary for intended
  functionality). [#5574]

skymatch
--------

- Fixed a bug due to which sky matching may fail under certain circumstances
  such as using 'mode' statistics on a single pixel (after sigma-clipping). [#5567]

stpipe
------

- Removed unused LinearPipeline class. [#5590]

wavecorr
--------
- Fixed bugs in wavecorr. [#5570]

0.18.0 (2020-12-21)
===================

ami
---
- Update code to use two new input parameters: psf_offset,rotation_search [#5548]

- Update code and unit tests to use new ami_analyze algorithms [#5390]

- Update ami_analyze to extract a SUB80 subarray from full-frame images [#5437]

assign_wcs
----------

- Add nrs_verify to the NIRSpec exposure list [#5403]

- Enable resample_spec for NIRSpec line lamp exposures [#5484]

- Added SIP approximation to WCS for imaging modes. FITS WCS keywords added to meta.wcsinfo. [#5507]

- Fix bug where subarray bounding boxes were 1 pixel too small. [#5543]

- Mark Nirspec slits which project on less than one pixel as invalid. [#5554]

associations
------------

- Asn_Lv2WFSS: Add segmentation map exposure to Level2 WFSS associations [#5532]

- Add new dither keyword subpxpts to constraints [#5525]

- Add further constraints to rule Asn_Lv2NRSLAMPSpectral such that associations
  are created only when LAMP is on and OPMODE indicates a valid optical path. [#5496]

- Restrict association creation based on optical path for NIRSpec Fixed-slit and IFU [#5504]

- Asn_Lv3SpecAux: Add optical element constraint [#5479]

- Add utility asn_gather [#5468]

- Do not allow target acqs to be considered TSO [#5385]

- Add NRS_VERIFY to the list of target acq/confirmation images [#5395]

cube_build
----------

- When making SINGLE type cubes for outlier detection or mrs_imatch data not in the
  appropriate channel/grating is skipped [#5347]

- If outlier detection has flagged all the data on a input file as DO_NOT_USE, then
  skip the file in creating an ifucube [*5347]

- Refactor DataTypes handling of ModelContainer. [#5409]

datamodels
----------

- Skip serializing `None` in datamodels to be compatible with `asdf>=2.8` [#5371]

- Implement full class deprecator decorator and use for MIRIRampModel [#5382]

- Add NRS_VERIFY to the core schema as an allowed EXP_TYPE [#5395]

- Remove logging from DataModel.close [#5413]

- Updated keyword schemas for EXP_TYPE and MODULE, to keep in sync with the
  JWST Keyword Dictionary [#5452]

- Added flatfield and photom correction arrays to slit data models [#5460]

- Move core ``jwst.datamodels`` code to ``stdatamodels`` package and add it as
  an install dependency [#5433]

- Update schemas to include new allowed SUBARRAY values for FGS ASIC tuning
  modes [#5531]

- Add meta.visit.pointing_engdb_quality entry to correspond to ENGQLPTG keyword [#5556]

- Update Moving Target CHEBY table extension schema for changes to column
  definitions in the JWSTKD and SDP [#5558]

- Update distortion reference file schema to have ``meta.instrument.channel``
  keyword [#5553]

extract_1d
----------

- Fixed bug involving the determination of source RA/Dec for resampled Slit
  data. [#5353]

- Updated to use an EXTRACT1D reference file for NIRCam TSGRISM exposures;
  added step param "bkg_fit" to allow for mean and median options in background
  computation, in addition to the existing polynomial fit; fixed bug in
  background computation that was preventing background subtraction from
  ever happening. [#5414]

- Fixed bug involving the processing of WFSS observations when there's only
  one spectrum instance for a given source. [#5439]

fits_generator
--------------

- Addressed deprecated get_children method of XML parser.  Changed type of PATTSIZE from
  float to string in templates. [#5536]

flatfield
---------

- Fixed bug in sending NIRSpec AUTOWAVE exposures to the spectroscopic
  processing branch. [#5356]

- Updated branch logic to handle NRS_LAMP exposures as spectroscopic. [#5370]

- Updated NIRSpec fixed-slit processing to compute and save correction
  values for both point and uniform sources in the primary slit when it
  contains a point source, in order to support master background corrections.
  [#5462]

jump
----

- Fixed bug in the minimum number of groups per integration for the jump
  detection step by changing it from 3 to 5. [#5376]

- Various rework to reduce memory usage and increase readability. [#5404]

master_background
-----------------

- Update the NIRSpec MOS master background logic to only proceed with processing
  after verifying that there are both background and source slits available in
  the input dataset. [#5370]

outlier_detection
-----------------

- Implement memory check in resample to prevent huge arrays [#5354]

photom
------

- Updated NIRSpec fixed-slit processing to compute and save correction
  values for both point and uniform sources in the primary slit when it
  contains a point source, in order to support master background corrections.
  [#5463]

pipeline
--------

- Update ``Image3Pipeline`` to allow sky subtraction when input contains
  only one image (group). [#5423]
- Enable resample_spec for NIRSpec line lamp exposures in Spec2Pipeline [#5484]

ramp_fitting
------------

- Update to store output as an `IFUImageModel` for NIRSpec AUTOWAVE exposures
  using the IFU mode. [#5356]

- Update to add 'DO_NOT_USE' DQ flag to pixels with all groups flagged as
  saturated. [#5367]

resample
--------

- Implement memory check in resample to prevent huge arrays [#5354]

- Add ``pixel_scale_ratio`` parameter to allow finer output grid. [#5389]
- Enable resample_spec for NIRSpec line lamp exposures [#5484]

reset
-----
- Turn the step back on for the calwebb_detector1 pipeline [#5485]

saturation
----------

- Set saturation threshold to A-to-D limit of 65535 for pixels flagged with
  NO_SAT_CHECK in the saturation reference file, instead of skipping any
  test of those pixels. [#5394]
- Flag groups values below A/D floor (0 DN) (#5422)

set_telescope_pointing
----------------------

- Add logging of the found quaternion information [#5495]
- Handle cases where engineering database's pointing mnemonics are all zero over the requested time range [#5540]
- Set value of keyword ENGQLPTG to CALCULATED or PLANNED depending on whether pointing telemetry was used to
  update the WCS [#5556]

skymatch
--------

- Fix a bug in ``skymatch`` that would result in a crash when ``skymethod``
  contains ``'global'`` and the *single image group*'s sky cannot be computed
  (e.g., because all pixels are flagged as "bad"). [#5440]

stpipe
------

- Implement utility function all_steps and fix crds reference file retrieval for non-datamodels [#5492]

tso_photometry
--------------

- Place aperture using header keywords XREF_SCI and YREF_SCI instead of
  CRPIX1 and CRPIX2 [#5533]

- Fixed the flux units in the output photometry catalog. [#5529]

tweakreg
--------

- Add support for the new ``fitgeom`` mode: ``'rshift'`` that can fit only
  for shifts and a rotation. [#5475]

wfs_combine
-----------

- Add checking for bad pixels by using DO_NOT_USE rather than DQ>0. [#5500, #5519]

white_light
-----------

- Add support for step parameters ``min_wavelength`` and ``max_wavelength`` to modify
  the wavelength region over which the flux integration is calculated. [#5501]

0.17.1 (2020-09-15)
===================

associations
------------

- Add product name override to the `IFUGratingBkg` class, to prevent the default
  "clear" suffix showing up in NIRSpec IFU product names. [#5326]

barshadow
---------

- Implement using a user-supplied correction which overrides all references. [#5302]

- Implement applying the inverse operation. [#5302]

blendmeta
---------

- Do not close files that were not opened by blendmodels [#5299]

cube_build
----------

- If every wavelength plane of the IFU cube contains 0 data, cube_build is skipped [#5294]

- Remove "clear" suffix from MIRI MRS product name templates [#5326]

flat_field
----------

- Update how the flat field reference dq mask is used for NIRSpec MOS data [#5284]

- Implement providing a user-supplied flat field which overrides all references. [#5302]

- Implement applying the inverse operation. [#5302]

master_background
-----------------

- Create new step `MasterBackgroundNrsSlits` step to handle NIRSpec MOS data in `Spec2Pipeline` [#5317]

- Implement option to save the 2d version of the calculated master background [#5317]

outlier_detection
-----------------

- Fix bug where background was being subtracted on the input data [#4858]

pathloss
--------

- Implement using a user-supplied correction which overrides all references. [#5302]

- Implement applying the inverse operation. [#5302]

photom
------

- Implement using a user-supplied correction which overrides all references. [#5302]

- Implement applying the inverse operation. [#5302]

pipeline
--------

- Spec3Pipeline check whether master background subtraction has already occurred. [#5308]

- Implement master background subtraction in Spec2Pipeline for NIRSpec MOS data. [#5302]

- Include the per-slit failure traceback in any RuntimeError raised in Spec2Pipeline. [#5315]

scripts
-------

- Add pointing analysis commands v1_calculate and pointing_summary. [#5311]

stpipe
------

- Do not attempt prefetch on pipelines that are set to not allow prefetch. [#5363]

ramp_fitting
------------

- Reinstate copying of INT_TIMES table to output rateints product for TSO exposures. [#5321]

tso_photometry
--------------

- Fix a bug in the computation of integration time stamps when the INT_TIMES
  table is not available. [#5318]

0.17.0 (2020-08-28)
===================

align_refs
----------

- Add bad pixel replacement for target and psf images [#4973]

assign_mtwcs
------------

- Skip the step if any input MT_RA/DEC keyword values are missing. [#5015]

assign_wcs
----------

- Enabled ``filteroffset`` correction for NIRISS and NIRCAM imaging modes. [#5018, #5027]

- Pass an optional ``input_frame`` parameter in ``assign_wcs.util.wcs_from_footprints``. [#5120]

- Improved calculation of bounding boxes in grism images. [#5122]

- Added two new optional parameters to ``utils.cerate_grism_bbox`` - ``wfss_extract_half_height``
  and ``wavelength_range``. [#5140]

- Shifted the bounding box of a resampled WCS by - 0.5 px to account for the
  center of the pixel. [#5241]

- Enable NIRSpec lamp processing in calspec2 pipeline. [#5267]

associations
------------

- Update diagrams in documentation to change sloper to detector1. [#4986]

- Update level-3 rules to exclude IFU exposures from ``calwebb_tso3`` associations. [#5202]

- Fix formatting error in Asn_IFUGrating product name construction. [#5231]

barshadow
---------

- Correct bar shadow parity bug for yslit. [#5095]

combine_1d
----------

- Skip spectra that are degenerate when combining [#5037]

cube_build
----------

- Changed default weighting to 'emsm'. [#5277]

- Fixed formatting of NIRSpec s3d output product names. [#5231]

- Modified NIRSpec blotting to the find min and max ra and dec for each slice and only
  invert those values on slice that fall in range [#5144]

- Changed default weighting back to 'msm' until NIRSPEC cube pars ref file contains emsm info [#5134]

- Added checks read from cube pars reference file that parameters have valid data [#5134]

- Change the name of default cube type from ``world`` to ``skyalign`` [#4974]

- Add ``ifualign`` cubes to be cubes rotated on sky to align with ifu instrument plane [#4974]

- Change the name of MIRI ``alpha-beta`` cube type to ``internal_cal`` [#4974]

- Add ability to make NIRSpec ``internal_cal`` ifu cubes aligned with slicer plane [#4974]

- Change default weighting from ``msm`` to ``emsm`` [#4974]

- NIRSpec IFU cubes built from all wavelengths rather than those defined in cube par ref file [#4974]

- Removed wavelength planes that contained only 0 data. These planes are edge cases [#4974]

datamodels
----------

- Add iscopy to ModelContainer init [#5256]

- Re-enable FITS-hash by default. [#5191]

- Add blend rule for keywords DETECTOR and MODULE. [#4998]

- Add methods ``Model.info`` and ``Model.search``. [#4660]

- Trim MT_RA, MT_DEC keyword comments to fit within FITS record. [#4994]

- Add enum list and default value of 'NONE' for ``meta.instrument.lamp_mode`` [#5022]

- Add TIMEUNIT keyword to schemas. [#5109]

- Split ``pathloss`` object into ``pathloss_ps`` and ``pathloss_un`` in schemas. [#5112]

- Add "PERSISTENCE" DQ flag definition. [#5137]

- Fix nonsensical premature closing of FITS file of a ``DataModel``. [#4930]

- Add a hash set/check to DataModel I/O to check whether schema traversal is necessary. [#5110]

- Update underlying MultiExposureModel from the SourceModelContainer models. [#5154]

- Add new MIRI LRS dither patterns to PATTTYPE enum list. [#5254]

extract_1d
----------

- Implement aperture corrections in the Extract1dStep. [#4902]

- Fix bug in creating a polynomial fit used in background extraction. [#4970]

- Recheck the input model container in run_extract1d to select the correct processing [#5076]

- Rework/refactor many functions for style and readability. [#5079]

- Checks subwcs and new_slit variables exist before trying to delete them. [#5093]

- Parameter ``mmag_extract`` is now propagated to the extraction routine. [#5122]

- Updated the logic for when and how to use the source position to offset the
  location of the extraction regions specified in the EXTRACT1D reference file. [#5157]

- Fixed the conversion of flux to surface brightness for IFU extended source case [#5201]

- Fixed bugs in aperture correction for NIRSpec multi-slit modes. [#5260]

extract_2d
----------

- Check that ``subwcs`` and ``new_slit`` variables exist before trying to delete them [#5093]

- Move NIRSpec wavecorr routines to the ``wavecorr`` step. [#5133]

- Added a new optional integer parameter to extract_2d (``wfss_extract_half_height``)
  which allows a user to specify the extraction height in the
  cross-dispersion direction for WFSS mode. [#5140]

flat_field
----------
- For NIRSpec BOTS and ALLSLITS add the slit start corner to the subarray start corner
  when determining what region of the flat_field reference files to extract. [#5269]

- Enable NIRSpec lamp processing in calspec2 pipeline. [#5267]

fringe
------

- Update the fringe step to handle 3D inputs for MIRI MRS TSO mode. [#5202]


master_background
-----------------

- Fix open files bug [#4995]

- Update to include pathloss corrections to NIRSpec IFU background [#5125]

mrs_imatch
----------

- MRSIMatchStep to create its ModelContainers with `iscopy=True` [#5256]

outlier_detection
-----------------

- Update median filter to use numpy's nanmedian. [#5114]

- Fix outlier_detection bug when saving intermediate results. [#5108]

- Update logic to correctly handle input ``CubeModel`` that have only
  1 integration. [#5211]

pathloss
--------

- Fix bug in NIRSpec IFU data that causes valid pixel dq flags to set to
  NON-SCIENCE in the region of an overlapping bounding box slice [#5047]

- Update to save both point source and uniform source 2D pathloss correction
  arrays to output. [#5112]

persistence
-----------

- Flag pixels with high persistence using "PERSISTENCE" DQ flag instead
  of "DO_NOT_USE". [#5137]

pipeline
--------

- Refactor the ``Image3Pipeline`` to use ``stpipe`` infrastructure. [#4990]

- Fix ``Coron3Pipeline`` to blend headers just from each input science model,
  not every integration. [#5007]

- Fix open files bug in ``get_config_from_reference`` class method, and in
  ``Spec2Pipeline``, ``Spec3Pipeline`` and ``tso3``. [#4995]

- Update ``calwebb_tso3`` to do more robust checking of input data type.
  [#5107]

- Update the ``Spec2Pipeline`` to include the new ``wavecorr`` step and put
  ``srctype`` before ``wavecorr``. [#5133]

- Update the ``Spec2Pipeline`` to skip ``extract_1d`` for IFU data that
  have not had a cube built (e.g. MIRI MRS TSO), and update the
  ``calwebb_tso-spec2.cfg`` configuration to turn on the ``fringe`` step
  and turn off ``cube_build`` for MIRI MRS TSO. [#5202]

- Update the ``Coron3Pipeline`` logic to correctly handle inputs that have
  only 1 integration. [#5211]

- Refactor Spec2Pipeline for execution logic and step flow isolation [#5214]

- Update ``Ami3Pipeline`` to only process psf and science members from the
  input ASN. [#5243]

- Enable NIRSpec lamp processing in calspec2 pipeline. [#5267]

photom
------

- Fix bug in NIRSpec IFU data that causes valid pixel dq flags to set to
  NON-SCIENCE in the region of an overlapping bounding box slice [#5047]

ramp_fitting
------------

- Add multi-processing capability. [#4815]

- Fix crash when DRPFRMS1 is not set [#5096]

- Update to always create the rateints product, even when NINTS=1. [#5211]

resample_spec
-------------

- Fix artifacts in resampled NIRSpec slit data caused by NaNs in the WCS [#5217]

source_catalog
--------------

- Use ``gwcs.WCS`` instead of FITS WCS. [#5120]

- Changed the type of column ``is_star`` from float to bool. [#5140]

- Implemented algorithm for determining whether a source is a star.
  [#5234]

stpipe
------

- Limit reference file prefetch to the first "science" exptype
  when a pipeline has an association as input. [#5031]

- Remove further sloper references. [#4989]

- Enable prefetch of pars reference files for associations. [#5249]

transforms
----------

- Wrap first spherical angle ("RA") at 360 degrees in the forward ``V23ToSky``
  transformation and to 180 degrees for the inverse transformation ("V2").
  This is now done using models defined in ``astropy`` and ``gwcs`` packages
  replacing ``V23ToSky`` model in JWST's WCS pipeline. [#5206]

wavecorr
--------

- Implemented the ``wavecorr`` step by pulling routines from the
  ``extract_2d`` step. [#5133]

0.16.2 (2020-06-10)
===================

- Fixed ``packaging`` dependency installation issue.  [#4977]


0.16.1 (2020-05-19)
===================

assign_wcs
----------

- Update keyword and attribute usage around SkyObject to reflect updated keywords. [#4943]

- Refactor PPS origin of NIRSpec MOS shutters from top left to bottom left. [#4959]

associations
------------

- Modify NIRSpec IFU level-3 ASN rules to include only one grating per association [#4926]

calwebb_coron3
--------------

- Update coron3 for new outlier detection application [#4968]

datamodels
----------

- Add ``to_container`` to ``CubeModel`` to convert a cube to a list of images [#4968]

- Add ``getarray_noinit`` to ``DataModel`` to access arrays without causing initialization [#4968]

- Limit looping over HDU's while resolving arrays in schema [#4951]

- Relax asdf requirement and use validator flag when asdf 2.6.x is installed [#4905]

- Updated core schema to include recent Keyword Dictionary changes
  (remove TIME-END; add TDB-BEG, TDB-MID, TDB-END, XPOSURE, TELAPSE)
  [#4925]

- Populate meta.asn.table_name when an association is loaded into a
  ``ModelContainer``. [#4873]

extract_1d
----------

- Add aperture correction in extract_1d processing. [#4902]

lib
---

- Update SkyObject keys. [#4943]

mrs_imatch
----------

- Fix ``mrs_imatch`` to avoid calls to ``sigma_clipped_stats`` with all-zero
  arrays. [#4944]

photom
------

- Fix flux units in photom for MultiSlit cases. [#4958]

pipeline
--------

- Updated calwebb_image3 pipeline to only load science and background member
  types from an input ASN. [#4937]

- Updated the calwebb_spec2 pipeline to only use the basename of the source
  catalog file when updating the source_catalogue keyword for WFSS inputs.
  [#4940]

rscd
----

- Fixed bug when the READPATT/SUBARRAY data is not found in RSCD reference file [#4934]

source_catalog
--------------

- Add more concentration indices and update step docs. [#4906, #4908]

- Added fallback background estimation method to make background
  estimation moare robust. [#4929]

- Fixed the nearest-neighbor code to handle the case of exactly one
  detected source. [#4929]

- Update abmag error calculation. [#4945]

- Exit gracefully if APCORR ref file is missing. [#4948]

tweakreg
--------

- Added align_to_gaia processing as an option [#4599]



0.16.0 (2020-05-04)
===================

ami
---

- Reorganized step documentation [#4697]

assign_wcs
----------

- Updated MIRI imaging distortion to use new filteroffset file format [#4776]

associations
------------

- Update asn_from_list to have default values in the asn header [#4720]

- Update rules so exclude dark files from associations [#4668]

- Update association rules so that nodded observations produce level 3 asn's [#4675]

cmdline
-------

- Re-enable exception tracebacks from strun for issues outside step processing [#4761]

coron
-----

- Reorganized step documentation [#4697]

datamodels
----------

- Update schemas to add moving_target_position and cheby tables to the level1b
  schema [#4760]

- Deprecate ``DrizProductModel`` and ``MultiProductModel`` and replace with
  updated versions of ``ImageModel`` and ``SlitModel`` that include "CON" and
  "WHT" arrays for resampled data. [#4552]

- Remove lev3_prod schema and move resample-related keywords to
  core schema. [#4552]

- Add data models for spectroscopic mode APCORR reference files. [#4770]

- Added ``pupil`` to the ``FilteroffsetModel`` to support NIRCAM and NIRISS WCS. [#4750]

- Removed old MIRI-specific filteroffset schema.  [#4776]

- Added FASTGRPAVG[8,16,32,64] to the READPATT keyword allowed values. [#4818]

- Added the SRCTYAPT keyword and moved SRCTYPE to the SCI extension header of
  all applicable data model schemas. [#4885]

exp_to_source
-------------

- Resulting MultiExposureModels are now updated with header information from the inputs. [#4771]

extract_1d
----------

- Updates for handling resampled input data as ``ImageModel``, ``SlitModel``,
  and ``MultiSlitModel``, instead of ``DrizProductModel`` and ``MultiProductModel``,
  which are deprecated. [#4552]

- Remove pixel-by-pixel calls to wcs; copy input keywords to output for
  more types of input data. [#4685]

- Updated to create a single ``x1d`` product per source for WFSS data, containing
  all extracted spectra for a given source, instead of multiple ``x1d`` files per
  source. [#4846]

extract_2d
----------

- Change the source type for NIRSpec MOS sources with stellarity = -1 from
  UNKNOWN to POINT. [#4686]

master_background
-----------------

- Updated step arguments in the documentation. [#4723]

- Fix issue with files left open at end of step [#4775]

mrs_imatch
----------

- Updated step to use EMSM cube weighting, and to perform iterative sigma
  rejection of sources prior to running the background solver.  [#4732]

outlier_detection
-----------------

- Updated step arguments in the documentation. [#4723]

- Change outlier and resample DQ bit usage.  [#4726]
  Default value of ``good_bits`` now includes all DQ flags except ``DO_NOT_USE``.
  Also, newly flagged outliers are flagged with ``DO_NOT_USE + OUTLIER``.

- Added a hardcoded declaration of a reasonable scale parameter for MIRI MRS as a stopgap
  measure until a parameter reference file can pass one more cleanly. [#4778]

pipeline
--------

- Update ``calwebb_detector1`` to reduce the memory used in processing. [#4643]

- Update ``calwebb_coron3`` to return ``ImageModel`` instead of ``DrizProductModel``,
  when necessary. [#4552]

- Fix issue with files left open at end of ``calwebb_spec2`` [#4775]

- Update ``calwebb_spec3`` to use suffix ``c1d`` for ``combine_1d`` products.
  [#4846]

- Update ``calwebb_spec3`` to update the ASNTABLE keyword in all output
  products, to reflect the name of the spec3 ASN used as input. [#4865]

resample
--------

- Update to return resampled data in an ``ImageModel``, instead of
  ``DrizProductModel``. [#4552]

- Updated documentation to include step arguments and reference file
  description. [#4723]

- Change outlier and resample DQ bit usage.  [#4726]
  The parameter ``good_bits`` has been removed in favor of allowing all
  DQ flags except for ``DO_NOT_USE``

- Updated to reject pixels with DQ flag NON_SCIENCE, in addition to
  DO_NOT_USE. [#4851]

resample_spec
-------------

- Update to return resampled data in a ``SlitModel`` or ``MultiSlitModel``,
  instead of ``DrizProductModel`` or ``MultiProductModel``. [#4552]

- Fix bug that was causing resampled MIRI LRS fixed-slit data to be all zero.
  [#4552]

- Enable model metadata blending [#4765]

rscd
----

- Added baseline algorithm that flags groups [#4669]

set_telescope_pointing
----------------------

- Update to add moving target coords to the header [#4760]

source_catalog
--------------

- Update to use ``ImageModel`` for resampled input data, instead of
  ``DrizProductModel``. [#4552]

- Updated step arguments in the documentation. [#4723]

- Updated to include aperture photometry and aperture corrections. [#4819]

- Rename AB-to-Vega reference file type to ABVEGAOFFSET. [#4872]

srctype
-------

- Change default source type for NRS_IFU from POINT to EXTENDED. Change the source
  type for NIRSpec MOS sources with stellarity = -1 from UNKNOWN to POINT. [#4686]

- Modified the step to use the SRCTYAPT keyword to get the user input value from
  the APT and store the derived source type in the SRCTYPE keyword. [#4885]

stpipe
------

- Unhide exceptions during CRDS steppars retrieval [#4691]

- Add command line and environmental options to not retrieve steppars references [#4676]

- Use only a single member of an association for CRDS STEPPARS checking [#4684]

- Fix handling of the boolean-like environmental variables PASS_INVALID_VALUES and STRICT_VALIDATION [#4842]

strun
-----

- Re-enable exception tracebacks from strun for issues outside step processing [#4761]

tweakreg
--------

- Updated step arguments in the documentation. [#4723]


wfs_combine
-----------

- Update the value of the ASNTABLE keyword in the output ``wfscmb`` product. [#4849]

0.15.1 (2020-03-10)
===================

assign_wcs
----------

- Fix NIRISS WFSS FWPOS angle bugs [#4653]

- Replaced FITS WCS transforms with GWCS transforms in computing bounding boxes of grisms slits. [#4665]

datamodels
----------

- Update schema-editor to match documentation and clarify execution [#4587]

- Remove the init file usage. Way too confusing [#4645]

mrs_imatch
----------

- If the background polynomial contains any Nan Values the mrs_imatch step is skipped [#4642]

stpipe
------

- Revert "JP-1090: Remove setLevel calls (#4621)" [#4667]


0.15.0 (2020-02-28)
===================

assign_wcs
----------

- A ``ValueError`` is now raised if input data is missing ``xref_sci`` or
  ``yref_sci`` keywords. [#4561]

associations
------------

- Cull Association tests [#4610]

- Correct PATTTYPE values in ASN level 3 rules [#4570]

- Update act_id format to allow base 36 values in product name [#4282]

- Refactor association logging configuration [#4510]

combine_1d
----------

- Check output pixel numbers for NaN [#4409]

datamodels
----------

- Update schema-editor to match documentation and clarify execution [#4578]

- Force data model type setting on save [#4318]

- Deprecate ``MIRIRampModel`` [#4328]

- Make ``memmap=False`` be the default in ``datamodels`` [#4445]

- Update schemas to add the ``id`` field and switch relative references
  from filesystem paths to URIs.  Make ``schema_url`` absolute to facilitate
  subclassing DataModel with schemas from other asdf extensions. [#4435]

- Update core.schema.yaml to include new allowed values for PATTTYPE
  [#4475, 4517, 4564]


- DataModel.update() now has ``extra_fits=False`` kwarg that controls whether
  an update happens from the ``extra_fits`` section of the datamodel.  Default
  is to stop doing this by default, i.e. ``False``. [#4593]

- Add units to filteroffset schema.  [#4595]

- Updated ``slitdata.schema.yaml`` to include ``SRCRA`` and ``SRCDEC`` for
  MOS slitlets to FITS SCI headers. These values are taken from the MOS
  metadata file. [#4613]

- Many keyword updates to bring us in-sync with KWD. [#4602, #4627]

- Update schemas to use transform-1.2.0. [#4604]

- Allow FileNotFoundError to be raised. [#4605]

extract_1d
----------

- Updated to work with the current output from photom [#4369]

- Fixed bug regarding background for NIRSpec or NIRISS (SOSS) point source
  spectra. [#4459]

extract_2d
----------

- For GRISM data, the variance arrays and INT_TIMES table are copied to output,
  and keywords SLTSTRT1 and SLTSTRT2 are set to the pixel location of the
  cutout in the input file. [#4504]

- A ``ValueError`` is now raised if the input data is missing ``xref_sci`` or
  ``yref_sci`` keywords. [#4561]

- Fix the WCS subarray offsets for NIRCam TSGRISM cutouts [#4573]

- Added ``source_ra`` and ``source_dec`` to MSA ``Slit`` with values
  from the MSA metadata file. [#4613]

master_background
-----------------

- Updated to fill the asn table and asn pool names. [#4240]

model_blender
-------------

- Do not overwrite rules with defaults. [#4521]

outlier_detection
-----------------

- Check for a zero array before sigma clipping [#4598]

- Fix bug and logic pertaining to detecting if the background has been
  subtracted or not. [#4523]

pipeline
--------

- Hardwire required pipeline outputs in the pipeline. [#4578]

- Added FGS_IMAGE to the exposure types to apply resampling in
  calwebb_image2.py [#4421]

- Make the naming and writing out of the resampled results to an `i2d` file
  in `Image2Pipeline` consistent between config and class invocations [#4333]

- Don't try to save the ``cube_build`` result if the step is skipped in the
  ``calwebb_spec2`` pipeline. [#4478]

- Use the `overwrite` option when saving the white-light photometry catalog in
  the ``calwebb_tso3`` pipeline. [#4493]

- Fixed error in formatting of example ASN file contents in the documents for
  the ``calwebb_coron3`` and ``calwebb_ami3`` pipelines. [#4496]

- Fixed the ``calwebb_tso3`` calculation of the number_of_integrations recorded
  in the photometric table product to avoid ``astropy.table`` merge conflicts.
  [#4502]

photom
------

- Added ``spectral_order`` to the fields matching the ``photom`` reference files
  for NIRCAM WFSS and TSGRISM modes. [#4538, 4558]

refpix
------

- Interchanged alpha and beta reference arrays; use the DQ extension [#4575]

- Fixed bugs in PR #4575; added unit tests [#4596]

- Changed the data type of columns OUTPUT and ODD_EVEN in the section of the
  schema for the DQ table in the NIRSpec IRS2 refpix reference file [#4618]

residual_fringe
---------------

- First implementation of step. Added third party package (BayesicFitting) to setup.cfg [#6211]

- Fixed suffix defined in spec [#6347]

- Fixed an error when a filename was give as the input to apply the residual fringe correction [#6349]

- Updated residual fringe reference data model to support new delivery of reference files [#6357]


set_telescope_pointing
----------------------

- Round S_REGION values in ``set_telescope_pointing`` [#4476]

source_catalog
--------------

- Remove directory path when populating SCATFILE keyword. [#4597]

srctype
-------

- Updated logic to populate SRCTYPE in all slit instances of slit-based
  data models. [#4541]

stpipe
------

- Fix sub-step nesting in parameter reference files [#4488]

transforms
----------

- Removed ``TPCorr`` WCS correction model as it is now defined in ``tweakwcs``
  as a compound model of elementary ``astropy`` and ``gwcs`` models. [#4790]

- Refactored the WFSS transforms to improve performance. [#4603]

- Added ``source_ra`` and ``source_dec`` to the ``Slit`` namedtuple
  with default values of 0.0. These are populated from the MSA metadata
  file. [#4613]

tweakreg
--------

- Improved code to be more resilient to the case when none of the
  image groups has valid sources that can be used for image alignment.
  Now the code will gracefully skip the ``tweakreg`` step altogether in such
  situations. [#4299]

wfs_combine
-----------

- Use float64 data types internally in ``wfs_combine`` so as not to cause an
  error in ``scipy.signal.convolve``. [#4432]

tso_photometry
--------------

- A ``ValueError`` is now raised if the input data for ``call`` is missing
  ``crpix1`` or ``crpix2`` keywords. [#4561]


0.14.2 (2019-11-18)
===================

associations
------------

- Refactor target acquisition handling [#4254]

emission
--------

- Removed the emission step, documentation, and tests from the jwst package.
  [#4253]

photom
------

- Fixed a bug so that the reference table column "PHOTMJ" is used for NIRSpec IFU
  exposures. [#4263]

- The pixel area is now gotten from the photom reference file. [#4270]

white_light
-----------

- Fixed bug which produces NaN results when only some input has NaN [#4256]


0.14.1 (2019-11-11)
===================

associations
------------

- Updated level 3 rules so that target acquisitions in the pool files are listed as
  exp_type = 'target_acquisition', not as science exposures. [#4223]

datamodels
----------

- Updated the list of allowed NIRCam CORONMSK values in model schemas. [#4234]

flat_field
----------
 - Updated handling of error arrays for FGS Guider data, which has not been run
   through ramp fitting [#4309]

lib
---

- Updated the EngDB web service url in ``engdb_tools``. [#4187]

photom
------

- Updated unit tests to use proper names for the MIRI LRS fixedslit
  subarray. [#4205]

pipeline
--------

- Updated ``calwebb_spec3`` to allow for processing of non-TSO
  NIRISS SOSS exposures. [#4194]

resample_spec
-------------

- Updated unit tests for new name of MIRI LRS slitless subarray
  ('SUBPRISM' -> 'SLITLESSPRISM'). [#4205]

rscd
----

- Updated to handle science data and reference files that use the old
  'SUBPRISM' name for the MIRI LRS slitless subarray and update the values
  to 'SLITLESSPRISM'. [#4205]

stpipe
------

- Only allow science members in step parameter reference call [#4236]

- get_pars returns all available parameters for a step and all sub-steps [#4215]

tests_nightly
-------------

- Added a ``set_telescope_pointing`` test for a NIRCam TSGRISM exposure.
  [#4187]

transforms
----------

- Updated all transforms to be consistent with astropy v 4.0.
  Transform classes define now two class variables - ``n_inputs``
  and `n_outputs``. The variables ``inputs`` and ``outputs`` are
  now instance variables (previously they were class variables). [#4216]


0.14.0 (2019-10-25)
===================

- Remove references to deprecated collections.* ABCs that will be removed in
  Python 3.8. [#3732]

- Remove ``jwpsf`` module. [#3791]

- Update dependencies ``python>=3.6`` and ``numpy>=1.16``. [#4134]


ami
---

- Unit tests were added for the ami_analyze pipeline. [#4176]

assign_wcs
----------

- This step populates keyword DISPAXIS. [#3799]

- For NIRISS WFSS data, the wavelengths were incorrect because the function
  for horizontally oriented spectra was called for GR150R, and the function
  for vertically oriented spectra was called for GR150C. [#3891]


associations
------------
- Update level 3 rules to create image3 associations for FGS_IMAGE exposures [#3920]

- Add mir_taconfirm to the target acquisition exp_types [#4135]

- Exclude mir_lrs-slitless calibration data from level 3 processing [#3990]

- Fix in load_as_asn for UTF-8 errors [#3942]

- Update association rules so that MIMF exposures are processed as WFS observations [#4034]

- asn_from_list fills the level2  member exptype correctly if the input is a tuple [#2942]

- Update rules to make level 3 associations for slitless LRS mode [#3940]

- Update rules so that nOPS5 observations with "ALONG-SLIT-NOD" dither
   pattern generates level 3 associations [#3912]

- Update rules to have NRS_IFU backgrounds in science associations [#3824]

- Return filename with extensions based on file type [#2671]

- Ensured that all target acqs are processed by Level 2 [#3765]

- Add a check that backgrounds are included in level 3 associations [#3678]

- Will not constrain on uniqueness of the MSACONFIG keyword [#3770]

- Process non-science exposures taken during WFS&C observations [#3947]

barshadow
---------

- Update barshadow position [#3897]

- Unit tests were added. [#3930]

combine_1d
----------

- Fixed the number of inputs to the spectral WCS - one expected, two were passed. [#3827]

calwebb_tso3
-------------

- Update to exclude target_acquisitions from processing in the calwebb_tso3 pipeline [#3759]

cube_build
----------

- Schema for the ``WAVE-TAB`` WCS no longer requires fixed-length arrays for
  the wavelength "coordinates". The ``'nelem'`` field therefore is no longer
  necessary and has been removed. [#3976]

- To support outlier detection the blotting from the sky back to the detector was
  improved [#4301]

datamodels
----------

- Update to prevent target_acquisitions from processing in the spec3 pipeline [#3777]

- Use public API of jsonschema to ease upgrade to 3.x. [#3705]

- Fixed corruption of FITS tables with unsigned int columns. [#3736]

- Fixed missing TUNITn keywords caused by changes for unsigned int columns. [#3753]

- Write ``siaf_xref_sci`` and ``siaf_yref_sci`` to FITS keywords ``XREF_SCI``
  and ``YREF_SCI`` for ``NRC_TSGRISM`` exposures. [#3766]

- Updated multiexposure.schema to just import slitdata.schema instead of explicitly
  specifying all of its attributes. [#3809]

- Improved ``properties._cast()`` to be able to handle structured arrays
  schemas without a specified (in schema) shape. In addition, ``ndim``
  can be used to constrain the dimensionality of data in structured array
  fields. [#3976]

- Fixed an issue with the fix from [#3976] that was affecting "casting" to
  data types defined by schema of structured arrays when input values are not
  native Python types (tuples). [#3995]

- Fixed an issue with the fix from [#3995] that was affecting "casting" to
  data types defined by schema of structured arrays when input values are
  already structured arrays. [#4030]

- Added "MIR_TACONFIRM" to the list of allowed EXP_TYPE values in the
  keyword schemas. [#4039]

- Added new imaging-specific photom reference file data models ``FgsImgPhotomModel``,
  ``MirImgPhotomModel``, ``NrcImgPhotomModel``, and ``NisImgPhotomModel``. [#4052]

- Add EXP_TYPE and P_EXP_TY keywords to new imaging photom reference file
  data model schemas. [#4068]

- Introduced a flag ``ignore_missing_extensions=True`` to the `DataModel` initializer
  which is propagated to the ``asdf.open`` function. It allows control over a warning
  asdf issues when opening files written with an extension version older than the
  extension version the file was written with. An example message is

  ``asdf/asdf.py:202: UserWarning: File was created with extension
  'astropy.io.misc.asdf.extension.AstropyAsdfExtension' from package astropy-4.0.dev24515,
  but older version astropy-3.2.1 is installed``. [#4070]

- Added new spectroscopic mode photom reference file data models. [#4096]

- Added new imaging mode aperture correction (apcorr) reference file data
  models ``FgsImgApcorrModel``, ``MirImgApcorrModel``, ``NrcImgApcorrModel``,
  and ``NisImgApcorrModel``. [#4168]

- Removed old photom reference file data models. [#4173]

- Add support for streaming reference files directly from S3. [#4170]

exp_to_source
-------------

- Updated the documentation and added some logging to the step. [#3803]

- Close input files after creating the new outputs. [#3828]

extract_1d
----------

- Parameters were added to ``ExtractBase.__init__``, and most of the initialization
  is done there rather than in the subclasses. [#3714]

- This step uses keyword DISPAXIS. [#3799]

- Fixed a bug in ``pixel_area`` when the input is a ``CubeModel``. [#3827]

- Computing the solid angle of a pixel is only done for the first integration
  of a multi-integration exposure, and it's not done at all for WFSS data
  [#3863]

extract_2d
----------

- For grism data, this step copies keyword DISPAXIS from input to output. [#3799]

- For NIRCam TSO data, wavelengths are computed and assigned to the
  wavelength attribute. [#3863]

- Improved the computation of ``S_REGION`` of a slit. [#4111]

flat_field
----------

- For NIRSpec spectroscopic data, the flat_field step needs the dispersion
  direction.  The step now gets that information from keyword DISPAXIS.
  [#3799, #3807]

- The test_flatfield_step_interface unit test in test_flatfield.py has been
  temporarily disabled. [#3997]

gain_scale
----------

- Updated to apply gain factor to variance arrays. [#3794]

group_scale
-----------

- Updates to documentation and log messages. [#3738]

ipc
---

Function is_irs2 has been removed from x_irs2.py.  The version of this function
that is now in lib/pipe_utils.py is used instead. [#4054]

lib
---

- A function to determine the dispersion direction has been added. [#3756]

- Function is_irs2 has been added to pipe_utils.py, and unit tests were
  added to tests/test_pipe_utils.py. [#4054]

master_background
-----------------

- Updated the documentation to include more details. [#3776]

photom
------

- Add unit tests [#4022]

- The code was modified to work with the new photom reference files. [#4118]

- Two bugs were fixed.  For NIRSpec IFU data the code was trying to access
  an attribute of a "slit", but there were no slits for this type of data.
  For NIRISS extended-source data, the code tried to divide by the pixel
  area, but the pixel area was undefined.  [#4174]

- NRS_BRIGHTOBJ data were incorrectly treated the same as fixed-slit, but
  the data models are actually not the same.  Also, the logic for pixel area
  for fixed-slit data was incorrect. [#4179]

refpix
------

- Call is_irs2 from lib/pipe_utils.py instead of using PATTTYPE keyword to
  check for IRS2 readout mode. [#4054]

resample_spec
-------------

- This step uses keyword DISPAXIS and also copies it to output. [#3799]

saturation
----------

Function is_irs2 has been removed from x_irs2.py.  The version of this function
that is now in lib/pipe_utils.py is used instead. [#4054]

stpipe
------

- Fix ``Step.print_configspec()`` method.  Add test.  [#3768]

- Integrate retrieval of Step parameters from CRDS. [#4090]

- Change properties ``Step.pars`` and ``Step.pars_model`` to methods. [#4117]

- Fix bug in ``Step.call()`` where a config file referencing another config
  file was not merged into the final spec properly. [#4161]

- Set ``Step.skip = True`` in ``Step.record_step_status()`` if
  ``success == False``. [#4165]

tests_nightly
-------------

- Some tests in general/nirspec/ were marked as "expected to fail" because
  the new reference files are not being selected. [#4180]

tso_photometry
--------------

- Unit tests were added to tso_photometry. [#3909]

tweakreg
--------

- Fixed a bug in a ``try-except`` block in the ``tweakreg`` step. [#4133]

- removed original ``jwst.tweakreg`` alignment code and changed step's code
  to call similar functionality from ``tweakwcs`` package. [#3689]

- Fix deprecated call to photutils.detect_threshold [#3982]


0.13.7 (2019-06-21)
===================

datamodels
----------

- Reverted #3680 and #3709. [#3717, #3718]

flatfield
---------

- Three new unit tests were added.  Two existing files were modified to
  split the tests into separate functions. [#3704]

0.13.6 (2019-06-20)
===================

associations
------------

- Fixed constraints on WFSC processing. [#3710]

datamodels
----------

- Fixed corruption of FITS tables with unsigned int columns. [#3680]


0.13.5 (2019-06-19)
===================

associations
------------

- Reverted over-restrictive constraints on WFSC processing. [#3691]

- Removed the rule creating associations for NIRSpec LAMP exposures in image modes. [#3693]


0.13.4 (2019-06-17)
===================

assign_wcs
----------

- A unique integer ``source_id`` is now assigned to all MOS background slitlets
  and NRS Fixed Slits. [#3584]

associations
------------

- MIRI MRS dedicated background exposures are now listed as science observations in
  a new association. [#3542]

- Generate will no longer merge Level2 associations by default [#3631]

- Prevent inclusion of data files with exp_type="NIS_EXTCAL" in the association files [#3611]

- Implemented Level 2 re-sequencing to prevent overwriting of associations [#3674]

- Implemented Level 2 background exposure reprocessing [#3675]

combine_1d
----------

The input DQ column is temporarily replaced by a zero-filled array of
the right data type. [#3666]

datamodels
----------

- Changed PATTSIZE keyword data type from float to string. [#3606]

- Added enumeration of allowed values of ``FXD_SLIT`` to the core schema. [#3584]

- Changed ``WHT_TYPE`` keyword to ``RESWHT``. [#3653]

- Add missing pattern/enum values to keyword_pband, keyword_pfilter, keyword_channel [#3653]

- New keywords [#3653]
   - ``DSETSTRT``
   - ``NUMDSETS``
   - ``DITHDIRC``
   - ``DITHOPFR``
   - ``DITHPNTS``
   - ``MRSPRCHN``
   - ``NDITHPTS``
   - ``DWTSCL``
   - ``DOUTUN``
   - ``RESBITS``
   - ``DFVAL``
   - ``DPIXFR``
   - ``DKERN``
   - ``SCIEXT``
   - ``CONEXT``
   - ``WHTEXT``

extract_1d
----------

- Checks for input from a SourceModelContainer. [#3649]

exp_to_source
-------------

- Changed `exp_to_source`` to use ``source_id`` to group exposures. [#3584]

- Removed the enum list for the SUBPXPAT keyword to allow validation of any value. [#3616]

extract_1d
----------

- Checks for input from a SourceModelContainer. [#3649]

extract_2d
----------

- Nircam ``TSGRISM`` extraction uses now ``wcsinfo.siaf_x(y)ref_sci`` as the source position
  on the detector. [#3646]

- For grism data, a wavelength array is computed and saved, and the variance
  arrays are extracted and copied to output. [#3664]

lib
---

- ``set_telescope_pointing`` now retrieves CRPIX1/2 from the SIAF for Nircam TSGRISM
  observations and saves the values as ``meta.wcsinfo.siaf_x(y)ref_sci``. These are used
  by ``extract_2d`` as the source position on the detector. [#3646]

outlier_detection
-----------------

- Changed default value of good_pixel from 4 to 6 [#3638]

- Don't use NaNs or masked values in weight image for blotting. [#3651]

- When calling cube_build for IFU data fixed selecting correct channels (MIRI) or
  correct grating (NIRSPEC) [#4301]

pipeline
--------

- ``calwebb_spec2`` was changed to allow processing of exposures
  with ``EXP_TYPE=NRS_LAMP.`` [#3603]

- ``calwebb_tso3`` was changed to allow processing of exposures
  with ``EXP_TYPE=MIR_IMAGE.`` [#3633]

- - ``calwebb_tso3`` was changed to allow tso photometry processing of exposures
  with (``EXP_TYPE=MIR_IMAGE`` and tsovisit = True) or  with (``EXP_TYPE=MIR_IMAGE``) [#3650]

- Changed the default value of good_pixel from 4 to 6 for all outlier
  detection steps and both resample steps [#3638]

resample
--------

- Changed default value of good_pixel from 4 to 6 [#3638]

wfs_combine
-----------

- Allow handling of non-science members in input associations [#3947]


0.13.3 (2019-06-04)
===================

ami
---

- Fixed indentation bug in ami_analyze, so now all results are sufficiently
  close to the results of the stand-alone prototype. Other modifications include
  minor tweaks to more closely match those in the prototype code: changed some of
  initial values of the estimation parameters, and the filtering routine
  arguments.  [#3487]

- Updated ami_analyze.cfg to use default value of zero for rotation. [#3520]

- ``ami_analyze`` now emits a RuntimeError if the input is _calints or if a
  throughput reference file cannot be found.  [#3567]

- Remove change to filtering routine arguments of #3487.  [#3612]

assign_wcs
----------

- Fix a one pixel off problem with the NIRSpec NRS2 WCS transforms. [#3473]

- Raise a ``ValueError`` if the FWCPOS keyword isn't found in input NIRISS
  WFSS images. [#3574]

associations
------------

- Added the fxd_slit keyword as the third optical component [#3607]

- Orders the elements in Level3 naming in alphabetical order [#3614]

- Ensured that higher-order candidates only exist for Level2 associations [#3629]

- Improve member checking and removed duplicate product names [#3647]

combine_1d
----------

- Unit tests were added to combine_1d.  [#3490]

datamodels
----------

- Datamodels schemas should now be referenced with
  ``http://stsci.edu/schemas/jwst_datamodel/image.schema`` instead of
  ``http://jwst.stsci.edu/schemas/image.schema.yaml``.  The datamodels
  ``BaseExtension`` is renamed internally to ``DataModelExtension``. [#3437]

- Added the new column "relresperror" to the "MiriImgPhotomModel" data
  model schema. [#3512]

- Added all ``SlitModel`` data arrays to ``MultiExposureModel``, so that all input
  arrays appear in the output of ``exp_to_source``. [#3572]

extract_1d
----------

- An indexing bug was fixed. [#3497]

- Pixels with wavelength = NaN are no longer used. [#3539]

flatfield
---------

- Remove flatfield step parameter `flat_suffix`.  Add boolean step parameter
  `save_interpolated_flat`.  Refactor flatfield internals. [#3493]

- Propagate uncertainty from FFLAT, SFLAT and DFLAT flat fields into science
  ERR array and VAR_FLAT array for NIRSpec spectroscopic modes.  [#3538]

jump
----

- Add multiprocessing capability to JumpStep [#3440]

extract_2d
----------

- Replaced a white space in the names of grism objects with an underscore. [#3517]

- Update WFSS slit names to use simple integer value, and add accompanying unit
  test for NIRCAM grism extract_2d [#3632].

master_background
-----------------

- Fix bug in master_background where the flux from the input x1d files
  was being combined instead of the background columns.  [#3468]

- Use the surf_bright column instead of flux in master_background.  [#3476]

model_blender
-------------

- Allow blendmodels to ignore attributes in asdf tree not in schema [#3480]
- Add new rules for dates and times [#3554]

photom
------

- Updated to zero-out pixels outside the wavelength range of flux calibration
  and set DQ=DO_NOT_USE. [#3475, #3489]

pipeline
--------

- ``calwebb_spec3`` was changed to allow processing of WFSS modes. [#3517]

- ``calwebb_image2`` was changed to prevent 3D data from being sent to
  ``resample``. [#3544]

- ``calwebb_spec2`` was changed to check for an error in ``assign_wcs`` processing
  before executing the ``background`` step. [#3574]

refpix
------

- Fixed a bug where pixeldq arrays were being transformed from DMS to detector
  coordinates for every group instead of just once

skymatch
--------

- Improved reliability when matching sky in images with very close sky
  footprints. [#3557]

stpipe
------

- Capability to define reference overrides using a ``DataModel`` instead of
  a file path was added.  [#3514]

tweakreg
--------

- Mask and do not use NON-SCIENCE regions in tweakreg source detection. [#3461]


0.13.2 (2019-05-14)
===================

assign_wcs
----------

- The MIRI LRS WCS was updated to include an nverse transform. [#3106, #3360]

- The MIRI LRS spectral distortion is implemented now using a spline model. [#3106]

- Both ``dither_point_index`` and ``metadata_id`` are used now to match rows
  into the MSA meta file. [#3448]

- ``MissingMSAFileError`` was renamed to ``MSAFileError`` [#3448]

- Added two parameters to ``assign_wcs``, ``slit_y_low`` and ``slit_y_high``,
  to allow changing the lower and upper limit of a Nirspec slit in the instrument
  model. [#2819]

background
----------

- Verify the exposures to be used as background have the same NIRSpec GWA
  tilt values as the science exposures. If the background and science
  exposures do not have matching GWA tilt values, then skip the background
  subtraction step in calspec2. [#3252]

barshadow
---------

- Updated to apply the correction to the science data arrays, in addition
  to attaching as an extension. [#3319]

- Updated to apply the square of the correction to VAR_FLAT [#3427]

calwebb_spec3
-------------

- Add the ``master_background`` subtraction step to the pipeline. [#3296]

combine_1d
----------

- Fix call to wcs.invert, and don't weight flux by sensitivity if the net
  column is all zeros. [#3274]

- Modified to use the same columns as now written by extract_1d.
  The background parameter has been removed, since dividing by npixels
  is now done in extract_1d. [#3412]

datamodels
----------

- Fix ``url_mapper`` for fits-schema to allow URLs with of the format
  http://stsci.edu/schemas/fits-schema/ to map to the correct location
  in the ``jwst`` package. [#3239]

- Change ``ModelContainer`` to load and instantiate datamodels from an
  association on init.  This reverts #1027. [#3264]

- Keyword updates to data model schemas, including OBSFOLDR, MIRNGRPS,
  MIRNFRMS, and new PATTTYPE values. [#3266]

- Keyword updates to remove GS_STATE and change GUIDESTA to string
  type. [#3314]

- Added BUNIT keyword to gain and readnoise reference file schemas.
  [#3322]

- Update ``dq_def.schema``, ``group.schema`` and ``int_times.schema`` to comply
  with ASDF standard.  Remove unused ``extract1d.schema``.  [#3386]

- Update schemas to add new READPATT and BAND allowed values. [#3463]

extract_1d
----------

- This step can now use a reference image for IFU data.  The reference
  image (for IFU) may be either 2-D or 3-D.  When using a reference image
  for non-IFU data, background smoothing is now done after scaling the
  background count rate. [#3258]

- Unit tests were added for IFU data. [#3285]

- The target coordinates are used (for some modes) to determine the
  extraction location, i.e. correcting for nod/dither offset.  For IFU,
  the areas of the source aperture and background annulus are computed
  differently. [#3362

- For IFU data for an extended source, the extraction parameters are
  assigned values so that the entire image will be extracted, with no
  background subtraction.  For non-IFU data, a try/except block was added
  to check for a WCS that does not have an inverse.  Some code (but not
  all) for the now-obsolete RELSENS extension has been deleted. [#3390]

- This now writes columns SURF_BRIGHT and SB_ERROR instead of NET and
  NERROR.  The BACKGROUND column is divided by NPIXELS, so the units will
  be surface brightness.  This step no longer looks for a RELSENS
  extension. [#3412]

- The keywords that describe the units for the FLUX and ERROR columns
  have been corrected; the units are now specified as "Jy". [#3447]

extract_2d
----------

- An attribute ``dither_point`` was added to each slit in a ``MultiSlitModel``
  for MOS observations. [#3448]

flatfield
---------

- Propagate uncertainty from flat field into science ERR array and new
  VAR_FLAT array which holds the variance due to the flat field.  [#3384]

master_background
-----------------

- Modified the unit tests for ``expand_to_2d``. [#3242]

- Modified ``MasterBackgroundStep`` to be skipped if ``BackgroundStep``
  was already run on the data.  A new ``force_subtract`` parameter is
  added to override this logic.  [#3263]

- ``MasterBackgroundStep`` now can handle BACKGROUND association members
  that come from nodded exposures of the source. [#3311]

- Updated the DQFlags of the background subtracted data to be DO_NOT_USE
  for the pixels that have wavelengths outside the master background [#3326]

- Modified ``expand_to_2d`` to loop over pixels for WFSS data. [#3408]

outlier_detection
-----------------

- Fixed a bug that was causing the step to crash when calling the
  ``cube_build`` step for MIRI MRS data. [#3296]

pathloss
--------

- Updated to apply the correction to the science data and err arrays. [#3323]

- Updated to apply the square of the correction to VAR_FLAT [#3427]

photom
------

- Updated to apply the flux calibration to the science data and err arrays.
  [#3359]

- Updated to compute a wavelength array for NIRISS SOSS exposures using
  spectral order 1. [#3387]

- Updated to apply the square of the correction to VAR_FLAT [#3427]

reffile_utils
-------------

- Improved error messages when problems are encountered in extracting
  subarrays from reference files. [#3268]

resample_spec
-------------

- Fixed an issue with the spatial component of the WCS where the inverse
  transform gave different results for negative ``RA`` and ``360 + RA``. [#3404]


set_telescope_pointing
----------------------

- Fix ``populate_model_from_siaf`` to convert SIAF pixel scale from
  arcsec to degrees for CDELTn keywords. [#3248]

- Updates to prevent crashes when SIAF values needed for crpix or
  cdelt keywords are missing. [#3316]

- Convert FSM correction values from arcsec to radians. [#3367]

srctype
-------

- Updated logic for background targets and nodded exposures. [#3310]


transforms
----------

- A field ``dither_point`` was added to the ``Slit`` structure. [#3448]


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

- Fixed NIRISS WFSS catalog naming and implement NIRCam WFSS [#3515]

- Fixed treating non-science as TSO [#3601]

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

- Added dq flagging [#3804]

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

- Updated SourceContainer to wrap each exposure of a MultiExposure in a
  SlitModel, allowing pipeline code to simply treat each as DataModel.
  [#3438]

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

 - Updated twopoint_difference.py to not use groups with groupdq set to DO_NOT_USE [#3495]

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

- The `LRSWavelength` model was removed as obsolete.
  Instead a spline is used for the wavelength solution. [#3106]

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
  of ``GRATING``. If ``GRATING=MIRROR`` imaging mode is chosen regardless of ``EXP_TYPE``. [#2761]

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

wiimatch
--------

0.11.0 (2018-09-10)
===================

The 0.11.0 release is highlighted by the inclusion of steps for resampling
spectral images and time series grism observations.   In addition, this
release had 39 issues closed and a number of pull requests to improve PEP8
compliance, improve performance, and enhance the testing.  The release also
included updated documentation for accessing CRDS when running the JWST
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

- Fixed bug in the ordering of cube footprint [#2371]

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
  NIRSPEC case there was a 40% improvement in the speed of creating IFUCubes.

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

- Updated fits_generator to ignore files beginning with '.' [#2333]

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

- Update IRS2 data model and add regression tests [#2295]


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

- An example has been added to the model_blender documentation for how to blend meta information [#2206]

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


white_light
-----------

wiimatch
--------
