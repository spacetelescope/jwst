.. _resample_B6_design_:

**************************************
Python Step Design: ResampleStep
**************************************

.. moduleauthor:: Warren Hack, Howard Bushouse <help@stsci.edu>

The full capabilities for resampling JWST data will require this step to do a lot more than simply apply a distortion model to an image, as implemented in Build 4.  This step will need to perform these operations in order to optimally  align and resample all JWST data.

* Interpret inputs (ASN table) to identify all input observations to be combined/resampled
* Read in WCS information for each input observation

  - Use input WCSs to build output WCS

* Perform background matching using skymatch to create background levelled/subtracted images
* For images:

  - Perform initial source identification/location using astropy.
  - Classify sources and select highest quality sources for alignment (point sources and/or high S/N compact extended sources)

* For spectral data:

  - Locate center of spectral orders or locate emission lines and treat as source(s)

* Cross-match sources and perform fit (use 'tweakreg' as model for these operations for auto-mosaic building)
* For each input:

  - Create transformation array from input pixels to output pixels using WCS transformations
  - Apply transformation array to resample input onto output using drizzle algorithm
  - Perform 'OUTLIER DETECTION' (bad-pixel/cosmic-ray ID as per astrodrizzle) and update input observation's DQ arrays with flags for identified pixels

* For images:

  - Perform source identification on final combined resampled output
  - Determine photometric quantities for each source
  - Match sources in resampled output with positions from each input observation
  - Create catalog with source positions in output frame, each input frame, and photometric quantities
