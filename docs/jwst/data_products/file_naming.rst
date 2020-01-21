.. _file_naming_schemes:

File Naming Schemes
-------------------

.. _exp_file_names:

Exposure file names
^^^^^^^^^^^^^^^^^^^
The names of the exposure level data (stage 0 to 2) are constructed with information from the science header of the exposure, allowing users to map it to the observation in their corresponding APT files. The FITS file naming scheme for Stage 0, 1, and 2 "exposure-based" products is:

 jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>(-<"seg"NNN>)_<detector>_<prodType>.fits

where

 - ppppp: program ID number
 - ooo: observation number
 - vvv: visit number
 - gg: visit group
 - s: parallel sequence ID (1=prime, 2-5=parallel)
 - aa: activity number (base 36)
 - eeeee: exposure number
 - segNNN: the text "seg" followed by a three-digit segment number (optional)
 - detector: detector name (e.g. 'nrca1', 'nrcblong', 'mirimage')
 - prodType: product type identifier (e.g. 'uncal', 'rate', 'cal')

An example Stage 1 product FITS file name is::

 jw93065002001_02101_00001_nrca1_rate.fits

.. _src_file_names:

Stage 3 file names
^^^^^^^^^^^^^^^^^^
In this stage, the calibration pipeline uses the association information to identify the relationship between exposures 
that are to be combined by design to form a single product. These data products result from the combination of multiple 
exposures like dithers or mosaics.

The format for the file name of a Stage 3 association product has all alphabetic characters in lower case, underscores 
are only used to delineate between major fields, and dashes can be used as separators for optional fields. 
Just as for Stage 2, the suffix distinguishes the different file products of Stage 3 of the calibration pipeline.

The FITS file naming scheme for Stage 3 "source-based" products is:

 jw<ppppp>-<AC_ID>_[<"t"TargID | "s"SourceID>](-<"epoch"X>)_<instr>_<optElements>(-<subarray>)_<prodType>(-<ACT_ID>).fits

where

 - ppppp: program ID number
 - AC_ID: association candidate ID
 - TargID: a 3-digit target ID (either TargID or SourceID must be present)
 - SourceID: a 5-digit source ID
 - epochX: the text "epoch" followed by a single digit epoch number (optional)
 - instr: science instrument name (e.g. 'nircam', 'miri')
 - optElements: a single or hyphen-separated list of optical elements (e.g. filter, grating)
 - subarray: optional indicator of subarray name
 - prodType: product type identifier (e.g. 'i2d', 's3d', 'x1d')
 - ACT_ID: a 2-digit activity ID


An example Stage 3 product FITS file name is:

 jw87600-a3001_t001_niriss_f480m-nrm_amiavg.fits

In this case, the EPOCH is an optional parameter that is used for observations of the same target that are taken at different epochs.

.. _segmented_files:

Segmented Products
^^^^^^^^^^^^^^^^^^
When the raw data volume for an individual exposure is determined to be large enough to result in
Stage 2 products greater than 2 GB in size, all Stage 0-2 products for the exposure are broken into
multiple segments, so as to keep total file sizes to a reasonable level. This is often the case with
Time Series Observations (TSO), where individual exposures can have thousands of integrations.
The detector data are broken along integration boundaries (never within an integration) and stored
in "segmented" products. The segments are identified by the "segNNN" field in exposure-based file
names, where NNN is 1-indexed and always includes any leading zeros.

Segmented products contain extra keywords located in their primary headers that help to identify
the segments and the contents of each segment. The following segment-related keywords are used:

 - EXSEGNUM: The segment number of the current product
 - EXSEGTOT: The total number of segments
 - INTSTART: The starting integration number of the data in this segment
 - INTEND: The ending integration number of the data in this segment

All of the Stage 1 and Stage 2 calibration pipelines will process each segment independently,
creating the full set of intermediate and calibrated products for each segment. The calibrated data
for all segments is then combined by one of the Stage 3 pipelines into a source-based Stage 3
product.

