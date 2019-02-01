File Naming Scheme
------------------
The FITS file naming scheme for Stage 0, 1, and 2 "exposure-based" products is::

 jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>_<detector>_<prodType>.fits

where

 - ppppp: program ID number
 - ooo: observation number
 - vvv: visit number
 - gg: visit group
 - s: parallel sequence ID (1=prime, 2-5=parallel)
 - aa: activity number (base 36)
 - eeeee: exposure number
 - detector: detector name (e.g. 'nrca1', 'nrcblong', 'mirimage')
 - prodType: product type identifier (e.g. 'uncal', 'rate', 'cal')

An example Stage 1 product FITS file name is::

 jw93065002001_02101_00001_nrca1_rate.fits

The FITS file naming scheme for Stage 3 "source-based" products is::

 jw<ppppp>-<AC_ID>_[<"t"TargID | "s"SourceID>](-<"epoch"X>)_<instr>_<optElements>(-<subarray>)_<prodType>(-<ACT_ID>).fits

where

 - ppppp: program ID number
 - AC_ID: association candidate ID
 - TargID: a 3-digit target ID (either TargID or SourceID must be present)
 - SourceID: a 5-digit source ID
 - epochX: the text "epoch" followed by a single digit epoch number
 - instr: science instrument name (e.g. 'nircam', 'miri')
 - optElements: a single or hyphen-separated list of optical elements (e.g. filter, grating)
 - subarray: optional indicator of subarray name
 - prodType: product type identifier (e.g. 'i2d', 's3d', 'x1d')
 - ACT_ID: a 2-digit activity ID

An example Stage 3 product FITS file name is::

 jw87600-a3001_t001_niriss_f480m-nrm_amiavg.fits

