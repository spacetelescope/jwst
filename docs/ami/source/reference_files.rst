Reference File
==============
The ami_analyze step uses a THROUGHPUT reference file, which contains
throughput data for the filter used in the input AMI image.

CRDS Selection Criteria
-----------------------
Throughput reference files are selected on the basis of INSTRUME and 
FILTER values for the input science data set.

Throughput Reference File Format
--------------------------------
Throughput reference files are FITS files with one BINTABLE
extension. The FITS primary data array is assumed to be empty. The 
table extension uses ``EXTNAME=THROUGHPUT`` and the data table has the
following characteristics:

===========  =========  ==========
Column name  Data type  Units
===========  =========  ==========
wavelength   float      Angstroms
throughput   float      (unitless)
===========  =========  ==========

The ami_average and ami_normalize steps do not use any reference files.
