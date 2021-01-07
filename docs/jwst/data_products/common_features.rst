Common Features
---------------

All JWST FITS data products have a few common features in their structure
and organization:

1. The FITS primary Header Data Unit (HDU) only contains header information,
   in the form of keyword records, with an empty data array, which is
   indicated by the occurence of NAXIS=0 in the primary header. Meta
   data that pertains to the entire product is stored in keywords in the
   primary header. Meta data related to specific extensions (see below)
   is stored in keywords in the headers of each extension.

2. All data related to the product are contained in one or more FITS
   IMAGE or BINTABLE extensions. The header of each extension may contain
   keywords that pertain uniquely to that extension.

