Description
============
The Source Type (``srctype``) step in the calibration pipeline attempts to
determine whether a spectroscopic source should be considered to be a point
or extended object, populating the "SRCTYPE" keyword with a value of either
"POINT" or "EXTENDED."
This information is then used in some subsequent spectroscopic processing
steps.

Depending on the JWST observing mode, the observer may have the option of
designating a source type in the APT template for the observation. They have
the choice of declaring whether or not the source should be considered
extended. If they don't know the character of the source, they can also
choose a value of "UNKNOWN." The observer's choice in the APT is passed along
to DMS processing, which sets the value of the "SRCTYPE" keyword in the
primary header of the uncalibrated (_uncal.fits) product that's used as input
to the calibration pipeline. The resulting "SRCTYPE" values can be
"POINT", "EXTENDED", or "UNKNOWN."

The ``srctype`` calibration step checks to see if the "SRCTYPE" keyword
has already been populated. If the observer did not provide a source
type value or the "SRCTYPE" keyword is set to "UNKNOWN", the ``srctype``
step will choose a suitable value based on the observing mode and
other characteristics of the exposure. The following choices are used, in
order of priority:

 - Background target exposures default to a source type of "EXTENDED."
   Background exposures are identified by the keyword "BKGDTARG" set
   to True.

 - Exposures that are part of a nodded dither pattern, which are assumed
   to only be used with point-like targets, default to a source type
   of "POINT." Nodded exposures are identified by the "PATTTYPE" keyword
   either being set to a value of "POINT-SOURCE" (MIRI MRS) or containing
   the sub-string "NOD" (NIRSpec IFU and Fixed Slit).

 - If none of the above conditions apply, and the user did not choose a
   value in the APT, the following table of defaults is used, based on
   the "EXP_TYPE" keyword value:

+-------------------+------------------------+----------+
| EXP_TYPE          | Exposure Type          | SRCTYPE  |
+===================+========================+==========+
| MIR_LRS-FIXEDSLIT | MIRI LRS fixed-slit    | POINT    |
+-------------------+------------------------+----------+
| MIR_LRS-SLITLESS  | MIRI LRS slitless      | POINT    |
+-------------------+------------------------+----------+
| MIR_MRS           | MIRI MRS (IFU)         | EXTENDED |
+-------------------+------------------------+----------+
| NIS_SOSS          | NIRISS SOSS            | POINT    |
+-------------------+------------------------+----------+
| NRS_FIXEDSLIT     | NIRSpec fixed-slit     | POINT    |
+-------------------+------------------------+----------+
| NRS_BRIGHTOBJ     | NIRSpec bright object  | POINT    |
+-------------------+------------------------+----------+
| NRS_IFU           | NIRSpec IFU            | POINT    |
+-------------------+------------------------+----------+

For NIRSpec MOS exposures (EXP_TYPE="NRS_MSASPEC"), there are multiple
sources per exposure and hence a single parameter can't be used in the
APT, nor a single keyword in the science product, to record the type of
each source. For these exposures, a stellarity value can be supplied by
the observer for each source used in the MSA Planning Tool (MPT). The
stellarity values are
in turn passed from the MPT to the MSA metadata (_msa.fits) file
created by DMS and used in the calibration pipeline. The stellarity
values from the MSA metadata file are loaded for each source/slitlet
by the ``assign_wcs`` step of the ``calwebb_spec2`` pipeline and then
evaluated by the ``srctype`` step to determine whether each source
should be treated as point or extended.

If the stellarity value is less than zero, the source type is set to
"UNKNOWN." If the stellarity value is between zero and 0.75, it is
set to "EXTENDED", and if the stellarity value is greater than 0.75,
it is set to "POINT." The resulting choice as stored in a "SRCTYPE"
keyword located in the header of the SCI extension associated with
each source/slitlet.

In the future, reference files will be used
to set more detailed threshold values for stellarity, based on the
particular filters, gratings, etc. of each exposure.
