Description
============

:Class: `jwst.srctype.SourceTypeStep`
:Alias: srctype

The Source Type (``srctype``) step in the calibration pipeline attempts to
determine whether a spectroscopic source should be considered to be a point
or extended object, populating the "SRCTYPE" keyword with a value of either
"POINT" or "EXTENDED."
This information is then used in some subsequent spectroscopic processing
steps to apply source-dependent corrections.

Single Source Observations
--------------------------
For JWST observing modes that use a single primary target (e.g. MIRI MRS
and LRS spectroscopy and NIRSpec IFU and Fixed-Slit spectroscopy), the observer
has the option of designating a source type in the APT template for the
observation. They have the choice of declaring whether or not the source
should be considered extended. If they don't know the character of the source,
they can also choose a value of "UNKNOWN." The observer's choice is passed along
to DMS processing, which sets the value of the "SRCTYAPT" keyword in the
primary header of the products used as input to the calibration pipeline.
If the user has selected a value in the APT, the "SRCTYAPT" keyword will be set
to "POINT", "EXTENDED", or "UNKNOWN." If the selection is not available for a
given observing mode or a choice wasn't made, the "SRCTYAPT" keyword will not
appear in the uncalibrated product header.

The ``srctype`` step sets a value for the "SRCTYPE" keyword that is stored in
the "SCI" extension header(s) of data products. The step sets the value of
"SRCTYPE" based on input from the user given in the "SRCTYAPT" keyword, as
well as other rules that can override the "SRCTYAPT" values.

The ``srctype`` step first checks to see if the "SRCTYAPT" keyword
is present and has already been populated. If "SRCTYAPT" is not present or
is set to "UNKNOWN", the step determines a suitable value based on the
observing mode, command line input, and other characteristics of the
exposure. The following choices are used, in order of priority:

 - The source type can be specified by the user on the command line.
   Exposure types for which this is permitted contain a single pre-defined
   target, i.e. MIR_LRS-FIXEDSLIT, MIR_LRS-SLITLESS, MIR_MRS,NRC_TSGRISM,
   NRS_FIXEDSLIT, NRS_BRIGHTOBJ, and NRS_IFU. Other EXP_TYPEs will be
   ignored.  For NRS_FIXEDSLIT exposures, a user-supplied value can replace
   the value for the target in the primary slit only, while the other slits
   will retain their default settings of "EXTENDED" (which is appropriate
   for sky background).

 - Background target exposures default to a source type of "EXTENDED."
   Background exposures are identified by the keyword "BKGDTARG" set
   to True.

 - TSO exposures default to a source type of "POINT." TSO exposures are
   identified by EXP_TYPE="NRC_TSGRISM" or "NRS_BRIGHTOBJ", or
   TSOVISIT=True.

 - Exposures that are part of a nodded dither pattern, which are assumed
   to only be used with point-like targets, default to a source type
   of "POINT." Nodded exposures are usually identified by the "PATTTYPE"
   keyword either being set to a value of "POINT-SOURCE" or containing the
   sub-string "NOD" (NIRSpec IFU and Fixed Slit). For MIRI MRS exposures
   the keyword "DITHOPFR" (DITHer pattern OPtimized FoR) is used instead of
   "PATTTYPE". If it has a value of "POINT-SOURCE", the source type is set
   to "POINT".

 - If none of the above conditions apply, and the user did not choose a
   value in the APT, the following table of defaults is used, based on
   the "EXP_TYPE" keyword value:

.. _srctype_table:

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
| NRS_IFU           | NIRSpec IFU            | EXTENDED |
+-------------------+------------------------+----------+

If the EXP_TYPE value of the input image is not in the above list,
SRCTYPE will be set to "UNKNOWN".

NOTE: NIRSpec fixed-slit (EXP_TYPE="NRS_FIXEDSLIT") exposures are
unique in that a single target is specified in the APT, yet data for
multiple slits can be contained within an exposure, depending on the
size of the readout used (e.g. SUBARRAY="ALLSLITS"). For this observing
mode, the source type selection resulting from the logic outlined above
is used to populate the SRCTYPE keyword associated with the data for
the primary slit instance in the pipeline data products. The primary slit
is determined from the value of the "FXD_SLIT" keyword. Any additional
slit instances contained within the data product will have their
SRCTYPE value set to "EXTENDED", as non-primary slits are expected to contain
background.

Multi-Source Observations
-------------------------

NIRSpec MOS
+++++++++++

For NIRSpec MOS exposures (EXP_TYPE="NRS_MSASPEC"), there are multiple sources
per exposure and hence a single user-selected parameter can't be used in the
APT, nor a single keyword in the science product, to record the type of each
source. For these exposures, a stellarity value can be supplied by the observer
for each source used in the MSA Planning Tool (MPT). The stellarity values are
in turn passed from the MPT to the MSA metadata (_msa.fits) file created by DMS
and used in the calibration pipeline. The stellarity values from the MSA
metadata file are loaded for each source/slitlet by the ``assign_wcs`` step of
the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline and then evaluated by the
``srctype`` step to determine whether each source should be treated as point or
extended.

If the stellarity value for a given source in the MSA metadata is less
than zero, the source type defaults to "POINT." If the stellarity value is
between zero and 0.75, it is set to "EXTENDED", and if the stellarity value
is greater than 0.75, it is set to "POINT." The resulting choice is stored in
the "SRCTYPE" keyword located in the header of the SCI extension associated with
each slitlet.

In the future, reference files will be used
to set more detailed threshold values for stellarity, based on the
particular filters, gratings, etc. of each exposure.

NIRCam and NIRISS WFSS
++++++++++++++++++++++
It is not possible to specify ahead of time the source types for spectra that
may show up in a Wide-Field Slitless Spectroscopy exposure. So for these modes
the ``srctype`` step simply sets the SRCTYPE keyword value to "UNKNOWN" and the
actual source sizes are derived from the catalog information generated
from direct images that are obtained as part of a WFSS observation.
