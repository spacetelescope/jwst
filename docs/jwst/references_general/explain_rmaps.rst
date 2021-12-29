CRDS Integration in CAL Code
============================

For JWST, the Calibration Reference Data System (CRDS) is directly integrated
with calibration steps and pipelines resulting in conventions for how Steps and
Pipelines should be written.

Step Attribute .reference_file_types
------------------------------------

Each calibration Step is required to define an attribute or property named
*reference_file_types* which defines the CRDS reference file types that are
required for running the calibration step.  Note that for some Steps the
reference file types actually used vary so the minimal list of required types
may not be known if no science data is defined.

Note that the ``Step`` parameter reference files do not need to be specified.
These are automatically requested in the ``Step`` architecture.

CRDS Prefetch
-------------

To ensure all reference files required by a pipeline are available prior to
processing, Pipelines perform a "pre-fetch" of all references required by any
Step in the Pipeline.  This generic Pipeline behavior is intended to prevent
processing which fails due to missing CRDS files after running for lengthy
periods.

When *CRDS_SERVER_URL* and *CRDS_PATH* are properly configured, CRDS will
download and locally cache the minimal set of references required to calibrate
specific data using a particular pipeline.  This configuration supports remote
processing for users with no access or inefficient access to the STScI local
network.

The pre-fetch also enables CRDS to report on all reference file assignment and
availability problems a pipeline will encounter in a single CAL run.  This is
required in I&T scenarios where the total number of pipeline runs is very
limited (often weekly) so solving as many reference file issues per run as
possible is needed. 

While the prefetch runs for onsite users,  since the default CRDS configuration
points to a complete CRDS cache,  no downloads will occur.

Step Method .get_reference_file()
---------------------------------

During processing individual Steps make secondary calls to CRDS via the
*get_reference_file(input_file, reference_file_type)* method to fetch the cache
paths of individual reference files.  If no file is applicable, CRDS returns
the value 'N/A' which is sometimes used to skip related Step processing.  If
there is an error determining a reference file, CRDS will raise an exception
stopping the calibration; typically these occur due to missing reference files
or incorrectly specified dataset parameters.

While get_reference_file() returns the absolute path of CRDS reference files,
reference file assignments are recorded in output products using a *crds://*
URI prefix which translates to roughly "the path of this file in the local
cache you've defined using CRDS_PATH,  or under /grp/crds/cache if you didn't
define CRDS_PATH".

Best Reference Matching
=======================

The Calibration Reference Data System (CRDS) assigns the best reference files
needed to process a data set based on the dataset's metadata (FITS headers) and
plain text CRDS rules files.

CRDS rules and references are organized into a 4 tiered hierarchical network of
versioned files consisting of:

* .pmap  - The overall context for the pipeline (i.e. all instruments)
* .imap  - The rules for all reference types of one instrument
* .rmap  - The rules for all reference files of one type of one instrument
* .fits,.asdf,.json - Individual reference files assigned by .rmaps

Based on dataset parameters, CRDS traverses the hierarchy of rules files,
generally starting from the .pmap and descending until a particular reference
file is assigned.

Visiting the JWST operational website here:

https://jwst-crds.stsci.edu/

and opening up instrument panes of the *Operational References* display can
rapidly give an idea about how reference files should be assigned.

CRDS Parameter Naming
---------------------

For the sake of brevity,  the CRDS website often refers to matching parameters
using truncated names intended to give the gist of a parameter.

Within CRDS rules for JWST, CRDS refers to parameters using jwst datamodels
attribute paths converted to capitalized strings analogous to FITS keywords.
For instance the datamodels attribute::

   meta.instrument.name

corresponds to CRDS rules parameter name::

   'META.INSTRUMENT.NAME'

and FITS keyword::

  'INSTRUME'

Using e.g. 'META.INSTRUMENT.NAME' permits consistent naming regardless of the
underlying file format (.fits vs. .asdf vs. .json).

When creating or accessing reference files, Python code uses the lower case
object path to populate an attribute corresponding to the upper case string.

Example .pmap contents
----------------------

Generally CRDS reference file lookups begin with a .pmap (context) file.

The .pmap's serial number describes the overall version of rules for a pipeline.

The contents of context jwst_0493.pmap are shown below::

  header = {
     'mapping' : 'PIPELINE',
     'observatory' : 'JWST',
     'name' : 'jwst_0493.pmap',
     'parkey' : ('META.INSTRUMENT.NAME',),
     ...
 }

 selector = {
     'FGS' : 'jwst_fgs_0073.imap',
     'MIRI' : 'jwst_miri_0158.imap',
     'NIRCAM' : 'jwst_nircam_0112.imap',
     'NIRISS' : 'jwst_niriss_0117.imap',
     'NIRSPEC' : 'jwst_nirspec_0173.imap',
     'SYSTEM' : 'jwst_system_0017.imap',
 }

Based on the parameter META.INSTRUMENT.NAME (INSTRUME) CRDS selects an
appropriate .imap for further searching.

In all CRDS rules files, the header's **parkey** field defines the parameter
names used to select a file.  These parkey names correspond to the values shown
in the selector's keys.

Conceptually all CRDS selectors consist of dictionaries which map parameter
values to either a file or a sub-selector.

If META.INSTRUMENT.NAME=NIRSPEC, then CRDS would choose *jwst_nirspec_0173.imap*
to continue it's search.

Example .imap contents
----------------------

A .imap file defines the appropriate version of .rmap to search for each
reference type supported by the corresponding instrument.   Below is an
example .imap taken from NIRSPEC::

  header = {
    'mapping' : 'INSTRUMENT',
    'instrument' : 'NIRSPEC',
    'name' : 'jwst_nirspec_0173.imap',
    'parkey' : ('REFTYPE',),
    ...
  }

  selector = {
    'AREA' : 'jwst_nirspec_area_0010.rmap',
    'BARSHADOW' : 'jwst_nirspec_barshadow_0002.rmap',
    'CAMERA' : 'jwst_nirspec_camera_0015.rmap',
    ...,
    'PATHLOSS' : 'jwst_nirspec_pathloss_0003.rmap',
    ...,
    'WAVECORR' : 'jwst_nirspec_wavecorr_0003.rmap',
    'WAVELENGTHRANGE' : 'jwst_nirspec_wavelengthrange_0015.rmap',
    'WCSREGIONS' : 'N/A',
    'WFSSBKG' : 'N/A',
  }

A value of N/A indicates that a particular reference type is not yet used by
this instrument and CRDS will return 'N/A' instead of a filename.

If the requested REFTYPE was PATHLOSS, CRDS would continue it's search with
*jwst_nirspec_pathloss_0003.rmap*.

Example .rmap contents
----------------------

Slightly modified contents of *jwst_nirspec_pathloss_0003.rmap* are shown
below::

 header = {
    'mapping' : 'REFERENCE',
    'observatory' : 'JWST',
    'instrument' : 'NIRSPEC',
    'filekind' : 'PATHLOSS',
    'name' : 'jwst_nirspec_pathloss_0003.rmap',
    'classes' : ('Match', 'UseAfter'),
    'parkey' : (('META.EXPOSURE.TYPE',), ('META.OBSERVATION.DATE', 'META.OBSERVATION.TIME')),
    ...
 }

 selector = Match({
    'NRS_AUTOWAVE' : 'N/A',
    'NRS_FIXEDSLIT|NRS_BRIGHTOBJ' : UseAfter({
        '1900-01-01 00:00:00' : 'jwst_nirspec_pathloss_0001.fits',
    }),
    'NRS_IFU' : UseAfter({
        '1900-01-01 00:00:00' : 'jwst_nirspec_pathloss_0003.fits',
    }),
    'NRS_MSASPEC' : UseAfter({
        '1900-01-01 00:00:00' : 'jwst_nirspec_pathloss_0002.fits',
        '2000-01-01 00:00:00' : 'jwst_nirspec_pathloss_0007.fits',
    }),
 })

Each class of CRDS rmap selector defines a search algorithm to be used at that
stage of the reference file lookup. 

Match Selector
++++++++++++++
 
In the example shown above, CRDS selects a nested UseAfter selector based on
the value of META.EXPOSURE.TYPE (EXP_TYPE).   The nested UseAfter is then
used for a secondary lookup to determine the assigned reference.

Parameters which contain or-bars, e.g.::
  
  'NRS_FIXEDSLIT|NRS_BRIGHTOBJ'

specify groups of values for which a file is equally applicable.

In this case the file *jwst_nirspec_pathloss_0001.fits* can be used to
calibrate either NRS_FIXEDSLIT or NRS_BRIGHTOBJ.

``Or'ed`` parameter combinations shown in rmaps are almost identical to the or'ed
parameter combinations taken from ``P\_`` pattern keywords; the only difference is
that rmaps do not specify the trailing or-bar required for ``P\_`` keyword values.

If a parameter combination maps to the value N/A,  then the reference type is
not applicable for that combination and CRDS returns the value N/A instead of
a filename.

UseAfter Selector
+++++++++++++++++

The UseAfter sub-selector applies a given reference file only to datasets which
occur at or after the specified date.  For cases where multiple references
occur prior to a dataset, CRDS chooses the most recent reference file as best.

Based on the dataset's values of::

   META.OBSERVATION.DATE (DATE-OBS) 
   META.OBSERVATION.TIME (TIME-OBS)

CRDS will choose the appropriate reference file by comparing them to the
date+time shown in the .rmap.  Conceptually, the date+time shown corresponds to
the value of::

   META.REFERENCE.USEAFTER (USEAFTER)

from each reference file with the USEAFTER's T replaced with a space.

* In the example above, if the dataset defines::

    EXP_TYPE=NRS_MSASPEC
    DATE-OBS=1999-01-01
    TIME-OBS=00:00:00

then CRDS will select *jwst_nirspec_pathloss_0002.fits* as best.

* In the example above, if the dataset defines::

    EXP_TYPE=NRS_MSASPEC
    DATE-OBS=2001-01-01
    TIME-OBS=00:00:00

then CRDS will select *jwst_nirspec_pathloss_0007.fits* as best.

* If the dataset defines e.g.::

    DATE-OBS=1864-01-01

then no reference match exists because the observation date precedes the
USEAFTER of all available reference files.

UseAfter selection is one of the rare cases where CRDS makes an
apples-to-oranges match and the dataset and reference file parameters being
correlated are not identical.  In fact,  not even the count of parameters
(DATE-OBS, TIME-OBS) vs. USEAFTER is identical.

Defining Reference File Applicability
-------------------------------------

Almost all reference files supply metadata which defines how CRDS should add
the file to its corresponding .rmap, i.e. each reference defines the science
data parameters for which it is *initially* applicable.

When creating reference files,  you will need to define a value for every
CRDS matching parameter and/or define a pattern using the ``P_`` version of the
matching parameter.

When CRDS adds a reference file to a .rmap, it uses literal matching between
the value defined in the reference file and the existing values shown in the
.rmap.  This enables CRDS to:

* add files to existing categories
* replace files in existing categories
* create new categories of files.

Because creating new categories is an unusual event which should be carefully
reviewed,  CRDS issues a warning when a reference file defines a new category.

Changing .rmaps to Reassign Reference Files
-------------------------------------------

While reference files generally specify their intended use, sometimes different
desired uses not specified in the reference file appear over time.  In CRDS it
is possible to alter only a .rmap to change the category or dates for which a
reference file applies.

This is a fundamental CRDS feature which enables changes to reference
assignment without forcing the re-delivery of an otherwise serviceable
reference file.  This feature is very commonly used, and the net consequence is
that **.rmap categories and dates do not have to match the contents of
reference files.**

It is better to view CRDS matching as a comparison between dataset parameters
and a .rmap.   Although references do state "initial intent",  reference file
metadata should not be viewed as definitive for how a file is assigned.

More Complex Matching
---------------------

CRDS matching supports more complex situations than shown in the example above.

Although reference files are generally constructed so that their metadata
defines the instrument modes for which they're applicable, conceptually, the
values shown in .rmaps correspond to values in the dataset.  Indeed, it is
possible to change the values shown in the rmap so that they differ from their
corresponding values in the reference file.  This makes it possible to reassign
reference files rather than redelivering them.

Match Parameter Combinations
++++++++++++++++++++++++++++

For matches using combinations of multiple parameters, the Match selector keys
will be shown as tuples, e.g.::

  ('NRS1|NRS2', 'ANY', 'GENERIC', '1', '1', '2048', '2048')

Because this match category matches either DETECTOR=NRS1 or NRS2, this single
rmap entry represents two discrete parameter combinations.  With multiple
pattern values (not shown here), a single match category can match many
different discrete combinations.

The *parkey* tuple from the NIRSPEC SUPERBIAS rmap which supplied the
example match case above looks like::

   (('META.INSTRUMENT.DETECTOR', 'META.EXPOSURE.READPATT',
   'META.SUBARRAY.NAME', 'META.SUBARRAY.XSTART', 'META.SUBARRAY.YSTART',
   'META.SUBARRAY.XSIZE', 'META.SUBARRAY.YSIZE'),
   ('META.OBSERVATION.DATE', 'META.OBSERVATION.TIME'))

The first sub-tuple corresponds to the Match cases,  and the second sub-tuple
corresponds to the nested UseAfters.

Weighted Matching
+++++++++++++++++

It's possible for CRDS to complete it's search without finding a unique match.
To help resolve these situations, the Match algorithm uses a weighting scheme.

Each parameter with an exact match contributes a value of 1 to the weighted
sum.   e.g. 'NRS1' matches 'NRS1|NRS2' exactly once patterns are accounted for.

An rmap value of ANY will match any dataset value and also has a weight of 1.

An rmap value of N/A or GENERIC will match any dataset value but have a weight
of 0, contributing nothing to the strength of the match.

Conceptually, the match with the highest weighting value is used.  It is
possible to create rmaps where ambiguity is not resolved by the weighting
scheme but it works fairly well when used sparingly and isolated to as few
parameters as possible.

Typically the value GENERIC corresponds to a full frame reference file which
can support the calibration of any SUBARRAY by performing a cut-out.

More Information
----------------

More information about CRDS can be found in the CRDS User's Guide maintained
on the CRDS server here:

https://jwst-crds.stsci.edu/static/users_guide/index.html

