Introduction
============

This document provides instructions on running the JWST Science Calibration
Pipeline (referred to as "the pipeline") and individual pipeline steps.

Pipeline modules are available for detector-level (stage 1) processing of
data from all observing modes, stage 2 processing for imaging and
spectroscopic modes, and stage 3 processing for imaging, spectroscopic,
coronagraphic, Aperture Masking Interferometry (AMI), and Time Series
Observations (TSO).

Stage 1 processing consists of detector-level
corrections that must be performed on a group-by-group basis
before ramp fitting is applied. The output of stage 1 processing
is a countrate image per exposure or per integration for some modes.
Details of this pipeline can be found at :ref:`stage1-flow`.

Stage 2 processing consists of additional corrections and
calibrations to produce fully calibrated exposures. The details
differ for imaging and spectroscopic exposures, and there are some
corrections that are unique to certain instruments or modes.
Details are at :ref:`stage2-imaging-flow`
and :ref:`stage2-spectroscopic-flow`.

Stage 3 processing consists of routines that work with multiple exposures
and in most cases produce some kind of combined product.
There are dedicated (and unique) pipeline modules for stage 3 processing of
imaging, spectroscopic, coronagraphic, AMI, and TSO observations. Details
of each are available at
:ref:`stage3-imaging-flow`,
:ref:`stage3-spectroscopic-flow`,
:ref:`stage3-coron-flow`,
:ref:`stage3-ami-flow`, and
:ref:`stage3-tso-flow`.

The remainder of this document discusses pipeline configuration files and
gives examples of running pipelines as a whole or in individual steps.

Reference Files
===============

Many pipeline steps rely on the use of a set of reference files essential to ensure the correct and accurate process of the data. The reference files are instrument-specific, and are periodically updated as the data process evolves and the understanding of the instruments improves. They are created, tested and validated by the JWST Instrument Teams. They ensure all the files are in the correct format and have all required header keywords. The files are then delivered to the Reference Data for Calibration and Tools (ReDCaT) Management Team. The result of this process is the files being ingested into CRDS (the JWST Calibration Reference Data System), and made available to the pipeline team and any other ground-subsystem that needs access to them.

Information about all the reference files used by the Calibration Pipeline can be found at
:ref:`reference-file-formats-documentation`
as well as in the documentation for the Calibration Step using them.
 
CRDS
====

CRDS reference file mappings are usually set by default to always give access
to the most recent reference file deliveries and selection rules. On
occasion it might be necessary or desirable to use one of the non-default
mappings in order to, for example, run different versions of the pipeline
software or use older versions of the reference files. This can be
accomplished by setting the environment variable `CRDS_CONTEXT` to the
desired project mapping version, e.g.
::

$ export CRDS_CONTEXT='jwst_0421.pmap'

The current storage location for all JWST CRDS reference files is:
::

/grp/crds/jwst/references/jwst/

Each pipeline step records the reference file that it used in the value of
a header keyword in the output data file. The keyword names use the syntax
"R_<ref>", where <ref> corresponds to a 6-character version of the reference
file type, such as `R_DARK`, `R_LINEAR`, and `R_PHOTOM`.

Running From the Command Line
=============================
Individual steps and pipelines (consisting of a series of steps) can be run
from the command line using the ``strun`` command:
::

    $ strun <class_name or cfg_file> <input_file>

The first argument to ``strun`` must be either the python class name of the
step or pipeline to be run, or the name of a configuration (.cfg) file for the
desired step or pipeline (see `Configuration Files`_ below for more details).
The second argument to `strun` is the name of the input data file to be processed.

For example, running the full stage 1 pipeline or an individual step by
referencing their class names is done as follows:
::

  $ strun jwst.pipeline.Detector1Pipeline jw00017001001_01101_00001_nrca1_uncal.fits
  $ strun jwst.dq_init.DQInitStep jw00017001001_01101_00001_nrca1_uncal.fits

When a pipeline or step is executed in this manner (i.e. by referencing the
class name), it will be run using all default parameter values. The same thing
can be accomplished by using the default configuration file corresponding to
each:
::

  $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
  $ strun dq_init.cfg jw00017001001_01101_00001_nrca1_uncal.fits

If you want to use non-default parameter values, you can specify them as
keyword arguments on the command line or set them in the appropriate
cfg file.
To specify parameter values for an individual step when running a pipeline
use the syntax `--steps.<step_name>.<parameter>=value`.
For example, to override the default selection of a dark current reference
file from CRDS when running a pipeline:
::

    $ strun jwst.pipeline.Detector1Pipeline jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.dark_current.override_dark='my_dark.fits'
    $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.dark_current.override_dark='my_dark.fits'

You can get a list of all the available arguments for a given pipeline or
step by using the '-h' (help) argument to strun:
::

    $ strun dq_init.cfg -h
    $ strun jwst.pipeline.Detector1Pipeline -h

If you want to consistently override the default values of certain arguments
and don't want to specify them on the command line every time you
execute ``strun``, you can specify them in the configuration (.cfg) file for
the pipeline or the individual step.
For example, to always run ``Detector1Pipeline`` using the override in the
previous example, you could edit your `calwebb_detector1.cfg` file to
contain the following:
::

 name = "Detector1Pipeline"
 class = "jwst.pipeline.Detector1Pipeline"

    [steps]
      [[dark_current]]
        override_dark = 'my_dark.fits'

Note that simply removing the entry for a step from a pipeline cfg file will
**NOT** cause that step to be skipped when you run the pipeline (it will simply
run the step with all default parameters). In order to skip a step you must
use the `skip = True` argument for that step (see `Skip`_ below).

Alternatively, you can specify arguments for individual steps within the
step's configuration file and then reference those step cfg files in the pipeline
cfg file, such as:
::

 name = "Detector1Pipeline"
 class = "jwst.pipeline.Detector1Pipeline"

    [steps]
      [[dark_current]]
        config_file = my_dark_current.cfg

where `my_dark_current.cfg` contains:
::

 name = "dark_current"
 class = "jwst.dark_current.DarkCurrentStep"
 override_dark = 'my_dark.fits'


Exit Status
-----------
``strun`` produces the following exit status codes:

- 0: Successful completion of the step/pipeline
- 1: General error occurred
- 64: No science data found

The "No science data found" condition is returned by the ``assign_wcs`` step of
``calwebb_spec2`` when, after successfully determining the WCS solution for a
file, the WCS indicates that no science data will be found. This condition is
most often found with NIRSpec's NRS2 detector. There are certain optical and MSA
configurations in which dispersion will not cross to the NRS2 detector.

Running From Within Python
==========================

You can execute a pipeline or a step from within python by using the
``call`` method of the class:
::

 from jwst.pipeline import Detector1Pipeline
 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits')

 from jwst.linearity import LinearityStep
 result = LinearityStep.call('jw00001001001_01101_00001_mirimage_uncal.fits')

The easiest way to use optional arguments when calling a pipeline from
within python is to set those parameters in the pipeline cfg file and
then supply the cfg file as a keyword argument:
::

 Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits', config_file='calwebb_detector1.cfg')


Universal Parameters
====================

.. _intro_output_directory:

Output Directory
----------------

By default, all pipeline and step outputs will drop into the current
working directory, i.e., the directory in which the process is
running. To change this, use the `output_dir` argument. For example, to
have all output from `calwebb_detector1`, including any saved
intermediate steps, appear in the sub-directory `calibrated`, use
::

    $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
        --output_dir=calibrated

`output_dir` can be specified at the step level, overriding what was
specified for the pipeline. From the example above, to change the name
and location of the `dark_current` step, use the following
::

    $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
        --output_dir=calibrated
        --steps.dark_current.output_file='dark_sub.fits'
        --steps.dark_current.output_dir='dark_calibrated'

.. _intro_output_file:

Output File
-----------

When running a pipeline, the ``stpipe`` infrastructure automatically passes the
output data model from one step to the input of the next step, without
saving any intermediate results to disk. If you want to save the results from
individual steps, you have two options:

  - Specify `save_results`

    This option will save the results of the step, using a filename
    created by the step.

  - Specify a file name using `output_file <filename>`

    This option will save the step results using the name specified.

For example, to save the result from the dark current step of
`calwebb_detector1` in a file named `dark_sub.fits`, use

::

    $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.dark_current.output_file='dark_sub.fits'

You can also specify a particular file name for saving the end result of
the entire pipeline using the `--output_file` argument also
::
   
    $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
        --output_file='detector1_processed.fits'

Output File and Associations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stage 2 pipelines can take an individual file or an
:ref:`association <associations>` as input. Nearly all Stage 3
pipelines require an associaiton as input. Normally, the output file
is defined in each association's `product_name`.

If there is need to produce multiple versions of a calibration based
on an association, it is highly suggested to use `output_dir` to place
the results in a different directory instead of using `output_file` to
rename the output files. Stage 2 pipelines do not allow the override
of the output using `output_file`. Stage 3 pipelines do. However,
since Stage 3 pipelines generally produce many files per association,
using different directories via `output_dir` will make file keeping
simpler.

Override Reference File
-----------------------

For any step that uses a calibration reference file you always have the
option to override the automatic selection of a reference file from CRDS and
specify your own file to use. Arguments for this are of the form
`--override_\<ref_type\>`, where `ref_type` is the name of the reference file
type, such as `mask`, `dark`, `gain`, or `linearity`. When in doubt as to
the correct name, just use the `-h` argument to ``strun`` to show you the list
of available override arguments.

To override the use of the default linearity file selection, for example,
you would use:
::

  $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.linearity.override_linearity='my_lin.fits'

Skip
----

Another argument available to all steps in a pipeline is `skip`.
If `skip=True` is set for any step, that step will be skipped, with the
output of the previous step being automatically passed directly to the input
of the step following the one that was skipped. For example, if you want to
skip the linearity correction step, edit the calwebb_detector1.cfg file to
contain:
::

   [steps]
      [[linearity]]
        skip = True
      ...

Alternatively you can specify the `skip` argument on the command line:
::

    $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.linearity.skip=True

Logging Configuration
---------------------

If there's no `stpipe-log.cfg` file in the working directory, which specifies
how to handle process log information, the default is to display log messages
to stdout. If you want log information saved to a file, you can specify the
name of a logging configuration file either on the command line or in the
pipeline cfg file.

For example:
::

    $ strun calwebb_detector1.cfg jw00017001001_01101_00001_nrca1_uncal.fits
        --logcfg=pipeline-log.cfg

and the file `pipeline-log.cfg` contains:
::

    [*]
    handler = file:pipeline.log
    level = INFO

In this example log information is written to a file called `pipeline.log`.
The `level` argument in the log cfg file can be set to one of the standard
logging level designations of `DEBUG`, `INFO`, `WARNING`, `ERROR`, and
`CRITICAL`. Only messages at or above the specified level
will be displayed.


Input Files
===========

There are two general types of input to any stage: references files
and data files.  The references files, unless explicitly
overridden, are provided through CRDS.

The input data files - the exposure FITS files, association JSON files
and input catalogs - are presumed to all be in the same directory as
the primary input file. Sometimes the primary input is an association
JSON file, and sometimes it is an exposure FITS file.

Output File Names
=================

File names for the outputs from pipelines and steps come from
three different sources:

- The name of the input file
- The product name defined in an association
- As specified by the `output_file` argument

Regardless of the source, each pipeline/step uses the name as a "base
name", on to which several different suffixes are appended, which
indicate the type of data in that particular file.

.. _pipeline_step_suffix_definitions:

Pipeline/Step Suffix Definitions
--------------------------------

However the file name is determined (see above), the various stage 1,
2, and 3 pipeline modules will use that file name, along with a set of
predetermined suffixes, to compose output file names. The output file
name suffix will always replace any existing suffix of the input file
name. Each pipeline module uses the appropriate suffix for the
product(s) it is creating. The list of suffixes is shown in the
following table.

=============================================  ========
Product                                        Suffix
=============================================  ========
Uncalibrated raw input                         uncal
Corrected ramp data                            ramp
Corrected countrate image                      rate
Corrected countrate per integration            rateints
Optional fitting results from ramp_fit step    fitopt
Background-subtracted image                    bsub
Per integration background-subtracted image    bsubints
Calibrated image                               cal
Calibrated per integration images              calints
CR-flagged image                               crf
CR-flagged per integration images              crfints
1D extracted spectrum                          x1d
1D extracted spectra per integration           x1dints
Resampled 2D image                             i2d
Resampled 2D spectrum                          s2d
Resampled 3D IFU cube                          s3d
Source catalog                                 cat
Time Series photometric catalog                phot
Time Series white-light catalog                whtlt
Coronagraphic PSF image stack                  psfstack
Coronagraphic PSF-aligned images               psfalign
Coronagraphic PSF-subtracted images            psfsub
AMI fringe and closure phases                  ami
AMI averaged fringe and closure phases         amiavg
AMI normalized fringe and closure phases       aminorm
=============================================  ========

Individual Step Outputs
-----------------------

If individual steps are executed without an output file name specified via
the `output_file` argument, the `stpipe` infrastructure
automatically uses the input file name as the root of the output file name
and appends the name of the step as an additional suffix to the input file
name. If the input file name already has a known suffix, that suffix
will be replaced. For example:
::

 $ strun dq_init.cfg jw00017001001_01101_00001_nrca1_uncal.fits

produces an output file named
`jw00017001001_01101_00001_nrca1_dq_init.fits`.

Configuration Files
===================

Configuration (.cfg) files can be used to specify parameter values
when running a pipeline or individual steps, as well as for
specifying logging options.

You can use the ``collect_pipeline_cfgs`` task to get copies of all the cfg
files currently in use by the jwst pipeline software. The task takes a single
argument, which is the name of the directory to which you want the cfg files
copied. Use '.' to specify the current working directory, e.g.
::

 $ collect_pipeline_cfgs .

Each step and pipeline has their own cfg file, which are used to specify
relevant parameter values. For each step in a pipeline, the pipeline cfg file
specifies either the step's arguments or the cfg file containing the step's
arguments.

The name of a file in which to save log information, as well as the desired
level of logging messages, can be specified in an optional configuration file
"stpipe-log.cfg". This file must be in the same directory in which you run the
pipeline in order for it to be used. If this file does not exist, the default
logging mechanism is STDOUT, with a level of INFO. An example of the contents
of the stpipe-log.cfg file is:
::

    [*]
    handler = file:pipeline.log
    level = INFO

which specifies that all log messages will be directed to a file called
"pipeline.log" and messages at a severity level of INFO and above will be
recorded.

For a given step, the step's cfg file specifies parameters and their default
values; it includes parameters that are typically not changed between runs.
Parameters that are usually reset for each run are not included in the cfg file,
but instead specified on the command line. An example of a cfg file for the
jump detection step is:
::

    name = "jump"
    class = "jwst.jump.JumpStep"
    rejection_threshold = 4.0

You can list all of the parameters for this step using:
::

 $ strun jump.cfg -h

which gives the usage, the positional arguments, and the optional arguments.
More information on configuration files can be found in the ``stpipe`` User's
Guide at :ref:`stpipe-user-steps`.

Available Pipelines
===================
There are many pre-defined pipeline modules for processing
data from different instrument observing modes through each of the 3 stages
of calibration. For all of the details see :ref:`pipelines`.


For More Information
====================
More information on logging and running pipelines can be found in the ``stpipe``
User's Guide at :ref:`stpipe-user-steps`.

More detailed information on writing pipelines can be found
in the ``stpipe`` Developer's Guide at :ref:`stpipe-devel-steps`.
