Introduction
============

This document provides instructions on running the JWST Science Calibration
Pipeline (referred to as "the pipeline") and individual pipeline steps.
 
The pipeline currently consists of ramps-to-slopes processing for all
observing modes, Level-2b processing for all imaging modes, and Level-2b
processing for some spectroscopic modes (mainly fixed slit and single
object slitless).

The ramps-to-slopes (Level-2a) processing consists of

* DQ initialization
* saturation detection
* Inter-Pixel Capacitance (IPC) correction
* superbias subtraction
* bias drift correction
* reset correction (MIRI only)
* last frame correction (MIRI only)
* linearity correction
* dark subtraction
* jump (CR) detection
* ramp fitting

Level-2b processing consists of

* WCS and distortion model assignment
* slit extraction (only for some spectroscopy modes)
* flat fielding
* persistence correction (currently a no-op)
* straylight correction (only for MIRI MRS)
* fringe correction (only for MIRI MRS)
* telescope emission correction (currently a no-op)  
* photometric calibration assignment

Details of the available pipeline configurations can be found at:

http://ssb.stsci.edu/doc/jwst_dev/jwst.pipeline.doc/html/index.html

This document contains setup instructions, discussion of the pipeline 
configuration files, and examples of running pipelines either as
a whole or in individual steps.

This document is a work in progress and will be updated frequently.  The most
recent version of this document is built nightly from the subversion source
code repository.  The most recent version of this document can be found on the
site:

http://ssb.stsci.edu/doc/jwst_dev/


Software Access
===============

The JWST calibration pipeline software is pre-installed, along with the
rest of the STScI SSB python-based software packages that are distributed
via the Ureka system, on Linux cluster systems at STScI, such as "jwcalibdev"
and "witserv."
The jwst pipeline can be installed on a private or non-cluster system
using the SSB Ureka installer and installing "jwst" as an add-on package.
Instructions for installing and configuring the Ureka environment are at:

http://ssb.stsci.edu/ssb_software.shtml

We suggest the following procedures for setting up the environment and
necessary working directories on your system:

* type 'ssbx' (on Linux cluster systems) or 'ur_setup common ssbx' (on Macs and other private systems)
* create a working directory in which to hold your test data files and pipeline configuration files, and 'cd' to it
* type 'collect_pipeline_cfgs .' to copy all the pipeline configuration files to your current working directory
* copy JWST data files into your working directory

Sample JWST data files that are in the correct Level-1b format to use as
input to the calibration pipeline can be found at 
\/grp\/jwst\/ssb\/sample_jwst_data\/level1b_wcs\/.


CRDS Reference Files
====================

CRDS reference file mappings are usually set by default to always give access
to the most recent reference file deliveries and selection rules. On
occasion it might be necessary or desirable to use one of the non-default
mappings in order to, for example, run different versions of the pipeline
software or use older versions of the reference files. This can be
accomplished by setting the environment variable 'CRDS_CONTEXT' to the
desired project mapping version, e.g.
::

$ setenv CRDS_CONTEXT jwst_0034.pmap

The current storage location for all JWST CRDS reference files is:
::

/grp/crds/jwst/references/jwst/

Each pipeline step records the reference file that it used in the value of
a header keyword in the output data file. The keyword names use the syntax
"R_<ref>", where <ref> corresponds to the first 6 characters of the reference
file type, such as "R_DARK", "R_LINEAR", and "R_PHOTOM".


Running From the Command Line
=============================

Individual steps and pipelines (consisting of a series of steps) can be run
from the command line using the `strun` command:
::

    $ strun <class_name or cfg_file> <input_file>

The first argument to `strun` must be either the python class name of the
step or pipeline to be run, or the name of a configuration (.cfg) file for the
desired step or pipeline (see `Configuration Files`_ below for more details).
The second argument to `strun` is the name of the input data file to be processed.

For example, running the full ramps-to-slopes pipeline or an individual step by
referencing their class names is done as follows:
::

  $ strun jwst.pipeline.SloperPipeline jw00017001001_01101_00001_NRCA1_uncal.fits
  $ strun jwst.dq_init.DQInitStep jw00017001001_01101_00001_NRCA1_uncal.fits

When a pipeline or step is executed in this manner (i.e. by referencing the 
class name), it will be run using all default parameter values. The same thing
can be accomplished by using the default configuration file corresponding to
each:
::

  $ strun calwebb_sloper.cfg jw00017001001_01101_00001_NRCA1_uncal.fits
  $ strun dq_init.cfg jw00017001001_01101_00001_NRCA1_uncal.fits

If you want to use non-default parameter values, you can specify them as
keyword arguments on the command line.
To specify parameter values for an individual step when running a pipeline
use the syntax `--steps.<step_name>.<parameter>=value`.
For example, to override the default selection of a dark current reference
file from CRDS when running a pipeline:
::

    $ strun jwst.pipeline.SloperPipeline jw00017001001_01101_00001_NRCA1_uncal.fits
          --steps.dark_current.override_dark='my_dark.fits'
    $ strun calwebb_sloper.cfg jw00017001001_01101_00001_NRCA1_uncal.fits
          --steps.dark_current.override_dark='my_dark.fits'

You can get a list of all the available arguments for a given pipeline or
step by using the '-h' (help) argument to strun:
::

    $ strun dq_init.cfg -h
    $ strun jwst.pipeline.SloperPipeline -h

If you want to consistently override the default values of certain arguments
and don't want to have to specify them on the command line every time you
execute `strun`, you can specify them in the configuration (.cfg) file for
either the pipeline or the individual step.
For example, to always run 'SloperPipeline' using the override in the
previous example, you could edit your 'calwebb_sloper.cfg' file to
contain the following:
::

 name = "SloperPipeline"
 class = "jwst.pipeline.SloperPipeline"

    [steps]
      [[dark_current]]
        override_dark = 'my_dark.fits'

Note that simply removing the entry for a step from a pipeline cfg file will
**NOT** cause that step to be skipped when you run the pipeline (it will simply
run the step with all default parameters). In order to skip a step you must
use the 'skip = True' argument for that step (see `Skip`_ below).

Alternatively, you can specify arguments for individual steps within the
step configuration file and then reference those step cfg files in the pipeline
cfg file, such as:
::

 name = "SloperPipeline"
 class = "jwst.pipeline.SloperPipeline"

    [steps]
      [[dark_current]]
        config_file = my_dark_current.cfg

where "my_dark_current.cfg" contains:
::

 name = "dark_current" 
 class = "jwst.dark_current.DarkCurrentStep"
 override_dark = 'my_dark.fits'


Running From Within Python
==========================

You can execute a pipeline or a step from within python by using the 
`call` method of the class:
::

 from jwst.pipeline import SloperPipeline
 SloperPipeline.call('jw00017001001_01101_00001_NRCA1_uncal.fits')

 from jwst.linearity import LinearityStep
 LinearityStep.call('jw00001001001_01101_00001_MIRIMAGE_uncal.fits')

The easiest way to use optional arguments when calling a pipeline from
within python is to set those parameters in the pipeline cfg file and
then supply the cfg file as a keyword argument:
::

 SloperPipeline.call('jw00017001001_01101_00001_NRCA1_uncal.fits', config_file='calwebb_sloper.cfg')


Universal Parameters
====================

Output File
-----------

When running a pipeline, the `stpipe` infrastructure automatically passes the 
output data model from one step to the input of the next step, without
saving any intermediate results to disk. If you want to save the results from
individual steps, you can use the `output_file` argument for each step.
For example, if you specify
::

    $ strun calwebb_sloper.cfg jw00017001001_01101_00001_NRCA1_uncal.fits
        --steps.dark_current.output_file='dark_sub.fits'

the results at the end of the dark current subtraction step would be saved
to the file `dark_sub.fits`.

You can also specify a particular file name for saving the end result of
the entire pipeline using the `--output_file` argument.

Override Reference File
-----------------------

For any step that uses a calibration reference file you always have the
option to override the automatic selection of a reference file from CRDS and
specify your own file to use. Arguments for this are of the form
`--override_\<ref_type\>`, where `ref_type` is the name of the reference file
type, such as `mask`, `dark`, `gain`, or `linearity`. When in doubt as to
the correct name, just use the `-h` argument to `strun` to show you the list
of available override arguments.

To override the use of the default linearity file selection, for example,
we would use:
::

  $ strun calwebb_sloper.cfg jw00017001001_01101_00001_NRCA1_uncal.fits
          --steps.linearity.override_linearity='my_lin.fits'

Skip
----

Another argument available to all steps in a pipeline is `skip`.
If 'skip=True' is set for any step, that step will be skipped, with the
output of the previous step being automatically passed directly to the input
of the step following the one that was skipped. For example, when
applying Level-2b spectroscopic processing to MIRI LRS data we do not want
to run the extract_2d step, which extracts a region of the image around the
spectrum. So the cfg appears as:
::

   [steps]
      [[extract_2d]]
        skip = True
      ...

Alternatively you can specify the `skip` argument on the command line:
::

    $ strun calwebb_sloper.cfg jw00017001001_01101_00001_NRCA1_uncal.fits
        --steps.dark_current.skip=True

Logging Configuration
---------------------

If there's no `stpipe-log.cfg` file in the working directory, which specifies
how to handle process log information, the default is to display log messages
to stdout. If you want log information saved to a file, you can specify the
name of a logging configuration file either on the command line or in the
pipeline cfg file. 

For example:
::

    $ strun calwebb_sloper.cfg jw00017001001_01101_00001_NRCA1_uncal.fits
        --logcfg=pipeline-log.cfg

and the file `pipeline-log.cfg` contains:
::

    [*]
    handler = file:pipeline.log
    level = INFO

In this example log information is written to a file called "pipeline.log." 
The `level` argument in the log cfg file can be set to one of the standard
logging level designations of `DEBUG`, `INFO`, `WARNING`, `ERROR`, and 
`CRITICAL`. Only messages at or above the specified level
will be displayed.


Output File Names
=================

Pipelines and steps will use default output file names or names provided by
the user via the `output_file` argument. In the absence of a user-specified
output file name, pipelines and steps use different schemes for setting a
default output name, which are explained below.

Pipeline Outputs
----------------

In the absence of a user-specified output file name, the various level-2a,
2b, and 3 pipeline modules will use the input root file name along with a set
of predetermined suffixes to compose output file names. The output file name
suffix will always replace the suffix of the input file name. Each pipeline
module uses the appropriate suffix for the product(s) it is creating. The
list of suffixes is shown in the following table.

=============================================  ========
Product                                        Suffix
=============================================  ========
Uncalibrated Level-1b input                    uncal
Partially calibrated Level-2a ramp data        ramp
Partially calibrated Level-2a countrate image  rate
Level-2a countrate per integration             rateints
Optional fitting results from ramp_fit step    fitopt
Level-2b background-subtracted image           bsub
Per integration background-subtracted image    bsubints
Fully-calibrated Level-2b image                cal
Fully-calibrated per integration images        calints
1D extracted spectrum                          spec
1D extracted spectra per integration           specints
Drizzled/combined image                        drz
3D IFU cube                                    cube
=============================================  ========

Individual Step Outputs
-----------------------

If individual steps are executed without an output file name specified via
the `output_file` argument, the `stpipe` infrastructure
automatically uses the input file name as the root of the output file name
and appends the name of the step as an additional suffix to the input file
name. For example:
::

 $ strun dq_init.cfg jw00017001001_01101_00001_NRCA1_uncal.fits

produces an output file named 
"jw00017001001_01101_00001_NRCA1_uncal_dq_init.fits."

Configuration Files
===================

Configuration (.cfg) files can be used to specify parameter values
when running a pipeline or individual steps, as well as for
specifying logging options.

You can use the "collect_pipeline_cfgs" task to get copies of all the cfg
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
    do_yintercept = False
    yint_threshold = 1.0

You can list all of the parameters for this step using:
::

 $ strun jump.cfg -h

which gives the usage, the positional arguments, and the optional arguments.
More information on configuration files can be found in the `stpipe` User's
Guide at:

http://ssb.stsci.edu/doc/jwst_dev/jwst.stpipe.doc/html/index.html

Available Pipelines
===================

There are currently several pre-defined pipelines available for processing 
the data from different instrument observing modes. For all of the details
see:

http://ssb.stsci.edu/doc/jwst_dev/jwst.pipeline.doc/html/index.html


For More Information
====================

More information on logging and running pipelines can be found in the `stpipe`
User's Guide at:

http://ssb.stsci.edu/doc/jwst_dev/jwst.stpipe.doc/html/user/index.html

More detailed information on writing pipelines can be found 
in the `stpipe` Developer's Guide at:

http://ssb.stsci.edu/doc/jwst_dev/jwst.stpipe.doc/html/devel/index.html


Another Source of JWST Test Data
================================

Sample JWST Level-1b data files can be found on the system
"jwcalibdev.stsci.edu".  The data conforms to the latest version of Daryl
Swade's Level 1 and 2 Data Product Design document (JWST-STScI-002111 Revision
A). They contain no ERR or DQ HDU's, as will be the state of Level-1b products
generated by SDP.

The data can be found in
::

    /grp/jwst/ssb/sample_jwst_data/level1b_wcs/dg000xx

where `xx` goes from 01 to 36.

These data should be treated as private and not shared outside of STScI.

