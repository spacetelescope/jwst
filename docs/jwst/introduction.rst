.. _introduction:

Introduction
============

This document provides instructions on running the JWST Science Calibration
Pipeline (referred to as "the pipeline") and individual pipeline steps.

Multiple pipeline modules are used for different stages of processing and for
different JWST observing modes. The modules are broken into 3 stages:

 - Stage 1: Detector-level corrections and ramp fitting for individual exposures
 - Stage 2: Instrument-mode calibrations for individual exposures
 - Stage 3: Combining data from multiple exposures within an observation

Stage 1 corrections are applied nearly universally for all instruments and modes.
Stage 2 is divided into separate modules for imaging and spectroscopic modes.
Stage 3 is divided into five separate modules for imaging, spectroscopic,
coronagraphic, Aperture Masking Interferometry (AMI), and Time Series
Observation (TSO) modes.

Details of all the pipeline modules can be found at :ref:`pipeline-modules`. The
remainder of this document discusses the general usage of the pipelines and
steps.

Reference Files
===============

Many pipeline steps rely on the use of reference files that contain different types of
calibration data or information necessary for processing the data. The reference files are
instrument-specific and are periodically updated as the data processing evolves and the
understanding of the instruments improves. They are created, tested, and validated by the
JWST Instrument Teams. They ensure all the files are in the correct format and have all
required header keywords. The files are then delivered to the Reference Data for Calibration
and Tools (ReDCaT) Management Team. The result of this process is the files being ingested
into the JWST Calibration Reference Data System (CRDS), and made available to the pipeline
team and any other ground subsystem that needs access to them.

Information about all the reference files used by the Calibration Pipeline can be found at
:ref:`reference_file_information`,
as well as in the documentation for each Calibration Step that uses a reference file.
Information on reference file types and their correspondence to calibration steps is
described within the table at
:ref:`reference_file_types`.


CRDS
====

CRDS reference file mappings are usually set by default to always give access
to the most recent reference file deliveries and selection rules. On
occasion it might be necessary or desirable to use one of the non-default
mappings in order to, for example, run different versions of the pipeline
software or use older versions of the reference files. This can be
accomplished by setting the environment variable ``CRDS_CONTEXT`` to the
desired project mapping version, e.g.
::

$ export CRDS_CONTEXT='jwst_0421.pmap'

For all information about CRDS, including context lists, see the JWST CRDS
website:

    `https://jwst-crds.stsci.edu/ <https://jwst-crds.stsci.edu/>`_

Each pipeline step records the reference file that it used in the value of
a header keyword in the output data file. The keyword names use the syntax
"R_<ref>", where <ref> corresponds to a 6-character version of the reference
file type, such as ``R_DARK``, ``R_LINEAR``, and ``R_PHOTOM``.

There are two servers available:

  - JWST OPS: https://jwst-crds.stsci.edu
  - JWST PUB: https://jwst-crds-pub.stsci.edu

JWST OPS supports the automatic processing pipeline at STScI. JWST PUB supports
the latest public release of the `jwst` package. Most often, the reference
contexts are one and the same. Regardless, if one wishes to calibrate using the
same exact information as the automatic processing, use JWST OPS. Otherwise, use
of JWST PUB is recommended.

Inside the STScI network, the pipeline defaults the CRDS setup to use JWST OPS with no modifications.
To run the pipeline outside the STScI network or to use a different server, CRDS must be configured by setting
two environment variables:

  - CRDS_PATH: Local folder where CRDS content will be cached.
  - CRDS_SERVER_URL: The server from which to pull reference information

To setup to use JWST OPS, use the following settings:

::

    export CRDS_PATH=$HOME/crds_cache/jwst_ops
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

To setup to use JWST PUB, use the following settings:

::

    export CRDS_PATH=$HOME/crds_cache/jwst_pub
    export CRDS_SERVER_URL=https://jwst-crds-pub.stsci.edu

.. _strun_command_line:

Running From the Command Line
=============================

.. note::

   For seasoned users who are familiar with using ``collect_pipeline_cfgs`` and
   running pipelines by the default configuration (CFG) files, please note that
   this functionality has been deprecated. Please read :ref:`CFG Usage
   Deprecation Notice<cfg_usage_deprecation_notice>`.

Individual steps and pipelines (consisting of a series of steps) can be run
from the command line using the ``strun`` command:
::

    $ strun <pipeline_name, class_name, or parameter_file> <input_file>

The first argument to ``strun`` must be one of either a pipeline name, python
class of the step or pipeline to be run, or the name of a parameter file for the
desired step or pipeline (see :ref:`parameter_files`). The second argument to
``strun`` is the name of the input data file to be processed.

For example, the Stage 1 pipeline is implemented by the class
:ref:`jwst.pipeline.Detector1Pipeline <calwebb_detector1>`. The command to run this pipeline is as
follows:
::

  $ strun jwst.pipeline.Detector1Pipeline jw00017001001_01101_00001_nrca1_uncal.fits

Pipeline classes also have a **pipeline name**, or **alias**, that can be used instead of the
full class specification. For example, ``jwst.pipeline.Detector1Pipeline`` has the
alias ``calwebb_detector1`` and can be run as
::

  $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits

A full list of pipeline aliases can be found in :ref:`Pipeline Stages <pipelines>`.

Exit Status
-----------
``strun`` produces the following exit status codes:

- 0: Successful completion of the step/pipeline
- 1: General error occurred
- 64: No science data found

The "No science data found" condition is returned by the ``assign_wcs`` step of
the ``calwebb_spec2`` pipeline when, after successfully determining the WCS
solution for a file, the WCS indicates that no science data will be found. This
condition most often occurs with NIRSpec's Multi-object Spectroscopy (MOS) mode:
There are certain optical and MSA configurations in which dispersion will not
cross one or the other of NIRSpec's detectors.

.. _intro_file_conventions:

Input and Output File Conventions
=================================

.. _intro_input_file_discussion:

Input Files
-----------

There are two general types of input to any step or pipeline: references files
and data files.  The references files, unless explicitly
overridden, are provided through CRDS.

Data files are the science input, such as exposure FITS files and association
files. All input files must be co-resident in the directory where the primary
input file is located. This is particularly important for associations: JWST
associations contain file names only. All files referred to by an association
are expected to be located in the directory in which the association file is located.

.. _intro_output_file_discussion:

Output Files
------------

Output files will be created either in the current working directory, or where
specified by the :ref:`output_dir <intro_output_directory>` parameter.

File names for the outputs from pipelines and steps come from
three different sources:

- The name of the input file
- The product name defined in an association
- As specified by the :ref:`output_file <intro_output_file>` parameter

Regardless of the source, each pipeline/step uses the name as a base
name, onto which several different suffixes are appended, which
indicate the type of data in that particular file. A list of the main suffixes
can be :ref:`found below <pipeline_step_suffix_definitions>`.

The pipelines do not manage versions. When re-running a pipeline, previous files
will be overwritten.


Output Files and Associations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stage 2 pipelines can take an individual file or an :ref:`association
<associations>` as input. Nearly all Stage 3 pipelines require an association as
input. Normally, the output file is defined in each association's "product name"
which defines the basename that will be used for output file naming.

Often, one may reprocess the same set of data multiple times, such as to change
reference files or parameters. When doing so, it is highly suggested to use
``output_dir`` to place the results in a different directory instead of using
``output_file`` to rename the output files. Most pipelines and steps create sets
of output files. Separating runs by directory may be much easier to manage.


Individual Step Outputs
^^^^^^^^^^^^^^^^^^^^^^^

If individual steps are executed without an output file name specified via
the ``output_file`` parameter, the ``stpipe`` infrastructure
automatically uses the input file name as the root of the output file name
and appends the name of the step as an additional suffix to the input file
name. If the input file name already has a known suffix, that suffix
will be replaced. For example:
::

   $ strun jwst.dq_init.DQInitStep jw00017001001_01101_00001_nrca1_uncal.fits

produces an output file named
``jw00017001001_01101_00001_nrca1_dq_init.fits``.

See :ref:`pipeline_step_suffix_definitions` for a list of the more common
suffixes used.

Parameters
==========

All pipelines and steps have **parameters** that can be set to change various
aspects of how they execute. To see what parameters are available for any given
pipeline or step, use the ``-h`` option on ``strun``. Some examples are:
::

   $ strun calwebb_detector1 -h
   $ strun jwst.dq_init.DQInitStep -h

To set a parameter, simply specify it on the command line. For example, to have
:ref:`calwebb_detector1 <calwebb_detector1>` save the calibrated ramp files, the
``strun`` command would be as follows:
::

   $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits --save_calibrated_ramp=true

To specify parameter values for an individual step when running a pipeline
use the syntax ``--steps.<step_name>.<parameter>=value``.
For example, to override the default selection of a dark current reference
file from CRDS when running a pipeline:
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.dark_current.override_dark='my_dark.fits'

If there is need to re-use a set of parameters often, parameters can be stored
in **parameter files**. See :ref:`parameter_files` for more information.

Universal Parameters
--------------------

The set of parameters that are common to all pipelines and steps are referred to
as **universal parameters** and are described below.

.. _intro_output_directory:

Output Directory
^^^^^^^^^^^^^^^^

By default, all pipeline and step outputs will drop into the current
working directory, i.e., the directory in which the process is
running. To change this, use the ``output_dir`` parameter. For example, to
have all output from ``calwebb_detector1``, including any saved
intermediate steps, appear in the sub-directory ``calibrated``, use
::
   
    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --output_dir=calibrated

``output_dir`` can be specified at the step level, overriding what was
specified for the pipeline. From the example above, to change the name
and location of the ``dark_current`` step, use the following
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --output_dir=calibrated
        --steps.dark_current.output_file='dark_sub.fits'
        --steps.dark_current.output_dir='dark_calibrated'

.. _intro_output_file:

Output File
^^^^^^^^^^^

When running a pipeline, the ``stpipe`` infrastructure automatically passes the
output data model from one step to the input of the next step, without
saving any intermediate results to disk. If you want to save the results from
individual steps, you have two options:

  - Specify ``save_results``.
    This option will save the results of the step, using a filename
    created by the step.

  - Specify a file name using ``output_file <basename>``.
    This option will save the step results using the name specified.

For example, to save the result from the dark current step of
``calwebb_detector1`` in a file named based on ``intermediate``, use

::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.dark_current.output_file='intermediate'

A file, ``intermediate_dark_current.fits``, will then be created. Note that the
suffix of the step is always appended to any given name.

You can also specify a particular file name for saving the end result of
the entire pipeline using the ``--output_file`` parameter also
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --output_file='stage1_processed'

In this situation, using the default configuration, three files are created:

  - ``stage1_processed_trapsfilled.fits``
  - ``stage1_processed_rate.fits``
  - ``stage1_processed_rateints.fits``


Override Reference File
^^^^^^^^^^^^^^^^^^^^^^^

For any step that uses a calibration reference file you always have the
option to override the automatic selection of a reference file from CRDS and
specify your own file to use. Parameters for this are of the form
``--override_<ref_type>``, where ``ref_type`` is the name of the reference file
type, such as ``mask``, ``dark``, ``gain``, or ``linearity``. When in doubt as to
the correct name, just use the ``-h`` argument to ``strun`` to show you the list
of available override parameters.

To override the use of the default linearity file selection, for example,
you would use:
::

  $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
          --steps.linearity.override_linearity='my_lin.fits'

Skip
^^^^

Another parameter available to all steps in a pipeline is ``skip``. If
``skip=True`` is set for any step, that step will be skipped, with the output of
the previous step being automatically passed directly to the input of the step
following the one that was skipped. For example, if you want to skip the
linearity correction step, one can specify the ``skip`` parameter for the
``strun`` command:
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --steps.linearity.skip=True

Alternatively, if using a :ref:`parameter file<parameter_files>`, edit the
file to add the following snippet:
::

  steps:
  - class: jwst.linearity.linearity_step.LinearityStep
    parameters:
      skip: true


Logging Configuration
=====================

The name of a file in which to save log information, as well as the desired
level of logging messages, can be specified in an optional configuration file.
Two options exist - if the configuration file should be used for all instances
of the pipeline, the configuration file should be named "stpipe-log.cfg".
This file must be in the same directory in which you run the pipeline in order
for it to be used.

If instead the configuration should be active only when specified,
you should name it something other than "stpipe-log.cfg"; this filename should be
specified using either the ``--logcfg`` parameter to the command line ``strun`` or
using the ``logcfg`` keyword to a .call() execution of either a Step or Pipeline
instance.

If this file does not exist, the default logging mechanism is STDOUT,
with a level of INFO. An example of the contents of the stpipe-log.cfg file is:

::

    [*]
    handler = file:pipeline.log
    level = INFO

If there's no ``stpipe-log.cfg`` file in the working directory, which specifies
how to handle process log information, the default is to display log messages
to stdout.

For example:
::

    $ strun calwebb_detector1 jw00017001001_01101_00001_nrca1_uncal.fits
        --logcfg=pipeline-log.cfg

Or in an interactive python environment:
::

    result = Detector1Pipeline.call("jw00017001001_01101_00001_nrca1_uncal.fits",
                                    logcfg="pipeline-log.cfg")

and the file ``pipeline-log.cfg`` contains:
::

    [*]
    handler = file:pipeline.log
    level = INFO

In this example log information is written to a file called ``pipeline.log``.
The ``level`` argument in the log cfg file can be set to one of the standard
logging level designations of ``DEBUG``, ``INFO``, ``WARNING``, ``ERROR``, and
``CRITICAL``. Only messages at or above the specified level
will be displayed.

.. note::

   Setting up ``stpipe-log.cfg`` can lead to confusion, especially if it is
   forgotten about. If one has not run a pipeline in awhile, and then sees no
   logging information, most likely it is because ``stpipe-log.cfg`` is
   present. Consider using a different name and specifying it explicitly on the
   command line.

.. _run_from_python:

Running From Within Python
==========================

You can execute a pipeline or a step from within python by using the
``call`` method of the class.

The ``call`` method creates a new instance of the class and runs the pipeline or
step. Optional parameter settings can be specified by via keyword arguments or
supplying a parameter file. Some examples are shown below. For more information,
see :ref:`Execute via call()<call_examples>`::

 from jwst.pipeline import Detector1Pipeline
 result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits')

 from jwst.linearity import LinearityStep
 result = LinearityStep.call('jw00001001001_01101_00001_mirimage_uncal.fits')

For more details on the different ways to run a pipeline step, see
the :ref:`Configuring a Step<configuring-a-step>` page.

CRDS Environment Variables
--------------------------

The CRDS environment variables need to be defined *before* importing anything
from `jwst` or `crds`. In general, any scripts should assume the environment
variables have been set before the scripts have run. If one needs to define the
CRDS environment variables within a script, the following code snippet is the
suggested method. These lines should be the first executable lines:

::

   import os
   os.environ['CRDS_PATH'] = 'path_to_local_cache'
   os.environ['CRDS_SERVER_URL'] = 'url-of-server-to-use'

   # Now import anything else needed
   import jwst


Available Pipelines
===================
There are many pre-defined pipeline modules for processing
data from different instrument observing modes through each of the 3 stages
of calibration. For all of the details see :ref:`pipelines`.

.. _pipeline_step_suffix_definitions:

Pipeline/Step Suffix Definitions
--------------------------------

However the output file name is determined (:ref:`see above
<intro_output_file_discussion>`), the various stage 1, 2, and 3 pipeline modules
will use that file name, along with a set of predetermined suffixes, to compose
output file names. The output file name suffix will always replace any known
suffix of the input file name. Each pipeline module uses the appropriate suffix
for the product(s) it is creating. The list of suffixes is shown in the
following table. Replacement occurs only if the suffix is one known to the
calibration code. Otherwise, the new suffix will simply be appended to the
basename of the file.

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
Resampled 2D image                             i2d
Resampled 2D spectrum                          s2d
Resampled 3D IFU cube                          s3d
1D extracted spectrum                          x1d
1D extracted spectra per integration           x1dints
1D combined spectrum                           c1d
Source catalog                                 cat
Segmentation map                               segm
Time Series photometric catalog                phot
Time Series white-light catalog                whtlt
Coronagraphic PSF image stack                  psfstack
Coronagraphic PSF-aligned images               psfalign
Coronagraphic PSF-subtracted images            psfsub
AMI fringe and closure phases                  ami
AMI averaged fringe and closure phases         amiavg
AMI normalized fringe and closure phases       aminorm
=============================================  ========


For More Information
====================
More information on logging and running pipelines can be found in the ``stpipe``
User's Guide at :ref:`stpipe-user-steps`.

More detailed information on writing pipelines can be found
in the ``stpipe`` Developer's Guide at :ref:`stpipe-devel-steps`.

If you have questions or concerns regarding the software, please open an issue
at https://github.com/spacetelescope/jwst/issues or contact
the `JWST Help Desk <https://jwsthelp.stsci.edu>`_.
