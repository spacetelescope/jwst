.. _introduction:

Introduction to the JWST Pipeline
=================================

Introduction
------------

The JWST Science Calibration Pipeline processes data from all JWST instruments
and observing modes by applying various science corrections sequentially,
producing both fully-calibrated individual exposures and high-level data
products (mosaics, extracted spectra, etc.). The pipeline is written in Python,
hosted open-source on Github, and can be run either via
:ref:`command line interface <run_from_strun>` (`strun`) or via
the :ref:`Python interface <run_from_python>`.

The full end-to-end 'pipeline' (from raw data to high-level data products)
is comprised of three seperate pipeline stages that are run individually
to produce output products at different calibration levels:

	:Stage 1: Detector-level corrections and ramp fitting for individual
			  exposures.
	:Stage 2: Instrument-mode calibrations for individual exposures.
	:Stage 3: Combining data from multiple exposures within an observation

As such, the term 'pipeline' may refer to a single pipeline stage or to the full
three-stage series. 

Because JWST has many different instruments and observing modes, there are
several different pipeline modules available for each stage. There is one single
pipeline for Stage 1 - corrections are applied nearly universally for all
instruments and modes. There are two pipeline modules for Stage 2: one for
imaging and one for spectroscopic modes. Stage 3 is divided into five separate
modules for imaging, spectroscopic, coronagraphic, Aperture Masking
Interferometry (AMI), and Time Series Observation (TSO) modes. Details of all
the available pipeline modules can be found at :ref:`pipeline-modules`.

Each pipeline stage consists of a series of sequential steps (e.g, saturation
correction, ramp fitting). Each full pipeline stage and every individual step
has a unique module name (i.e Detector1Pipeline, or DarkCurrentCorrection).
Steps can also be run individually on data to apply a single correction. The
output of each pipeline stage is the input to the next, and within a pipeline
stage the output of each step is the input to the next.

The pipeline relies on three components to direct processing: input data,
step parameters, and reference files. The inputs to the pipeline modules are
individual exposures (`.fits` files) or associations of multiple exposures
(`asn.json` files). The parameters for each pipeline step are determined
hierarchically from the parameter defaults, parameter reference files, and any
specified overrides at run time. Finally, reference files provide data for each
calibration step that is specific to the dataset being processed. These files
may depend on things like instrument, observing mode, and date. In both the
command line and Python interface, a pipeline or step module may be configured
before running. Reference files can be overridden from those chosen by CRDS,
steps in a pipeline can be skipped, step parameters can be changed, and the
output and intermediate output products can be controlled. 

A pipeline (or individual step) outputs corrected data either by writing an output
file on disk  or returning an in-memory datamodel object. The output file suffix
(i.e `cal.fits`, `rate.fits`) depends on level of calibration - each full pipeline
stage as well as each individual step have a unique file suffix so that outputs
may be obtained at any level of calibration. Other pipeline outputs include
photometry catalogs and alignment catalogs (at stage 3).


Overview of Pipeline Code
-------------------------

The following is a brief overview of how the pipeline code in `jwst` is
organized.


**Pipeline and Step Classes** 

The JWST pipeline is organized into two main classes - `pipeline` classes and
`step` classes. Pipelines are made up of sequential `step` classes chained together,
the output of one step being piped to the next, but both pipelines and steps
are represented as objects that can be configured and run on input data.

::

	Detector1Pipeline  # an example of a pipeline class
	DarkCurrentStep    # an example of a step class

Each pipeline or step has a unique module name, which is the identifier used to
invoke the correct pipeline/step when using either the Python or the Command
Line Interface.

**Package Structure** 

Within the `jwst` repository, there are separate modules for each pipeline step. 
There is also a `pipeline` module, where the `pipeline` classes, consisting of
`step` classes called in sequence, are defined. 

::

	jwst/
		assign_wcs/
			assign_wcs_step.py  # contains AssignWcsStep
			...
		dark_current/
			dark_current_step.py  # contains DarkCurrent Step
			...
		pipeline/
			calwebb_detector1.py  # contains Detector1Pipeline
			calwebb_image2.py  # contains Image2Pipeline
		...

**Dependencies** 

The `jwst` package has several dependencies (see the `setup.cfg` file in the
top-level directory of `jwst` for a full list). Some notable dependencies
include:

**asdf**

`ASDF <https://asdf.readthedocs.io/en/latest/>`_, the Advanced Scientific Data
Format is the file format the JWST uses to encode world coordinate system (WCS)
information. 

**gwcs**

`GWCS <https://gwcs.readthedocs.io/en/latest/>`_, Generalized World Coordinate
System - is an generalized alternative to FITS WCS which makes use of astropy
models to describle the translation between detector and sky coordinates. In
JWST data, WCS information is encoded in an ASDF extension in the FITS file that
contains GWCS object. In contrast, FITS WCS is limited because it stores the WCS
transformation as header keywords, which is not sufficient to describe many of
the transformations JWST uses.

**stpipe**

`STPIPE <https://github.com/spacetelescope/stpipe>`_ contains base classes for
`pipeline` and `step`, and command line tools that are shared between the JWST
and `Nancy Grace Roman Telescope <https://roman-pipeline.readthedocs.io/en/latest/>`_
(Roman) pipelines.

**stcal**

The `stcal` package contains step code that is common to both JWST and the Roman 
telescope, to avoid redundancy. All step classes for the JWST
pipeline are still defined in `jwst`, but some of the underlying code for these
steps lives in `stcal` if the algorithm is shared by Roman (for example, ramp
fitting, saturation).
