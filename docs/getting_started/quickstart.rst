.. _quickstart:

=================
Quickstart Guide 
=================

The following is a quickstart guide to installing and running and the
latest version of `jwst`.

In short, the only setup required to run the JWST pipeline is to `pip` install
the `jwst` package into a `conda` environment, and then to set correct
environment variables for accessing reference files through CRDS. From there,
the JWST pipeline can be :ref:`run in a Python session<run_from_python>` or with
the :ref:`command line interface<run_from_strun>`, and
.. comment out until stdatamodels is released
.. ref JWST datamodels<data-models>
JWST datamodels
and other pipeline utilities can be imported
and used in a Python session.

**1. Create a conda environment.**

Python environments allow you to install different versions of packages and
their dependencies and keep them isolated from one another. While there are
several possible ways to achieve this (e.g `venv`), we will use `conda` in this
example.

If you don't already have `conda`, please follow the
`install instructions <https://docs.conda.io/en/latest/miniconda.html>`_.

To create a conda environment specifically for the latest stable release of
`jwst` (in this example, called jwst_latest):

::

	conda create --name jwst_latest python=3.10

This will create a new, (nearly) empty Python 3.10 environment in which you can
install the `jwst` package.

**2. Install jwst from PyPi**

Once you have created your conda environment, make sure it is active by doing:
::

	conda activate jwst_latest

To install the last stable release of `jwst`, and all its basic dependencies
(e.g numpy, stcal):

::

	pip install jwst

For detailed installation instructions, including how to install the development
version of `jwst` from Github or how to install a previous released version, see
the :ref:`installation` page.

**3. Set environment variables for Calibration References Data System (CRDS)**

CRDS is the system that manages the reference files needed to run the
pipeline. Inside the STScI network, the pipeline works with default CRDS
setup with no modifications. To run the pipeline outside the STScI
network, CRDS must be configured by setting two environment variables:
::

	export CRDS_PATH=$HOME/crds_cache
	export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

The `CRDS_PATH` is the directory on your filesystem that contains your local
CRDS cache, where reference files are accessed by the pipeline. The
`CRDS_SERVER_URL` variable specifies from which CRDS server reference files should
be obtained. For more information, see :ref:`reference_files_crds`.

**4. Running the Pipeline**

With `jwst` installed and CRDS configured for JWST, you can now run the pipeline
and use JWST `datamodels`.

For information on how to run the pipeline using the Python interface, see
:ref:`Running the JWST pipeline: Python Interface<run_from_python>`.

For information on how to run the pipeline using the command line interface, see
:ref:`Running the JWST pipeline: Command Line Interface<run_from_strun>`.

For information on how to read and write data files with JWST `datamodels`, see
.. comment out until stdatamodels is released
.. ref JWST datamodels<data-models>
JWST datamodels.
