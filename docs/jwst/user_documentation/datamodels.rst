.. _jwst-data-models:

===============
JWST Datamodels
===============

The `stdatamodels` package contains the interface for JWST Datamodels. ``Datamodels``
are the recommended way of reading and writing JWST data files (.fits) and
reference files (.fits and .asdf). JWST data are encoded in FITS files, and reference
files consist of a mix of FITS and ASDF - datamodels were designed to
abstract away these intricacies and provide a simple interface to the data. They
represent the data in FITS extensions and meta data in FITS headers in a Python object
with a tree-like structure. The following section gives a brief overview of
``Datamodels`` as they pertain to the pipeline - see :ref:`stdatamodels:data-models`
for more detailed documentation on Datamodels.

Datamodels and the JWST pipeline
================================

When :ref:`running the pipeline in python <run_from_python>`, the inputs and 
outputs of running a pipeline or a step are JWST ``Datamodels``. 

The input to a pipeline/step can be a ``Datamodel``, created from an input
file on disk. E.g:

::

	# running a single pipeline step, input is datamodel object
	from jwst.linearity import LinearityStep
	import stdatamodels.jwst.datamodels as dm
	input_model = dm.open('jw00001001001_01101_00001_mirimage_uncal.fits')
	result = LinearityStep.call(input_model)

If a string path to a file on disk is passed in, a ``DataModel`` object will be
created internally when the pipeline/step is run.

By default, when running in Python, the corrected data will be returned in-memory
as a ``DataModel`` instead of being written as an output file.
See :ref:`controlling output file behavior`<python_outputs>` for instructions on
how to write the returned ``DataModel`` to an output file.

Groups of Datamodels
====================

Many of the JWST calibration steps and pipelines expect an
:ref:`Association <associations>` file as input. When opened with
:meth:`~jwst.datamodels.open`, a
:class:`~jwst.datamodels.ModelContainer` is returned. `ModelContainer`
is a list-like object where each element is the
`DataModel` of each member of the association. The `asn_table` attribute is
populated with the association data structure, allowing direct access
to the association itself.  The association file, as well as the files
listed in the association file, must be in the input directory.

Data structures that handle groups of datamodels are stored in the `jwst`
repository instead of inside `stdatamodels`. The API for these data structures
is documented here:

.. automodapi:: jwst.datamodels.container
	:no-inheritance-diagram:

.. automodapi:: jwst.datamodels.library
	:no-inheritance-diagram:

.. automodapi:: jwst.datamodels.source_container
	:no-inheritance-diagram:
