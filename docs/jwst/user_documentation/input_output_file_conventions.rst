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
files. All files are assumed to be co-resident in the directory where the primary
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

The pipelines do not file manage versions. When re-running a pipeline, previous files
will be overwritten.


Output Files and Associations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stage 2 pipelines can take an individual file or an :ref:`association
<associations>` as input. Nearly all Stage 3 pipelines require an association as
input. Normally, the output file is defined in each association's "product name",
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
