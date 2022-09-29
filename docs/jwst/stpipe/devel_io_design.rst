.. _step_io_design:

===============
Step I/O Design
===============

API Summary
===========

`Step` command-line options
---------------------------


    - `--output_dir`: :ref:`Directory <intro_output_directory>` where all output will go.
    - `--output_file`: :ref:`File name <intro_output_file>` upon which
      output files will be based.

`Step` configuration options
----------------------------

    - `output_dir`: :ref:`Directory <intro_output_directory>` where all output will go.
    - `output_file`: :ref:`File name <intro_output_file>` upon which
      output files will be based.
    - `suffix`: :ref:`Suffix <pipeline_step_suffix_definitions>` defining the output of this step.
    - `save_results`: True to create output files. :ref:`[more] <devel_io_when_files_are_created>`
    - `search_output_file`: True to retrieve the `output_file` from a
      parent `Step` or `Pipeline`. :ref:`[more]<devel_io_substeps_and_output>`
    - `output_use_model`: True to always base output file names on the
      `DataModel.meta.filename` of the `DataModel` being saved.
    - `input_dir`: Generally defined by the location of the primary
      input file unless otherwise specified.  All input files must be
      in this directory.

Classes, Methods, Functions
---------------------------

    - :meth:`Step.open_model <jwst.stpipe.step.Step.open_model>`: Open
      a `DataModel`
    - :meth:`Step.load_as_level2_asn`: Open a list or file as Level2 association.
    - :meth:`Step.load_as_level3_asn`: Open a list or file as Level3 association.
    - :meth:`Step.make_input_path
      <jwst.stpipe.step.Step.make_input_path>`: Create a file name to
      be used as input
    - :meth:`Step.save_model <jwst.stpipe.step.Step.save_model>`: Save a `DataModel` immediately
    - :attr:`Step.make_output_path
      <jwst.stpipe.step.Step._make_output_path>`: Create a file name
      to be used as output

Design
======

The :class:`~jwst.stpipe.step.Step` architecture is designed such that
a `Step`'s intended sole responsibility is to perform the calculation
required. Any input/output operations are handled by the surrounding
`Step` architecture. This is to help facilitate the use of `Step`'s
from both a command-line environment, and from an interactive Python
environment, such as Jupyter notebooks or `ipython`.

For command-line usage, all inputs and outputs are designed to come
from and save to files.

For interactive Python use, inputs and outputs are expected to be
Python objects, negating the need to save and reload data after every
`Step` call. This allows users to write Python scripts without having
to worry about doing I/O at every step, unless, of course, if the user
wants to do so.

The high-level overview of the input/output design is given in
:ref:`writing-a-step`. The following discusses the I/O API and
best practices.

To facilitate this design, a basic `Step` is suggested to have the
following structure::

  class MyStep(jwst.stpipe.step.Step):

      spec = ''  # Desired configuration parameters

      def process(self, input):

          with jwst.datamodels.open(input) as input_model:

              # Do awesome processing with final result
              # in `result`
              result = final_calculation(input_model)

          return (result)

When run from the command line::

  strun MyStep input_data.fits

the result will be saved in a file called::

  input_data_mystep.fits

Similarly, the same code can be used in a Python script or interactive
environment as follows:

.. doctest-skip::

  >>> import jwst
  >>> input = jwst.datamodels.open('input_data.fits')
  >>> result = MyStep.call(input)
      # `result` contains the resulting data
      # which can then be used by further `Steps`'s or
      # other functions.
      #
      # when done, the data can be saved with the `DataModel.save`
      # method
  >>> result.save('my_final_results.fits')


Input and JWST Conventions
==========================

A `Step` gets its input from two sources:

    - Configuration parameters
    - Arguments to the `Step.process` method

The definition and use of parameters is documented in :ref:`writing-a-step`.

When using the `Step.process` arguments, the code must at least expect
strings. When invoked from the command line using `strun`, how many
arguments to expect are the same number of arguments defined by
`Step.process`. Similarly, the arguments themselves are passed to
`Step.process` as strings.

However, to facilitate code development and interactive usage, code
is expected to accept other object types as well.

A `Step`'s primary argument is expected to be either a string containing
the file path to a data file, or a JWST
:class:`~jwst.datamodels.DataModel` object. The method
:meth:`~jwst.stpipe.step.Step.open_model` handles either type of
input, returning a `DataModel` from the specified file or a shallow
copy of the `DataModel` that was originally passed to it. A typical
pattern for handling input arguments is::

  class MyStep(jwst.stpipe.step.Step):

      def process(self, input_argument):

          input_model = self.open_model(input_argument)

          ...

`input_argument` can either be a string containing a path to a data
file, such as `FITS` file, or a `DataModel` directly.

:meth:`~jwst.stpipe.step.Step.open_model` handles `Step`-specific
issues, such ensuring consistency of input directory handling.

If some other file type is to be opened, the lower level method
:meth:`~jwst.stpipe.step.Step.make_input_path` can be used to specify
the input directory location.

Input and Associations
----------------------

Many of the JWST calibration steps and pipelines expect an
:ref:`Association <associations>` file as input. When opened with
:meth:`~jwst.stpipe.step.Step.open_model`, a
:class:`~jwst.datamodels.ModelContainer` is returned. `ModelContainer`
is, among other features, a list-like object where each element is the
`DataModel` of each member of the association. The `meta.asn_table` is
populated with the association data structure, allowing direct access
to the association itself.  The association file, as well as the files
listed in the association file, must be in the input directory.

To read in a list of files, or an association file, as an association,
use the `load_as_level2_asn` or `load_as_level3_asn` methods.

Input Source
------------

All input files, except for references files provided by CRDS,
are expected to be co-resident in the same directory. That directory
is determined by the directory in which the primary input file
resides. For programmatic use, this directory is available in the
`Step.input_dir` attribute.

Output
======

.. _devel_io_when_files_are_created:

When Files are Created
----------------------

Whether a `Step` produces an output file or not is ultimately
determined by the built-in parameter option `save_results`. If
`True`, output files will be created. `save_results` is set under a
number of conditions:

    - Explicitly through a parameter file or as a command-line option.
    - Implicitly when a step is called by `strun`.

Output File Naming
------------------

File names are constructed based on three components: basename,
suffix, and extension::

  basename_suffix.extension

The extension will often be the same as the primary input file. This
will not be the case if the data format of the output needs to be
something different, such as a text table with `.ecsv` extension.

Similarly, the basename will usually be derived from the primary input
file. However, there are some :ref:`caveats <basename_determination>`
discussed below.

Ultimately, the suffix is what `Step`'s use to identify their output.
The most common suffixes are listed in the
:ref:`pipeline_step_suffix_definitions`.

A `Step`'s suffix is defined in a couple of different ways:

    - By the `Step.name` attribute. This is the default.
    - By the `suffix` parameter.
    - Explicitly in the code. Often this is done in ``Pipelines`` where
      a single pipeline creates numerous different output files.

.. _basename_determination:

Basename Determination
``````````````````````

Most often, the output file basename is determined through any of the
following, given from higher precedence to lower:

    - The `--output_file` command-line option.
    - The `output_file` parameter option.
    - Primary input file name.
    - If the output is a `DataModel`, from the `DataModel.meta.filename`.

In all cases, if the originating file name has a known suffix on it,
that suffix is removed and replaced by the `Step`'s own suffix.

In very rare cases, when there is no other source for the basename, a
basename of `step_\<step_name\>` is used.  This can happen when a
`Step` is being programmatically used and only the `save_results`
parameter option is given.

.. _devel_io_substeps_and_output:

Sub-Steps and Output
````````````````````
Normally, the value of a parameter option is completely local to
the `Step`: A `Step`, called from another `Step` or `Pipeline`, can
only access its own parameters. Hence, options such as
`save_results` do not affect a called `Step`.

The exceptions to this are the parameters `output_file` and
`output_dir`. If either of these parameters are queried by a `Step`,
but are not defined for that `Step`, values will be retrieved up
through the parent. The reason is to provide consistency in output
from `Step` and `Pipelines`. All file names will have the same
basename and will all appear in the same directory.

As expected, if either parameter is specified for the `Step` in
question, the local value will override the parent value.

Also, for `output_file`, there is another option,
`search_output_file`, that can also control this behavior. If set to
`False`, a `Step` will never query its parent for its value.

Basenames, Associations, and Stage 3 Pipelines
``````````````````````````````````````````````

Stage 3 pipelines, such as :ref:`calwebb_image3<calwebb_image3>`
or :ref:`calwebb_spec3<calwebb_spec3>`, take associations
as their primary input. In general, the association defines what the
output basename should be. A typical pattern used to handle
associations is::

  class MyStep(jwst.stpipe.step.Step):

      spec = ''  # Desired configuration parameters

      def process(self, input):

          with jwst.datamodels.open(input) as input_model:

              # If not already specified, retrieve the output
              # file name from the association.
              if self.save_results and self.output_file is None:
                  try:
                     self.output_file = input_model.meta.asn_table.products[0].name

                  except AttributeError:
                      pass

              # Do awesome processing with final result
              # in `result`
              result = final_calculation(input_model)

          return (result)

Some pipelines, such as :ref:`calwebb_spec3 <calwebb_spec3>`, call steps which
are supposed to save their results, but whose basenames should not be based on
the association product name. An example is the
`~jwst.outlier_detection.OutlierDetectionStep` step. For such steps, one can
prevent using the `Pipeline.output_file` specification by setting the parameter
`search_output_file=False`. When such steps then save their output, they will go
through the standard basename search. If nothing else is specified, the basename
will be based on `DataModel.meta.filename` that step's result, creating
appropriate names for that step.

Output API: When More Control Is Needed
---------------------------------------

In summary, the standard output API, as described so far, is basically "set a
few parameters, and let the `Step` framework handle the rest". However, there
are always the exceptions that require finer control, such as saving
intermediate files or multiple files of different formats. This section
discusses the method API and conventions to use in these situations.

Save That Model: Step.save_model
````````````````````````````````

If a `Step` needs to save a `DataModel` before the step completes, use
of :meth:`Step.save_model <jwst.stpipe.step.Step.save_model>` is the recommended over
directly calling :meth:`DataModel.save <jwst.datamodels.DataModel.save>`.
`Step.save_model` uses the `Step` framework and hence will honor the
following:

    - If `Step.save_results` is `False`, nothing will happen.
    - Will ensure that `Step.output_dir` is used.
    - Will use `Step.suffix` if not otherwise specified.
    - Will determine the output basename through the `Step`
      framework, if not otherwise specified.

The basic usage, in which nothing is overridden, is::

  class MyStep(Step):

      def process(self, input):
          ...
          result = some_DataModel
          self.save_model(result)

The most common use case, however, is for saving some intermediate
results that would have a different suffix::

  self.save_model(intermediate_result_datamodel, suffix='intermediate')

See :meth:`jwst.stpipe.step.Step.save_model` for further information.

Make That Filename: Step.make_output_path
`````````````````````````````````````````

For the situations when a filename is needed to be constructed before
saving, either to know what the filename will be or for data that is
not a `DataModel`, use :meth:`Step.make_output_path
<jwst.stpipe.step.Step.make_output_path>`. By default, calling
`make_output_path` without any arguments will return what the default
output file name will be::

  output_path = self.make_output_path()

This method encapsulates the following `Step` framework functions:

    - Will ensure that `Step.output_dir` is used.
    - Will use `Step.suffix` if not otherwise specified.
    - Will determine the output basename through the `Step`
      framework, if not otherwise specified.

A typical use case is when a `Step` needs to save data that is not a
`DataModel`. The current `Step` architecture does not know how to
handle these, so saving needs to be done explicitly. The pattern of
usage would be::

  # A table need be saved and needs a different
  # suffix than what the Step defines.
  table = some_astropy_table_data
  if self.save_results:
      table_path = self.make_output_path(suffix='cat', ext='ecsv')
      table.save(table_path, format='ascii.ecsv', overwrite=True)
