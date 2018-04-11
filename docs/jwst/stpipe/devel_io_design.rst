===============
Step I/O Design
===============

The `Step` architecture is designed such that a `Step`'s intended sole
responsibility is to perform the calculation required. Any
input/output operations are handled by the surrounding `Step`
architecture. This is to help facilitate the use of `Step`'s from both
a command-line environment, and from an interactive Python
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
environment as follows::

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

    - Conifguration parameters
    - Arguments to the `Step.process` method

The definition and use of the configuration parameters is
documented in :ref:`writing-a-step`.

When using the `Step.process` arguments, the code must at least expect
strings. When invoked from the command line using `strun`, how many
arguments to expect are the same number of arguments defined by
`Step.process`. Similarly, the arguments themselves are passed to
`Step.process` as strings.

However, to facilitate code development and interactive usage, code
is expected to accept other object types as well.

For JWST-produced code, nearly all `Step`'s primary argument is
expected to be either a string containing the file path to a data
file, or a JWST :ref:`jwst.datamodels.DataModel` object. The utility
function :ref:`jwst.datamodels.open` handles either type of input,
returning a `DataModel` from the specified file or a shallow copy of
the `DataModel` that was originally passed to it. The code to
accomplish this looks like::

  class MyStep(jwst.stpipe.step.Step):

      def process(self, input_argument):

          input_model = jwst.datamodels.open(input_argument)

          ...

`input_argument` can either be a string containing a path to a data
file, such as `FITS` file, or a `DataModel` directly.

Input and Associations
----------------------

Many of the JWST calibration steps and pipelines expect an
:ref:`Association` file as input. When opened with
:ref:`jwst.datamodels.open`, a :ref:`jwst.datamodels.ModelContainer`
is returned. `ModelContainer` is, among other features, a list-like
object where each element is the `DataModel` of each member of the
association. The `meta.asn_table` is populated with the association
data structure, allowing direct access to the association itself.

Output
======

When Files are Created
----------------------

Whether a `Step` produces an output file or not is ultimately
determined by the built-in configuration option `save_results`. If
`True`, output files will be created. `save_results` is set under a
number of conditions:

    - Explicitly through the `cfg` file or as a command-line option.
    - Implicitly when a step is called by `strun`.
    - Implicitly when the configuration option `output_file` is given
      a value.

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
    - By the `suffix` configuration parameter.
    - Explicitly in the code. Often this is done in `Pipeline`s where
      a single pipeline creates numerous different output files.

.. _basename_determination:

Basename Determination
``````````````````````

Most often, the output file basename is determined from three sources:

    - Primary input file name.
    - The `--output_file` command-line option.
    - The `output_file` configuration option.

In all cases, if the originating file name has a known suffix on it,
that suffix is removed and replaced by the `Step`'s own suffix.

In very rare cases, when there is no other source for the basename, a
basename of `step_\<step_name\>` is used.  This can happen when a
`Step` is being programmatically used and only the `save_results`
configuration option is given.

Basenames, Associations, and Stage 3 Pipelines
```````````````````````````````````````````````
