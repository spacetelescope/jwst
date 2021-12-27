.. _config_asdf_files:

ASDF Parameter Files
====================

ASDF is the format of choice for parameter files. `ASDF
<https://asdf-standard.readthedocs.io/>`_ stands for "Advanced Scientific Data
Format", a general purpose, non-proprietary, and system-agnostic format for the
dissemination of data. Built on `YAML <https://yaml.org/>`_, the most basic file
is text-based requiring minimal formatting.

ASDF replaces the original :ref:`CFG <config_cfg_files>` format for step
configuration. Using ASDF allows the configurations to be stored and retrieved
from CRDS, selecting the best parameter file for a given set of criteria, such
as instrument and observation mode.

.. _asdf_minimal_file:

To create a parameter file, the most direct way is to choose the Pipeline
class, Step class, or already existing .asdf or .cfg file, and run that step
using the ``--save-parameters`` option. For example, to get the parameters for
the ``Spec2Pipeline`` pipeline, do the following: ::

   $ strun jwst.pipeline.Spec2Pipeline jw00017001001_01101_00001_nrs1_uncal.fits --save-parameters my_spec2.asdf

Once created and modified as necessary, the file can now be used by ``strun`` to run the step/pipeline with the desired parameters:
::

   $ strun my_spec2.asdf jw00017001001_01101_00001_nrs1_uncal.fits

The remaining sections will describe the file format and contents.

File Contents
-------------

To describe the contents of an ASDF file, the configuration for the step
``CubeBuildStep`` will be used as the example:

.. code-block::

   #ASDF 1.0.0
   #ASDF_STANDARD 1.5.0
   %YAML 1.1
   %TAG ! tag:stsci.edu:asdf/
   --- !core/asdf-1.1.0
   asdf_library: !core/software-1.0.0 {author: Space Telescope Science Institute, homepage: 'http://github.com/spacetelescope/asdf',
     name: asdf, version: 2.7.3}
   history:
     extensions:
     - !core/extension_metadata-1.0.0
       extension_class: asdf.extension.BuiltinExtension
       software: !core/software-1.0.0 {name: asdf, version: 2.7.3}
   class: jwst.cube_build.cube_build_step.CubeBuildStep
   name: CubeBuildStep
   parameters:
     band: all
     channel: all
     coord_system: skyalign
     filter: all
     grating: all
     input_dir: ''
     output_ext: .fits
     output_type: band
     output_use_index: true
     output_use_model: true
     post_hooks: []
     pre_hooks: []
     rois: 0.0
     roiw: 0.0
     save_results: false
     scale1: 0.0
     scale2: 0.0
     scalew: 0.0
     search_output_file: false
     single: false
     skip: false
     skip_dqflagging: false
     weight_power: 2.0
     weighting: emsm
     ...

Required Components
~~~~~~~~~~~~~~~~~~~

Preamble
++++++++

The first 5 lines, up to and including the "---" line, define the file as an
ASDF file. The rest of the file is formatted as one would format YAML data.
Being YAML, the last line, containing the three ``...`` is essential.

class and name
++++++++++++++

There are two required keys at the top level: ``class`` and ``parameters``.
``parameters`` is discussed below.

``class`` specifies the Python class to run.  It should be a
fully-qualified Python path to the class.  Step classes can ship with
``stpipe`` itself, they may be part of other Python packages, or they
exist in freestanding modules alongside the configuration file.  For
example, to use the ``SystemCall`` step included with ``stpipe``, set
``class`` to ``stpipe.subprocess.SystemCall``.  To use a class called
``Custom`` defined in a file ``mysteps.py`` in the same directory as
the configuration file, set ``class`` to ``mysteps.Custom``.

``name`` defines the name of the step.  This is distinct from the
class of the step, since the same class of Step may be configured in
different ways, and it is useful to be able to have a way of
distinguishing between them.  For example, when Steps are combined
into :ref:`stpipe-user-pipelines`, a Pipeline may use the same Step class
multiple times, each with different configuration parameters.

Parameters
++++++++++

``parameters`` contains all the parameters to pass onto the step. The order of
the parameters does not matter. It is not necessary to specify all parameters
either. If not defined, the default, as defined in the code or values from CRDS
parameter references, will be used.

Formatting
**********

YAML has two ways of formatting a list of key/value pairs. In the above example,
each key/value pair is on separate line. The other way is using a form that is similar to a Python ``dict``.
For example, the ``parameters`` block above could also have been formatted as:

.. code-block::

    parameters: {band: all, channel: all, coord_system: world, filter: all,
      grating: all, output_type: band, output_use_model: true, rois: 0.0,
      roiw: 0.0, scale1: 0.0, scale2: 0.0, scalew: 0.0, search_output_file: false,
      single: false, skip_dqflagging: false, weight_power: 2.0, weighting: msm}

Optional Components
~~~~~~~~~~~~~~~~~~~

The ``asdf_library`` and ``history`` blocks are necessary only when a parameter
file is to be used as a parameter reference file in CRDS. See `Parameter Files
as Reference Files`_ below.

.. _`Completeness`:

Completeness
~~~~~~~~~~~~

For any parameter file, it is not necessary to specify all step/pipeline
parameters. Any parameter left unspecified will get, at least, the default value
define in the step's code. If a parameter is defined without a default value,
and the parameter is never assigned a value, an error will be produced when the
step is executed.

Remember that parameter values can come from numerous sources. Refer to
:ref:`Parameter Precedence` for a full listing of how parameters can be set.

From the ``CubeBuildStep`` example, if all that needed to change is the
``weight_power`` parameter with a setting of ``4.0``, the ``parameters`` block
need only contain the following:

.. code-block::

    parameters:
      weight_power: 4.0


Pipeline Configuration
~~~~~~~~~~~~~~~~~~~~~~

Pipelines are essentially steps that refer to sub-steps. As in the original cfg
format, parameters for sub-steps can also be specified. All sub-step parameters
appear in a key called `steps`. Sub-step parameters are specified by using the
sub-step name as the key, then underneath and indented, the parameters to change
for that sub-step. For example, to define the ``weight_power`` of the
``cube_build`` step in a ``Spec2Pipeline`` parameter file, the parameter
block would look as follows:

.. code-block::

   class: jwst.pipeline.Spec2Pipeline
   parameters: {}
   steps:
   - class: jwst.cube_build.cube_build_step.CubeBuildStep
     parameters:
       weight_power: 4.0

As with step parameter files, not all sub-steps need to be specified. If left
unspecified, the sub-steps will be run with their default parameter sets. For
the example above, the other steps of ``Spec2Pipeline``, such as ``assign_wcs``
and ``photom`` would still be executed.

Similarly, to skip a particular step, one would specify ``skip: true`` for that
substep. Continuing from the above example, to skip the ``msa_flagging`` step,
the parameter file would look like:

.. code-block::

   class: jwst.pipeline.Spec2Pipeline
   parameters: {}
   steps:
   - class: jwst.msaflagopen.msaflagopen_step.MSAFlagOpenStep
     parameters:
       skip: true
   - class: jwst.cube_build.cube_build_step.CubeBuildStep
     parameters:
       weight_power: 4.0

.. note::

   In the previous examples, one may have noted the line ``parameters: {}``. In
   neither example, and is a common situation when defining pipeline
   configurations, there is no need to set any of the parameters for the
   pipeline itself. However, the keyword ``parameters`` is required. As such,
   the value for ``parameters`` is defined as an empty dictionary, ``{}``.

Python API
----------

There are a number of ways to create an ASDF parameter file. From the
command line utility ``strun``, the option ``--save-parameters`` can be used.

Within a Python script, the method ``Step.export_config(filename: str)`` can be
used. For example, to create a parameter file for ``CubeBuildStep``, use the
following:

.. doctest-skip::

   >>> from jwst.cube_build import CubeBuildStep
   >>> step = CubeBuildStep()
   >>> step.export_config('cube_build.asdf')

Parameter Files as Reference Files
----------------------------------

ASDF-formatted parameter files are the basis for the parameter reference
reftypes in CRDS. There are two more keys that are needed to be added which CRDS
requires: ``meta`` and ``history``.

The direct way of creating a parameter reference file is through the
``Step.export_config`` method, just as one would to get a basic parameter file.
The only addition is the argument ``include_meta=True``. For example, to get a
reference-file ready version of the ``CubeBuildStep``, use the following Python
code:

.. doctest-skip::

   >>> from jwst.cube_build import CubeBuildStep
   >>> step = CubeBuildStep()
   >>> step.export_config('pars-cubebuildstep.asdf', include_meta=True)


The explanations for the ``meta`` and ``history`` blocks are given below.

META Block
~~~~~~~~~~

When a parameter file is to be ingested into CRDS, there is another key
required, ``meta``, which defines the information needed by CRDS parameter file
selection. A basic reference parameter file will look as follows:

.. code-block:: yaml

   #ASDF 1.0.0
   #ASDF_STANDARD 1.3.0
   %YAML 1.1
   %TAG ! tag:stsci.edu:asdf/
   --- !core/asdf-1.1.0
   history:
     entries:
     - !core/history_entry-1.0.0 {description: Base values, time: !!timestamp '2019-10-29
         21:20:50'}
     extensions:
     - !core/extension_metadata-1.0.0
       extension_class: asdf.extension.BuiltinExtension
       software: {name: asdf, version: 2.4.2}
   meta:
      author: Alfred E. Neuman
      date: '2019-07-17T10:56:23.456'
      description: MakeListStep parameters
      instrument: {name: MIRI}
      pedigree: GROUND
      reftype: pars-spec2pipeline
      telescope: JWST
      title: Spec2Pipeline default parameters
      useafter: '1990-04-24T00:00:00'
   class: jwst.pipeline.calwebb_spec2.Spec2Pipeline
   parameters: {}
   ...

All of the keys under ``meta`` are required, most of which are
self-explanatory. For more information, refer to the `CRDS documentation
<https://jwst-crds.stsci.edu/static/users_guide/>`_.

The one keyword to explain further is ``reftype``. This is what CRDS uses to
determine which parameter file is being sought after. This has the format
``pars-<step_name>`` where ``<step_name>`` is the Python class name, in
lowercase.


History
~~~~~~~

Parameter reference files also require at least one history entry. This can be found in the ``history`` block under ``entries``:

.. code-block::

    history:
      entries:
      - !core/history_entry-1.0.0 {description: Base values, time: !!timestamp '2019-10-29
          21:20:50'}

It is highly suggested to use the ASDF API to add history entries:

.. doctest-skip::

   >>> import asdf
   >>> cfg = asdf.open('config.asdf')
       #
       # Modify `parameters` and `meta` as necessary.
       #
   >>> cfg.add_history_entry('Parameters modified for some reason')
   >>> cfg.write_to('config_modified.asdf')

JWST, Parameters and Parameter References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In general, the default parameters for any pipeline or step are valid for nearly
all instruments and observing modes. This means that when a pipeline or step is
run without any explicit parameter setting, that pipeline or step will usually
do the desired operation. Hence, most of the time there is no need for a
parameter reference to be available in CRDS, or provided by the user. Only for a
small set of observing mode/step combinations, will there be need to create a
parameter reference. Even then, nearly all cases will involve changing a subset
of a pipeline or step parameters.

Keeping this sparse-population philosophy in mind, for most parameter
references, only those parameters that are explicitly changed should be
specified in the reference. If adhered to, when a pipeline/step default value
for a particular parameter needs to change, the change will be immediately
available. Otherwise, all references that mistakenly set said parameter will
need to be updated. See :ref:`Completeness` for more information.

Furthermore, every pipeline/step have a common set of parameters, listed
below. These parameters generally affect the infrastructure operation of
pipelines/steps, and should not be included in a parameter reference.

- input_dir
- output_ext
- output_use_index
- output_use_model
- post_hooks
- pre_hooks
- save_results
- search_output_file
