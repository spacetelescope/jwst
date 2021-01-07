.. _config_asdf_files:

ASDF Configuration Files
========================

The format of choice to use for step configuration files. `ASDF <https://asdf-standard.readthedocs.io/>`_ stands for "Advanced Scientific Data Format", a general purpose, non-proprietary, and system-agnostic format for the dissemination of data. Built on `YAML <https://yaml.org/>`_, the most basic file is text-based requiring minimal formatting.

ASDF replaces the original :ref:`CFG <config_cfg_files>` format for step
configuration. Using ASDF allows the configurations to be stored and retrieved
from CRDS, selecting the best configuration for a given set of criteria, such as
instrument and observation mode.

.. _asdf_minimal_file:

To create a configuration file, the most direct way is to choose the Pipeline class, Step class, or already existing .asdf or .cfg file, and run that step using the ``--save-parameters`` option. For example, to get the parameters for the ``Spec2Pipeline`` pipeline, do the following:
::

   $ strun jwst.pipeline.Spec2Pipeline jw00017001001_01101_00001_nrs1_uncal.fits --save-parameters my_spec2.asdf

Once created and modified as necessary, the file can now be used by ``strun`` to run the step/pipeline with the desired parameters:
::

   $ strun my_spec2.asdf jw00017001001_01101_00001_nrs1_uncal.fits

The remaining sections will describe the file format and contents.

File Contents
-------------

To describe the contents of an ASDF file, the configuration for the step ``CubeBuildStep`` will be used as the example:

.. code-block::

    #ASDF 1.0.0
    #ASDF_STANDARD 1.3.0
    %YAML 1.1
    %TAG ! tag:stsci.edu:asdf/
    --- !core/asdf-1.1.0
    asdf_library: !core/software-1.0.0 {author: Space Telescope Science Institute, homepage: 'http://github.com/spacetelescope/asdf',
      name: asdf, version: 2.4.2}
    history:
      entries:
      - !core/history_entry-1.0.0 {description: Base values, time: !!timestamp '2019-10-29
          21:20:50'}
      extensions:
      - !core/extension_metadata-1.0.0
        extension_class: asdf.extension.BuiltinExtension
        software: {name: asdf, version: 2.4.2}
    meta:
      author: SPECIFY AUTHOR
      date: '2019-10-29T16:02:59.377'
      description: Parameters for calibration step SPECIFY
      filename: pars-cubebuildstep.asdf
      instrument: {name: SPECIFY}
      model_type: StepParsModel
      origin: STScI
      pedigree: SPECIFY PEDIGREE
      reftype: pars-cubebuildstep
      telescope: JWST
      useafter: SPECIFY
    parameters: {band: all, channel: all, class: jwst.cube_build.cube_build_step.CubeBuildStep,
      coord_system: world, filter: all, grating: all, name: CubeBuildStep, output_type: band,
      output_use_model: true, rois: 0.0, roiw: 0.0, scale1: 0.0, scale2: 0.0, scalew: 0.0,
      search_output_file: false, single: false, skip_dqflagging: false, weight_power: 2.0,
      weighting: msm}
    ...

Required Components
~~~~~~~~~~~~~~~~~~~

Preamble
++++++++

The first 5 lines, up to and including the "---" line, define the file as an
ASDF file. The rest of the file is formatted as one would format YAML data.
Being YAML, the last line, containing the three ``...`` is essential.

Parameters
++++++++++

A step configuration requires one key, called ``parameters``, which
contains all the parameters to pass onto the step.

The only key required in the ``parameters`` block is ``class``. This defines
which ``Step`` or ``Pipeline`` the parameters belong to. This key is used by
``strun`` to determine which class to actually execute when given a
configuration file.

Another key that will often be seen is ``name``. This defines an alias to use
for the class referenced by the ``class`` key. Pipelines use this alias to refer
to their sub-steps.

All other keys are the parameters and their values to be used when the
step/pipeline is run. The order of the parameters does not matter. Except for
the ``class`` key, no other parameter needs to be specified. If not defined, the
default, as defined in the code, will be used.

Formatting
**********

YAML has two ways of formatting a list of key/value pairs. In the above example, the formatting is very similar to how a Python ``dict`` would be defined. The other way is by simply splitting out all the key/value pairs on separate lines. For example, the ``parameters`` block above could also have been formatted as:

.. code-block::

    parameters:
      band: all
      channel: all
      class: jwst.cube_build.cube_build_step.CubeBuildStep
      coord_system: world
      filter: all
      grating: all
      name: CubeBuildStep
      output_type: band
      output_use_model: true
      rois: 0.0
      roiw: 0.0
      scale1: 0.0
      scale2: 0.0
      scalew: 0.0
      search_output_file: false
      single: false
      skip_dqflagging: false
      weight_power: 2.0
      weighting: msm

Optional Components
~~~~~~~~~~~~~~~~~~~

The ``meta`` and ``history`` blocks are necessary only when the configuration
file is to be used as a parameter reference file in CRDS. See `Configuration as
Reference File`_ below.

Completeness
~~~~~~~~~~~~

For any configuration file, it is not necessary to specify all step/pipeline
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
      class: jwst.cube_build.cube_build_step.CubeBuildStep
      weight_power: 4.0


Pipeline Configuration
~~~~~~~~~~~~~~~~~~~~~~

Pipelines are essentially steps that refer to sub-steps. As in the original cfg
format, parameters for sub-steps can also be specified. All sub-step parameters
appear in a key called `steps`. Sub-step parameters are specified by using the
sub-step name as the key, then underneath and indented, the parameters to change
for that sub-step. For example, to define the ``weight_power`` of the
``cube_build`` step in a ``Spec2Pipeline`` configuration file, the parameter
block would look as follows:

.. code-block::

   parameters:
       class: jwst.pipeline.Spec2Pipeline
       name: calwebb_spec2
       steps:
           cube_build:
               weight_power: 4.0

As with step configuration files, not all sub-steps need to be specified. If
left unspecified, the sub-steps will be run with their default parameter sets.
For the example above, the other steps of ``Spec2Pipeline``, such as
``assign_wcs`` and ``photom`` would still be executed.

Similarly, to skip a particular step, one would specify ``skip: true`` for that
substep. Continuing from the above example, to skip the ``msa_flagging`` step,
the configuration file would look like:

.. code-block::

   parameters:
       class: jwst.pipeline.Spec2Pipeline
       name: calwebb_spec2
       steps:
           msa_flagging:
               skip: true
           cube_build:
               weight_power: 4.0

Python API
----------

Configuration files can be created and modified through the use of the `asdf` package.
The example below demonstrates how to create, from scratch, the ``Spec2Pipeline`` configuration
file shown in the above example:

.. code-block::

   import asdf

   parameters = {'class': 'jwst.pipeline.Spec2Pipeline',
                 'name': 'calwebb_spec2', 'steps': {
                 'msa_flagging': {'skip': True},
                 'cube_build': {'weight_power': 4.0}}
                }

   cfg = asdf.AsdfFile({'parameters': parameters})
   cfg.write_to('my_spec2_weight_40.asdf')

The following example show modifying the ``weight_power`` value to ``8.0``:

.. code-block::

   import asdf

   cfg = asdf.open('my_spec2_weight_40.asdf')

   cfg['parameters']['steps']['cube_build']['weight_power'] = 8.0

   cfg.write_to('my_spec2_weight_80.asdf')
              

Configuration as Reference File
-------------------------------

META Block
~~~~~~~~~~

When a configuration file is to be ingested into CRDS, there is another key
required, ``meta``, which defines the information needed by CRDS to select a
configuration file. A basic reference configuration will look as follows:

.. code-block::

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
      instrument: {name: GENERIC}
      pedigree: GROUND
      reftype: pars-makeliststep
      telescope: JWST
      title: MakeListStep default parameters
      useafter: '1990-04-24T00:00:00'
   parameters:
      class: jwst.stpipe.tests.steps.MakeListStep
   ...

All of the keys under ``meta`` are required, most of which are
self-explanatory. For more information, refer to the `CRDS documentation
<https://jwst-crds.stsci.edu/static/users_guide/>`_.

The one keyword to explain further is ``reftype``. This is what CRDS uses to
determine which reference file is being sought after. For step configurations,
this has the format ``pars-<step_name>`` where ``<step_name>`` will be the Python
class name, in lowercase.


History
~~~~~~~

Parameter reference files also require at least one history entry. This can be found in the ``history`` block under ``entries``:

.. code-block::

    history:
      entries:
      - !core/history_entry-1.0.0 {description: Base values, time: !!timestamp '2019-10-29
          21:20:50'}

It is highly suggested to use the Python API to add history entries:

.. doctest-skip::

   >>> import asdf
   >>> cfg = asdf.open('config.asdf')
       #
       # Modify cfg['parameters'] as necessary
       #
   >>> cfg.add_history_entry('Parameters modified for some reason')
   >>> cfg.write_to('config_modified.asdf')
