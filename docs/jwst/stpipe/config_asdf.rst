.. _config_asdf_files:

ASDF Configuration Files
========================

The format of choice to use for step configuration files. `ASDF <https://asdf-standard.readthedocs.io/>`_ stands for "Advanced Scientific Data Format", a general purpose, non-proprietary, and system-agnostic format for the dissemination of data. Built on `YAML <https://yaml.org/>`_, the most basic file is text-based requiring minimal formatting.

ASDF replaces the original :ref:`CFG <config_cfg_files>` format for step
configuration. Using ASDF allows the configurations to be stored and retrieved
from CRDS, selecting the best configuration for a given set of criteria, such as
instrument and observation mode.

.. _asdf_minimal_file:

Minimal File
------------

All configuration files must have at least the following:

.. code-block::

   #ASDF 1.0.0
   #ASDF_STANDARD 1.3.0
   %YAML 1.1
   %TAG ! tag:stsci.edu:asdf/
   --- !core/asdf-1.1.0
   parameters:
     class: jwst.stpipe.tests.steps.MakeListStep
   ...

The first 5 lines define the file as an ASDF file. The rest of the file is
formatted as one would format YAML data. Being YAML, the last line, containing
the three ``...`` is essential.

A step configuration requires one key, called ``parameters``, which
contains all the parameters to pass onto the step. All parameters should be
indented. The amount of indentation does not matter, as long as they are all
indented equally.

See :ref:`Running a Step from a configuration file<running_a_step_from_a_configuration_file>` for a full discussion of required and optional parameters that will be found in the ``parameters`` block.

Configuration as Reference File
-------------------------------

When a configuration file is to be ingested into CRDS, there is another key
required, ``meta``, which defines the information needed by CRDS to select a
configuration file. A basic reference configuration will look as follows:

.. code-block::

   #ASDF 1.0.0
   #ASDF_STANDARD 1.3.0
   %YAML 1.1
   %TAG ! tag:stsci.edu:asdf/
   --- !core/asdf-1.1.0
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
