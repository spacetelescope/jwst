.. _config_asdf_files:

ASDF Configuration Files
========================

The format of choice to use to configure steps.

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
