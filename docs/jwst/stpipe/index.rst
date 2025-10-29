======
STPIPE
======

.. _stpipe-user-steps:

For Users
=========
.. toctree::
   :maxdepth: 2

   user_step.rst
   user_pipeline.rst
   config_asdf.rst
   config_cfg.rst
   call_via_call.rst
   call_via_run.rst
   parameter_files.rst
   cfg_deprecation.rst


.. _stpipe-devel-steps:

For Developers
==============
.. toctree::
   :maxdepth: 2

   devel_step.rst
   devel_pipeline.rst
   devel_logging.rst
   devel_io_design.rst

.. automodapi:: jwst.stpipe
   :no-inheritance-diagram:

Base Classes
------------

.. autoclass:: jwst.stpipe.core::JwstStep
   :members: load_as_level2_asn, load_as_level3_asn

.. autoclass:: jwst.stpipe.core::JwstPipeline
