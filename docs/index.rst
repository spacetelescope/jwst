.. jwst documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:ref:`genindex`  |  :ref:`modindex`

.. image:: _static/webb_logo.png
   :width: 400
   :align: center

.. warning::

   As of November 10, 2022, the process of deprecating the CRDS PUB Server will start.

   For details, refer to the :ref:`pub-deprecation` page.

Welcome to the documentation for `jwst`. This package contains the Python
software suite for the James Webb Space Telescope (JWST) calibration pipeline,
which processes data from all JWST instruments by applying various corrections to
produce science-ready, calibrated output products including fully calibrated
individual exposures as well as high-level data products (mosaics, extracted
spectra, etc.). The tools in this package allow users to run and configure the
pipeline to custom process their JWST data. Additionally, the `jwst` package
contains the interface to JWST datamodels, the recommended method of reading and
writing JWST data files in Python.


--------------------------------

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting_started/quickstart.rst
   getting_started/install.rst
   getting_started/contributing.rst

.. _user-docs:

.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   jwst/user_documentation/introduction.rst
   jwst/user_documentation/reference_files_crds.rst
   jwst/user_documentation/pub_deprecation.rst
   jwst/user_documentation/parameters.rst
   jwst/user_documentation/running_pipeline_python.rst
   jwst/user_documentation/running_pipeline_command_line.rst
   jwst/user_documentation/available_pipelines.rst
   jwst/user_documentation/input_output_file_conventions.rst
   jwst/user_documentation/logging_configuration.rst
   jwst/user_documentation/datamodels.rst

.. toctree::
   :maxdepth: 2
   :caption: Data Products Documentation

   jwst/data_products/index.rst


.. toctree::
   :maxdepth: 2
   :caption: Error Propagation Documentation

   jwst/error_propagation/index.rst

--------------------------------

=====================
Package Documentation
=====================

.. toctree::
   :maxdepth: 2
   :caption: Package Documentation

   jwst/package_index.rst
