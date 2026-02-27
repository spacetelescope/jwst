.. jwst documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:ref:`genindex`  |  :ref:`modindex`

.. image:: _static/webb_logo.png
   :width: 400
   :align: center
   :target: https://stsci.edu


**Version**: |release|

This package processes uncalibrated data from both imagers and spectrographs onboard the `James Webb Space Telescope (JWST) <https://science.nasa.gov/mission/webb/>`_, an orbiting infrared observatory stationed at Earth-Sun L :subscript:`2`.
The pipeline performs a series of calibration steps that result in standard data products,
applying various corrections to produce science-ready, calibrated output products including
individual exposures and high-level data products (mosaics, extracted spectra, etc.).

This package allows users to run and configure the calibration pipeline themselves for custom processing of JWST data, 
either :ref:`from the command line <run_from_strun>` with ``strun``
or :ref:`from a Python script via the public API <run_from_python>`.
Additionally, this package provides :ref:`JWST datamodel classes <jwst-data-models>`, the recommended method for reading and writing JWST data files in Python.

See `README.md <https://github.com/spacetelescope/jwst>`_ for installation and usage instructions.

.. note::

   If you have trouble installing this package, have encountered a bug while running the pipeline, or wish to request a new feature,
   please `open an issue on GitHub <https://github.com/spacetelescope/jwst/issues>`_ or `contact the JWST Help Desk <https://jwsthelp.stsci.edu>`_.

Detailed explanations of specific calibration stages, reference files, and pipeline builds can be found on `JDox <https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline>`_.

--------------------------------

.. _user-docs:

.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   jwst/user_documentation/background_subtraction_methods/index.rst
   jwst/user_documentation/introduction.rst
   jwst/user_documentation/reference_files_crds.rst
   jwst/user_documentation/parameters.rst
   jwst/user_documentation/running_pipeline_python.rst
   jwst/user_documentation/running_pipeline_command_line.rst
   jwst/user_documentation/available_pipelines.rst
   jwst/user_documentation/input_output_file_conventions.rst
   jwst/user_documentation/logging.rst
   jwst/user_documentation/datamodels.rst
   jwst/user_documentation/more_information.rst

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
   jwst/changes.rst

============
Contributing
============

``jwst`` is an open source package written in Python.
The source code is `available on GitHub <https://github.com/spacetelescope/jwst>`_.
New contributions and contributors are very welcome!

Please read `CONTRIBUTING.md <https://github.com/spacetelescope/jwst/blob/main/CONTRIBUTING.md>`_,
the :ref:`public API definition <jwst-public-vs-private-api>`,
and the :ref:`public API deprecation policy <jwst-deprecation-policy>`.

We strive to provide a welcoming community by abiding with our `CODE_OF_CONDUCT.md <https://github.com/spacetelescope/jwst/blob/main/CODE_OF_CONDUCT.md>`_.

