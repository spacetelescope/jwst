.. _quickstart:

================
Quickstart Guide
================

To run the JWST pipeline:

1. ``pip install jwst`` into an environment 
2. activate that environment 
3. set environment variables for accessing CRDS reference files
4. run the JWST pipeline :ref:`from the command line with strun <run_from_strun>` or :ref:`via the Python API <run_from_python>`

Additionally, the :ref:`JWST datamodels <jwst-data-models>` and other pipeline utilities are accessible for import and use in Python.

1. Install the JWST pipeline
============================

See `the installation instructions in README.md <https://github.com/spacetelescope/jwst/blob/main/README.md#installation>`_.

2. Activate the environment
===========================

.. tab:: ``mamba`` / ``conda``
    
    .. code-block:: shell

        mamba activate -n jwst_env

    > [!NOTE]
    > ``mamba`` is the recommended drop-in replacement for the ``conda`` command. If you only have ``conda``, just replace ``mamba`` with ``conda`` in the above commands.

.. tab:: ``virtualenv``
    
    .. code-block:: shell

        source ~/venvs/jwst_env/bin/activate


3. Set CRDS environment variables
=================================

The `Calibration References Data System (CRDS) <https://jwst-crds.stsci.edu/static/users_guide/index.html>`_ serves and manages reference files for several telescope calibration pipelines, including ``jwst``.

Set ``CRDS_SERVER_URL`` and ``CRDS_PATH`` to run the calibration pipeline with access to reference files from CRDS:

.. code-block:: shell

    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
    export CRDS_PATH=$HOME/data/crds_cache/

The pipeline will automatically download individual reference files and cache them in the ``CRDS_PATH`` directory.
Expect to use upwards of 200 gigabytes of disk space for reference files.

.. tip::

    If you are inside the STScI network, physically or via VPN, you do not need to set the ``CRDS_PATH`` environment variable (it defaults to shared network storage).

To use a specific CRDS context other than that `automatically associated with a given pipeline version <https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/crds-migration-to-quarterly-calibration-updates>`_, explicitly set the ``CRDS_CONTEXT`` environment variable:

.. code-block:: shell

    export CRDS_CONTEXT=jwst_1179.pmap

.. warning::
    
    The CRDS PUB Server (``https://jwst-crds-pub.stsci.edu``) is decommissioned as of March 2023.
    To use historical files from the PUB server, contact the [JWST Help Desk](https://jwsthelp.stsci.edu).

For more information, see :ref:`reference_files_crds`.

4. Run the Pipeline
===================

For information on how to run the pipeline using the Python interface, see
:ref:`Running the JWST pipeline: Python Interface <run_from_python>`.

For information on how to run the pipeline using the command line interface, see
:ref:`Running the JWST pipeline: Command Line Interface <run_from_strun>`.

For information on how to read and write data files with JWST ``datamodels``, see
:ref:`JWST datamodels <jwst-data-models>`.
