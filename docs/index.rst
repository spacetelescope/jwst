.. jwst documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:ref:`genindex`  |  :ref:`modindex`

============
Installation
============


Stable releases of the ``jwst`` package are registered at
`PyPI <https://pypi.org/project/jwst/>`_. The latest released version can be
installed into a fresh virtualenv or conda environment using pip:

::

   pip install jwst

Installation details (via conda)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``jwst`` package should be installed into a virtualenv or conda
environment via ``pip``. We recommend that for each installation you
start by creating a fresh environment that only has Python installed and
then install the ``jwst`` package into that bare environment.

If using conda environments, first make sure you have a
recent version of Anaconda or Miniconda installed.

Installation is generally a 3-step process:

-  Create a conda environment
-  Activate that environment
-  Install the ``jwst`` package into that environment

In a bash-compatible shell:

::

   conda create -n <env_name> python
   conda activate <env_name>
   pip install jwst

For more detailed instructions and alternate installation methods see
`our Github README <https://github.com/spacetelescope/jwst>`_.


===========
User Manual
===========

.. toctree::
   :maxdepth: 2

   jwst/introduction.rst


===========================
Data products Documentation
===========================

.. toctree::
   :maxdepth: 1

   jwst/data_products/index.rst


===============================
Error Propagation Documentation
===============================

.. toctree::
   :maxdepth: 1

   jwst/error_propagation/index.rst


=====================
Package Documentation
=====================

.. toctree::
   :maxdepth: 1

   jwst/package_index.rst


