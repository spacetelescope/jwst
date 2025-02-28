.. _installation:

============
Installation
============

Stable releases of the ``jwst`` package are registered at
`PyPI <https://pypi.org/project/jwst/>`_. The development version of `jwst` is
installable from the
`Github repository <https://github.com/spacetelescope/jwst>`_.

``jwst`` is also available as part of
`stenv <https://stenv.readthedocs.io/en/latest/>`_ (Space Telescope Environment).

Detailed Installation Instructions
==================================

The `jwst` package can be installed into a virtualenv or conda environment via
`pip`. We recommend that for each installation you start by creating a fresh
environment that only has Python installed and then install the `jwst` package
and its dependencies into that bare environment. If using conda environments,
first make sure you have a recent version of Anaconda or Miniconda
`installed <https://docs.conda.io/en/latest/miniconda.html>`_. If desired, you
can create multiple environments to allow for switching between different
versions of the `jwst` package (e.g. a released version versus the current
development version).

In all cases, the installation is generally a 3-step process

#. Create a conda environment
#. Activate that environment
#. Install the desired version of the `jwst` package into that environment

Details are given below on how to do this for different types of installations,
including tagged releases, DMS builds used in operations, and development
versions. Remember that all conda operations must be done from within a bash/zsh
shell.

.. warning::

    JWST requires a C compiler for dependencies.

.. warning::
    Installation on MacOS Mojave 10.14 will fail due to lack of a stable build for dependency ``opencv-python``.

Installing Latest Release
-------------------------

You can install the latest released version via `pip`.  From a bash/zsh shell:

    | >> conda create -n <env_name> python=3.12
    | >> conda activate <env_name>
    | >> pip install jwst

.. _installing_previous_release:

Installing Previous Releases
----------------------------

You can also install a specific version (from `jwst 0.17.0` onward):

    | >> conda create -n <env_name> python=3.12
    | >> conda activate <env_name>
    | >> pip install jwst==1.12.5

.. _installing_dev:

Installing the Development Version from Github
----------------------------------------------

You can install the latest development version (not as well tested) from the
Github main branch:

    | >> conda create -n <env_name> python=3.12
    | >> conda activate <env_name>
    | >> pip install git+https://github.com/spacetelescope/jwst

.. _upgrade_install:

Upgrading Installed Version
---------------------------

.. Important:: Do NOT use `pip install jwst --upgrade` to upgrade your
    installation. This does not check if dependencies are upgraded and will cause
    issues. Instead, use the method detailed below.

If you have previously installed `jwst` and you would like to upgrade to keep your
install up-to-date, we recommend that you first uninstall the package in your
environment of choice and then reinstall:

    | >> pip uninstall jwst
    | >> pip install jwst

This will ensure that all dependency packages are also upgraded. This also
applies when using the development version of jwst - to upgrade and grab recent
changes, uninstall and re-install the main branch from Github:

    | >> pip uninstall jwst
    | >> pip install git+https://github.com/spacetelescope/jwst

Installing with ``stenv``
-------------------------

``jwst`` is available as part of ``stenv``, a set of installable Conda
environments that bundle software for astronomical data analysis with JWST, HST,
and other observatories. See `the stenv documentation <https://stenv.readthedocs.io/en/latest/>`_
for more information.


**For more install instructions, including how to install jwst for development**
**or how to install a DMS operational build, see** `the Github README <https://github.com/spacetelescope/jwst>`_.
