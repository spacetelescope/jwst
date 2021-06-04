JWST Calibration Pipeline
=========================

|Build Status| |codecov| |Documentation Status| |Powered by STScI Badge|
|Powered by Astropy Badge|

.. |Build Status| image:: https://github.com/spacetelescope/jwst/workflows/CI/badge.svg?branch=master
   :target: https://github.com/spacetelescope/jwst/actions
.. |codecov| image:: https://codecov.io/gh/spacetelescope/jwst/branch/master/graph/badge.svg?token=Utf5Zs9g7z
   :target: https://codecov.io/gh/spacetelescope/jwst
.. |Documentation Status| image:: https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest
   :target: http://jwst-pipeline.readthedocs.io/en/latest/?badge=latest
.. |Powered by STScI Badge| image:: https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat
   :target: http://www.stsci.edu
.. |Powered by Astropy Badge| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
   :target: http://www.astropy.org/


.. image:: docs/_static/stsci_logo.png
   :alt: STScI Logo

**JWST requires Python 3.7 or above and a C compiler for dependencies.**

**Linux and MacOS platforms are tested and supported. Windows is not
currently supported.**

Installation
------------

The easiest way to install the latest ``jwst`` release into a fresh
virtualenv or conda environment is

::

   pip install jwst


Calibration References Data System (CRDS) Setup
-----------------------------------------------

CRDS is the system that manages the reference files needed to run the
pipeline. Inside the STScI network, the pipeline works with default CRDS
setup with no modifications. To run the pipeline outside the STScI
network, CRDS must be configured by setting two environment variables:

::

   export CRDS_PATH=$HOME/crds_cache
   export CRDS_SERVER_URL=https://jwst-crds.stsci.edu


Contributions and Feedback
--------------------------

We welcome contributions and feedback on the project. Please follow the
`contributing guidelines`_ to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by
abiding with the `Code of Conduct`_.

If you have questions or concerns regarding the software, please open an
issue at https://github.com/spacetelescope/jwst/issues or contact the
`JWST Help Desk`_.

.. _contributing guidelines: CONTRIBUTING.md
.. _Code of Conduct: CODE_OF_CONDUCT.md
.. _JWST Help Desk: https://jwsthelp.stsci.edu


