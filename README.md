JWST Calibration Pipeline
=========================
[![Documentation Status](https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest)](http://jwst-pipeline.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/spacetelescope/jwst.svg?branch=master)](https://travis-ci.org/spacetelescope/jwst)
[![codecov](https://codecov.io/gh/spacetelescope/jwst/branch/master/graph/badge.svg)](https://codecov.io/gh/spacetelescope/jwst)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

![STScI Logo](docs/_static/stsci_logo.png)

**JWST requires Python 3.5 or above and a C compiler for dependencies**.

Installation
------------

The ``jwst`` package can be installed into a virtualenv or conda environment via ``pip``.  We recommend creating a fresh environment with only python installed.  Via conda:

    conda create -n jwst_env python
    conda activate jwst_env

### Installing for end-users ###

To install a released (tagged) version, you can install directly from Github.  To install tagged release ``jwst 0.13.8``:

    pip install numpy
    pip install git+https://github.com/spacetelescope/jwst@0.13.8

The latest development version (from ``master``) can also be installed from Github:

    pip install git+https://github.com/spacetelescope/jwst

As can a particular commit hash:

    pip install git+https://github.com/spacetelescope/jwst@3f03323c

### Installing a DMS release ###

We still package our releases to DMS via a conda spec file that lists the exact versions of all packages to be installed.
To create a new environment with a specific release, use:

    conda create -n jwst_env --file <URL>
    conda activate jwst_env

where `<URL>` is of the form:

    Linux: http://ssb.stsci.edu/releases/jwstdp/0.13.7/latest-linux
    OS X: http://ssb.stsci.edu/releases/jwstdp/0.13.7/latest-osx

Other particular versions can be installed by choosing a different version tag in place of "0.13.7" in the URL path.
See the "Software vs DMS build version map" table below for a list of tags corresponding to particular releases.

### Installing for developers ###

Fork and clone the repo:

    git clone https://github.com/spacetelescope/jwst
    cd jwst

*Note: `setup.py install` and `setup.py develop` commands no longer works! Use `pip install` as stated below.*

Install from your local checked out copy as an "editable" install:

    pip install numpy
    pip install -e .

If you want to run the tests and/or build the docs, you can make sure those dependencies are installed too:

    pip install -e .[test]
    pip install -e .[docs]
    pip install -e .[test,docs]

Need other useful packages in your development environment?

    pip install ipython flake8

### CRDS Setup ###

Inside the STScI network, the pipeline works with default CRDS setup with no modifications.  To run the pipeline outside the STScI network, CRDS must be configured by setting two environment variables:

    export CRDS_PATH=$HOME/crds_cache
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu


Documentation
-------------

Documentation (built daily from `master`) is available at:

https://jwst-pipeline.readthedocs.io/en/latest/

To build the docs yourself, clone this repository and build the documentation with:

    pip install -e .[docs]
    cd docs
    make html
    make latexpdf


Contributions and Feedback
--------------------------
We welcome contributions and feedback on the project. Please follow the [contributing guidelines](CONTRIBUTING.md) to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by abiding to the [Code of Conduct](CODE_OF_CONDUCT.md).

If you have questions or concerns regarding the software, please open an issue
at https://github.com/spacetelescope/jwst/issues or
contact the [JWST Help Desk](https://jwsthelp.stsci.edu).

Software vs DMS build version map
---------------------------------

| jwst tag | DMS build | CRDS_CONTEXT |   Date     |          Notes                           |
| -------- | --------- | ------------ | ---------- | -----------------------------------------|
|  0.13.7  | B7.3      | 0535         | 06/21/2019 | Final release candidate for Build 7.3    |
|  0.13.6  | B7.3rc4   | 0534         | 06/20/2019 | Fourth release candidate for Build 7.3   |
|  0.13.5  | B7.3rc3   | 0534         | 06/19/2019 | Third release candidate for Build 7.3    |
|  0.13.4  | B7.3rc2   | 0534         | 06/18/2019 | Second release candidate for Build 7.3   |
|  0.13.3  | B7.3rc1   | 0532         | 06/04/2019 | First release candidate for Build 7.3    |
|  0.13.2  |           | 0500*        | 05/14/2019 | DMS test, no delivery to I&T             |
|  0.13.1  |           | 0500*        | 03/08/2019 | DMS test, no delivery to I&T             |
|  0.13.0  |           | 0500         | 02/15/2019 | DMS test, no delivery to I&T             |
|  0.12.3  | B7.2.1    | 0500         | 01/15/2019 | DMS Build 7.2.1 patch release            |
|  0.12.2  | B7.2      | 0495         | 11/07/2018 | Final release candidate for Build 7.2    |
|  0.12.1  | B7.2rc2   | 0495         | 11/01/2018 | Second release candidate for Build 7.2   |
|  0.12.0  | B7.2rc1   | 0493*        | 10/09/2018 | First release candidate for Build 7.2    |
|  0.11.0  |           | 0482*        | 09/10/2018 | DMS test, no delivery to I&T             |
|  0.10.0  |           | 0477*        | 07/31/2018 | DMS test, no delivery to I&T             |
|  0.9.6   | B7.1.3    | 0468         | 06/08/2018 | Final release candidate for Build 7.1.3  |
|  0.9.5   | B7.1.3rc3 | 0468         | 06/06/2018 | Third release candidate for Build 7.1.3  |
|  0.9.4   | B7.1.3rc2 | 0463*        | 05/29/2018 | Second release candidate for Build 7.1.3 |
|  0.9.3   | B7.1.3rc1 | 0457*        | 05/11/2018 | First release candidate for Build 7.1.3  |
|  0.9.2   |           | 0441*        | 03/28/2018 | DMS test, no delivery to I&T             |
|  0.9.1   |           | 0432*        | 02/16/2018 | DMS test, no delivery to I&T             |
|  0.9.0   | B7.1.2    | 0422         | 12/22/2017 | DMS patch release to I&T 02/15/2018      |
|  0.8.0   | B7.1.1    | 0422         |            | DMS patch release to I&T 01/17/2018      |
|  0.8.0   | B7.1      | 0422         | 11/14/2017 | Final, delivered to I&T 11/17/2017       |
|  0.7.0rc7| B7.0      | 0303         | 12/13/2016 | Final, delivered to I&T                  |

Note: CRDS_CONTEXT values flagged with an asterisk in the above table are estimates
(formal CONTEXT deliveries are only provided with final builds).

Unit Tests
----------

Unit tests can be run via `pytest`.  Within the top level of your local `jwst` repo checkout:

    pip install -e .[test]
    pytest

Need to parallelize your test runs over 8 cores?

    pip install pytest-xdist
    pytest -n 8


Regression Tests
----------------

Latest regression test results can be found here (STScI staff only):

https://plwishmaster.stsci.edu:8081/job/RT/job/JWST/

The test builds start at 6pm local Baltimore time Monday through Saturday on `jwcalibdev`.

To run the regression tests on your local machine, get the test dependencies and set the environment variable TEST_BIGDATA to our Artifactory server (STSci staff members only)

    pip install -e .[test]
    export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory

When you run the tests, control where the test results are written with the `--basetemp` arg to `pytest`.  So to run all the regression tests:

    pytest --bigdata --basetemp=<PATH> jwst/tests_nightly

If you would like to run a specific test, find its name or ID and use the `-k` option:

    pytest --bigdata --basetemp=<PATH> jwst/tests_nightly -k image3_pipeline

If developers need to update the truth files in our nightly regression tests, there are instructions in the repository wiki.

https://github.com/spacetelescope/jwst/wiki/Maintaining-Regression-Tests

JupyterHub Access
-----------------

**Note:** This is currently still in research-and-development stage and is subject to change.

To run a pre-installed pipeline in JupyterHub:

* Click on https://dev.science.stsci.edu/hub/spawn?image=793754315137.dkr.ecr.us-east-1.amazonaws.com/datb-tc-pipeline-nb:jwstdp-snapshot and sign in.
* Click "Terminal" to:
    * Run `pip freeze` to see what is installed.
    * Grab your notebooks (e.g., using `git clone`) and install any optional software (e.g., using `pip install`).
* Launch your notebook to run the JWST pipeline.

Latest release of any packages is not guaranteed in this environment. Amazon Web Services charges may apply.
