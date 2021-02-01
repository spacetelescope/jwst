# JWST Calibration Pipeline

[![Build Status](https://github.com/spacetelescope/jwst/workflows/CI/badge.svg?branch=master)](https://github.com/spacetelescope/jwst/actions)
[![codecov](https://codecov.io/gh/spacetelescope/jwst/branch/master/graph/badge.svg?token=Utf5Zs9g7z)](https://codecov.io/gh/spacetelescope/jwst)
[![Documentation Status](https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest)](http://jwst-pipeline.readthedocs.io/en/latest/?badge=latest)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

![STScI Logo](docs/_static/stsci_logo.png)

**JWST requires Python 3.7 or above and a C compiler for dependencies.**

**Linux and MacOS platforms are tested and supported.  Windows is not currently supported.**


## Installation

The easiest way to install the latest `jwst` release into a fresh virtualenv or conda environment is

    pip install jwst

### Detailed Installation

The `jwst` package can be installed into a virtualenv or conda environment via `pip`.
We recommend that for each installation you start by creating a fresh
environment that only has Python installed and then install the `jwst` package and
its dependencies into that bare environment.
If using conda environments, first make sure you have a recent version of Anaconda
or Miniconda installed.
If desired, you can create multiple environments to allow for switching between different
versions of the `jwst` package (e.g. a released version versus the current development version).

In all cases, the installation is generally a 3-step process:
* Create a conda environment
* Activate that environment
* Install the desired version of the `jwst` package into that environment

Details are given below on how to do this for different types of installations,
including tagged releases, DMS builds used in operations, and development versions.
Remember that all conda operations must be done from within a bash shell.


### Installing latest releases

You can install the latest released version via `pip`.  From a bash shell:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install jwst

You can also install a specific version (from `jwst 0.17.0` onward):

    conda create -n <env_name> python
    conda activate <env_name>
    pip install jwst==0.18.3

Installing specific versions before `jwst 0.17.0` need to be installed from Github:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install git+https://github.com/spacetelescope/jwst@0.16.2


### Installing the development version from Github

You can install the latest development version (not as well tested) from the
Github master branch:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install git+https://github.com/spacetelescope/jwst
 

### Installing a DMS Operational Build

There may be occasions where an exact copy of an operational DMS build is
desired (e.g. for validation testing or debugging operational issues).
We package releases for DMS builds via environment snapshots that specify the
exact versions of all packages to be installed. This method may result in more
stable processing than what was outlined above for installing a particular
tagged release, because that method installs the latest versions of dependency
packages, while this method installs dependencies pinned to particular versions
that have been well tested.

To install a particular DMS build, consult the
[Software vs DMS build version map](https://github.com/spacetelescope/jwst#software-vs-dms-build-version-map)
table shown below to determine the correct jwst tag. For example, to install the
version of `jwst` used in DMS build 7.5, use jwst tag 0.16.1. The overall
procedure is similar to the 3-step process outlined in the previous section, but the
details of each command vary, due to the use of environment snapshot files that specify
all of the particular packages to install. Also note that different snapshot files are
used for Linux and Mac OS systems.

Linux:

    conda create -n <env_name> --file https://ssb.stsci.edu/releases/jwstdp/0.16.1/conda_python_stable-deps.txt
    conda activate <env_name>
    pip install -r https://ssb.stsci.edu/releases/jwstdp/0.16.1/reqs_stable-deps.txt

MacOS:

    conda create -n <env_name> --file https://ssb.stsci.edu/releases/jwstdp/0.16.1/conda_python_macos-stable-deps.txt
    conda activate <env_name>
    pip install -r https://ssb.stsci.edu/releases/jwstdp/0.16.1/reqs_macos-stable-deps.txt

Each DMS delivery has its own installation instructions, which may be found in
the corresponding release documentation linked from this page:
https://github.com/astroconda/astroconda-releases/tree/master/jwstdp
The installation procedures may change from time to time, so consulting the
documentation page for the specific version in question is the best way to get
that version installed.


### Installing for Developers

If you want to be able to work on and test the source code with the `jwst` package,
the high-level procedure to do this is to first create a conda environment using
the same procedures outlined above, but then install your personal copy of the
code overtop of the original code in that environment. Again, this should be done
in a separate conda environment from any existing environments that you may have
already installed with released versions of the `jwst` package.

As usual, the first two steps are to create and activate an environment:

    conda create -n <env_name> python
    conda activate <env_name>

To install your own copy of the code into that environment, you first need to
fork and clone the `jwst` repo:

    cd <where you want to put the repo>
    git clone https://github.com/spacetelescope/jwst
    cd jwst

*Note: `python setup.py install` and `python setup.py develop` commands do not work.*

Install from your local checked-out copy as an "editable" install:

    pip install -e .

If you want to run the unit or regression tests and/or build the docs, you can make
sure those dependencies are installed too:

    pip install -e .[test]
    pip install -e .[docs]
    pip install -e .[test,docs]

Need other useful packages in your development environment?

    pip install ipython pytest-xdist


## Calibration References Data System (CRDS) Setup

CRDS is the system that manages the reference files needed to run the pipeline.
Inside the STScI network, the pipeline works with default CRDS setup with no modifications.
To run the pipeline outside the STScI network, CRDS must be configured by setting
two environment variables:

    export CRDS_PATH=$HOME/crds_cache
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu


## Documentation

Documentation (built daily from the Github `master` branch) is available at:

https://jwst-pipeline.readthedocs.io/en/latest/

To build the docs yourself, clone this repository and build the documentation with:

    pip install -e .[docs]
    cd docs
    make html
    make latexpdf


## Contributions and Feedback

We welcome contributions and feedback on the project. Please follow the
[contributing guidelines](CONTRIBUTING.md) to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by abiding with
the [Code of Conduct](CODE_OF_CONDUCT.md).

If you have questions or concerns regarding the software, please open an issue
at https://github.com/spacetelescope/jwst/issues or
contact the [JWST Help Desk](https://jwsthelp.stsci.edu).


## Software vs DMS build version map

| jwst tag | DMS build | CRDS_CONTEXT |   Date     |          Notes                                |
| -------- | --------- | ------------ | ---------- | ----------------------------------------------|
|  0.18.3  | B7.7      | 0670         | 01/25/2021 | Final release candidate for B7.7              |
|  0.18.2  | B7.7rc3   | 0668         | 01/19/2021 | Third release candidate for B7.7              |
|  0.18.1  | B7.7rc2   | 0664         | 01/08/2021 | Second release candidate for B7.7             |
|  0.18.0  | B7.7rc1   | 0645         | 12/21/2020 | First release candidate for B7.7              |
|  0.17.1  | B7.6      | 0641         | 09/15/2020 | Final release candidate for B7.6              |
|  0.17.0  | B7.6rc1   | 0637         | 08/28/2020 | First release candidate for B7.6              |
|  0.16.2  | B7.5      | 0619         | 06/10/2020 | Same as 0.16.1, but with installation bug fix |
|  0.16.1  | B7.5      | 0619         | 05/19/2020 | Final release candidate for B7.5              |
|  0.16.0  | B7.5rc1   | 0614         | 05/04/2020 | First release candidate for B7.5              |
|  0.15.1  | B7.4.2    | 0586         | 03/10/2020 | Final release candidate for B7.4.2            |
|  0.15.0  | B7.4.2rc1 | 0585         | 02/28/2020 | First release candidate for B7.4.2            |
|  0.14.2  | B7.4      | 0570         | 11/18/2019 | Final release candidate for B7.4              |
|  0.14.1  | B7.4rc2   | 0568         | 11/11/2019 | Second release candidate for B7.4             |
|  0.14.0  | B7.4rc1   | 0563         | 10/25/2019 | First release candidate for B7.4              |
|  0.13.8  | B7.3.1    | 0541         | 09/05/2019 | Patch for Build 7.3 released as Build 7.3.1   |
|  0.13.7  | B7.3      | 0535         | 06/21/2019 | Final release candidate for Build 7.3         |
|  0.13.6  | B7.3rc4   | 0534         | 06/20/2019 | Fourth release candidate for Build 7.3        |
|  0.13.5  | B7.3rc3   | 0534         | 06/19/2019 | Third release candidate for Build 7.3         |
|  0.13.4  | B7.3rc2   | 0534         | 06/18/2019 | Second release candidate for Build 7.3        |
|  0.13.3  | B7.3rc1   | 0532         | 06/04/2019 | First release candidate for Build 7.3         |
|  0.13.2  |           | 0500*        | 05/14/2019 | DMS test, no delivery to I&T                  |
|  0.13.1  |           | 0500*        | 03/08/2019 | DMS test, no delivery to I&T                  |
|  0.13.0  |           | 0500         | 02/15/2019 | DMS test, no delivery to I&T                  |
|  0.12.3  | B7.2.1    | 0500         | 01/15/2019 | DMS Build 7.2.1 patch release                 |
|  0.12.2  | B7.2      | 0495         | 11/07/2018 | Final release candidate for Build 7.2         |
|  0.12.1  | B7.2rc2   | 0495         | 11/01/2018 | Second release candidate for Build 7.2        |
|  0.12.0  | B7.2rc1   | 0493*        | 10/09/2018 | First release candidate for Build 7.2         |
|  0.11.0  |           | 0482*        | 09/10/2018 | DMS test, no delivery to I&T                  |
|  0.10.0  |           | 0477*        | 07/31/2018 | DMS test, no delivery to I&T                  |
|  0.9.6   | B7.1.3    | 0468         | 06/08/2018 | Final release candidate for Build 7.1.3       |
|  0.9.5   | B7.1.3rc3 | 0468         | 06/06/2018 | Third release candidate for Build 7.1.3       |
|  0.9.4   | B7.1.3rc2 | 0463*        | 05/29/2018 | Second release candidate for Build 7.1.3      |
|  0.9.3   | B7.1.3rc1 | 0457*        | 05/11/2018 | First release candidate for Build 7.1.3       |
|  0.9.2   |           | 0441*        | 03/28/2018 | DMS test, no delivery to I&T                  |
|  0.9.1   |           | 0432*        | 02/16/2018 | DMS test, no delivery to I&T                  |
|  0.9.0   | B7.1.2    | 0422         | 12/22/2017 | DMS patch release to I&T 02/15/2018           |
|  0.8.0   | B7.1.1    | 0422         |            | DMS patch release to I&T 01/17/2018           |
|  0.8.0   | B7.1      | 0422         | 11/14/2017 | Final, delivered to I&T 11/17/2017            |
|  0.7.0rc7| B7.0      | 0303         | 12/13/2016 | Final, delivered to I&T                       |

Note: CRDS_CONTEXT values flagged with an asterisk in the above table are estimates
(formal CONTEXT deliveries are only provided with final builds).


## Unit Tests

Unit tests can be run via `pytest`.  Within the top level of your local `jwst` repo checkout:

    pip install -e .[test]
    pytest

Need to parallelize your test runs over 8 cores?

    pip install pytest-xdist
    pytest -n 8


## Regression Tests

Latest regression test results can be found here (STScI staff only):

https://plwishmaster.stsci.edu:8081/job/RT/job/JWST/

The test builds start at 6pm local Baltimore time Monday through Saturday on `jwcalibdev`.

To run the regression tests on your local machine, get the test dependencies
and set the environment variable TEST_BIGDATA to our Artifactory server
(STSci staff members only):

    pip install -e .[test]
    export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory

To run all the regression tests:

    pytest --bigdata jwst/regtest

You can control where the test results are written with the
`--basetemp=<PATH>` arg to `pytest`.  _NOTE that `pytest` will wipe this directory clean
for each test session, so make sure it is a scratch area._

If you would like to run a specific test, find its name or ID and use the `-k` option:

    pytest --bigdata jwst/regtest -k nirspec

If developers need to update the truth files in our nightly regression tests,
there are instructions in the repository wiki.

https://github.com/spacetelescope/jwst/wiki/Maintaining-Regression-Tests


## JupyterHub Access

**NOTE:** This is currently still in research-and-development stage and is subject to change.

To run a pre-installed pipeline in JupyterHub:

* Click on https://dev.science.stsci.edu/hub/spawn?image=793754315137.dkr.ecr.us-east-1.amazonaws.com/datb-tc-pipeline-nb:jwstdp-snapshot and sign in.
* Click "Terminal" to:
    * Run `pip freeze` to see what is installed.
    * Grab your notebooks (e.g., using `git clone`) and install any optional software (e.g., using `pip install`).
* Launch your notebook to run the JWST pipeline.

Latest release of any packages is not guaranteed in this environment.
Amazon Web Services charges may apply.
