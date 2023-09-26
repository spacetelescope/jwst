# JWST Calibration Pipeline

[![Build Status](https://github.com/spacetelescope/jwst/workflows/CI/badge.svg?branch=master)](https://github.com/spacetelescope/jwst/actions)
[![codecov](https://codecov.io/gh/spacetelescope/jwst/branch/master/graph/badge.svg?token=Utf5Zs9g7z)](https://codecov.io/gh/spacetelescope/jwst)
[![Documentation Status](https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest)](http://jwst-pipeline.readthedocs.io/en/latest/?badge=latest)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![DOI](https://zenodo.org/badge/60551519.svg)](https://zenodo.org/badge/latestdoi/60551519)

![STScI Logo](docs/_static/stsci_logo.png)

**JWST requires Python 3.9 or above and a C compiler for dependencies.**

**Linux and MacOS platforms are tested and supported.  Windows is not currently supported.**

**If installing on MacOS Mojave 10.14, you must install
  into an environment with python 3.9. Installation will fail on python 3.10 due
  to lack of a stable build for dependency ``opencv-python``.**

## Installation

Please contact the [JWST Help Desk](https://jwsthelp.stsci.edu) for installation issues.

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
Remember that all conda operations must be done from within a bash/zsh shell.


### Installing latest releases

You can install the latest released version via `pip`.  From a bash/zsh shell:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install jwst

You can also install a specific version:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install jwst==1.9.4

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
exact versions of all packages to be installed.

To install a particular DMS build, consult the
[Software vs DMS build version map](https://github.com/spacetelescope/jwst#software-vs-dms-build-version-map)
table shown below to determine the correct jwst tag. For example, to install the
version of `jwst` used in DMS build 9.0, use jwst tag 1.8.2. The overall
procedure is similar to the 3-step process outlined in the previous section, but the
details of each command vary, due to the use of environment snapshot files that specify
all of the particular packages to install. Also note that different snapshot files are
used for Linux and Mac OS systems.

Linux:

    conda create -n jwstdp-1.8.2 --file https://ssb.stsci.edu/releases/jwstdp/1.8.2/conda_python_stable-deps.txt
    conda activate jwstdp-1.8.2
    pip install -r https://ssb.stsci.edu/releases/jwstdp/1.8.2/reqs_stable-deps.txt

MacOS:

    conda create -n jwstdp-1.8.2 --file https://ssb.stsci.edu/releases/jwstdp/1.8.2/conda_python_macos-stable-deps.txt
    conda activate jwstdp-1.8.2
    pip install -r https://ssb.stsci.edu/releases/jwstdp/1.8.2/reqs_macos-stable-deps.txt

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
    git clone https://github.com/<your_github_username>/jwst.git
    cd jwst

*Note: `python setup.py install` and `python setup.py develop` commands do not work.*

Install from your local checked-out copy as an "editable" install:

    pip install -e .

If you want to run the unit or regression tests and/or build the docs, you can make
sure those dependencies are installed too:

    pip install -e ".[test]"
    pip install -e ".[docs]"
    pip install -e ".[test,docs]"

Need other useful packages in your development environment?

    pip install ipython jupyter matplotlib pylint


## Calibration References Data System (CRDS) Setup

**Note: As of November 10, 2022, the process of deprecating the CRDS PUB Server will start.
For details, refer to the [CRDS PUB Server Freeze
and Deprecation page](https://jwst-pipeline.readthedocs.io/en/stable/jwst/pub_deprecation.html#pub-deprecation)**


CRDS is the system that manages the reference files needed to run the pipeline.
For details about CRDS, see the [User's
Guide](https://jwst-crds.stsci.edu/static/users_guide/index.html)

The JWST CRDS server is available at  https://jwst-crds.stsci.edu

It supports the automatic processing pipeline at STScI.
Inside the STScI network, the same server is used by the pipeline by default with no modifications.
To run the pipeline outside the STScI network, CRDS must be configured by setting
two environment variables:

    export CRDS_PATH=<locally-accessable-path>/crds_cache/jwst_ops
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu


``<locally-accessable-path>`` can be any the user has permissions to use, such as `$HOME`.
Expect to use upwards of 200GB of disk space to cache the latest couple of contexts.

## Documentation

Documentation (built daily from the Github `master` branch) is available at:

https://jwst-pipeline.readthedocs.io/en/latest/

To build the docs yourself, clone this repository and build the documentation with:

    pip install -e ".[docs]"
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

The table below provides information on each release of the `jwst` package
and its relationship to software builds used in the STScI JWST DMS operations
environment. The `Released` column gives the date on which the `jwst` tag
was released on PyPi and the `Ops Install` column gives the date on which
the build incorporating that release was installed in DMS operations.
Note that the `CRDS_CONTEXT` listed is a minimum context that can be used with
that release. A release should work with any contexts between
the specified context and less than the context for the next release.

| jwst tag            | DMS build | SDP_VER  | CRDS_CONTEXT | Released   | Ops Install | Notes                                         |
|---------------------|-----------|----------|--------------|------------|-------------|-----------------------------------------------|
| 1.12.1              | B10.0rc2  | 2023.3.x | 1132         | 2023-09-26 |             | Second release candidate for B10.0            |
| 1.12.0              | B10.0rc1  | 2023.3.x | 1130         | 2023-09-18 |             | First release candidate for B10.0             |
| 1.11.4              | B9.3.1    | 2023.2.1 | 1107         | 2023-08-14 |             | Final release for B9.3.1 patch                |
| 1.11.3              | B9.3      | 2023.2.0 | 1097         | 2023-07-17 |             | Final release candidate for B9.3              |
| 1.11.2              | B9.3rc3   |          | 1097         | 2023-07-12 |             | Third release candidate for B9.3              |
| 1.11.1              | B9.3rc2   |          | 1094         | 2023-06-29 |             | Second release candidate for B9.3             |
| 1.11.0              | B9.3rc1   |          | 1094         | 2023-06-21 |             | First release candidate for B9.3              |
| 1.10.2              |           |          | 1077         | 2023-04-14 |             | Pinning dependencies for external users       |
| 1.10.1              | B9.2.x    | 2023.1.1 | 1077         | 2023-04-13 | 2023-05-23  | Final release candidate for B9.2              |
| 1.10.0              | B9.2rc1   |          | 1075         | 2023-03-31 |             | First release candidate for B9.2              |
| 1.9.6               | B9.1.2    | 2022.5.2 | 1068         | 2023-03-09 | 2023-03-15  | Final release candidate for B9.1.2            |
| 1.9.5               |           |          | 1061         | 2023-03-02 |             | First release candidate for B9.1.2            |
| 1.9.4               | B9.1.1    | 2022.5.1 | 1041         | 2023-01-27 | 2023-02-28  | Final release candidate for B9.1.1            |
| 1.9.3               | B9.1      | 2022.5.0 | 1030         | 2023-01-12 | 2023-02-28  | Final release candidate for B9.1              |
| 1.9.2               | B9.1rc2   |          |              | 2023-01-04 |             | Second release candidate for B9.1 (hotfix)    |
| 1.9.1               | B9.1rc2   |          |              | 2023-01-03 |             | Second release candidate for B9.1             |
| 1.9.0               | B9.1rc1   |          |              | 2022-12-27 |             | First release candidate for B9.1              |
| 1.8.5               | B9.0      |          | 1019         | 2022-12-12 |             | Documentation patch release for B9.0          |
| 1.8.4               | B9.0      |          |              | 2022-11-16 |             | Documentation patch release for B9.0          |
| 1.8.3               | B9.0      |          |              | 2022-11-11 |             | Documentation patch release for B9.0          |
| 1.8.2               | B9.0      | 2022.4.0 | 1017         | 2022-10-19 | 2022-11-17  | Final release candidate for B9.0              |
| 1.8.1               | B9.0rc2   |          |              | 2022-10-17 |             | Second release candidate for B9.0             |
| 1.8.0               | B9.0rc1   |          |              | 2022-10-10 |             | First release candidate for B9.0              |
| 1.7.2               | B8.1.2    | 2022.3.1 | 0984         | 2022-09-12 | 2022-09-21  | Final release candidate for B8.1.2            |
| 1.7.1               | B8.1.2rc2 |          |              | 2022-09-07 |             | Second release candidate for B8.1.2           |
| 1.7.0               | B8.1.2rc1 |          |              | 2022-09-01 |             | First release candidate for B8.1.2            |
| 1.6.2               | B8.1      | 2022.3.0 | 0953         | 2022-07-19 | 2022-08-19  | Final release candidate for B8.1              |
| 1.6.1               | B8.1rc2   |          |              | 2022-07-15 |             | Second release candidate for B8.1             |
| 1.6.0               | B8.1rc1   |          |              | 2022-07-11 |             | First release candidate for B8.1              |
| 1.5.3               | B8.0.1    | 2022.2.1 | 0913         | 2022-06-20 | 2022-06-30  | Patch release B8.0.1                          |
| 1.5.2               | B8.0      | 2022.2.0 | 0874         | 2022-05-20 | 2022-06-16  | Final release candidate for B8.0              |
| 1.5.1               | B8.0rc2   |          |              | 2022-05-17 |             | Second release candidate for B8.0             |
| 1.5.0               | B8.0rc1   |          |              | 2022-05-05 |             | First release candidate for B8.0              |
| 1.4.6               | B7.9.3    | 2022.1.2 | 0800         | 2022-03-25 |             | Final release candidate for B7.9.3            |
| 1.4.5               | B7.9.3rc2 |          |              | 2022-03-23 |             | Second release candidate for B7.9.3           |
| 1.4.4               | B7.9.3rc1 |          |              | 2022-03-16 |             | First release candidate for B7.9.3            |
| 1.4.3               | B7.9.1    | 2022.1.1 | 0800         | 2022-02-03 |             | Final B7.9.1                                  |
| 1.4.2               | B7.9      | 2022.1.0 | 0797         | 2022-01-20 |             | Final release candidate for B7.9              |
| 1.4.1               | B7.9rc2   |          |              | 2022-01-15 |             | Second release candidate for B7.9             |
| 1.4.0               | B7.9rc1   |          |              | 2022-01-10 |             | First release candidate for B7.9              |
| Pre-launch releases |           |          |              |            |             |                                               |
| 1.3.3               | B7.8.2    | 2021.4.0 | 0764         | 2021-10-05 |             | Same as 1.3.2, but with installation bug fix  |
| 1.3.2               | B7.8.2    | 2021.4.0 | 0764         | 2021-09-03 |             | Final release candidate for B7.8.2            |
| 1.3.1               | B7.8.1    | 2021.3.0 | 0742         | 2021-08-09 |             | Final release candidate for B7.8.1            |
| 1.3.0               | B7.8.1rc1 |          | 0741         | 2021-08-02 |             | First release candidate for B7.8.1            |
| 1.2.3               | B7.8      | 2021.2.0 | 0732         | 2021-06-08 |             | Final release candidate for B7.8              |
| 1.2.2               | B7.8rc3   |          |              | 2021-06-08 |             | Third release candidate for B7.8              |
| 1.2.1               | B7.8rc2   |          |              | 2021-06-07 |             | Second release candidate for B7.8             |
| 1.2.0               | B7.8rc1   |          | 0723         | 2021-05-24 |             | First release candidate for B7.8              |
| 1.1.0               | B7.7.1    | 2021.1.0 | 0682         | 2021-02-26 |             | Final release candidate for B7.7.1            |
| 1.0.0               | B7.7.1rc1 |          | 0678         | 2021-02-22 |             | First release candidate for B7.7.1            |
| 0.18.3              | B7.7      | 2020.4.0 | 0670         | 2021-01-25 |             | Final release candidate for B7.7              |
| 0.18.2              | B7.7rc3   |          | 0668         | 2021-01-19 |             | Third release candidate for B7.7              |
| 0.18.1              | B7.7rc2   |          | 0664         | 2021-01-08 |             | Second release candidate for B7.7             |
| 0.18.0              | B7.7rc1   |          | 0645         | 2020-12-21 |             | First release candidate for B7.7              |
| 0.17.1              | B7.6      | 2020.3.0 | 0641         | 2020-09-15 |             | Final release candidate for B7.6              |
| 0.17.0              | B7.6rc1   |          | 0637         | 2020-08-28 |             | First release candidate for B7.6              |
| 0.16.2              | B7.5      | 2020.2.0 | 0619         | 2020-06-10 |             | Same as 0.16.1, but with installation bug fix |
| 0.16.1              | B7.5      | 2020.2.0 | 0619         | 2020-05-19 |             | Final release candidate for B7.5              |
| 0.16.0              | B7.5rc1   |          | 0614         | 2020-05-04 |             | First release candidate for B7.5              |
| 0.15.1              | B7.4.2    | 2020.1.0 | 0586         | 2020-03-10 |             | Final release candidate for B7.4.2            |
| 0.15.0              | B7.4.2rc1 |          | 0585         | 2020-02-28 |             | First release candidate for B7.4.2            |
| 0.14.2              | B7.4      | 2019.3.0 | 0570         | 2019-11-18 |             | Final release candidate for B7.4              |
| 0.14.1              | B7.4rc2   |          | 0568         | 2019-11-11 |             | Second release candidate for B7.4             |
| 0.14.0              | B7.4rc1   |          | 0563         | 2019-10-25 |             | First release candidate for B7.4              |
| 0.13.8              | B7.3.1    | 2019.2.0 | 0541         | 2019-09-05 |             | Patch for Build 7.3 released as Build 7.3.1   |
| 0.13.7              | B7.3      | 2019.1.0 | 0535         | 2019-06-21 |             | Final release candidate for Build 7.3         |
| 0.13.6              | B7.3rc4   |          | 0534         | 2019-06-20 |             | Fourth release candidate for Build 7.3        |
| 0.13.5              | B7.3rc3   |          | 0534         | 2019-06-19 |             | Third release candidate for Build 7.3         |
| 0.13.4              | B7.3rc2   |          | 0534         | 2019-06-18 |             | Second release candidate for Build 7.3        |
| 0.13.3              | B7.3rc1   |          | 0532         | 2019-06-04 |             | First release candidate for Build 7.3         |
| 0.13.2              |           |          | 0500         | 2019-05-14 |             | DMS test, no delivery to I&T                  |
| 0.13.1              |           |          | 0500         | 2019-03-08 |             | DMS test, no delivery to I&T                  |
| 0.13.0              |           |          | 0500         | 2019-02-15 |             | DMS test, no delivery to I&T                  |
| 0.12.3              | B7.2.1    |          | 0500         | 2019-01-15 |             | DMS Build 7.2.1 patch release                 |
| 0.12.2              | B7.2      | 2018_2   | 0495         | 2018-11-07 |             | Final release candidate for Build 7.2         |
| 0.12.1              | B7.2rc2   |          | 0495         | 2018-11-01 |             | Second release candidate for Build 7.2        |
| 0.12.0              | B7.2rc1   |          | 0493         | 2018-10-09 |             | First release candidate for Build 7.2         |
| 0.11.0              |           |          | 0482         | 2018-09-10 |             | DMS test, no delivery to I&T                  |
| 0.10.0              |           |          | 0477         | 2018-07-31 |             | DMS test, no delivery to I&T                  |
| 0.9.6               | B7.1.3    | 2018_1   | 0468         | 2018-06-08 |             | Final release candidate for Build 7.1.3       |
| 0.9.5               | B7.1.3rc3 |          | 0468         | 2018-06-06 |             | Third release candidate for Build 7.1.3       |
| 0.9.4               | B7.1.3rc2 |          | 0463         | 2018-05-29 |             | Second release candidate for Build 7.1.3      |
| 0.9.3               | B7.1.3rc1 |          | 0457         | 2018-05-11 |             | First release candidate for Build 7.1.3       |
| 0.9.2               |           |          | 0441         | 2018-03-28 |             | DMS test, no delivery to I&T                  |
| 0.9.1               |           |          | 0432         | 2018-02-16 |             | DMS test, no delivery to I&T                  |
| 0.9.0               | B7.1.2    |          | 0422         | 2017-12-22 |             | DMS patch release to I&T 2018-02-15           |
| 0.8.0               | B7.1.1    |          | 0422         | 2017-11-06 |             | DMS patch release to I&T 2018-01-17           |
| 0.8.0               | B7.1      | 2017_1   | 0422         | 2017-11-06 |             | Final release for Build 7.1                   |
| 0.7.7               | B7.0      | 2016_2   | 0303         | 2016-12-13 |             | Final release for Build 7.0                   |


## Unit Tests

Unit tests can be run via `pytest`.  Within the top level of your local `jwst` repo checkout:

    pip install -e ".[test]"
    pytest

Need to parallelize your test runs over all available cores?

    pip install pytest-xdist
    pytest -n auto


## Regression Tests

Latest regression test results can be found here (STScI staff only):

https://plwishmaster.stsci.edu:8081/job/RT/job/JWST/

The test builds start at 6pm local Baltimore time Monday through Saturday on `jwcalibdev`.

To run the regression tests on your local machine, get the test dependencies
and set the environment variable TEST_BIGDATA to our Artifactory server
(STSci staff members only):

    pip install -e ".[test]"
    export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory

To run all the regression tests (except the very slow ones):

    pytest --bigdata jwst/regtest

You can control where the test results are written with the
`--basetemp=<PATH>` arg to `pytest`.  _NOTE that `pytest` will wipe this directory clean
for each test session, so make sure it is a scratch area._

If you would like to run a specific test, find its name or ID and use the `-k` option:

    pytest --bigdata jwst/regtest -k nirspec

If developers need to update the truth files in our nightly regression tests,
there are instructions in the repository wiki.

https://github.com/spacetelescope/jwst/wiki/Maintaining-Regression-Tests
