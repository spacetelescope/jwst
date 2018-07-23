JWST Calibration Pipeline
=========================
[![Documentation Status](https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest)](http://jwst-pipeline.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/spacetelescope/jwst.svg?branch=master)](https://travis-ci.org/spacetelescope/jwst)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

![STScI Logo](docs/_static/stsci_logo.png)

Note
----
Beginning with version 0.9.0, **JWST requires Python 3.5 or above**.

Installing
----------

To install a particular released version of the package, and all dependencies, we recommend using
[conda](https://conda.io/docs/index.html) and a spec file that lists the exact versions of all packages to be installed.
To create a new environment, use:

    conda create -n jwst --file <URL>
    source activate jwst

where `<URL>` is of the form:

    Linux: http://ssb.stsci.edu/releases/jwstdp/0.9.6/latest-linux
    OS X: http://ssb.stsci.edu/releases/jwstdp/0.9.6/latest-osx

Other particular versions can be installed by choosing a different version tag in place of "0.9.6" in the URL path.
See the "Software vs DMS build version map" table below for a list of tags corresponding to particular releases.

To update to the latest nightly build:

    conda update -n jwst --override-channels -c http://ssb.stsci.edu/astroconda-dev -c defaults --all

To install the development version of the repository, we recommend creating a new
environment, using the [astroconda](https://astroconda.readthedocs.io) channel
to install the dependencies, and then installing from the github repository:

    conda create -n jwst_dev --only-deps --override-channels -c http://ssb.stsci.edu/astroconda-dev -c defaults python=3.6 jwst
    source activate jwst_dev
    git clone https://github.com/spacetelescope/jwst.git
    cd jwst
    python setup.py develop

Once installed, the software can be updated to the lastest development version by updating the dependencies,
pulling the latest version of `master` from the Github repository inside the `jwst` directory:

    conda update -n jwst_dev --override-channels -c http://ssb.stsci.edu/astroconda-dev -c defaults --all
    git pull origin master
    python setup.py develop


Documentation
-------------

Documentation (built daily from `master`) is available here:

https://jwst-pipeline.readthedocs.io/en/latest/


Contributions and Feedback
--------------------------
We welcome contributions and feedback on the project. Please follow the [contributing guidelines](CONTRIBUTING.md) to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by abiding to the [Code of Conduct](CODE_OF_CONDUCT.md).


Software vs DMS build version map
---------------------------------

| jwst tag | DMS build |    Date    |          Notes                           |
| -------- | --------- | ---------- | -----------------------------------------|
|  0.9.6   | B7.1.3    | 06/08/2018 | Final release candidate for Build 7.1.3  |
|  0.9.5   | B7.1.3rc3 | 06/06/2018 | Third release candidate for Build 7.1.3  |
|  0.9.4   | B7.1.3rc2 | 05/29/2018 | Second release candidate for Build 7.1.3 |
|  0.9.3   | B7.1.3rc1 | 05/11/2018 | First release candidate for Build 7.1.3  |
|  0.9.2   |           | 03/28/2018 | DMS test, no delivery to I&T             |
|  0.9.1   |           | 02/16/2018 | DMS test, no delivery to I&T             |
|  0.9.0   |   B7.1.2  | 12/22/2017 | DMS patch release to I&T 02/15/2018      |
|  0.8.0   |   B7.1.1  |            | DMS patch release to I&T 01/17/2018      |
|  0.8.0   |   B7.1    | 11/14/2017 | Final, delivered to I&T 11/17/2017       |
|  0.7.0rc7|   B7.0    | 12/13/2016 | Final, delivered to I&T                  |


Unit Tests
----------

Unit tests can be run via `pytest`.  We recommend using `pytest-xdist` so you can run them in parallel.  Install `pytest-xdist` and run pytest in the top level of the repository

    conda install pytest-xdist
    pytest -n <cores>

where `cores` is the number of cores you'd like to use on your machine for the tests.

Regression Tests
----------------

Regression tests - both the data and the result reports - are currently only accessible to STScI staff members. If you do need information about this, please open an issue.

Latest regression test results can be found here:

https://boyle.stsci.edu:8081/job/RT/job/JWST/

The test builds start at 11am and 6pm local Baltimore time every day on jwcalibdev.

To run the regression tests on your local machine, `rsync` or `scp` the input and comparison data locally

    rsync -av <username>@jwcalibdev:/data4/jwst_test_data /my/local/path/

set the environment variable TEST_BIGDATA to this location

    export TEST_BIGDATA=/my/local/path/jwst_test_data/

and then run the tests in the repository

    cd /path/to/jwst/jwst/tests_nightly/general
    pytest --bigdata .


