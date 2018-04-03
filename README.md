JWST Calibration Pipeline
=========================
[![Documentation Status](https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest)](http://jwst-pipeline.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/STScI-JWST/jwst.svg?branch=master)](https://travis-ci.org/STScI-JWST/jwst)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

![STScI Logo](docs/_static/stsci_logo.png)

Note
----
Beginning with the next version (0.9.0), **JWST requires Python 3.5 and above**.

Installing
----------

To install the latest stable version of the library, we recommend using [conda](https://conda.io/docs/index.html) and
the [astroconda](https://astroconda.readthedocs.io) channel:

    % conda config --add channels http://ssb.stsci.edu/astroconda-dev
    % conda install jwst

To install the development version of the repository, we recommend using the [astroconda](https://astroconda.readthedocs.io) channel
to install the dependencies, and then installing from the github repository:

    % conda config --add channels http://ssb.stsci.edu/astroconda-dev

    % conda install jwst

    % git clone https://github.com/STScI-JWST/jwst.git

    % cd jwst

    % python setup.py install

    or

    % python setup.py develop

Once installed, the software can be updated to the lastest development version by running the following command inside the `jwst` 
repository:

    % git pull origin master


Contributing Code, Documentation or Feedback
--------------------------------------------
We welcome feedback and contributions to the project. Please follow the [contributing guidelines](CONTRIBUTING.md) to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by abiding to the [Code of Conduct](CODE_OF_CONDUCT.md).

Using
-----

Documentation (latest off `master`) is available here:

https://jwst-pipeline.readthedocs.io/en/latest/


Software vs DMS build version map
---------------------------------

| jwst tag | DMS build |    Date    |          Notes               |
| -------- | --------- | ---------- | ---------------------------- |
|  0.9.2   |   B7.1x   | 03/28/2018 | DMS test, no delivery to I&T |
|  0.9.1   |   B7.1x   | 02/16/2018 | DMS test, no delivery to I&T |
|  0.9.0   |   B7.1x   | 12/22/2017 | DMS test, no delivery to I&T |
|  0.8.0   |   B7.1    | 11/14/2017 | Final, delivered to I&T      |
|  0.7.0rc7|   B7.0    | 12/13/2016 | Final, delivered to I&T      |
