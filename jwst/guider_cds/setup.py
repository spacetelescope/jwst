#!/usr/bin/env python
 
from __future__ import division # confidence high


######################################################################
# PyFITS
try:
    import pyfits
except ImportError:
    print "WARNING: PyFITS must be installed to use guider_cds."
    print "         Since this is not a build-time dependency, the "
    print "         build will proceed."


try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

setup(
    setup_requires=['d2to1>=0.2.11', 'stsci.distutils>=0.3.7'],
    namespace_packages=['jwst_pipeline'], packages=['jwst_pipeline'],
    d2to1=True
)
