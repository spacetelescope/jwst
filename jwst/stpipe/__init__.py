from __future__ import absolute_import, division, print_function

import sys
if sys.version_info[0] >= 3:
    __builtins__['unicode'] = str
    __builtins__['basestring'] = str

from .step import Step
from .pipeline import Pipeline
from .linear_pipeline import LinearPipeline

__version__ = '0.8.0'

#def test(verbose=False):
    #import nose

    ## get the pandokia plugin if it is available (it will only
    ## do anything if we are run from pandokia).
    #try:
        #import pandokia.helpers.nose_plugin as nose_plugin
    #except ImportError:
        #nose_plugin = None

    #if nose_plugin:
        #addplugins = [nose_plugin.Pdk()]
    #else:
        #addplugins = None

    ## get the name of the test package
    #argv = ['nosetests', '--exe', __name__ + '.tests']

    ## run nose
    #return nose.main(argv=argv, addplugins=addplugins)
