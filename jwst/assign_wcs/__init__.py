from __future__ import absolute_import
from .assign_wcs_step import AssignWcsStep
#from .assign_wcs import nrs_set_inputs


def test( verbose=False ) :
    import os, pytest

    # get the pandokia plugin if it is available (it will only
    # do anything if we are run from pandokia).
    try :
        import pandokia.helpers.pytest_plugin as pytest_plugin
    except ImportError :
        pytest_plugin = None

    if pytest_plugin :
        addplugins = [ pytest_plugin]
    else :
        addplugins = None

    from . import tests
    dir = os.path.dirname( tests.__file__ )

    args=[dir]

    return pytest.main( args)
