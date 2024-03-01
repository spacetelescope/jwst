import pytest
import os


@pytest.fixture(autouse=True)
def _docdir(request):

    # Trigger ONLY for doctestplus.
    doctest_plugin = request.config.pluginmanager.getplugin("doctestplus")
    if isinstance(request.node.parent, doctest_plugin._doctest_textfile_item_cls):
        tmp_path = request.getfixturevalue('tmp_path')
        old_cwd = os.getcwd()
        os.chdir(tmp_path)
        yield
        os.chdir(old_cwd)
    else:
        yield
