import pytest

@pytest.fixture(autouse=True)
def _docdir(request):

    # Trigger ONLY for doctestplus.
    doctest_plugin = request.config.pluginmanager.getplugin("doctestplus")
    if isinstance(request.node.parent, doctest_plugin._doctest_textfile_item_cls):
        tmpdir = request.getfixturevalue('tmpdir')
        with tmpdir.as_cwd():
            yield
    else:
        yield
