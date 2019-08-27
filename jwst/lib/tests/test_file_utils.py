"""Test file utilities"""
import os
import pytest

from jwst.lib.file_utils import pushdir


def test_pushdir():
    """Test for successful change"""
    # Retrieve what the /tmp folder really is
    current = os.getcwd()
    try:
        os.chdir('/tmp')
    except Exception:
        pytest.xfail('Cannot change to the tmp directory. Test cannot run')
    tmp_dir = os.getcwd()
    os.chdir(current)

    # Now on with the test.
    with pushdir(tmp_dir):
        assert tmp_dir == os.getcwd(), 'New directory is not what was specified.'
    assert current == os.getcwd(), 'Directory was not restored'


def test_pusdir_fail():
    """Test for failing changing."""
    current = os.getcwd()
    with pytest.raises(Exception):
        with pushdir('Really_doesNOT-exist'):
            # Nothing should happen here. The assert should never be checked.
            assert False
    assert current == os.getcwd()
