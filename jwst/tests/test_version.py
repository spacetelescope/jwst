def test_version_is_string():
    from jwst import __version__
    assert isinstance(__version__, str)


def test_version_commit_is_string():
    from jwst import __version_commit__
    assert isinstance(__version_commit__, str)
