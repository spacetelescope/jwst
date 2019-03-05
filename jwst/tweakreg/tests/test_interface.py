def test_tweakreg_chelp_importable():
    """Test that C extension for arrxyzero can be imported"""
    from jwst.tweakreg import chelp

    assert chelp.arrxyzero
