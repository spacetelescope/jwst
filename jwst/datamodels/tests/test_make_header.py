from jwst.datamodels import make_header


def test_make_header():
    # This is mostly a "smoke test" to determine the code works

    instrument = 'nircam'
    for mode in ('image', 'spectrum'):
        for level in (1, 2, 3):
            im = make_header.main(instrument, mode, level)
            assert im.meta.model_type == 'ImageModel'
            assert im.meta.instrument.filter == 'F115W'
            assert im.meta.instrument.name == instrument.upper()
