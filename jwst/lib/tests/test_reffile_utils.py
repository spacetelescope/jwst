from jwst.lib.reffile_utils import find_row


def test_find_row():
    filters = [
        {'column_offset': 1.0, 'filter': 'F277W', 'pupil': 'FLAT', 'row_offset': 2.0},
        {'column_offset': 0.0, 'filter': 'F356W', 'pupil': 'FLAT', 'row_offset': 0.0},
    ]
    match_keys = {'filter': 'F277W', 'pupil': 'FLAT'}
    missing_key = {'filter': 'F277H', 'pupil': 'FLAT'}

    result = find_row(filters, match_keys)
    assert result == {'column_offset': 1.0, 'filter': 'F277W', 'pupil': 'FLAT', 'row_offset': 2.0}

    result = find_row(filters, missing_key)
    assert result is None
