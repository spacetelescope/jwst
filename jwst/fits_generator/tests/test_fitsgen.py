import os

from jwst.fits_generator import create_data

test_directory = os.path.dirname(__file__)

dotfile_directory = test_directory + '/dotfile/'
okfile_directory = test_directory + '/okfile/'

def test_get_proposals():
    """Test the get_proposals function. """ 
    dotfile_list = create_data.get_proposals(dotfile_directory)
    assert dotfile_list == []

    okfile_list = create_data.get_proposals(okfile_directory)
    assert okfile_list == [(okfile_directory, '05181.prop')]

    return

