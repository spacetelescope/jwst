from __future__ import absolute_import, unicode_literals, division, print_function

import os.path
from astropy.extern import six

def check(init):
    """
    Determine the type of a file and return it as a string
    
    Parameters
    ----------

    init : file path or file object

    Returns
    -------
    file_type: a string with the file type ("asdf", "asn", or "fits")
    """
    
    if isinstance(init, six.string_types):
        fd = open(init, "rb")
        magic = fd.read(5)
        fd.close()
        
    elif hasattr(init, "read") and hasattr(init, "seek"):
        magic = init.read(5)
        init.seek(0, 0)
        
    else:
        magic = None

    if magic is None or len(magic) < 5:
        raise ValueError("Cannot get file type of " + str(init))

    if magic == b'#ASDF':
        file_type = "asdf"
    elif magic == b'SIMPL':
        file_type = "fits"
    else:
        file_type = "asn"

    return file_type
