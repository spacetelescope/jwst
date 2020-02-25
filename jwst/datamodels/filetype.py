import os


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
    supported = ('asdf', 'fits', 'json')

    if isinstance(init, str):
        path, ext = os.path.splitext(init)
        ext = ext.strip('.')

        if not ext:
            raise ValueError(f'Input file path does not have an extension: {init}')

        if ext not in supported:  # Could be the file is zipped; try splitting again
            path, ext = os.path.splitext(path)
            ext = ext.strip('.')

            if ext not in supported:
                raise ValueError(f'Unrecognized file type for: {init}')

        if ext == 'json':
            return 'asn'

        return ext

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
