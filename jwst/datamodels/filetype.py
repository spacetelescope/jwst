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
    if isinstance(init, str):
        exts = init.split(os.extsep)  # Will support multiple extensions, for example in cases of zipped files.
        exts.pop(0)  # Don't need the filename

        if not len(exts):
            raise ValueError(f'Input file path does not have an extension: {init}')

        if exts[0] not in ('fits', 'json', 'asdf'):
            raise ValueError(f'Unrecognized file type: {exts[0]}')

        if exts[0] == 'json':
            return 'asn'

        return exts[0]

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
