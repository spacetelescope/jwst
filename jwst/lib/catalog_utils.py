"""
Utilities for naming source catalogs.
"""

import re
from os.path import split, splitext, join, abspath, expanduser


def replace_suffix_ext(filepath, old_suffix_list, new_suffix,
                       output_ext='escv', output_dir=None):
    """
    Replace the suffix and extension of a filepath.
    """

    path, filename = split(filepath)
    name, ext = splitext(filename)
    remove_suffix = '^(.+?)(_(' + '|'.join(old_suffix_list) + '))?$'
    match = re.match(remove_suffix, name)
    name = match.group(1)

    output_path = '{0}_{1}.{2}'.format(name, new_suffix, output_ext)
    if output_dir is not None:
        output_path = abspath(expanduser(join(output_dir, output_path)))

    return output_path
