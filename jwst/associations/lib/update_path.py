"""Update path of members in an association."""

from pathlib import Path


def update_path(asn, file_path, target="expname"):
    """
    Update path of members in an association.

    Parameters
    ----------
    asn : Association
        An association. The association is modified in-place.
    file_path : str
        New path to prepend to each member
    target : str
        Key to replace
    """
    update_key_value(asn, target, (file_path,), mod_func=replace_path)


def update_key_value(obj, target, func_args, mod_func=None):
    """
    Update all instances of a key using a modifier.

    Parameters
    ----------
    obj : object
        Hashable object to modify
    target : str
        The target key to modify
    func_args : (arg(, ...))
        Arguments to pass to the modification function
    mod_func : function
        Function to modify the target key with. If `None`,
        the key will be replaced by the arg list.

    Notes
    -----
    The first argument to `mode_func` will always be the value of the
    target key. Any other arguments given will then be passed to the
    the function.
    """
    if mod_func is None:

        def mod_func(_value, args):
            return args

    if hasattr(obj, "items"):
        for key, value in obj.items():
            if key == target:
                obj[key] = mod_func(value, *func_args)
            if isinstance(value, dict):
                update_key_value(value, target, func_args, mod_func=mod_func)
            elif isinstance(value, list):
                for item in value:
                    update_key_value(item, target, func_args, mod_func=mod_func)


def replace_path(old_path, new_path):
    """
    Replace the path prefix of the basename.

    Parameters
    ----------
    old_path : str or Path
        A string with a file path.
    new_path : str or Path
        The new path to prepend to the basename.

    Returns
    -------
    str
        The basename of the path with the new path prepended.
    """
    return str(Path(new_path) / Path(old_path).name)
