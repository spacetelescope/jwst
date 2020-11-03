"""asn_gather: Copy or Move data that is listed in an association"""
from pathlib import Path

__all__ = ['asn_gather']


def asn_gather(source_asn_path, destination=None, exp_types=None, copy=True, recurse=False, member_root=None):
    """Copy/Move members of an association from one location to another

    The association is copied into the destination, re-written such that the member
    list points to the new location of the members. If `copy` is `False`, the member
    list will simply point back to the original location of the members.

    If members cannot be copied, warnings will be generated and the member information
    in the copied association will be left untouched.

    Parameters
    ----------
    source_asn_path : str, pathlib.Path, Association, or dict-like
        The association to gather.

    destination : str, pathlib.Path, or None
        The folder to place the association and its members.
        If None, the current working directory is used.

    exp_types : [str[,...]] or None
        List of exposure types to gather.
        If None, all are copied.

    copy : bool
        Copy the members to the destination. Otherwise, just copy
        the association itself with the members pointing to the original
        location.

    recurse : bool
        If members appear to come from other associations, attempt
        to gather.

    member_root : str, pathlib.Path, or None
        The folder where the members are found.
        If None, the folder path of the requested association is used.

    Returns
    -------
    dest_asn : pathlib.Path
        The association.
    """
    from .load_as_asn import LoadAsAssociation

    if destination is None:
        dest_folder = Path('./')
    else:
        dest_folder = Path(destination)

    # Open the source association
    dest_asn = LoadAsAssociation.load(source_asn_path)

    # Copy the members

    # Create the new association.

    # Save new association.
    dest_path = dest_folder / source_asn_path.name
    _, serialized = dest_asn.dump()
    with open(dest_path, 'w') as fh:
        fh.write(serialized)

    # That's all folks
    return dest_path


def from_cmdline():
    """Collect asn_gather arguments from the commandline"""
